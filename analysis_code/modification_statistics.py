import pathlib
import re
from collections import defaultdict

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.axes import Axes
from tqdm import tqdm

from .utils import get_residue_name


def combine_spectra_in_peptide_lists(list_directory: pathlib.Path, save_directory: pathlib.Path) -> None:
    """
    Combine the peptides in the peptide list and save it in a new directory.

    :param list_directory: The directory with the peptide lists.
    :param save_directory: The directory where the peptide list should be saved.
    :return:
    """
    files_list = [(file, pathlib.Path(save_directory, f"{file.name}")) for file in list_directory.glob("*.xlsx")
                  if file.is_file()]

    for input_file, output_file in tqdm(files_list):
        # Read the file
        df: pd.DataFrame = pd.read_excel(io=input_file, usecols=["from", "to", "seq", "modifs", "#"])
        df.columns = ["Start", "End", "Sequence", "Modification", "Spectra"]
        # Remove Q/E loss
        for i, row in df.iterrows():
            if row["Modification"] == "-":
                continue
            modification_list = []
            # Loop over the modification
            for mod in row["Modification"].split(" "):
                residue_id = int(mod.split("@")[0])
                residue = str(row["Sequence"])[residue_id - row["Start"]]
                if residue not in ["Q", "E"]:
                    modification_list.append(mod)

            if len(modification_list) == 0:
                df.at[i, "Modification"] = "-"
            else:
                df.at[i, "Modification"] = ";".join(modification_list)
        # Sum the spectra (Assume that spectra with only Q/E modifications are the same as non-modified spectra)
        df["Spectra"] = df.groupby(["Sequence", "Modification"])["Spectra"].transform("sum")
        df = df.drop_duplicates(subset=["Sequence", "Modification"]).reset_index(drop=True)
        df.to_excel(output_file)


def calculate_modification_percentages(peptide_list_directory: pathlib.Path, sequence: str,
                                       residue_str: str, mod_mass: float) -> dict:
    """
    Calculate the modification percentages for each of the given lists in the directory.

    :param peptide_list_directory: The directory with the peptide lists.
    :param sequence: The protein sequence.
    :param residue_str: The residues to calculate modifications for as a string.
    :param mod_mass: The modification mass.
    :return: The modification percentage and counts for each of the given lists.
    """
    directory = pathlib.Path(peptide_list_directory)
    files = [(file.stem, file) for file in directory.glob("*.xlsx") if file.is_file()]
    residues = [match.start() + 1 for match in re.finditer(f"[{residue_str}]", sequence.upper())]
    mod_dict: defaultdict = defaultdict(lambda: {"Percentage": 0, "ModifiedSpectra": 0, "TotalModSpectra": 0})

    for condition_name, file in files:
        oxidized_spectra_count: int = 0
        total_oxidation_spectra_count: int = 0

        df: pd.DataFrame = pd.read_excel(file, index_col=0)
        for _, row in df.iterrows():
            if any([res for res in residue_str if res in row["Sequence"]]):
                total_oxidation_spectra_count += row["Spectra"]
                if not not re.search(f"[{'|'.join([str(res) for res in residues])}]@{mod_mass}", row["Modification"]):
                    oxidized_spectra_count += row["Spectra"]

        if total_oxidation_spectra_count == 0:
            mod_dict[condition_name]["Percentage"] = 0
        else:
            mod_dict[condition_name]["Percentage"] = round((oxidized_spectra_count / total_oxidation_spectra_count)
                                                           * 100, 2)
        mod_dict[condition_name]["ModifiedSpectra"] = oxidized_spectra_count
        mod_dict[condition_name]["TotalModSpectra"] = total_oxidation_spectra_count

    return mod_dict


def create_modification_barplot(data: pd.DataFrame, mod_name: str, residues: str, condition_title: str) -> None:
    """
    Create the modification bar plot.

    :param data: The data.
    :param mod_name: The modification name.
    :param residues: The residues to show.
    :param condition_title: The condition information to show in the plot title.
    """
    chart: Axes = sns.barplot(data=data, x=data.index, y="Percentage")
    chart.set_title(f"Total {mod_name} of {' and '.join(get_residue_name(res) for res in residues)}")
    chart.set_xlabel("Condition")
    chart.set_ylabel("Percentage\n(Spectra count/Total spectra count)")

    for idx, p in enumerate(chart.patches):
        chart.annotate(f"{p.get_height()} %\n"
                       f"({int(data.iloc[idx]['ModifiedSpectra'])}/{int(data.iloc[idx]['TotalModSpectra'])})",
                       (p.get_x() + p.get_width() / 2., p.get_height()),
                       ha='center', va='center', fontsize=12, color='black', xytext=(0, 12),
                       textcoords='offset points')

    plt.show()
