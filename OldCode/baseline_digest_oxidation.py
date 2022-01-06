import re
from typing import List, Tuple, Dict

from matplotlib.axes import Axes
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from analysis_code.utils import get_residue_positions


def read_peptide_list(list_file_name: str) -> pd.DataFrame:
    # Read the file
    df: pd.DataFrame = pd.read_excel(io=list_file_name, usecols=["from", "to", "seq", "modifs", "#"])
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
                modification_list.append(f"{residue}{mod}")

        if len(modification_list) == 0:
            df.at[i, "Modification"] = "-"
        else:
            df.at[i, "Modification"] = ";".join(modification_list)
    # Sum the spectra (Assume that spectra with only Q/E modifications are the same as non-modified spectra)
    df["Spectra"] = df.groupby(["Sequence", "Modification"])["Spectra"].transform("sum")
    df = df.drop_duplicates(subset=["Sequence", "Modification"]).reset_index(drop=True)

    return df


def get_positions_and_total_spectra(peptide_list: pd.DataFrame, residues: str) \
        -> Tuple[List[Tuple[int, str]], Dict[str, int]]:
    """
    Get the relevant positions and the spectra count per position.

    :param peptide_list: The peptide list.
    :param residues: The residues
    :return: The list of positions with number and residue and the dictionary with the position counts.
    """
    residues: List[Tuple[int, str]] = get_residue_positions(residues=residues)
    residues = [(res_num, res) for res_num, res in residues if res_num > 17]
    total_position_spectra: Dict[str, int] = {r: 0 for _, r in residues}

    for idx, row in peptide_list.iterrows():
        start_position: int = row["Start"]
        end_position: int = row["End"]
        spectra: int = row["Spectra"]

        for res_num, res in residues:
            if start_position <= res_num <= end_position:
                total_position_spectra[res] += spectra

    return residues, total_position_spectra


def analyse_methionine_proline_cysteine_oxidation(peptide_list: pd.DataFrame):
    """
    Analyse the methionine, proline, cysteine (All modifications) oxidation.

    :param peptide_list: The peptide list.
    :return:
    """
    residues, total_position_spectra = get_positions_and_total_spectra(peptide_list=peptide_list, residues="MPC")
    oxidation_spectra_count: Dict[str, int] = {k: 0 for k in total_position_spectra.keys()}

    for idx, row in peptide_list.iterrows():
        for match in re.finditer("([MPC][0-9]{1,3})", row["Modification"]):
            oxidation_spectra_count[match.group(1)] += row["Spectra"]

    percentage_spectra_count: Dict[str, float] = {
        k: [(round(oxidation_spectra_count[k] / float(total_position_spectra[k]) * 100,
                   2) if total_position_spectra[k] != 0 else 0.0)] for k
        in total_position_spectra.keys()}

    data: pd.DataFrame = pd.DataFrame.from_dict(percentage_spectra_count, orient="index", columns=["Percentage"])

    chart: Axes = sns.barplot(data=data, x=data.index, y="Percentage")
    chart.set_title("Baseline modification/oxidiation of methionine, cysteine, and proline")
    chart.text(x=26, y=129, s="*) Cysteines are a summation of multiple modifications")
    chart.set_xlabel("Position")
    chart.set_ylabel("Percentage\n(Spectra count/Total spectra count)")


    for idx, p in enumerate(chart.patches):
        cys_warning = list(oxidation_spectra_count.keys())[idx].startswith("C")
        ox_spec_count = list(oxidation_spectra_count.values())[idx]
        total_spec_count = list(total_position_spectra.values())[idx]
        chart.annotate(f"{p.get_height()} %\n({ox_spec_count}/{total_spec_count}){'*' if  cys_warning else ''}",
                       (p.get_x() + p.get_width() / 2., p.get_height()),
                       ha='center', va='center', fontsize=10, color='black', xytext=(0, 10),
                       textcoords='offset points')
    plt.show()


def main():
    peptide_list_file: str = \
        r"C:\Users\spec-makie17\Documents\Experiments\211027_Baseline_Digest_Oxidation_PH2120_ONGO\PeptideLists\EXP3_1000271_NKH_rCRT_JC_0_auto.xlsx"

    peptide_list: pd.DataFrame = read_peptide_list(peptide_list_file)
    analyse_methionine_proline_cysteine_oxidation(peptide_list=peptide_list)

    # peptide_list.to_excel(r"C:\Users\spec-makie17\Documents\Experiments\211027_Baseline_Digest_Oxidation_PH2120_ONGO\Data.xlsx")


if __name__ == "__main__":
    main()
