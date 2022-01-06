import json
from typing import List, Tuple

import pandas as pd
import os

from tqdm import tqdm

from analysis_code import utils


def read_and_save_peptide_lists(list_directory: str, save_directory: str):
    """
    Read the peptide lists.

    :param list_directory: The directory with the peptide lists.
    :param save_directory: The directory where the peptide list should be saved.
    :return:
    """
    files_list = [(os.path.join(list_directory, file), os.path.join(save_directory, file.split("_rCRT_")[1]))
                  for file in os.listdir(list_directory) if file.endswith(".xlsx")]

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


def combine_peptide_lists(directory: str):
    full_files_list: List[tuple, pd.DataFrame] = [(file.split(".xlsx")[0].split("_"),
                                                   pd.read_excel(os.path.join(directory, file), index_col=0))
                                                  for file in os.listdir(directory) if file.endswith(".xlsx")]
    dfs = []
    for info, df in full_files_list:
        print(info, info[2:], len(info[2:]))
        df = df.assign(crt_type=info[0], Time=info[1], Digestion="Autodigestion" if info[2] == "auto" else "Trypsin",
                       mod_oxy="18O" if len(info[2:]) == 3 else "16O")
        df.rename(columns={"crt_type": "CRT Type", "mod_oxy": "Oxygen modification"}, inplace=True)
        dfs.append(df)

    full_df: pd.DataFrame = pd.concat(dfs)
    full_df.to_excel(os.path.join(directory, "Combined.xlsx"))
    return full_df


def get_total_modification_statistics(df: pd.DataFrame, residues: str):
    residues: List[Tuple[int, str]] = utils.get_residue_positions(residues=residues)
    modifications_16O: dict = {
        -33.988: "Dehydro",
        15.995: "Oxidation",
        31.990: "Sulfinic",
        63.962: "SulfDiOx",
        47.967: "SulfOx",
        47.985: "Sulfonic",
        79.957: "SSulfonic",
        91.957: "SO3"
    }
    modifications_18O: dict = {
        -33.988: "Dehydro",
        17.999: "Oxidation",
        35.998: "Sulfinic",
        67.970: "SulfDiOx",
        49.971: "SulfOx",
        53.997: "Sulfonic",
        85.970: "SSulfonic",
        97.970: "SO3"
    }

    result_dict: dict = {}

    for i, row in df.iterrows():
        if row["CRT Type"] not in result_dict:
            result_dict[row["CRT Type"]] = {}
        if row["Time"] not in result_dict[row["CRT Type"]]:
            result_dict[row["CRT Type"]][row["Time"]] = {}
        if row["Digestion"] not in result_dict[row["CRT Type"]][row["Time"]]:
            result_dict[row["CRT Type"]][row["Time"]][row["Digestion"]] = {}
        if row["Oxygen modification"] not in result_dict[row["CRT Type"]][row["Time"]][row["Digestion"]]:
            result_dict[row["CRT Type"]][row["Time"]][row["Digestion"]][row["Oxygen modification"]] = {}

    print(json.dumps(result_dict, indent=4))


def main():
    peptide_list_directory = \
        r"C:\Users\spec-makie17\Documents\Experiments\211027_Baseline_Digest_Oxidation_PH2120_ONGO\PeptideLists"
    parsed_peptide_list_directory = \
        r"C:\Users\spec-makie17\Documents\Experiments\211027_Baseline_Digest_Oxidation_PH2120_ONGO\PeptideLists_Graphs"

    read_and_save_peptide_lists(list_directory=peptide_list_directory, save_directory=parsed_peptide_list_directory)
    combined_df: pd.DataFrame = combine_peptide_lists(directory=parsed_peptide_list_directory)

    get_total_modification_statistics(df=combined_df, residues="C")


if __name__ == '__main__':
    main()
