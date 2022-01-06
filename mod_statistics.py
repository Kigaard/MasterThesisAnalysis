import pathlib

import pandas as pd

from analysis_code import calculate_modification_percentages, create_modification_barplot, \
    combine_spectra_in_peptide_lists


def perform_analysis(peptide_list_directory: pathlib.Path, modifications: list[tuple[str, str, float]],
                     condition_title: str):
    for mod_name, res, mod_mass in modifications:
        print("*" * 5, f"{mod_name} ({res}@{mod_mass})", "*" * 5)
        mod_dict: dict = calculate_modification_percentages(peptide_list_directory=peptide_list_directory,
                                                            residue_str=res,
                                                            mod_mass=mod_mass)
        mod_df: pd.DataFrame = pd.DataFrame.from_dict(mod_dict, orient="index")
        # Rename and reindex columns
        mod_df.index = ["Native CRT", "Native Lacto", "Native Ribo"]

        print(mod_df)
        print()

        # create_modification_barplot(data=mod_df, mod_name=mod_name, residues=res, condition_title=condition_title)


def main():
    base_directory = pathlib.Path(r"C:\Users\spec-makie17\Documents\Experiments\211215_Bait_Take3_PH2121_ONGO\MGF")
    raw_peptide_list_folder = pathlib.Path("List")
    peptide_list_folder = pathlib.Path(r"Graph")

    if True:
        list_dir = base_directory / raw_peptide_list_folder
        save_dir = base_directory / peptide_list_folder
        combine_spectra_in_peptide_lists(list_directory=list_dir, save_directory=save_dir)

    modifications_cbm_cysteine: list[tuple[str, str, float]] = [
        ("Oxidation", "MP", 15.995),
        ("Carbamidomethyl", "C", 57.022),
        ("Dehydroalanine", "C", -87.986),
        ("Oxidation", "C", -41.027),
        ("Sulfinic", "C", -25.032),
        ("SulfDiOx", "C", 6.940),
        ("SulfOx", "C", -9.054),
        ("Sulfonic", "C", -9.037),
        ("SSulfonic", "C", 22.935),
        ("SO3", "C", 34.935)]

    modifications_free_cysteine: list[tuple[str, str, float]] = [
        ("Carbamidomethyl", "C", 57.022),
        ("Dehydroalanine", "C", 0),
        ("Oxidation", "C", 0),
        ("Sulfinic", "C", 0),
        ("SulfDiOx", "C", 0),
        ("SulfOx", "C", 0),
        ("Sulfonic", "C", 0),
        ("SSulfonic", "C", 0),
        ("SO3", "C", 0)]

    modifications_cbm_cysteine_oxygen18: list[tuple[str, str, float]] = [
        ("Oxidation", "MP", 17.999),
        ("Carbamidomethyl", "C", 57.022),
        ("Dehydroalanine", "C", -87.986),
        ("Oxidation", "C", -39.935),
        ("Sulfinic", "C", -21.023),
        ("SulfDiOx", "C", 10.949),
        ("SulfOx", "C", -7.050),
        ("Sulfonic", "C", -3.024),
        ("SSulfonic", "C", 28.948),
        ("SO3", "C", 40.948)
    ]

    modification_oxygen18: list[tuple[str, str, float]] = [("Oxidation", "MP", 17.999),
                                                           ("Oxidation", "C", -39.935),
                                                           ("Sulfinic", "C", -21.023),
                                                           ("SulfDiOx", "C", 10.949),
                                                           ("SulfOx", "C", -7.050),
                                                           ("Sulfonic", "C", -3.024),
                                                           ("SSulfonic", "C", 28.948),
                                                           ("SO3", "C", 40.948)]

    perform_analysis(peptide_list_directory=base_directory / peptide_list_folder,
                     modifications=modifications_cbm_cysteine, condition_title="")


if __name__ == "__main__":
    main()
