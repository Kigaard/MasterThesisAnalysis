from typing import List

import pandas as pd

from analysis_code import calculate_modification_percentages, create_modification_barplot,\
    combine_spectra_in_peptide_lists


def perform_analysis(peptide_list_directory: str, modifications: list[tuple[str, str, float]], condition_title: str):
    for mod_name, res, mod_mass in modifications:
        print("*" * 5, f"{mod_name} ({res}@{mod_mass})", "*" * 5)
        mod_dict: dict = calculate_modification_percentages(peptide_list_directory=peptide_list_directory, residues=res,
                                                            mod_mass=mod_mass)
        mod_df: pd.DataFrame = pd.DataFrame.from_dict(mod_dict, orient="index")
        # Rename and reindex columns
        mod_df.index = ["1 day 37 °C (Autodigestion)", "1 day 37 °C (Trypsin)", "1 day 42 °C (Autodigestion)",
                        "1 day 42 °C (Trypsin)", "4 day 37 °C (Autodigestion)", "4 day 37 °C (Trypsin)",
                        "4 day 42 °C (Autodigestion)", "4 day 42 °C (Trypsin)", "Day 0 (Autodigestion)",
                        "Day 0 (Trypsin)"]
        mod_df = mod_df.reindex(["Day 0 (Autodigestion)", "Day 0 (Trypsin)", "1 day 37 °C (Autodigestion)",
                                 "1 day 37 °C (Trypsin)", "1 day 42 °C (Autodigestion)", "1 day 42 °C (Trypsin)",
                                 "4 day 37 °C (Autodigestion)", "4 day 37 °C (Trypsin)", "4 day 42 °C (Autodigestion)",
                                 "4 day 42 °C (Trypsin)"])
        print(mod_df)
        print()

        create_modification_barplot(data=mod_df, mod_name=mod_name, residues=res, condition_title=condition_title)


def main():
    raw_peptide_list_dir = \
        r"C:\Users\spec-makie17\Documents\Experiments\211210_MDH_Oxidation_PH2118_ONGO\PeptideLists"
    peptide_list_dir = \
        r"C:\Users\spec-makie17\Documents\Experiments\211210_MDH_Oxidation_PH2118_ONGO\PeptideLists_Graphs\CRT-MDH_MDH"

    combine_spectra_in_peptide_lists()

    mods: List[tuple[str, str, float]] = [('Oxidation', 'MP', 15.995)]
    #perform_analysis(peptide_list_directory=parsed_peptide_list_directory, modifications=mods,
    #                 condition_title="MDH in combined CRT + MDH")

    # # read_and_save_peptide_lists(list_directory=peptide_list_directory, save_directory=parsed_peptide_list_directory)
    # oxidation_dict: dict = calculate_oxidation_percentages(peptide_list_directory=parsed_peptide_list_directory,
    #                                                        residues="MP", mod_mass=15.995)
    # data: pd.DataFrame = pd.DataFrame.from_dict(oxidation_dict, orient="index")
    # print(data)
    #
    # # data = data.reindex(["Control", "Acetyl-L-Cysteine", "Hydroquinone", "Hydroxyurea"])
    #
    # create_plot(data)


if __name__ == '__main__':
    main()
