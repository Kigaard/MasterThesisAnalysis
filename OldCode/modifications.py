import os
from typing import List

import pandas as pd
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt


def _read_data(file: str, pass_percentage: int) -> pd.DataFrame:
    """
    Read and filter data.

    :param file: The input file.
    :param pass_percentage: The passing modification percentage.
    :return: The filtered dataframe.
    """
    modifications_df: pd.DataFrame = pd.read_csv(filepath_or_buffer=file,
                                                 usecols=["Protein Position", "Modifications",
                                                          "Crtl_37_Tryp modified", "Crtl_37_Tryp unmodified",
                                                          "Crtl_37_TrypChym modified",
                                                          "Crtl_37_TrypChym unmodified", "Crtl_42_Tryp modified",
                                                          "Crtl_42_Tryp unmodified", "Crtl_42_TrypChym modified",
                                                          "Crtl_42_TrypChym unmodified",
                                                          "CompleteEDTA_37_Tryp modified",
                                                          "CompleteEDTA_37_Tryp unmodified",
                                                          "CompleteEDTA_37_TrypChym modified",
                                                          "CompleteEDTA_37_TrypChym unmodified",
                                                          "CompleteEDTA_42_Tryp modified",
                                                          'CompleteEDTA_42_Tryp unmodified',
                                                          "CompleteEDTA_42_TrypChym modified",
                                                          "CompleteEDTA_42_TrypChym unmodified",
                                                          "Complete_37_Tryp modified",
                                                          "Complete_37_Tryp unmodified",
                                                          "Complete_37_TrypChym modified",
                                                          "Complete_37_TrypChym unmodified",
                                                          "Complete_42_Tryp modified",
                                                          "Complete_42_Tryp unmodified",
                                                          "Complete_42_TrypChym modified",
                                                          "Complete_42_TrypChym unmodified"])
    print("Raw length:", len(modifications_df))
    # Filter and calculate the percentages
    modifications_df["Crtl_37_Tryp"], modifications_df["Crtl_37_TrypChym"], modifications_df["Crtl_42_Tryp"], \
        modifications_df["Crtl_42_TrypChym"], modifications_df["CompleteEDTA_37_Tryp"], \
        modifications_df["CompleteEDTA_37_TrypChym"], modifications_df["CompleteEDTA_42_Tryp"], \
        modifications_df["CompleteEDTA_42_TrypChym"], modifications_df["Complete_37_Tryp"], \
        modifications_df["Complete_37_TrypChym"], modifications_df["Complete_42_Tryp"], \
        modifications_df["Complete_42_TrypChym"] \
        = zip(*modifications_df.apply(_filter_rows, axis=1, args=(pass_percentage,)))
    modifications_df = modifications_df.dropna()

    print("Filtered length (After removing rows where none of the conditions has at least",
          pass_percentage, "% modification):", len(modifications_df))

    modifications_df = pd.melt(modifications_df, var_name="Condition", id_vars=["Protein Position", "Modifications"],
                               value_name="Percentage",
                               value_vars=["Crtl_37_Tryp", "Crtl_37_TrypChym", "Crtl_42_Tryp", "Crtl_42_TrypChym",
                                           "CompleteEDTA_37_Tryp", "CompleteEDTA_37_TrypChym", "CompleteEDTA_42_Tryp",
                                           "CompleteEDTA_42_TrypChym", "Complete_37_Tryp", "Complete_37_TrypChym",
                                           "Complete_42_Tryp", "Complete_42_TrypChym"])
    return modifications_df


def _calculate_percentages(row):
    """
    Calculate the modification percentages for each row.

    :param row: The row.
    :return: The tuple with the percentages.
    """
    crtl_37_tryp_percentage = row["Crtl_37_Tryp modified"] / (
            row["Crtl_37_Tryp modified"] + row["Crtl_37_Tryp unmodified"]) * 100 \
        if row["Crtl_37_Tryp modified"] + row["Crtl_37_Tryp unmodified"] != 0 else 0
    crtl_37_trypchym_percentage = row["Crtl_37_TrypChym modified"] / (
            row["Crtl_37_TrypChym modified"] + row["Crtl_37_TrypChym unmodified"]) * 100 \
        if row["Crtl_37_TrypChym modified"] + row["Crtl_37_TrypChym unmodified"] != 0 else 0
    crtl_42_tryp_percentage = row["Crtl_42_Tryp modified"] / (
            row["Crtl_42_Tryp modified"] + row["Crtl_42_Tryp unmodified"]) * 100 \
        if row["Crtl_42_Tryp modified"] + row["Crtl_42_Tryp unmodified"] != 0 else 0
    crtl_42_trypchym_percentage = row["Crtl_42_TrypChym modified"] / (
            row["Crtl_42_TrypChym modified"] + row["Crtl_42_TrypChym unmodified"]) * 100 \
        if row["Crtl_42_TrypChym modified"] + row["Crtl_42_TrypChym unmodified"] != 0 else 0
    complete_edta_37_tryp_percentage = row["CompleteEDTA_37_Tryp modified"] / (
            row["CompleteEDTA_37_Tryp modified"] + row["CompleteEDTA_37_Tryp unmodified"]) * 100 \
        if row["CompleteEDTA_37_Tryp modified"] + row["CompleteEDTA_37_Tryp unmodified"] != 0 else 0
    complete_edta_37_trypchym_percentage = row["CompleteEDTA_37_TrypChym modified"] / (
            row["CompleteEDTA_37_TrypChym modified"] + row["CompleteEDTA_37_TrypChym unmodified"]) * 100 \
        if row["CompleteEDTA_37_TrypChym modified"] + row["CompleteEDTA_37_TrypChym unmodified"] != 0 else 0
    complete_edta_42_tryp_percentage = row["CompleteEDTA_42_Tryp modified"] / (
            row["CompleteEDTA_42_Tryp modified"] + row["CompleteEDTA_42_Tryp unmodified"]) * 100 \
        if row["CompleteEDTA_42_Tryp modified"] + row["CompleteEDTA_42_Tryp unmodified"] != 0 else 0
    complete_edta_42_trypchym_percentage = row["CompleteEDTA_42_TrypChym modified"] / (
            row["CompleteEDTA_42_TrypChym modified"] + row["CompleteEDTA_42_TrypChym unmodified"]) * 100 \
        if row["CompleteEDTA_42_TrypChym modified"] + row["CompleteEDTA_42_TrypChym unmodified"] != 0 else 0
    complete_37_tryp_percentage = row["Complete_37_Tryp modified"] / (
            row["Complete_37_Tryp modified"] + row["Complete_37_Tryp unmodified"]) * 100 \
        if row["Complete_37_Tryp modified"] + row["Complete_37_Tryp unmodified"] != 0 else 0
    complete_37_trypchym_percentage = row["Complete_37_TrypChym modified"] / (
            row["Complete_37_TrypChym modified"] + row["Complete_37_TrypChym unmodified"]) * 100 \
        if row["Complete_37_TrypChym modified"] + row["Complete_37_TrypChym unmodified"] != 0 else 0
    complete_42_tryp_percentage = row["Complete_42_Tryp modified"] / (
            row["Complete_42_Tryp modified"] + row["Complete_42_Tryp unmodified"]) * 100 \
        if row["Complete_42_Tryp modified"] + row["Complete_42_Tryp unmodified"] != 0 else 0
    complete_42_trypchym_percentage = row["Complete_42_TrypChym modified"] / (
            row["Complete_42_TrypChym modified"] + row["Complete_42_TrypChym unmodified"]) * 100 \
        if row["Complete_42_TrypChym modified"] + row["Complete_42_TrypChym unmodified"] != 0 else 0

    return crtl_37_tryp_percentage, crtl_37_trypchym_percentage, crtl_42_tryp_percentage, crtl_42_trypchym_percentage, \
        complete_edta_37_tryp_percentage, complete_edta_37_trypchym_percentage, complete_edta_42_tryp_percentage, \
        complete_edta_42_trypchym_percentage, complete_37_tryp_percentage, complete_37_trypchym_percentage, \
        complete_42_tryp_percentage, complete_42_trypchym_percentage


def _filter_rows(row, pass_percentage):
    """
    Filter the rows to require at least one sample above or equal to the pass percentage.

    :param row: The row.
    :param pass_percentage: The pass percentage.
    :return: True, if row passed the requirement; Otherwise. False.
    """
    crtl_37_tryp_percentage, crtl_37_trypchym_percentage, crtl_42_tryp_percentage, crtl_42_trypchym_percentage, \
        complete_edta_37_tryp_percentage, complete_edta_37_trypchym_percentage, complete_edta_42_tryp_percentage, \
        complete_edta_42_trypchym_percentage, complete_37_tryp_percentage, complete_37_trypchym_percentage, \
        complete_42_tryp_percentage, complete_42_trypchym_percentage = _calculate_percentages(row=row)

    valid_entry: bool = crtl_37_tryp_percentage >= pass_percentage or \
        crtl_37_trypchym_percentage >= pass_percentage or \
        crtl_42_tryp_percentage >= pass_percentage or \
        crtl_42_trypchym_percentage >= pass_percentage or \
        complete_edta_37_tryp_percentage >= pass_percentage or \
        complete_edta_37_trypchym_percentage >= pass_percentage or \
        complete_edta_42_tryp_percentage >= pass_percentage or \
        complete_edta_42_trypchym_percentage >= pass_percentage or \
        complete_37_tryp_percentage >= pass_percentage or \
        complete_37_trypchym_percentage >= pass_percentage or \
        complete_42_tryp_percentage >= pass_percentage or \
        complete_42_trypchym_percentage >= pass_percentage

    if valid_entry:
        return crtl_37_tryp_percentage, crtl_37_trypchym_percentage, crtl_42_tryp_percentage, \
               crtl_42_trypchym_percentage, complete_edta_37_tryp_percentage, complete_edta_37_trypchym_percentage, \
               complete_edta_42_tryp_percentage, complete_edta_42_trypchym_percentage, complete_37_tryp_percentage, \
               complete_37_trypchym_percentage, complete_42_tryp_percentage, complete_42_trypchym_percentage
    else:
        return [None] * 12


def _create_plot(modification_data: pd.DataFrame, mod_name: str, save_file_name: str):
    sns.set(style="darkgrid")
    sns.set_palette(palette="Paired", n_colors=12)
    g = sns.catplot(x="Protein Position", y="Percentage", hue="Condition", data=modification_data,
                    height=10, aspect=1.5, s=7.5)
    plt.xlabel("Position")
    plt.ylabel("Modified peptides (%)")
    plt.title(mod_name)  # TODO: Find good title
    #plt.show()
    plt.close()
    g.savefig(save_file_name)


if __name__ == "__main__":
    file_name: str = r"C:\Users\spec-makie17\Documents\Experiments\210903_Complete_Incubation\ptmprofile.csv"
    output_directory: str = \
        r"C:\Users\spec-makie17\Documents\Experiments\210903_Complete_Incubation\Open_search_modifications1"
    valid_entry_modification_percentage = 1
    modification_df: pd.DataFrame = _read_data(file=file_name, pass_percentage=valid_entry_modification_percentage)

    modifications: List[str] = modification_df['Modifications'].unique()
    print("Found", len(modifications), "modifications:", ", ".join(modifications))

    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    for mod in tqdm(modifications):
        mod_file_name: str = os.path.join(output_directory,
                                          f"{mod.lower().replace(' ', '_').replace('(', '').replace(')', '')}.png")

        mod_df: pd.DataFrame = modification_df[modification_df["Modifications"] == mod]
        _create_plot(modification_data=mod_df, mod_name=mod, save_file_name=mod_file_name)
