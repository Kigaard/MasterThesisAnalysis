import os
from typing import List, Tuple

import pandas as pd
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.axes


def create_plots(bait_statistics: str):
    alpha_bait: pd.DataFrame = pd.read_excel(io=bait_statistics, sheet_name="Alpha")
    bsa_bait: pd.DataFrame = pd.read_excel(io=bait_statistics, sheet_name="BSA")
    lacto_bait: pd.DataFrame = pd.read_excel(io=bait_statistics, sheet_name="Lacto")
    ribonuclease_bait: pd.DataFrame = pd.read_excel(io=bait_statistics, sheet_name="Ribonuclease")

    alpha_bait = pd.melt(alpha_bait, id_vars=["Type", "Protein", "Condition"], value_name="Value",
                         value_vars=["Immediate", "1 day", "2 days", "4 days"], var_name="Time")
    bsa_bait = pd.melt(bsa_bait, id_vars=["Type", "Protein", "Condition"], value_name="Value",
                       value_vars=["Immediate", "1 day", "2 days", "4 days"], var_name="Time")
    lacto_bait = pd.melt(lacto_bait, id_vars=["Type", "Protein", "Condition"], value_name="Value",
                         value_vars=["Immediate", "1 day", "2 days", "4 days"], var_name="Time")
    ribonuclease_bait = pd.melt(ribonuclease_bait, id_vars=["Type", "Protein", "Condition"], value_name="Value",
                                value_vars=["Immediate", "1 day", "2 days", "4 days"], var_name="Time")

    fig, ax = plt.subplots(nrows=4, ncols=4)
    _create_bait_plot(bait_df=alpha_bait, bait_name="Alpha-1-acid glycoprotein 1",
                      coverage_crt_plot=ax[0, 0], coverage_bait_plot=ax[1, 0],
                      hits_crt_plot=ax[2, 0], hits_bait_plot=ax[3, 0])

    _create_bait_plot(bait_df=bsa_bait, bait_name="BSA",
                      coverage_crt_plot=ax[0, 1], coverage_bait_plot=ax[1, 1],
                      hits_crt_plot=ax[2, 1], hits_bait_plot=ax[3, 1])
    _create_bait_plot(bait_df=lacto_bait, bait_name="Beta-lactoglobulin",
                      coverage_crt_plot=ax[0, 2], coverage_bait_plot=ax[1, 2],
                      hits_crt_plot=ax[2, 2], hits_bait_plot=ax[3, 2])
    _create_bait_plot(bait_df=ribonuclease_bait, bait_name="Ribonuclease B",
                      coverage_crt_plot=ax[0, 3], coverage_bait_plot=ax[1, 3],
                      hits_crt_plot=ax[2, 3], hits_bait_plot=ax[3, 3])

    #fig.suptitle("Coverage and hits for the CRT and reduced and alkylated bait")
    #fig.legend(labels=["Autodigest", "Trypsin 0x18O - 16O", "Trypsin 1x18O - 16O", "Trypsin 2x18O - 16O",
    #                   "Trypsin 0x18O - 18O", "Trypsin 1x18O - 18O", "Trypsin 2x18O - 18O"], loc=8, ncol=7)
    plt.show()


def _create_bait_plot(bait_df: pd.DataFrame, bait_name: str,
                      coverage_crt_plot: matplotlib.axes.Axes, coverage_bait_plot: matplotlib.axes.Axes,
                      hits_crt_plot: matplotlib.axes.Axes, hits_bait_plot: matplotlib.axes.Axes):
    """
    Create the plot for the bait.

    :param bait_df: The dataframe containing the bait information-
    :param bait_name: The name of the bait.
    :param coverage_plot: The coverage plot.
    :param hit_plot: The hit plot.
    :return: None
    """
    coverage_crt_df: pd.DataFrame = bait_df[(bait_df["Type"] == "Coverage") & (bait_df["Protein"] == "CRT")]
    coverage_bait_df: pd.DataFrame = bait_df[(bait_df["Type"] == "Coverage") & (bait_df["Protein"] == "Bait")]
    hits_crt_df: pd.DataFrame = bait_df[(bait_df["Type"] == "Hits") & (bait_df["Protein"] == "CRT")]
    hits_bait_df: pd.DataFrame = bait_df[(bait_df["Type"] == "Hits") & (bait_df["Protein"] == "Bait")]
    # TODO: Change to seaborn
    # Create the coverage plot
    sns.lineplot(data=coverage_crt_df, x="Time", y="Value", hue="Condition",
                 ci=None, ax=coverage_crt_plot, legend=None)
    coverage_crt_plot.set_title(f"{bait_name} - CRT Coverage")
    coverage_crt_plot.set_xlabel("")
    coverage_crt_plot.set_ylabel("Coverage (%)")

    sns.lineplot(data=coverage_bait_df, x="Time", y="Value", hue="Condition",
                 ci=None, ax=coverage_bait_plot, legend=None)
    coverage_bait_plot.set_title(f"{bait_name} - Bait Coverage")
    coverage_bait_plot.set_xlabel("")
    coverage_bait_plot.set_ylabel("Coverage (%)")

    # Create the hit plot
    sns.lineplot(data=hits_crt_df, x="Time", y="Value", hue="Condition",
                 ci=None, ax=hits_crt_plot, legend=None)
    hits_crt_plot.set_title(f"{bait_name} - CRT Hits")
    hits_crt_plot.set_xlabel("")
    hits_crt_plot.set_ylabel("Total hits")

    sns.lineplot(data=hits_bait_df, x="Time", y="Value", hue="Condition",
                 ci=None, ax=hits_bait_plot, legend=None)
    hits_bait_plot.set_title(f"{bait_name} - Bait Hits")
    hits_bait_plot.set_xlabel("")
    hits_bait_plot.set_ylabel("Total hits")


if __name__ == "__main__":
    create_plots(bait_statistics=r"C:\Users\spec-makie17\Documents\Experiments\211012_Bait_PH2121\Bait statistics.xlsx")
