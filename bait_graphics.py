import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


def main():
    filename: str = r"C:\Users\spec-makie17\Documents\Experiments\FigureGeneration\Bait\HitStats.xlsx"
    # lacto_t1_df: pd.DataFrame = pd.read_excel(filename, sheet_name="Lacto_T1")
    # ribo_t1_df: pd.DataFrame = pd.read_excel(filename, sheet_name="Ribo_T1")
    # lacto_t2_df: pd.DataFrame = pd.read_excel(filename, sheet_name="Lacto_T2")
    # ribo_t2_df: pd.DataFrame = pd.read_excel(filename, sheet_name="Ribo_T2")
    # data_t3_df: pd.DataFrame = pd.read_excel(filename, sheet_name="Take3")
    #
    # #create_barplot_take_12(lacto_t1_df, ribo_t1_df)
    # #create_barplot_take_12(lacto_t2_df, ribo_t2_df)
    # create_barplot_take3(data_t3_df)
    df: pd.DataFrame = pd.read_excel(filename, sheet_name="Take3_CBM")
    create_barplot_take3(df)


def create_barplot_take_12(lacto_data: pd.DataFrame, ribo_data: pd.DataFrame):
    def create_chart(data_f: pd.DataFrame, bait_name: str, condition: str, axis):
        sns.barplot(data=data_f, x="Time", y="Values", hue="Condition", ci=None, ax=axis)
        for container in axis.containers:
            axis.bar_label(container, fontsize=17, padding=-5)
        axis.get_legend().remove()
        axis.set_title(f"Hit counts for bait experiment with CRT and {bait_name} - {condition}", fontsize=19)
        axis.set_xlabel("", fontsize=17)
        axis.set_ylabel("Hit count", fontsize=17)
        axis.tick_params(axis='both', which='major', labelsize=16)

    lacto_crt_df = lacto_data[lacto_data["Type"] == "CRT"]
    lacto_bait_df = lacto_data[lacto_data["Type"] == "Bait"]

    ribo_crt_df = ribo_data[ribo_data["Type"] == "CRT"]
    ribo_bait_df = ribo_data[ribo_data["Type"] == "Bait"]

    fig, ax = plt.subplots(nrows=4, ncols=1)

    create_chart(lacto_crt_df, "Beta-lactoglobulin", "CRT", ax[0])
    create_chart(lacto_bait_df, "Beta-lactoglobulin", "Bait", ax[1])
    create_chart(ribo_crt_df, "Ribonuclease B", "CRT", ax[2])
    create_chart(ribo_bait_df, "Ribonuclease B", "Bait", ax[3])

    ax[3].set_xlabel("Incubation time (Days)", fontsize=17)
    plt.legend(loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.4), fontsize=15)
    plt.subplots_adjust(left=0.04, bottom=0.08, right=0.97, top=0.95)

    plt.show()


def create_barplot_take3(data: pd.DataFrame):
    def create_subplot(hits_df: pd.DataFrame, condition: str, axis):
        sns.barplot(data=hits_df, x="Type", y="UniqueHits", hue="Condition", ax=axis, ci=None)
        for container in axis.containers:
            axis.bar_label(container, fontsize=17)
        axis.get_legend().remove()
        axis.set_title(f"Unique hit counts for the {condition}", fontsize=19)
        axis.set_xlabel("Proteins", fontsize=17)
        axis.set_ylabel("Unique hit count", fontsize=17)
        axis.tick_params(axis='both', which='major', labelsize=17)

    fig, ax = plt.subplots(nrows=2, ncols=1)
    create_subplot(data[data["Time"] == 0], "the native proteins at 0 hours", ax[0])
    create_subplot(data[data["Time"] == 72], "bait experiment at 72 hours", ax[1])

    plt.legend(loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.2), fontsize=13)
    plt.subplots_adjust(left=0.04, bottom=0.08, right=0.97, top=0.95)

    plt.show()


if __name__ == "__main__":
    main()
