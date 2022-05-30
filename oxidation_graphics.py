import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def main():
    data_take: pd.DataFrame = pd.read_excel(
        r"C:\Users\spec-makie17\Documents\Experiments\FigureGeneration\OxidationCysMod\ModsT3.xlsx",
        sheet_name="CysOxidation", index_col=0)
    create_modification_plot_take3(data_take)


def create_modification_plot(data: pd.DataFrame):
    axis = sns.barplot(data=data, x=data.index, y="Percentage")
    #axis = sns.barplot(data=data, x=data.index, y="Percentage", hue="Portion")
    # axis.set_title("Total oxidation of methionine, proline and histine after 72 hour incubation", fontsize=19)
    axis.set_title("Total oxidation of methionine, proline and histine after incubation with ROS inhibitors", fontsize=19)
    #axis.set_title("Total dehydroalanine conversion of cysteine", fontsize=19)
    axis.set_xlabel("Conditions", fontsize=17)
    axis.set_ylabel("Percentage\n(Modified hit count/Modifiable hit count)", fontsize=17)
    for idx, p in enumerate(axis.patches):
        axis.annotate(
            f"{p.get_height()} %\n({int(data.iloc[idx]['ModifiedSpectra'])}/{int(data.iloc[idx]['TotalModSpectra'])})",
            (p.get_x() + p.get_width() / 2., p.get_height()),
            ha='center', va='center', fontsize=17, color='black', xytext=(0, 11),
            textcoords='offset points')

    plt.subplots_adjust(left=0.04, bottom=0.08, right=0.97, top=0.95)
    plt.tick_params(axis='both', which='major', labelsize=14)
    #plt.legend(loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.085), fontsize=13)
    plt.show()


def create_modification_plot_take3(data: pd.DataFrame):
    def create_subplots(df: pd.DataFrame, time: str, axis):
        print(df)
        sns.barplot(data=df, x=df.index, y="Percentage", hue="Condition", ax=axis)
        #axis.set_title(f"Total oxidation of methionine, proline and histine {time} hours", fontsize=19)
        #axis.set_title(f"Total dehydroalanine conversion of cysteine {time} hours", fontsize=19)
        axis.set_title(f"Total oxidation of cysteine {time} hours", fontsize=19)
        axis.set_xlabel("Condition", fontsize=17)
        axis.set_ylabel("Percentage\n(Modified hit count/Modifiable hit count)", fontsize=17)
        for idx, p in enumerate(axis.patches):
            axis.annotate(
                f"{p.get_height()} %\n({int(df.iloc[idx]['ModifiedSpectra'])}/{int(df.iloc[idx]['TotalModSpectra'])})",
                (p.get_x() + p.get_width() / 2., p.get_height()),
                ha='center', va='center', fontsize=17, color='black', xytext=(0, 11),
                textcoords='offset points')

        axis.tick_params(axis='both', which='major', labelsize=14)

    fig, ax = plt.subplots(nrows=2, ncols=1)
    create_subplots(data[data["Time"] == 0], "at 0", ax[0])
    create_subplots(data[data["Time"] == 72], "after 72", ax[1])

    ax[0].legend().remove()
    ax[1].legend(loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.2), fontsize=13)
    plt.subplots_adjust(left=0.04, bottom=0.08, right=0.97, top=0.95)

    plt.show()


if __name__ == "__main__":
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    pd.set_option('display.max_colwidth', None)
    main()
