import pathlib
import re
from collections import defaultdict

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def has_oxidation(mod_str: str, pos: int) -> bool:
    if mod_str == "-":
        return False

    mod_list: list[str] = []

    if ';' in mod_str:
        for part in mod_str.split(';'):
            mod_list.append(part)
    else:
        mod_list.append(mod_str)

    return f"{pos}@{MASS_CHANGE}" in mod_list


def create_dataframe(file_path: pathlib.Path) -> pd.DataFrame:
    df: pd.DataFrame = pd.read_excel(file_path, index_col=0)

    raw_pos_mod: defaultdict = defaultdict(lambda: {"Percentage": 0, "ModifiedSpectra": 0, "TotalModSpectra": 0})

    for _, row in df.iterrows():
        for pos in RESIDUE_POSITIONS:
            if int(row["Start"]) <= int(pos) <= int(row["End"]):
                raw_pos_mod[int(pos)]["TotalModSpectra"] += 1
                if has_oxidation(mod_str=row["Modification"], pos=int(pos)):
                    raw_pos_mod[int(pos)]["ModifiedSpectra"] += 1

    for pos in raw_pos_mod:
        if raw_pos_mod[int(pos)]["TotalModSpectra"] == 0:
            raw_pos_mod[int(pos)]["Percentage"] = 0
        else:
            raw_pos_mod[int(pos)]["Percentage"] = round(
                (raw_pos_mod[int(pos)]["ModifiedSpectra"] / raw_pos_mod[int(pos)]["TotalModSpectra"]) * 100, 2)

    pos_mod: defaultdict = defaultdict(lambda: {"Percentage": 0, "ModifiedSpectra": 0, "TotalModSpectra": 0})

    for pos in raw_pos_mod:
        if raw_pos_mod[pos]["Percentage"] != 0:
            new_id: str = f"{SEQUENCE[pos - 1]}{pos + 17}"
            pos_mod[new_id] = raw_pos_mod[pos]

    percentage_df: pd.DataFrame = pd.DataFrame.from_dict(pos_mod, orient='index')
    return percentage_df


def create_modification_plot(data: pd.DataFrame, condition: str):
    axis = sns.barplot(data=data, x=data.index, y="Percentage", color="tab:blue")

    #axis.set_title(f"Oxidation of methionine, proline, and histidine per position - {condition}", fontsize=19)
    axis.set_title(f"Dehydroalanine conversion of cysteine per position - {condition}", fontsize=19)
    axis.set_xlabel("Position", fontsize=17)
    axis.set_ylabel("Percentage\n(Modified hit count/Modifiable hit count)", fontsize=17)
    for idx, p in enumerate(axis.patches):
        axis.annotate(
            f"{p.get_height()} %\n({int(data.iloc[idx]['ModifiedSpectra'])}/{int(data.iloc[idx]['TotalModSpectra'])})",
            (p.get_x() + p.get_width() / 2., p.get_height()),
            ha='center', va='center', fontsize=17, color='black', xytext=(0, 11),
            textcoords='offset points')

    plt.subplots_adjust(left=0.04, bottom=0.08, right=0.97, top=0.95)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.show()


def main():
    file_directory: pathlib.Path = pathlib.Path(
        r"C:\Users\spec-makie17\Documents\Experiments\220425_Batches_PH22006\PeptideLists_Combined")

    files: list[tuple[str, pathlib.Path]] = [

        ("Batch 3.5 - Retained", pathlib.Path("CRT35_Retained.xlsx")),
        ("Batch 2.5 - Retained", pathlib.Path("CRT25_Retained.xlsx")),
        ("Batch 19/03F - Retained", pathlib.Path("CRT1903_Retained.xlsx")),
        ("Batch 18/11F - Retained", pathlib.Path("CRT1811_Retained.xlsx")),
        ("Batch 18/06F - Retained", pathlib.Path("CRT1806_Retained.xlsx")),
        ("Batch 3.5 - Flowthrough", pathlib.Path("CRT35_Flowthrough.xlsx")),
        ("Batch 2.5 - Flowthrough", pathlib.Path("CRT25_Flowthrough.xlsx")),
        ("Batch 19/03F - Flowthrough", pathlib.Path("CRT1903_Flowthrough.xlsx")),
        ("Batch 18/11F - Flowthrough", pathlib.Path("CRT1811_Flowthrough.xlsx")),
        ("Batch 18/06F - Flowthrough", pathlib.Path("CRT1806_Flowthrough.xlsx"))
    ]

    with pd.ExcelWriter(file_directory / f"../Oxidations.xlsx") as writer:
        for cond, file in files:
            df: pd.DataFrame = create_dataframe(file_path=file_directory / file)
            if len(df) < 1:
                continue
            df.to_excel(writer, cond.replace("/", ""))
            create_modification_plot(df, cond)


if __name__ == "__main__":
    SEQUENCE: str = "EPAVYFKEQFLDGDGWTSRWIESKHKSDFGKFVLSSGKFYGDEEKDKGLQTSQDARFYALSASFEPFSNKGQTLVVQFTVKHEQNIDCGGGYVKLFP" + \
                    "NSLDQTDMHGDSEYNIMFGPDICGPGTKKVHVIFNYKGKNVLINKDIRCKDDEFTHLYTLIVRPDNTYEVKIDNSQVESGSLEDDWDFLPPKKIKDP" + \
                    "DASKPEDWDERAKIDDPTDSKPEDWDKPEHIPDPDAKKPEDWDEEMDGEWEPPVIQNPEYKGEWKPRQIDNPDYKGTWIHPEIDNPEYSPDPSIYAY" + \
                    "DNFGVLGLDLWQVKSGTIFDNFLITNDEAYAEEFGNETWGVTKAAEKQMKDKQDEEQRLKEEEEDKKRKEEEEAEDKEDDEDKDEDEEDEEDKEEDEEEDVPGQAKDEL"
    #RESIDUES: str = "MPH"
    RESIDUES: str = "C"
    #MASS_CHANGE: str = "15.995"
    MASS_CHANGE: str = "-87.986"
    RESIDUE_POSITIONS = [match.start() + 1 for match in re.finditer(f"[{RESIDUES}]", SEQUENCE)]
    print(RESIDUE_POSITIONS)
    main()
