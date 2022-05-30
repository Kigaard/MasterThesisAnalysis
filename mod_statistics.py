import pathlib

import pandas as pd

from analysis_code import calculate_modification_percentages, combine_spectra_in_peptide_lists


def perform_analysis(peptide_list_directory: pathlib.Path, sequence: str, modifications: list[tuple[str, str, float]]):
    with pd.ExcelWriter(peptide_list_directory / f"../CRTModStats.xlsx") as writer:
        for mod_name, res, mod_mass in modifications:
            print("*" * 5, f"{mod_name} ({res}@{mod_mass})", "*" * 5)
            mod_dict: dict = calculate_modification_percentages(peptide_list_directory=peptide_list_directory,
                                                                sequence=sequence,
                                                                residue_str=res,
                                                                mod_mass=mod_mass)
            mod_df: pd.DataFrame = pd.DataFrame.from_dict(mod_dict, orient="index")
            mod_df = mod_df.reindex(["Nat_Crt_0", "Nat_lacto_0", "Nat_Ribo_0", "Nat_Crt_72", "RedAlk_Crt_72", "Nat_LactoCrt_72_Crt", "RedAlk_LactoCrt_72_Crt", "Nat_LactoCrt_72_Bait", "RedAlk_LactoCrt_72_Bait", "Nat_RiboCrt_72_Crt", "RedAlk_RiboCrt_72_Crt", "Nat_RiboCrt_72_Bait", "RedAlk_RiboCrt_72_Bait"])

            print(mod_df)
            print()
            mod_df.to_excel(writer, sheet_name=mod_name)
            


def main():
    base_directory = pathlib.Path(
        r"C:\Users\spec-makie17\Documents\Experiments\FigureGeneration\Bait_Take3")
    raw_peptide_list_folder = pathlib.Path("Lists_CBM")
    peptide_list_folder = pathlib.Path("Lists_Combined_CBM")
    bsa_sequence: str = "DTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLASSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA"
    crt_sequence: str = "EPAVYFKEQFLDGDGWTSRWIESKHKSDFGKFVLSSGKFYGDEEKDKGLQTSQDARFYALSASFEPFSNKGQTLVVQFTVKHEQNIDCGGGYVKLFPNSLDQTDMHGDSEYNIMFGPDICGPGTKKVHVIFNYKGKNVLINKDIRCKDDEFTHLYTLIVRPDNTYEVKIDNSQVESGSLEDDWDFLPPKKIKDPDASKPEDWDERAKIDDPTDSKPEDWDKPEHIPDPDAKKPEDWDEEMDGEWEPPVIQNPEYKGEWKPRQIDNPDYKGTWIHPEIDNPEYSPDPSIYAYDNFGVLGLDLWQVKSGTIFDNFLITNDEAYAEEFGNETWGVTKAAEKQMKDKQDEEQRLKEEEEDKKRKEEEEAEDKEDDEDKDEDEEDEEDKEEDEEEDVPGQAKDEL"
    lacto_sequence: str = "LIVTQTMKGLDIQKVAGTWYSLAMAASDISLLDAQSAPLRVYVEELKPTPEGDLEILLQKWENGECAQKKIIAEKTKIPAVFKIDALNENKVLVLDTDYKKYLLFCMENSAEPEQSLACQCLVRTPEVDDEALEKFDKALKALPMHIRLSFNPTQLEEQCHI"
    ribo_sequence: str = "KETAAAKFERQHMDSSTSAASSSNYCNQMMKSRNLTKDRCKPVNTFVHESLADVQAVCSQKNVACKNGQTNCYQSYSTMSITDCRETGSSKYPNCAYKTTQANKHIIVACEGNPYVPVHFDASV"
    if False:
        list_dir = base_directory / raw_peptide_list_folder
        save_dir = base_directory / peptide_list_folder
        combine_spectra_in_peptide_lists(list_directory=list_dir, save_directory=save_dir)
        exit()

    ox_mod: list[tuple[str, str, float]] = [("Oxidation", "MPH", 15.995)]

    modifications_cbm_cysteine: list[tuple[str, str, float]] = [
        ("Oxidation", "MPH", 15.995),
        ("Carbamidomethyl", "C", 57.022),
        ("Dehydroalanine", "C", -87.986),
        ("Cys_Oxidation", "C", -41.027),
        ("Sulfinic", "C", -25.032),
        ("SulfDiOx", "C", 6.940),
        ("SulfOx", "C", -9.054),
        ("Sulfonic", "C", -9.037),
        ("SSulfonic", "C", 22.935),
        ("SO3", "C", 34.935)]

    modifications_cysteine: list[tuple[str, str, float]] = [
        ("Oxidation", "MPH", 15.995),
        ("Carbamidomethyl", "C", 57.022),
        ("Dehydroalanine", "C", -33.988),
        ("Cys_Oxidation", "C", 15.995),
        ("Sulfinic", "C", 31.990),
        ("SulfDiOx", "C", 63.962),
        ("SulfOx", "C", 47.967),
        ("Sulfonic", "C", 47.985),
        ("SSulfonic", "C", 79.957),
        ("SO3", "C", 91.957)]


    perform_analysis(peptide_list_directory=base_directory / peptide_list_folder,
                     modifications=modifications_cysteine, sequence=crt_sequence)


if __name__ == "__main__":
    main()
