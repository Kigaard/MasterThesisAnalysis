from typing import List, Tuple, Dict

import re


def get_residue_positions(residues: str) -> List[Tuple[int, str]]:
    """
    Get residues in CRT.

    :param residues: The residues to search for.
    :return: The tuple containing the position and the residue and position. Fx (163, 'C163').
    """

    sequence = \
        "MLLSVPLLLGLLGLAVAEPAVYFKEQFLDGDGWTSRWIESKHKSDFGKFVLSSGKFYGDE" + \
        "EKDKGLQTSQDARFYALSASFEPFSNKGQTLVVQFTVKHEQNIDCGGGYVKLFPNSLDQT" + \
        "DMHGDSEYNIMFGPDICGPGTKKVHVIFNYKGKNVLINKDIRCKDDEFTHLYTLIVRPDN" + \
        "TYEVKIDNSQVESGSLEDDWDFLPPKKIKDPDASKPEDWDERAKIDDPTDSKPEDWDKPE" + \
        "HIPDPDAKKPEDWDEEMDGEWEPPVIQNPEYKGEWKPRQIDNPDYKGTWIHPEIDNPEYS" + \
        "PDPSIYAYDNFGVLGLDLWQVKSGTIFDNFLITNDEAYAEEFGNETWGVTKAAEKQMKDK" + \
        "QDEEQRLKEEEEDKKRKEEEEAEDKEDDEDKDEDEEDEEDKEEDEEEDVPGQAKDEL"
    result = [(match.start() + 1, f"{match.group(0)}{match.start() + 1}") for match in
              re.finditer(f"[{residues}]", sequence.upper())]

    return result


def get_residue_name(three_letter: str) -> str:
    """
    Get the residue name from the three letter code.

    :param three_letter: The tree letter code.
    :return: The residue name.
    """
    residues: Dict[str, str] = {"M": "Methionine", "P": "Proline", "C": "Cysteine"}
    return residues[three_letter.upper()]
