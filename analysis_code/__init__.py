"""
The package containing the final analysis code for the master thesis.
"""
from .convert_peptide_list import convert_peptide_list_files

from .utils import get_residue_positions, get_residue_name

from .modification_statistics import combine_spectra_in_peptide_lists, calculate_modification_percentages,\
    create_modification_barplot
