"""
Filter data for protein domain identification.

This module contains constants that are used for the filtering.
"""

# All InterPro IDs that have the word "armadillo" in the description
arm_accs_ipr = [
    'IPR013636',
    'IPR024574',
    'IPR041322',
    'IPR055241',
    'IPR056252',
    'IPR016617',
    'IPR031524',
    'IPR038739',
    'IPR038905',
    'IPR039868',
    'IPR040268',
    'IPR042462',
    'IPR042834',
    'IPR043379',
    'IPR044282',
    'IPR051303',
    'IPR052441',
    'IPR011989',
    'IPR016024',
    'IPR000225',
    'IPR041209',
    'IPR049152',
    'IPR006911'
]

# All Pfam accessions that have "armadillo" in the description (with version numbers)
arm_accs_pfam = [
    'PF00514.29',
    'PF17822.6',
    'PF08427.15',
    'PF15767.10',
    'PF22915.1',
    'PF04826.19',
    'PF23295.1',
    'PF16629.10',
    'PF18770.7',
    'PF21052.3',
    'PF11841.14',
    'PF14726.11',
    'PF18581.7'
]

# Definitions for Disordered Regions IUPred3
IUPRED3_THRESHOLD = 0.5
MIN_LENGTH_DISORDERED_REGION = 20
IUPRED_CACHE_DIR = '/home/markus/MPI_local/production1/IUPred3'
IUPRED3_PATH = '../../iupred3'

# Limit for AF3 tokens
AF_TOKEN_LIMIT = 5120

AF_TRAINING_CUTOFF = '2021-09-30'
