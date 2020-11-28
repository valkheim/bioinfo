from dataclasses import dataclass
from typing import List


@dataclass
class Base:
    name: str
    letter: str
    index: int


@dataclass
class AminoAcid:
    full_name: str
    short_name: str
    letter: str
    codons: List[str]


CODON_LENGTH = 3


START = [
    'AUG', # methionine
    'GUG', # valine
    'UUG'  # leucine
]

STOP = [
    'UAA', # ochre
    'UAG', # ambre
    'UGA'  # opale
]


AMINO_ACIDS = [
    AminoAcid('Alanine',       'Ala', 'A', [ 'GCU', 'GCC', 'GCA', 'GCG' ]),
    AminoAcid('Arginine',      'Arg', 'R', [ 'CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG' ]),
    AminoAcid('Asparagine',    'Asn', 'N', [ 'AAU', 'AAC' ]),
    AminoAcid('Aspartic acid', 'Asp', 'D', [ 'GAU', 'GAC' ]),
    AminoAcid('Cysteine',      'Cys', 'C', [ 'UGU', 'UGC' ]),
    AminoAcid('Glutamine',     'Gln', 'Q', [ 'CAA', 'CAG' ]),
    AminoAcid('Glutamic acid', 'Glu', 'E', [ 'GAA', 'GAG' ]),
    AminoAcid('Glycine',       'Gly', 'G', [ 'GGU', 'GGC', 'GGA', 'GGU' ]),
    AminoAcid('Histidine',     'His', 'H', [ 'CAU', 'CAC' ]),
    AminoAcid('Isoleucine',    'Ile', 'I', [ 'AUU', 'AUC', 'AUA' ]),
    AminoAcid('Leucine',       'Leu', 'L', [ 'UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG' ]),
    AminoAcid('Lysine',        'Lys', 'K', [ 'AAA', 'AAG' ]),
    AminoAcid('Methionine',    'Met', 'M', [ 'AUG' ]),
    AminoAcid('Phenylalanine', 'Phe', 'F', [ 'UUU', 'UUC' ]),
    AminoAcid('Proline',       'Pro', 'P', [ 'CCU', 'CCC', 'CCA', 'CCG' ]),
    AminoAcid('Serine',        'Ser', 'S', [ 'UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC' ]),
    AminoAcid('Threonine',     'Thr', 'T', [ 'ACU', 'ACC', 'ACA', 'ACG' ]),
    AminoAcid('Tryptophan',    'Trp', 'W', [ 'UGG' ]),
    AminoAcid('Tyrosine',      'Tyr', 'Y', [ 'UAU', 'UAC' ]),
    AminoAcid('Valine',        'Val', 'V', [ 'GUU', 'GUC', 'GUA', 'GUG' ]),
]
