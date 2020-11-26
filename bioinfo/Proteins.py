from dataclasses import dataclass
from typing import List

from . import Sequence
from . import Base


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


@dataclass
class AminoAcid:
    full_name: str
    short_name: str
    letter: str
    codons: List[str]


AMINO_ACIDS = [  # Expresssed with latin alphabet without B, J, O, U, X, Z
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


def get_amino_acid_from_base(base: Base.Base) -> AminoAcid:
    for amino_acid in AMINO_ACIDS:
        if amino_acid.letter == base.letter:
            return amino_acid


def get_amino_acid_from_codon(codon):
    for amino_acid in AMINO_ACIDS:
        if codon in amino_acid.codons:
            return amino_acid


def get_codons(rna_strand):
    for i in range(0, len(rna_strand), 3):
        yield rna_strand[i:i + 3]


def get_protein(rna: Sequence.Sequence):
    amino_acids = []
    for codon in get_codons(rna.strand):
        if codon in STOP:
            return amino_acids

        amino_acid = get_amino_acid_from_codon(codon)
        if amino_acid:
            amino_acids.append(amino_acid)

    return amino_acids
