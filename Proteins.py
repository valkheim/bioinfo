from dataclasses import dataclass
from typing import List

import Sequence


@dataclass
class AminoAcid:
    full_name: str
    short_name: str
    codons: List[str]


AMINO_ACIDS = [
    AminoAcid('Alanine',       'Ala', [ 'GCU', 'GCC', 'GCA', 'GCG' ]),
    AminoAcid('Arginine',      'Arg', [ 'CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG' ]),
    AminoAcid('Asparagine',    'Asn', [ 'AAU', 'AAC' ]),
    AminoAcid('Aspartic acid', 'Asp', [ 'GAU', 'GAC' ]),
    AminoAcid('Cysteine',      'Cys', [ 'UGU', 'UGC' ]),
    AminoAcid('Glutamine',     'Gln', [ 'CAA', 'CAG' ]),
    AminoAcid('Glutamic acid', 'Glu', [ 'GAA', 'GAG' ]),
    AminoAcid('Glycine',       'Gly', [ 'GGU', 'GGC', 'GGA', 'GGU' ]),
    AminoAcid('Histidine',     'His', [ 'CAU', 'CAC' ]),
    AminoAcid('Isoleucine',    'Ile', [ 'AUU', 'AUC', 'AUA' ]),
    AminoAcid('Leucine',       'Leu', [ 'UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG' ]),
    AminoAcid('Lysine',        'Lys', [ 'AAA', 'AAG' ]),
    AminoAcid('Methionine',    'Met', [ 'AUG' ]),
    AminoAcid('Phenylalanine', 'Phe', [ 'UUU', 'UUC' ]),
    AminoAcid('Proline',       'Pro', [ 'CCU', 'CCC', 'CCA', 'CCG' ]),
    AminoAcid('Serine',        'Ser', [ 'UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC' ]),
    AminoAcid('Threonine',     'Thr', [ 'ACU', 'ACC', 'ACA', 'ACG' ]),
    AminoAcid('Tryptophane',   'Trp', [ 'UGG' ]),
    AminoAcid('Tyrosine',      'Tyr', [ 'UAU', 'UAC' ]),
    AminoAcid('Valine',        'Val', [ 'GUU', 'GUC', 'GUA', 'GUG' ]),
]


def find_amino_acid(codon):
    for amino_acid in AMINO_ACIDS:
        if codon in amino_acid.codons:
            return amino_acid


def get_codons(rna_strand):
    for i in range(0, len(rna_strand), 3):
        yield rna_strand[i:i + 3]


def get_protein(rna: Sequence.Sequence):
    STOP = [ 'UAA', 'UAG', 'UGA' ]
    amino_acids = []
    for codon in get_codons(rna.strand):
        if codon in STOP:
            return amino_acids

        amino_acid = find_amino_acid(codon)
        if amino_acid:
            amino_acids.append(amino_acid)

    return amino_acids
