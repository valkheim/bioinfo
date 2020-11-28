from . import Sequence
from . import Base
from . import Data


def get_amino_acid_from_base(base: Data.Base) -> Data.AminoAcid:
    for amino_acid in Data.AMINO_ACIDS:
        if amino_acid.letter == base.letter:
            return amino_acid


def get_amino_acid_from_codon(codon):
    for amino_acid in Data.AMINO_ACIDS:
        if codon in amino_acid.codons:
            return amino_acid


def get_codons(rna_strand):
    for i in range(0, len(rna_strand), 3):
        yield rna_strand[i:i + 3]


def get_protein(rna: Sequence.Sequence):
    amino_acids = []
    for codon in get_codons(rna.strand):
        if codon in Data.STOP:
            return amino_acids

        amino_acid = get_amino_acid_from_codon(codon)
        if amino_acid:
            amino_acids.append(amino_acid)

    return amino_acids
