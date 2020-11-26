import bioinfo.Sequence as Sequence
import bioinfo.Proteins as Proteins


def _get_amino_acid_names(protein):
    return [ aa.full_name for aa in protein ]


def test_dna_to_protein():
    sequence = Sequence.Sequence('GAAAAG')
    protein = Proteins.get_protein(sequence)
    assert _get_amino_acid_names(protein) == [ 'Glutamic acid', 'Lysine' ]


def test_dna_to_protein_with_stop_triplet():
    sequence = Sequence.Sequence('GCCUGA')
    protein = Proteins.get_protein(sequence)
    assert _get_amino_acid_names(protein) == [ 'Alanine' ]


def test_get_codons():
    sequence = Sequence.Sequence('GAAATTCTT')
    assert list(Proteins.get_codons(sequence.strand)) == ['GAA', 'ATT', 'CTT']


def test_get_amino_acid_from_codon():
    assert Proteins.get_amino_acid_from_codon('AAA') == Proteins.AminoAcid(full_name='Lysine', short_name='Lys', letter='K', codons=['AAA', 'AAG'])
    assert Proteins.get_amino_acid_from_codon('UAG') == None  # STOP ambre
    assert Proteins.get_amino_acid_from_codon('AUG') == Proteins.AminoAcid(full_name='Methionine', short_name='Met', letter='M', codons=['AUG'])


def test_get_amino_acid_map():
    sequence = Sequence.Sequence('ARNDCQEGHILKMFPSTWYV')
    amino_acids = []
    sequence.walk_bases(lambda b: amino_acids.append([ b.letter, Proteins.get_amino_acid_from_base(b).full_name ]))
    expected = [
        [ 'A', 'Alanine' ],
        [ 'R', 'Arginine' ],
        [ 'N', 'Asparagine' ],
        [ 'D', 'Aspartic acid' ],
        [ 'C', 'Cysteine' ],
        [ 'Q', 'Glutamine' ],
        [ 'E', 'Glutamic acid' ],
        [ 'G', 'Glycine' ],
        [ 'H', 'Histidine' ],
        [ 'I', 'Isoleucine' ],
        [ 'L', 'Leucine' ],
        [ 'K', 'Lysine' ],
        [ 'M', 'Methionine' ],
        [ 'F', 'Phenylalanine' ],
        [ 'P', 'Proline' ],
        [ 'S', 'Serine' ],
        [ 'T', 'Threonine' ],
        [ 'W', 'Tryptophan' ],
        [ 'Y', 'Tyrosine' ],
        [ 'V', 'Valine' ]
    ]
    assert amino_acids == expected

