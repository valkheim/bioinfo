import bioinfo.Sequence as Sequence


def test_dna_to_rna():
    def _test(dna_strand, expected):
        sequence = Sequence.Sequence(dna_strand)
        assert sequence.to_rna() == expected
        assert sequence.rna_strand == expected

    _test('CAGCTGACTTTACTTCAGTA', 'CAGCUGACUUUACUUCAGUA')
    _test('CCTCATAAGC', 'CCUCAUAAGC')