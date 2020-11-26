import bioinfo.Sequence as Sequence


def test_dna_to_rna():
    def _test(dna_strand, expected):
        sequence = Sequence.Sequence(dna_strand)
        assert sequence.to_rna() == expected
        assert sequence.rna_strand == expected

    _test('CAGCTGACTTTACTTCAGTA', 'CAGCUGACUUUACUUCAGUA')
    _test('CCTCATAAGC', 'CCUCAUAAGC')


def test_complement():
    def _test(dna_strand, expected):
        sequence = Sequence.Sequence(dna_strand)
        assert sequence.complement().strand == expected

    _test('CAATCGTTACGTTA', 'GTTAGCAATGCAAT')


def test_inverse():
    def _test(dna_strand, expected):
        sequence = Sequence.Sequence(dna_strand)
        assert sequence.invert().strand == expected

    _test('', '')
    _test('ACGT', 'TGCA')
