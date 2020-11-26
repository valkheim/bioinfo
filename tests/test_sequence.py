import bioinfo.Sequence as Sequence


def test_dna_to_rna():
    def _case(dna_strand, expected):
        sequence = Sequence.Sequence(dna_strand)
        assert sequence.to_rna() == expected
        assert sequence.rna_strand == expected

    _case('CAGCTGACTTTACTTCAGTA', 'CAGCUGACUUUACUUCAGUA')
    _case('CCTCATAAGC', 'CCUCAUAAGC')


def test_complement():
    def _case(dna_strand, expected):
        sequence = Sequence.Sequence(dna_strand)
        assert sequence.complement().strand == expected

    _case('CAATCGTTACGTTA', 'GTTAGCAATGCAAT')


def test_inverse():
    def _case(dna_strand, expected):
        sequence = Sequence.Sequence(dna_strand)
        assert sequence.invert().strand == expected

    _case('', '')
    _case('ACGT', 'TGCA')


def test_inverse_complement():
    def _case(dna_strand, expected):
        sequence = Sequence.Sequence(dna_strand)
        assert sequence.complement().invert().strand == expected

    _case('TTGCATACGTAAGCA', 'TGCTTACGTATGCAA')
