import bioinfo.Sequence as Sequence
import bioinfo.data as data


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


def test_coding_regions():
    def _case(dna_strand, expected):
        sequence = Sequence.Sequence(dna_strand)
        assert sequence.coding_regions() == expected

    _case(data.STOP[0], [])
    _case(data.STOP[1], [])
    _case(data.STOP[2], [])
    _case('AAA' + data.STOP[0], [])
    _case(data.START[0] + 'AA' + data.STOP[0], [])
    _case(data.START[0] + 'AAA' + data.STOP[0], ['AUGAAAUAA'])
    _case(data.START[0] + data.STOP[0], [])
    _case(data.START[0] + data.START[0] + data.STOP[0], [])
    _case(data.START[0] + data.STOP[0] + data.STOP[0], [])
    _case(data.START[0] + 'AAA' + data.STOP[1] + data.START[0] + 'BBB' + data.STOP[2], ['AUGAAAUAG', 'AUGBBBUGA'])
    _case(data.START[0] + 'AAA' + data.STOP[0] + 'XXX' + data.START[0] + 'BBB' + data.STOP[1], ['AUGAAAUAA', 'AUGBBBUAG'])
