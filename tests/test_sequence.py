import bioinfo.Sequence as Sequence
import bioinfo.Data as Data


def test_strand():
    def _test(sequence, dna_strand, rna_strand, length):
        assert sequence.strand == dna_strand
        assert sequence.length == length

    _test(Sequence.Sequence('CAA'), 'CAA', 'CAA', 3)
    _test(Sequence.Sequence('ATG'), 'ATG', 'AUG', 3)


def test_strand_concat():
    def _test(sequence, strand, bases):
        assert sequence.strand == strand
        assert sequence.bases == bases

    sequence = Sequence.Sequence('A')
    _test(sequence, 'A', [ Data.Base(name='A', letter='A', index=0) ])
    sequence.append('T')
    _test(sequence, 'AT', [ Data.Base(name='A', letter='A', index=0), Data.Base(name='T', letter='T', index=1) ])


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

    _case(Data.STOP[0], [])
    _case(Data.STOP[1], [])
    _case(Data.STOP[2], [])
    _case('AAA' + Data.STOP[0], [])
    _case(Data.START[0] + 'AA' + Data.STOP[0], [])
    _case(Data.START[0] + 'AAA' + Data.STOP[0], ['AUGAAAUAA'])
    _case(Data.START[0] + Data.STOP[0], [])
    _case(Data.START[0] + Data.START[0] + Data.STOP[0], [])
    _case(Data.START[0] + Data.STOP[0] + Data.STOP[0], [])
    _case(Data.START[0] + 'AAA' + Data.STOP[1] + Data.START[0] + 'BBB' + Data.STOP[2], ['AUGAAAUAG', 'AUGBBBUGA'])
    _case(Data.START[0] + 'AAA' + Data.STOP[0] + 'XXX' + Data.START[0] + 'BBB' + Data.STOP[1], ['AUGAAAUAA', 'AUGBBBUAG'])
