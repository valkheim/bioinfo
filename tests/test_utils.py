import bioinfo.Sequence as Sequence
import bioinfo.utils as utils

def test_prediction_sensibility():
    # 4100 genes, 3500 found but 1200 errors
    assert round(utils.prediction_sensitivity(3500, 4100-3500), 2) == 0.85


def test_prediction_precision():
    # 4100 genes, 3500 found but 1200 errors
    assert round(utils.prediction_precision(3500, 1200), 2) == 0.74


def test_hamming_distance():
    def _case(strand_a, strand_b, expected_distance):
        assert utils.hamming(strand_a, strand_b) == expected_distance

    _case('', '', 0)
    _case('A', 'A', 0)
    # From lesson:
    _case(
            #   #          #
        'ACCTCTGTATCTATTCGGCATCATCAT',
        'ACCCCTGAATCTATTCGGGATCATCAT',
        3
    )
    _case(
                #      #         #
        'ACCTCTGTATCTATTCGGGATCATCAT',
        'ACCTCTGAATCTATCCGGGATCATGAT',
        3
    )
    # From exercices:
    _case(
        'TTGCATTGCTTAGGCATA',
        'TTGCGTTGCTTAGCCATA',
        2
    )
    _case(
        'TTGCATTGCTTAGGCATA',
        'TTGCAGTCCTTAGGCATT',
        3
    )
    _case(
        'TTGCGTTGCTTAGCCATA',
        'TTGCAGTCCTTAGGCATT',
        5
    )


def test_substitution_cost():
    s1 = Sequence.Sequence('TTCACTATGTA')
    s2 = Sequence.Sequence('TTAGCTATGTA')
    assert utils.substitution_cost(s1[2].letter, s2[2].letter) == 1
    assert utils.substitution_cost(s1[3].letter, s2[3].letter) == 0.5
    assert utils.substitution_cost(s1[1].letter, s2[6].letter) == 1
    assert utils.substitution_cost(s1[5].letter, s2[7].letter) == 0


def test_cost():
    s1 = Sequence.Sequence('TTCACTATGTA')
    s2 = Sequence.Sequence('TTAGCTATGTA')
    assert s1.length == 11
    assert s1.length == s2.length
    #assert utils.compute_cost(s1, s2, s1.length, s2.length) == 1.5  #  too much complexity here
