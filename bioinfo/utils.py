def prediction_sensitivity(true_positives, false_negatives):
    return true_positives / (true_positives + false_negatives)


def prediction_precision(true_positives, false_positives):
    return true_positives / (true_positives + false_positives)


def hamming(strand_a: str, strand_b: str) -> int:
    """ Hamming distance """
    if len(strand_a) != len(strand_b):
        raise ValueError("len mismatch")

    return sum(a != b for a, b in zip(strand_a, strand_b))


def substitution_cost(i, j):
    scale = "ACGT"
    #              A    C    G    T
    matrix = [ [   0,   1, 0.5,   1 ],  # A
             [     1,   0,   1, 0.5 ],  # C
             [   0.5,   1,   0,   1 ],  # G
             [     1, 0.5,   1,   0 ] ] # T

    return matrix[scale.index(i)][scale.index(j)]


def insertion_cost():
    return 1


def compute_cost(sequence_a, sequence_b, i, j) -> int:
    # Needleman-Wunch algorithm
    # version iterative plus efficace:
    # calculer tous les couts puis backtrack pour trouver le meilleur chemin
    if i == 0 and j == 0:
        return 0

    if i == 0 and j > 0:
        return j * insertion_cost()

    if i > 0 and j == 0:
        return i * insertion_cost()

    return min(
        compute_cost(sequence_a, sequence_b, i - 1, j - 1) + substitution_cost(sequence_a[i - 1].letter, sequence_b[j - 1].letter),
        compute_cost(sequence_a, sequence_b, i, j - 1) + insertion_cost(),
        compute_cost(sequence_a, sequence_b, i - 1, j) + insertion_cost()
    )
