def prediction_sensitivity(true_positives, false_negatives):
    return true_positives / (true_positives + false_negatives)


def prediction_precision(true_positives, false_positives):
    return true_positives / (true_positives + false_positives)


def hamming(strand_a: str, strand_b: str) -> int:
    """ Hamming distance """
    if len(strand_a) != len(strand_b):
        raise ValueError("len mismatch")

    return sum(a != b for a, b in zip(strand_a, strand_b))

