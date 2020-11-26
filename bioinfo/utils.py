def prediction_sensitivity(true_positives, false_negatives):
    return true_positives / (true_positives + false_negatives)

def prediction_precision(true_positives, false_positives):
    return true_positives / (true_positives + false_positives)
