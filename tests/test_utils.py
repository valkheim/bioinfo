import bioinfo.utils as utils

def test_prediction_sensibility():
    # 4100 genes, 3500 found but 1200 errors
    assert round(utils.prediction_sensitivity(3500, 4100-3500), 2) == 0.85


def test_prediction_precision():
    # 4100 genes, 3500 found but 1200 errors
    assert round(utils.prediction_precision(3500, 1200), 2) == 0.74

