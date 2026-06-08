import numpy as np


def calc_stats(vec):
    n = len(vec)
    if n == 1:
        return np.mean(vec), 0
    return np.mean(vec), np.sqrt(np.var(vec) / (n - 1))