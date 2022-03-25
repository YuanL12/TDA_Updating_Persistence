import bats
import numpy as np

def kendall_tau_distance(vals1, vals2):
    """
    compute kendall-tau distance between sort permutations for vals1, vals2
    """
    kd = 0
    nmax = 0
    for v1, v2 in zip(vals1, vals2):
        perm1 = np.argsort(v1, kind='stable')
        perm2 = np.argsort(v2, kind='stable')

        nperm = len(perm1)
        maxswaps = (nperm * (nperm-1)) // 2

        kt = bats.kendall_tau(perm1, perm2)
        kd += kt
        nmax += maxswaps

    return kd, nmax


def normalized_kt(vals1, vals2):
    """
    normalized kendal-tau distance
    """
    kd, nmax = kendall_tau_distance(vals1, vals2)
    return kd / nmax
