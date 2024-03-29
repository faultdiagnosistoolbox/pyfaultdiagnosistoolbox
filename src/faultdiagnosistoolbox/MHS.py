"""Minimal Hitting Set functionality."""

import numpy as np


def MHS(conflist) -> list[np.ndarray]:
    """Return Minimal Hitting Set for a set of sets.

    Parameters
    ----------
    conflist : list
        list of lists of integers

    Returns
    -------
    mmhs : list
        list of minimal hitting sets

    """
    mmhs = [np.array([], dtype=np.int64)]
    for c in conflist:
        j = 0
        while j < len(mmhs):
            if len(np.intersect1d(c, mmhs[j])) == 0:
                tmp = mmhs.pop(j)
                j = j - 1
                for ci in c:
                    newCand = np.concatenate((tmp, [ci]))
                    candmin = True
                    k = 0
                    while candmin and k < len(mmhs):
                        if len([x for x in mmhs[k] if x in newCand]) == len(mmhs[k]):
                            candmin = False
                            break
                        k = k + 1
                    if candmin:
                        mmhs.append(newCand)
            j = j + 1
    return mmhs
