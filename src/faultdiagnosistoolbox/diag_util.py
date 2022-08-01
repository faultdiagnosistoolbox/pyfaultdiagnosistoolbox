"""Some utility functions for the diagnosis toolbox."""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


def IsolabilityMatrix(fsm: np.ndarray) -> np.ndarray:
    """Compute isolability matrix based on a fault signature matrix"""
    nf = fsm.shape[1]
    im = np.ones((nf, nf), dtype=np.int)
    for ri in fsm:
        im[np.ix_(ri > 0, ri == 0)] = 0
    return im


def DiagnosesAndConfusionMatrix(data, residx=None) -> tuple[np.ndarray, np.ndarray]:
    """Compute consistency based diagnoses and corresponding confusion matrix
       based on a dataset."""
    if isinstance(residx, type(None)):
        dx = SingleFaultIsolability(data['res'], data['fsm'])
    else:
        dx = SingleFaultIsolability(data['res'][:, residx],
                                    data['fsm'][residx, :])

    nf = data['fsm'].shape[1]
    C = np.zeros((nf, nf))
    for fi in range(nf):
        for fj in range(nf):
            fjIdx = data['mode'] == fj
            C[fj, fi] = np.sum(dx[fjIdx, fi]) / np.sum(fjIdx)

    return dx, C


def SingleFaultIsolability(res, fsm) -> np.ndarray:
    """Compute single fault consistency based diagnoses."""
    N = res.shape[0]
    nf = fsm.shape[1]
    dx = np.zeros((N, nf))
    for k, rk in enumerate(res):
        alarm = (rk > 0.5)
        if np.any(alarm):
            dx[k] = np.all(fsm[alarm], axis=0)
        else:
            dx[k] = np.hstack((True, np.zeros(nf - 1)))
    return dx


def PlotConfusionMatrix(C, ax=None):
    """Plot a confusion matrix in a suitable colormap."""

    if not ax:
        ax = plt.gca()

    nf = C.shape[0]
    ax.imshow(C, cmap=summer_cmap)
    for fi in range(nf):
        for fj in range(nf):
            ax.text(fi, fj, f'{(C[fj, fi] * 100):.1f}',
                    ha='center', va='center', color='k')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_xticks(np.arange(C.shape[1] + 1) - .5, minor=True)
    ax.set_yticks(np.arange(C.shape[0] + 1) - .5, minor=True)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)


def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)


summer_cmap = make_colormap([
    (0.000, 0.500, 0.400),
    (0.016, 0.508, 0.400),
    (0.032, 0.516, 0.400),
    (0.048, 0.524, 0.400),
    (0.063, 0.532, 0.400),
    (0.079, 0.540, 0.400),
    (0.095, 0.548, 0.400),
    (0.111, 0.556, 0.400),
    (0.127, 0.563, 0.400),
    (0.143, 0.571, 0.400),
    (0.159, 0.579, 0.400),
    (0.175, 0.587, 0.400),
    (0.190, 0.595, 0.400),
    (0.206, 0.603, 0.400),
    (0.222, 0.611, 0.400),
    (0.238, 0.619, 0.400),
    (0.254, 0.627, 0.400),
    (0.270, 0.635, 0.400),
    (0.286, 0.643, 0.400),
    (0.302, 0.651, 0.400),
    (0.317, 0.659, 0.400),
    (0.333, 0.667, 0.400),
    (0.349, 0.675, 0.400),
    (0.365, 0.683, 0.400),
    (0.381, 0.690, 0.400),
    (0.397, 0.698, 0.400),
    (0.413, 0.706, 0.400),
    (0.429, 0.714, 0.400),
    (0.444, 0.722, 0.400),
    (0.460, 0.730, 0.400),
    (0.476, 0.738, 0.400),
    (0.492, 0.746, 0.400),
    (0.508, 0.754, 0.400),
    (0.524, 0.762, 0.400),
    (0.540, 0.770, 0.400),
    (0.556, 0.778, 0.400),
    (0.571, 0.786, 0.400),
    (0.587, 0.794, 0.400),
    (0.603, 0.802, 0.400),
    (0.619, 0.810, 0.400),
    (0.635, 0.817, 0.400),
    (0.651, 0.825, 0.400),
    (0.667, 0.833, 0.400),
    (0.683, 0.841, 0.400),
    (0.698, 0.849, 0.400),
    (0.714, 0.857, 0.400),
    (0.730, 0.865, 0.400),
    (0.746, 0.873, 0.400),
    (0.762, 0.881, 0.400),
    (0.778, 0.889, 0.400),
    (0.794, 0.897, 0.400),
    (0.810, 0.905, 0.400),
    (0.825, 0.913, 0.400),
    (0.841, 0.921, 0.400),
    (0.857, 0.929, 0.400),
    (0.873, 0.937, 0.400),
    (0.889, 0.944, 0.400),
    (0.905, 0.952, 0.400),
    (0.921, 0.960, 0.400),
    (0.937, 0.968, 0.400),
    (0.952, 0.976, 0.400),
    (0.968, 0.984, 0.400),
    (0.984, 0.992, 0.400),
    (1.000, 1.000, 0.400)])
