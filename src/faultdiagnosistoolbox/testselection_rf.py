"""Test selection using random forest machine learning models."""

import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix
from .diag_util import IsolabilityMatrix, DiagnosesAndConfusionMatrix


def RandomForestTestSelection(data, n_estimators=100):
    """Select residuals for fault isolation using Random Forests.

    Implementation of the residual selection approach described in the
    paper

      Frisk, Erik, and Mattias Krysander. "Residual Selection for
      Consistency Based Diagnosis Using Machine Learning Models."
      IFAC Safeprocess, Warszaw, Poland (2018)
      (https://doi.org/10.1016/j.ifacol.2018.09.547)

    The method computes a sorted list of residuals, where the most important
    residuals are first in the list. The method also computes performance
    indices that can be used to determine how many in the sorted list of
    residuals that should be selected.

    Parameters
    ----------
    data : dict
        Dictionary with data and supporting information needed by the algorithm. The dictionary must have the following keys
        ``modes``: an array of size m with names of no-fault and fault modes where the no-fault mode is assumed the first element in the list.
        ``res``: an (N x n) numpy array with N samples of the n residuals.
        ``mode``: an (N x 1) numpy array with the fault mode at each sample where a value of 0 corresponds to fault free operation.
        ``fsm``: a fault signature matrix of size (n x m)
    n_estimators : int
        Number of decision trees in the ensemble model. Defaults to 100

    Returns
    -------
    result : dict
        A dictionary containing the resulting test selection.
        The dictionary contains the keys
        ``sortidx``: A sorted list of indices to residuals with the most important residuals first.
        ``pfa``  : Performance curve for the false alarm probability
        ``pmd``  : Performace curve for missed detection
        ``pfi``  : Performance curve for fault isolation
        ``pmfi`` : Maximal isolation performance per fault mode
        ``residualimportance`` : List of relative residual importance. The list is not sorted.
    C : arraylike
        A confusion matrix when computing consistency based diagnoses using all residuals.
    rf : RandomForestClassifier
        The random forest object
    Crf : arraylike
        The confusion matrix for the random forest classifier.

    Example
    -------
    Let data be thresholded residual data

    >>> ts, C, rf, Crf = RandomForestTestSelection(data, n_estimators=100)

    Then plot performance as a function of selected number of residuals using
    >>> plt.plot(range(1, n), ts['pfa'], 'r', label='False alarm probability')
    >>> plt.plot(range(1, n), ts['pmd'], 'b',
    >>>          label='Missed detection probability')
    >>> plt.plot(range(1, n), ts['pfi'], 'y',
    >>>          label='False isolation probability')

    based on the plots, determine ntest, the number of tests to select, and
    then plot the corresponding fault isolation performance

    >>> _, C_select = DiagnosesAndConfusionMatrix(data, residx=ts['sortidx'][0:ntests])
    >>> PlotConfusionMatrix(C_select)

    that could be compared to

    >>> PlotConfusionMatrix(C)

    """

    nf = len(data["modes"])
    nr = data["res"].shape[1]
    im = IsolabilityMatrix(data["fsm"])

    rf = RandomForestClassifier(n_estimators=n_estimators)
    rf.fit(data["res"], data["mode"])
    sortIdx = np.argsort(rf.feature_importances_)[::-1]
    residualImportance = rf.feature_importances_

    # Compute classifiers confusion matrix (needed?)
    s = np.diag([1 / sum(data["mode"] == mi) for mi in range(nf)])
    Crf = s @ confusion_matrix(data["mode"], rf.predict(data["res"]))

    pfa = np.zeros(nr - 1)
    pmd = np.zeros(nr - 1)
    pfi = np.zeros(nr - 1)
    pmfi = np.zeros((nr - 1, nf))

    # Make sure isolability matrix for NF corresponds to diagnosis
    # statement computed by DiagnosesAndConfusionMatrix
    imk = IsolabilityMatrix(data["fsm"])
    imk[0] = np.zeros(nf)
    imk[0, 0] = 1
    C = None
    for k in range(1, nr):
        dx, C = DiagnosesAndConfusionMatrix(data, residx=sortIdx[0:k])
        pfa[k - 1] = 1 - C[0, 0]
        pmd[k - 1] = np.mean(C[1:, 0])
        pfi[k - 1] = np.mean(np.diag(1 - C[1:, 0]) @ np.abs(C[1:, 1:] - im[1:, 1:]))

        for fi in range(nf):
            pmfi[k - 1, fi] = np.mean(np.all(dx[data["mode"] == fi] == imk[fi, :], axis=1))
    return (
        {
            "sortidx": sortIdx,
            "pfa": pfa,
            "pmd": pmd,
            "pfi": pfi,
            "pmfi": pmfi,
            "residualimportance": residualImportance,
        },
        C,
        rf,
        Crf,
    )
