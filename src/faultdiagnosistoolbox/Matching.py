"""Compute matchings."""
import numpy as np
import faultdiagnosistoolbox.dmperm as dmperm
from dataclasses import dataclass, field


# class HallComponentMatching:
#     """Base class for Hall Component matching."""
#
#     matchType = ''
#     row = []
#     col = []
#     intState = []
#     derState = []
#
#     def __init__(self, t, r, c, derState=None, intState=None):
#         """Initialize matching."""
#         if intState is None:
#             intState = []
#         if derState is None:
#             derState = []
#         self.matchType = t
#         self.row = r
#         self.col = c
#         self.intState = intState
#         self.derState = derState
@dataclass
class HallComponentMatching:
    """Base class for Hall Component matching."""

    matchType: str
    row: list
    col: list
    intState: list = field(default_factory=lambda: [])
    derState: list = field(default_factory=lambda: [])


class Matching:
    """Base class for matching."""

    matchType = ''
    matching = []

    def __init__(self, Xmodel, eqs_in):
        """Initialize matching."""
        X = Xmodel.copy()

        eqs = np.array(eqs_in).copy()
        X = X[eqs, :]
        dm = dmperm.GetDMParts(X)
        if len(dm.M0eqs) != len(eqs):
            print("Error, matchings can only be computed for " +
                  "exactly determined systems")
        self.matching = []
        derCausal = False
        intCausal = False
        for hc in reversed(dm.M0):
            eqi = np.flipud(hc.row)  # Flip to go to computational order
            vi = np.flipud(hc.col)  # from matching order
            intEdge = np.any(X[eqi, :][:, vi] == 2)
            derEdge = np.any(X[eqi, :][:, vi] == 3)
            if not intEdge and not derEdge:  # Algebraic component
                self.matching.append(
                    HallComponentMatching('algebraic', eqs[eqi], vi))
            elif len(eqi) == 1 and derEdge:  # Trivial derivative component
                self.matching.append(
                    HallComponentMatching('der', eqs[eqi], vi, derState=vi))
                derCausal = True
            elif len(eqi) == 1 and intEdge:  # Trivial integral component
                self.matching.append(
                    HallComponentMatching('int', eqs[eqi], vi, intState=vi))
                intCausal = True
            else:  # Non-trivial component with dynamic constraints
                Gamma_mixed = MixedCausalityMatching(hc, X)
                if Gamma_mixed.matchType == 'int':
                    intCausal = True
                elif Gamma_mixed.matchType == 'der':
                    derCausal = True
                elif Gamma_mixed.matchType == 'mixed':
                    derCausal = True
                    intCausal = True
                self.matching.append(
                    HallComponentMatching(Gamma_mixed.matchType,
                                          eqs[Gamma_mixed.row],
                                          Gamma_mixed.col,
                                          derState=Gamma_mixed.derState,
                                          intState=Gamma_mixed.intState))

        if not derCausal and not intCausal:
            self.matchType = 'algebraic'
        elif derCausal and not intCausal:
            self.matchType = 'der'
        elif intCausal and not derCausal:
            self.matchType = 'int'
        else:
            self.matchType = 'mixed'


def MixedCausalityMatching(hallComponent, X):
    """Compute mixed causality matching."""
    # Create graph representation
    G = (hallComponent.row,
         hallComponent.col,
         X[hallComponent.row, :][:, hallComponent.col])

    Gamma = dmperm.MixedCausalityMatching(G)

    gamma_type = 'algebraic'
    derCausal = False
    intCausal = False
    intState = np.array([], dtype=np.int64)
    derState = np.array([], dtype=np.int64)

    for g in Gamma:
        if X[g[0], g[1]] == 2:  # Integrate
            intCausal = True
            intState = np.concatenate((intState, [g[1]]))
        elif X[g[0], g[1]] == 3:  # Differentiate
            derCausal = True
            derState = np.concatenate((derState, [g[1]]))
    if derCausal and intCausal:
        gamma_type = 'mixed'
    elif derCausal:
        gamma_type = 'der'
    elif intCausal:
        gamma_type = 'int'
    return HallComponentMatching(gamma_type,
                                 Gamma[:, 0],
                                 Gamma[:, 1],
                                 derState=derState,
                                 intState=intState)
