"""Main class for DiagnosisModel."""
import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
import copy
import faultdiagnosistoolbox.dmperm as dmperm
import faultdiagnosistoolbox.Matching as match
import faultdiagnosistoolbox.StructurePlotting as smplot
import faultdiagnosistoolbox.CodeGeneration as codegen
import faultdiagnosistoolbox.SensorPlacement as sensplace
import faultdiagnosistoolbox.TestSelection as testselection

from faultdiagnosistoolbox.VarIdGen import VarIdGen
import sys


class DiagnosisModel(object):
    """Base class for diagnosis models."""

    def __init__(self, modeldef, name=''):
        """Initialize DiagnosisModel object."""
        self.X = np.zeros((0, 0), dtype='float64')
        self.F = np.zeros((0, 0), dtype='float64')
        self.Z = np.zeros((0, 0), dtype='float64')
        self.x = []
        self.f = []
        self.z = []
        self.name = name
        self.parameters = []

        self.vGen = VarIdGen()
        self.modelType = 'Structural'
        if modeldef['type'] is 'VarStruc' or modeldef['type'] is 'Symbolic':
            self.X = _ModelStructure(modeldef['rels'], modeldef['x'])
            self.x = modeldef['x']
            self.F = _ModelStructure(modeldef['rels'], modeldef['f'])
            self.f = modeldef['f']
            self.Z = _ModelStructure(modeldef['rels'], modeldef['z'])
            self.z = modeldef['z']
            self.e = list(map(lambda x: self.vGen.NewE(), np.arange(0, self.ne())))
            if modeldef['type'] is 'Symbolic':
                self.modelType = 'Symbolic'

            if 'parameters' in modeldef:
                self.parameters = modeldef['parameters']
        elif modeldef['type'] is 'MatrixStruc':
            self.X = np.array(modeldef['X'], dtype=np.int64)
            ne = self.X.shape[0]
            if len(modeldef['F']) > 0:
                self.F = np.array(modeldef['F'], dtype=np.int64)
            else:
                self.F = np.zeros((ne, 0), dtype=np.int64)

            if len(modeldef['Z']) > 0:
                self.Z = np.array(modeldef['Z'], dtype=np.int64)
            else:
                self.Z = np.zeros((ne, 0), dtype=np.int64)

            if 'x' in modeldef:
                self.x = modeldef['x']
            else:
                self.x = list(map(lambda x: "x" + np.str(x + 1),
                                  np.arange(0, self.X.shape[1])))
            if 'f' in modeldef:
                self.f = modeldef['f']
            else:
                self.f = list(map(lambda x: "f" + np.str(x + 1),
                                  np.arange(0, self.F.shape[1])))
            if 'z' in modeldef:
                self.z = modeldef['z']
            else:
                self.z = list(map(lambda x: "z" + np.str(x + 1),
                                  np.arange(0, self.Z.shape[1])))

            self.e = list(map(lambda x: "e" + np.str(x + 1),
                              np.arange(0, self.ne())))
        else:
            print('Model definition type ' + modeldef['type'] +
                  ' is not supported (yet)')

        if modeldef['type'] is 'Symbolic':
            self.syme = np.array(_ToEquations(modeldef['rels']))

        if np.any(np.sum(self.X > 1, 0) > 1):
            print("The model has higher order derivatives, " +
                  "please rewrite as a set of first order " +
                  "differential equations")

        self.P = np.arange(0, len(self.x))
        self.Pfault = []

    def copy(self):
        """Return a new copy of the model object."""
        return copy.deepcopy(self)

    def ne(self):
        """Return number of equations in model."""
        return self.X.shape[0]

    def nx(self):
        """Return number of unknown variables in model."""
        return self.X.shape[1]

    def nf(self):
        """Return number of fault variables in model."""
        return self.F.shape[1]

    def nz(self):
        """Return number of known variables in model."""
        return self.Z.shape[1]

    def GetDMParts(self):
        """Return Dulmage-Mendelsohn decomposition of structure for unknown variables."""
        return dmperm.GetDMParts(self.X)

    def MSO(self):
        """Return the set of MSO sets.

        For details of the algorithm see the journal publication
        Krysander, Mattias, Jan Aslund, and Mattias Nyberg.
        "An efficient algorithm for finding minimal overconstrained
        subsystems for model-based diagnosis."
        Systems, Man and Cybernetics, Part A: Systems and Humans,
        IEEE Transactions on 38.1 (2008): 197-206.
        """
        return dmperm.MSO(self.X)

    def MTES(self):
        """Return the set of MTES sets.

        For details of the algorithm see the paper
        A Structural Algorithm for Finding Testable Sub-models and
        Multiple Fault Isolability Analysis. Mattias Krysander, Jan Åslund,
        and Erik Frisk (2010). 21st International Workshop on Principles of
        Diagnosis (DX-10). Portland, Oregon, USA.
        """
        return dmperm.MTES(self)

    def TestSelection(self, arr, method='aminc', isolabilitymatrix=''):
        """A minimal hitting set based test selection.

        Find sets of tests, based on a set of equations or a fault sensitivity
        matrix, FSM, that achieves isolability performance specifications

        Simple test selection strategy that finds sets of tests that, ideally,
        fulfills specified isolability performance specifications. Note that
        this is a purely structural method, no noise or model uncertainty
        considerations are used.

        Input
        -----
          arr               : An array of arrays, interpreted as a set of equations
                              sets used to design residuals. For eample the ouput of
                              the MSO or MTES functions
          isolabilitymatrix : Matrix specifying required isolability
                              performance. A 0 in position (i,j) represents that
                              fault i is isolable from fault j, a 1 indicates
                              that fault i is not isolable from fault j. A fault
                              can not be isolable from itself and therefore must
                              the diagonal always be 1.
          method            : Choice of test seleciton method.
                              'aminc' -  Searches for a subset minimal sets of tests
                                         that fulfills requirements (default)
                                         Uses aminc, an approximative minimal cardinality
                                         hitting set approach from:
                                         Cormen, L., Leiserson, C. E., and
                                         Ronald, L. (1990). Rivest,
                                         "Introduction to Algorithms.", 1990.

                                        Information also in De Kleer, Johan. "Hitting set
                                        algorithms for model-based diagnosis."
                                        22th International Workshop on Principles
                                        of Diagnosis, DX, 2011.

                              'full' - Finds all subset minimal sets of tests
                                       that fulfills requirements. Warning,
                                       might easily lead to computationally
                                       intractable problems.
        """
        return testselection.TestSelection(self, arr, method, isolabilitymatrix)

    def srank(self):
        """Return the structural rank of the incidence matrix for the unknown variables."""
        return dmperm.srank(self.X)

    def IsPSO(self, *args):
        """Return True if the model is a PSO set.

        Input
        -----
          eq (optional) : Set of equations to analyze. (default: all equations)
        """
        if len(args) > 0:
            eq = args[0]
        else:
            eq = np.arange(0, self.X.shape[0])

        dm = dmperm.GetDMParts(self.X[eq, :])
        return (len(dm.Mm.row) == 0) and (len(dm.M0) == 0)

    def Structural(self):
        """Convert the model to a structural model."""
        self.modelType = 'Structural'
        self.syme = []

    def LumpDynamics(self):
        """Lump the dynamics variables.

        If model is symbolic, it will be converted into a structural model
        """
        if self.modelType is 'Symbolic':
            self.Structural()

        X = self.X
        F = self.F
        Z = self.Z

        dvars = np.where(np.any(X == 3, axis=0))[0]
        diffEq = self.DifferentialConstraints()
        for e in diffEq:
            iv = np.where(X[e, :] == 2)[0][0]
            dv = np.where(X[e, :] == 3)[0][0]
            X[:, iv] = np.logical_or(X[:, iv] > 0, X[:, dv] > 0).astype(np.int64)

        self.X = np.delete(np.delete(X, dvars, axis=1), diffEq, axis=0)
        self.x = [self.x[xi] for xi in np.arange(0, len(self.x)) if not (xi in dvars)]
        self.F = np.delete(F, diffEq, axis=0)
        self.Z = np.delete(Z, diffEq, axis=0)
        self.e = [self.e[ei] for ei in np.arange(0, len(self.e)) if not (ei in diffEq)]
        self.parameters = []
        self.syme = []
        self.Pfault = []
        self.P = np.arange(0, len(self.x))

    def IsObservable(self, eq=[]):
        """Return true if the model is observable.

        Input
        -----
          eq : Set of equations to analyze. (default: all equations)
        """
        return dmperm.IsObservable(self.X, eq)

    def Pantelides(self):
        """Compute structural index and differentiation vector for exactly determined models

        idx, nu = model.Pantelides()

        Returns
        -----
          idx : Structural index of model
          nu  : Differentiation vector as defined in
                Pantelides, Constantinos C. "The consistent initialization of
                differential-algebraic systems." SIAM Journal on Scientific and
                Statistical Computing 9.2 (1988): 213-231.
        """
        X = self.X.copy()
        ne, nx = X.shape
        if dmperm.srank(X) < nx:
            print("Error: Pantelides can only be called for square, " +
                    "exactly determined models")
            return None
        der = np.zeros((ne, 1), dtype=np.int)
        X[X == 3] = 1
        X = X - 1
        cmatching = False
        while not cmatching:
            Xder = X + der@np.ones((1, nx))
            Xder[X < 0] = -1

            hd = np.max(Xder, axis=0)
            hod = (Xder == (np.ones((ne, 1))@hd.reshape((1, -1))))
            if dmperm.srank(hod) == nx:
                cmatching = True
            else:
                dm = dmperm.GetDMParts(hod)
                der[dm.Mp.row] = der[dm.Mp.row] + 1
        idx = np.max(der) + np.any(hd == 0)
        return (idx, der.reshape(-1))

    def IsHighIndex(self, eq=[]):
        """Return true if the model is high structural differential index.

        Input
        -----
          eq : Set of equations to analyze. (default: all equations)
        """
        if len(eq) == 0:
            eq = np.arange(0, self.X.shape[0])
        return dmperm.IsHighIndex(self.X, eq)

    def IsLowIndex(self, eq=[]):
        """Return true if the model is low structural differential index.

        Input
        -----
          eq : Set of equations to analyze. (default: all equations)
        """
        if len(eq) == 0:
            eq = np.arange(0, self.X.shape[0])
        return dmperm.IsLowIndex(self.X, eq)

    def IsUnderdetermined(self):
        """Return true if the model has underdetermined parts."""
        dm = dmperm.GetDMParts(self.X)
        return len(dm.Mm.row) > 0

    def DifferentialConstraints(self):
        """Return indices to the differential constraints."""
        return np.where(np.any(self.X == 2, axis=1))[0]

    def DynamicVariables(model):
        """Return variables and index to dynamic variables in model."""
        idx = np.where(np.any(model.X == 2, axis=0))[0]
        c = list(np.array(model.x)[idx])
        return (c, idx)

    def AlgebraicVariables(model):
        """Return variables and index to algebraic variables in model."""
        dyn_idx = np.where(np.any(model.X == 2, axis=0))[0]
        idx = [x for x in range(model.nx()) if not(x in dyn_idx)]
        c = list(np.array(model.x)[idx])
        return (c, idx)

    def FSM(self, eqs_sets, plot=False):
        """Return the fault signature matrix for a set of equation sets.

        Input
        -----
          eqs_sets : list of sets of equations, e.g., MSO sets
          Plot     : If True, plot the fault signature matrix (dafault: False)
        """
        r = np.zeros((len(eqs_sets), self.F.shape[1]), dtype=np.int64)
        for idx, eqs in enumerate(eqs_sets):
            r[idx, :] = np.any(self.F[eqs, :], axis=0)

        if plot:
            plt.spy(r, markersize=10, marker='o')
            plt.xticks(np.arange(0, self.nf()), self.f)
            plt.yticks(np.arange(0, len(eqs_sets)),
                       ["eq. set " + str(k + 1) for k in np.arange(0, len(eqs_sets))])
            plt.gca().xaxis.tick_bottom()

        return r

    def Matching(self, eqs):
        """Return a matching for a set of equations.

        Input
        -----
          eqs : Set of equations
        """
        return match.Matching(self.X, np.array(eqs))

    def MSOCausalitySweep(self, mso, diffres='int', causality=''):
        """Determine casuality for MSO set residual generator.

        Determine causality for sequential residual generator for
        each n residual equations for a given MSO set.

        Input
        -----
          mso       : list of equations
          diffres   : Can be 'int' or 'der' (default 'int'). Determines how
                      to treat differential constraints when used as a
                      residual equation.
          causality : Can be 'int' or 'der'. When causality is specified,
                      the call returns a boolean vector indicating if it is
                      possible to realize the residual generator in
                      derivative or integral causality respectively with the
                      corresponding equations as residual equation. If
                      this option is given, the diffres key have no
                      effect.
        """
        def IsDifferentialStructure(X, eq):
            return np.any(X[eq, :] == 3)

        if (causality is 'der') and (not (diffres is 'der')):
            diffres = 'der'
        if (causality is 'int') and (not (diffres is 'int')):
            diffres = 'int'
        res = []
        X = self.X
        for red in mso:
            m0 = np.sort([e for e in mso if e != red])
            Gamma = self.Matching(m0)

            if not IsDifferentialStructure(X, red):
                res.append(Gamma.matchType)
            elif diffres is 'int':
                if Gamma.matchType is 'der':
                    res.append('mixed')
                elif Gamma.matchType is 'int' or Gamma.matchType is 'mixed':
                    res.append(Gamma.matchType)
                elif Gamma.matchType is 'algebraic':
                    res.append('int')
            elif diffres is 'der':
                if Gamma.matchType is 'int':
                    res.append('mixed')
                elif Gamma.matchType is 'der' or Gamma.matchType is 'mixed':
                    res.append(Gamma.matchType)
                elif Gamma.matchType is 'algebraic':
                    res.append('der')
        if len(causality) > 0:
            res = np.array(list(map(lambda c: c is causality or c is 'algebraic', res)))
        return res

    def MeasurementEquations(self, m):
        """Return indices to measurement equations.

        Input
        -----
          m : List of names for measurement variables (strings)
        """
        mIdx = np.array([self.z.index(zi) for zi in m if zi in self.z])
        if len(mIdx) > 0:
            return np.where(np.any(self.Z[:, mIdx], axis=1))[0]
        else:
            return []

    def SubModel(self, m, v, clear=True, remove=False, verbose=False):
        """Extract submodel.

        Input
        -----
          eqs     : Set of indices to or logicals for equations to keep/remove
          vars    : Set of indices to or logicals for variables to keep/remove
          clear   : If true, non used varaibles in the submodel will be
                    eliminated (default: true)
          verbose : If true, verbose output (default: false)
          remove  : If true, supplied equations are removed instead of kept (default: false)
        """
        if remove:
            eqs = np.array([ei for ei in np.arange(0, self.ne()) if not (ei in m)])
            var = np.array([vi for vi in np.arange(0, self.nx()) if not (vi in v)])
        else:
            eqs = m
            var = np.arange(0, self.nx())

        if clear:
            xIdx = var[np.any(self.X[eqs, :][:, var], axis=0)]
            fIdx = np.argwhere(np.any(self.F[eqs, :], axis=0)).flatten()
            zIdx = np.argwhere(np.any(self.Z[eqs, :], axis=0)).flatten()
        else:
            xIdx = var
            fIdx = np.arange(0, self.nf())
            zIdx = np.arange(0, self.nz())

        if verbose and clear:
            cxIdx = [xi for xi in var if not (xi in xIdx)]
            cfIdx = [fi for fi in np.arange(0, self.nf()) if not (fi in fIdx)]
            czIdx = [zi for zi in np.arange(0, self.nf()) if not (zi in zIdx)]
            if len(cxIdx) > 0:
                sys.stdout.write('Removing unknown variables: ')
                for xi in cxIdx:
                    sys.stdout.write(self.x[xi] + '  ')
                sys.stdout.write('\n')
            if len(cfIdx) > 0:
                sys.stdout.write('Removing fault variables: ')
                for fi in cfIdx:
                    sys.stdout.write(self.f[fi] + '  ')
                sys.stdout.write('\n')
            if len(czIdx) > 0:
                sys.stdout.write('Removing known variables: ')
                for zi in czIdx:
                    sys.stdout.write(self.z[zi] + '  ')
                sys.stdout.write('\n')

        self.X = self.X[eqs, :][:, xIdx]
        self.F = self.F[eqs, :][:, fIdx]
        self.Z = self.Z[eqs, :][:, zIdx]
        self.e = [self.e[ei] for ei in eqs]
        self.x = [self.x[xi] for xi in xIdx]
        self.f = [self.f[fi] for fi in fIdx]
        self.z = [self.z[zi] for zi in zIdx]

        if len(self.syme) > 0:
            self.syme = self.syme[eqs]

        self.P = np.arange(0, len(xIdx))
        self.Pfault = []

    def IsDynamic(self):
        """Return True if model is dynamic."""
        return np.any(self.X == 2)

    def IsStatic(self):
        """Return True if model is static."""
        return not self.IsDynamic()

    def Redundancy(self, eqs=[]):
        """Return the redundancy of the model."""
        if len(eqs) == 0:
            eqs = np.arange(0, self.ne())

        dm = dmperm.GetDMParts(self.X[eqs, :])
        return len(dm.Mp.row) - len(dm.Mp.col)

    def MTESRedundancy(self):
        """Return the redundancy of an MTES for the model."""
        eqs = np.argwhere(np.any(self.F, axis=1))[:, 0]
        return self.Redundancy(eqs) + 1

    def PlotDM(self, **options):
        """Plot Dulmage-Mendelsohn decomposition of model structure.

        Plots a Dulmage-Mendelsohn decomposition, originally described in
        Dulmage, A. and Mendelsohn, N. "Coverings of bipartite graphs."
        Canadian Journal of Mathematics 10.4 (1958): 516-534.

        Input
        -----
          eqclass : If True, perform canonical decomposition of M+ and
                    plot equivalence classes
          fault   : If true, indicates fault equations in canonical
                    decomposition of M+
          verbose : If True, use full variable names in plot. If not specified, heuristic
                    decides based on model size of variable names should be printed.
        """
        labelVars = False
        if 'verbose' in options:
            labelVars = options['verbose']
        elif self.X.shape[0] < 40:
            labelVars = True
        if 'eqclass' in options:
            eqclass = options['eqclass']
        else:
            eqclass = False
        if 'fault' in options:
            fault = options['fault']
        else:
            fault = False
        smplot.PlotDM(self, verbose=labelVars, eqclass=eqclass, fault=fault)

    def PlotModel(self, **options):
        """Plot the model structure."""
        labelVars = False
        if 'verbose' in options:
            labelVars = options['verbose']
        elif (self.nx() + self.nf() + self.nz()) < 40:
            labelVars = True

        smplot.PlotModel(self, verbose=labelVars)

    def PlotMatching(self, Gamma, **options):
        """Plot a matchingself.

        Plot the structure in an upper-triangular
        incidence matrix with the matched variables in the diagonal.

        Input
        -----
          Gamma : A matching computed by the Matching class method
        """
        q = np.array([], dtype=np.int64)
        for g in Gamma.matching:
            q = np.concatenate((q, g.col))

        # Determine if axis should be labeled
        labelVars = False
        if 'verbose' in options:
            labelVars = options['verbose']
        elif len(q) < 40:
            labelVars = True
        smplot.PlotMatching(self, Gamma, verbose=labelVars)

    def PossibleSensorLocations(self, x=-1):
        """Set possible sensor locations.

        Input
        -----
          x : Specification of possible sensor locations
              The sensor positions x can be given either as
              indices into model.x or variable names (default: all positions)
        """
        if x == -1:
            self.P = np.arange(0, len(self.x))  # Assume all possible sensor locations
        else:
            if issubclass(type(x[0]), str):  # list of strings
                self.P = np.array([self.x.index(xi) for xi in x if xi in self.x])
            else:
                self.P = x.copy()

    def SensorLocationsWithFaults(self, x=[]):
        """Set possible sensor locations that has faults in new sensors.

        Input
        -----
          x : Index to those sensor locations that can become faulty.
              If no input argument is given, no sensors may
              become faulty. The sensor positions x can be
              given either as indices into model.x or variable
              names
        """
        if len(x) > 0:
            if issubclass(type(x[0]), str):  # list of strings
                self.Pfault = np.array([self.x.index(xi) for xi in x if xi in self.x])
            else:
                self.Pfault = x.copy()
        else:
            self.Pfault = []

    def SensorPlacementIsolability(self, isolabilityspecification=-1):
        """Compute all minimal sensor sets that achieves maximal fault isolability of the faults in the model.

        Krysander, Mattias, and Erik Frisk, "Sensor placement for fault
        diagnosis." Systems, Man and Cybernetics, Part A: Systems and Humans,
        IEEE Transactions on 38.6 (2008): 1398-1410.

        Input
        -----
          isolabilityspecification : Isolability specification, a 0/1-matrix with
                                     a 0 in position (i,j) if fault fi should be
                                     isolable from fault fj; 1 otherwise. Structural
                                     isolability is a symmetric relation; and if
                                     the specification is not symmetric; the
                                     specification is made symmetric. Defaults to
                                     the identity matrix, i.e., full isolability
                                     among faults.

        Outputs
        -------
        res - list of all minimal sensor sets, represented by strings
        idx - same as res, but represented with indicices into model.f
        """
        if isolabilityspecification == -1:
            Ispec = np.eye(self.nf())

        return sensplace.SensorPlacementIsolability(self, Ispec)

    def SensorPlacementDetectability(self):
        """Compute all minimal sensor sets that achieves maximal fault detectability of the faults in the model.

        Krysander, Mattias, and Erik Frisk, "Sensor placement for fault
        diagnosis." Systems, Man and Cybernetics, Part A: Systems and Humans,
        IEEE Transactions on 38.6 (2008): 1398-1410.

        Outputs
        -------
        res - list of all minimal sensor sets, represented by strings
        idx - same as res, but represented with indicices into model.f
        """
        return sensplace.SensorPlacementDetectability(self)

    def AddSensors(self, sensors, name=[], fault=[]):
        """Add sensor equations to a model.

        Input
        -----
          sensors : Description of sensors to add s can be a list of strings
                    with the names of sensors to add or indices into the known
                    variables (model.x) which sensors to add.

                    It is only possible to add sensors measuring single
                    variables in x. If functions of variables in x are
                    measured, extend the model with a new variable first.

                    If the corresponding variable can be faulty, see
                    class method SensorLocationsWithFaults, a fault will be added
                    automatically.

          name    : List with names of new sensor variables
          fault   : List with names of fault variables for new sensors
        """
        if issubclass(type(sensors[0]), str):  # list of strings, convert to indices into self.x
            s = np.array([self.x.index(xi) for xi in sensors if xi in self.x])
        else:
            s = sensors

        ns = len(s)
        nx = self.X.shape[1]
        nz = self.Z.shape[1]
        nf = self.F.shape[1]
        ne = self.X.shape[0]
        Xs = np.zeros((ns, nx), np.int64)
        Fs = np.zeros((ns, nf), np.int64)
        Zs = np.zeros((ns, nz + ns), np.int64)

        fs = np.zeros(ns).astype(np.bool)
        for sIdx, si in enumerate(s):
            Xs[sIdx, si] = 1
            Zs[sIdx, sIdx + nz] = 1
            if si in self.Pfault:
                nf = nf + 1
                Fs = np.hstack((Fs, np.zeros((ns, 1), dtype=np.int64)))  # Add column for new fault
                Fs[sIdx, -1] = 1  # Add fault
                fs[sIdx] = True
            else:
                Fs[sIdx, :] = np.zeros((1, nf))

        self.X = np.vstack((self.X, Xs))
        self.Z = np.hstack((self.Z, np.zeros((ne, ns), dtype=np.int64)))
        self.Z = np.vstack((self.Z, Zs))
        if np.sum(fs) > 0:
            self.F = np.hstack((self.F, np.zeros((ne, np.sum(fs)), dtype=np.int64)))
        self.F = np.vstack((self.F, Fs))

        self.e = self.e + list(map(lambda x: self.vGen.NewE(), s))

        for idx, zi in enumerate(s):
            if len(name) == 0:
                znum = np.sum(np.array(s)[0:idx] == s[idx]) + 1
                if znum > 1:
                    zName = "z" + str(znum) + self.x[zi]
                else:
                    zName = "z" + self.x[zi]
            else:
                zName = name[idx]
            self.z = self.z + [zName]

            if fs[idx]:
                if len(fault) == 0:
                    fName = 'f' + zName
                else:
                    fName = fault[idx]
                self.f = self.f + [fName]

            if self.modelType is 'Symbolic':
                if fs[idx]:
                    rel = sym.Eq(sym.symbols(zName), sym.symbols(self.x[zi]) +
                                 sym.symbols(fName))
                else:
                    rel = sym.Eq(sym.symbols(zName), sym.symbols(self.x[zi]))
                self.syme = np.concatenate((self.syme, [rel]))

    def DetectabilityAnalysis(model, causality='mixed'):
        """Perform a structural detectability analysis.

        Inputs
        ------
          causality : Can be 'mixed' (default), 'int', or 'der' for mixed,
                      integral, or derivative causality analysis respectively.
                      For details, see

                      Frisk, E., Bregon, A., Aaslund, J., Krysander, M.,
                      Pulido, B., Biswas, G., "Diagnosability analysis
                      considering causal interpretations for differential
                      constraints", IEEE Transactions on Systems, Man and
                      Cybernetics, Part A: Systems and Humans, 2012, 42(5),
                      1216-1229.

        Outputs
        -------
          df  : Detectable faults
          ndf : Non-detectable faults
        """
        MplusCausal = lambda X: dmperm.Mplus(X, causality=causality)

        dm = MplusCausal(model.X)

        df = [model.f[fidx] for fidx in np.arange(0, model.F.shape[1])
              if np.argwhere(model.F[:, fidx])[:, 0] in dm]
        ndf = [model.f[fidx] for fidx in np.arange(0, model.F.shape[1])
               if np.argwhere(model.F[:, fidx])[:, 0] not in dm]
        return (df, ndf)

    def IsolabilityAnalysis(self, plot=False, permute=True, causality='mixed'):
        """Perform structural single fault isolability analysis of model.

        Inputs
        ------
          plot      : If True, plot the isolability matrix
          permute   : If True, permute the fault variables such that the
                      isolability matrix gets a block structure for easier
                      interpretation when plotted. Does not affect the output
                      argument im, only the plot (default True)

          causality : Can be 'mixed' (default), 'int', or 'der' for mixed,
                      integral, or derivative causality analysis respectively.
                      For details, see

                      Frisk, E., Bregon, A., Aaslund, J., Krysander, M.,
                      Pulido, B., Biswas, G., "Diagnosability analysis
                      considering causal interpretations for differential
                      constraints", IEEE Transactions on Systems, Man and
                      Cybernetics, Part A: Systems and Humans, 2012, 42(5),
                      1216-1229.
        Outputs
        -------
          im        : Isolability matrix, im(i,j)=1 if fault i can be isolated
                      from fault j, 0 otherwise
        """
        MplusCausal = lambda X: dmperm.Mplus(X, causality=causality)
#        ne = self.X.shape[0]
        nf = len(self.f)

        plusRow = MplusCausal(self.X)

        # Determine equations for each fault
#        feq = list(map(lambda fi: np.argwhere(self.F[:, fi])[0][0],
#                       np.arange(0, nf)))

        # Determine non-detectable faults
#        ndrows = [x for x in np.arange(0, ne) if x not in plusRow]
#        ndf = [self.f[fi] for fi in np.arange(0, len(self.f)) if feq[fi] in ndrows]
#        df = [self.f[fi] for fi in np.arange(0, len(self.f)) if not feq[fi] in ndrows]

        im = np.ones((nf, nf), dtype=np.int64)
        for fi in np.arange(0, nf):
            # Decouple fi
            fieqs = [x[0] for x in np.argwhere(self.F[:, fi] == 0)]
            X = self.X[fieqs, :]
            plusRow = MplusCausal(X)
            fisolrows = [fieqs[ei] for ei in plusRow]
            idx = [fj for fj in np.arange(0, nf) if np.any(self.F[fisolrows, :],
                                                           axis=0)[fj]]
            im[idx, fi] = 0

        if plot:
            if permute:
                p, q, _, _, _, _ = dmperm.dmperm(im)
            else:
                p = np.arange(0, nf)
                q = p
            plt.spy(im[p, :][:, q], markersize=8, marker="o", color="b")
            plt.xticks(np.arange(0, nf), [self.f[fi] for fi in p])
            plt.yticks(np.arange(0, nf), [self.f[fi] for fi in p])
            if len(self.name) > 0:
                titleString = "Isolability matrix for '%s'" % self.name
            else:
                titleString = "Isolability matrix"
            if causality is 'der':
                titleString = "%s (derivative causality)" % titleString
            elif causality is 'int':
                titleString = "%s (integral causality)" % titleString
            plt.title(titleString)
            plt.gca().xaxis.tick_bottom()
        return im

    def IsolabilityAnalysisArrs(self, arrs, permute=True, plot=False):
        """Perform structural single fault isolability analysis of a set of ARRs.

        Inputs
        ------
          arrs      : List of ARRs, e.g., list of MSO sets
          plot      : If True, plot the isolability matrix
          permute   : If True, permute the fault variables such that the
                      isolability matrix gets a block structure for easier
                      interpretation when plotted. Does not affect the output
                      argument im, only the plot (default True)

        Outputs
        -------
          im        : Isolability matrix, im(i,j)=1 if fault i can be isolated
                      from fault j, 0 otherwise
        """
        FSM = self.FSM(arrs)
        return self.IsolabilityAnalysisFSM(FSM, permute=permute, plot=plot)

    def IsolabilityAnalysisFSM(self, FSM, permute=True, plot=False):
        """Perform structural single fault isolability analysis of a fault sensitivity matrix (FSM).

        Inputs
        ------
          FSM       : Fault sensitivity matrix
          plot      : If True, plot the isolability matrix
          permute   : If True, permute the fault variables such that the
                      isolability matrix gets a block structure for easier
                      interpretation when plotted. Does not affect the output
                      argument im, only the plot (default True)

        Outputs
        -------
          im        : Isolability matrix, im(i,j)=1 if fault i can be isolated
                      from fault j, 0 otherwise
        """
        nf = FSM.shape[1]
#        nr = FSM.shape[0]
        im = np.ones((nf, nf), dtype=np.int64)
        for f in FSM:
            zIdx = np.array([[x0, y0] for x0 in np.where(f > 0)[0] for y0 in np.where(f == 0)[0]])
            if len(zIdx) > 0:
                im[zIdx[:, 0], zIdx[:, 1]] = 0

        if plot:
            if permute:
                p, q, _, _, _, _ = dmperm.dmperm(im)
            else:
                p = np.arange(0, nf)
                q = p
            plt.spy(im[p, :][:, q], markersize=10, marker='o')
            plt.xticks(np.arange(0, self.nf()), self.f)
            plt.yticks(np.arange(0, self.nf()), self.f)
            plt.gca().xaxis.tick_bottom()
            if len(self.name) > 0:
                plt.title("Isolability matrix for a given FSM in '" +
                          self.name + "'")
            else:
                plt.title('Isolability matrix for a given FSM')
        return im

    def SeqResGen(self, Gamma, resEq, name, diffres='int', language='Python',
                  batch=False, api='Python', user_functions={},
                  external_src=[], external_headers=[]):
        """(Experimental) Generate Python/C code for sequential residual generator.

        Given a matching and a residual equation, generate code implementing the
        residual generator. Generates Python/C file. How to call/compile  the
        generated file is described in the generated file (or the user manual).

        Inputs
        ------
          Gamma            : Matching
          resEq            : Index to equation to use as residual equation
          name             : Name of residual equation. Will also be used as a basis for
                             filename.
          diffres          : Can be 'int' or 'der' (default 'int'). Determines how
                             to treat differential constraints when used as a
                             residual equation.
          language         : Defaults to Matlab but also C code can be generated

          batch            : Generate a batch mode residual generator, only
                             applicable when generating C code. Instead of
                             computing the residual for one data samnple,
                             batch mode runs the residual generator for a whole
                             data set. This can significantly decrease
                             computational time.
          user_functions   : Translation dictionary, from user function name
                             to generated code
          external_src     : List of Python files with external functions
                             used in residual generator
          external_headers : List of external header files to include in
                             generated C source
        """
        codegen.SeqResGen(self, Gamma, resEq, name, diffres=diffres,
                          language=language, batch=batch,
                          user_functions=user_functions,
                          external_src=external_src,
                          external_headers=external_headers)

    def Lint(self):
        """Print model information and checks for inconsistencies in model definition."""
        war = False
        err = False

#        dm = dmperm.GetDMParts(self.X)

        if len(self.name) > 0:
            print('Model: ' + self.name)
        else:
            print("Model information")

        sys.stdout.write("\n  Type:" + self.modelType)

        nd = np.sum(self.X == 3)
        if nd > 0:
            sys.stdout.write(", dynamic\n")
        else:
            sys.stdout.write(", static\n")

        print('\n  Variables and equations')
        print('    ' + str(self.nx()) + ' unknown variables')
        print('    ' + str(self.nz()) + ' known variables')
        print('    ' + str(self.nf()) + ' fault variables')
        print('    ' + str(self.ne()) + ' equations, including ' + str(nd) + ' differential constraints')
        print('\n  Degree of redundancy: ' + str(self.Redundancy()))

        if self.Redundancy() > 0 and len(self.f) > 0:
            print('  Degree of redundancy of MTES set: ' +
                  str(self.MTESRedundancy()))
        print('\n')

        if self.ne() != self.F.shape[0] or self.ne() != self.Z.shape[0]:
            print('Error! Inconsistent numnber of rows in incidence matrices')
            err = err + 1

        if self.nx() != len(self.x):
            print('Error! Inconsistent number of unknown variables')
            err = err + 1

        if self.nz() != len(self.z):
            print('Error! Inconsistent number of known variables')
            err = err + 1

        if self.nf() != len(self.f):
            print('Error! Inconsistent number of fault variables')
            err = err + 1

        if self.ne() != len(self.e):
            print('Error! Inconsistent number of equations')
            err = err + 1

        if len([v for v in self.P if not (v in np.arange(0, self.nx()))]) > 0:
            print('Error! Possible sensor locations outside set of unknown variables')
            err = err + 1

        if len([v for v in self.Pfault if not (v in np.arange(0, self.nx()))]) > 0:
            print('Error! Possible sensor locations with faults outside set of unknown variables')
            err = err + 1

        if np.any(np.sum(self.F > 0, axis=0) > 1):
            print('Error! Fault variables can only appear in 1 equation, rewrite model with intermediate variables')
            err = err + 1

        xIdx = np.where(np.all(self.X == 0, axis=0))[0]
        for ii in xIdx:
            print('Warning! Variable ' + self.x[ii] + ' is not included in model')
            war = war + 1

        zIdx = np.where(np.all(self.Z == 0, axis=0))[0]
        for ii in zIdx:
            print('Warning! Variable ' + self.z[ii] + ' is not included in model')
            war = war + 1

        fIdx = np.where(np.all(self.F == 0, axis=0))[0]
        for ii in fIdx:
            print('Warning! Variable ' + self.f[ii] + ' is not included in model')
            war = war + 1
        if self.IsUnderdetermined():
            print('Warning! Model is underdetermined')
            war = war + 1

        print('  Model validation finished with %d errors and %d warnings.' % (err, war))


def DiffConstraint(dvar, ivar):
    """Define a differential constraint."""
    return [dvar, ivar, "diff"]


def _ModelStructure(rels, x):
    """Compute incidence matrix."""
    ne = len(rels)
    nx = len(x)

    X = np.zeros((ne, nx), dtype='int64')
    for k, rel in enumerate(rels):
        if IsDifferentialConstraint(rel):
            if (rel[0] in x) and (rel[1] in x):
                dvIdx = x.index(rel[0])
                ivIdx = x.index(rel[1])
                X[k, dvIdx] = 3
                X[k, ivIdx] = 2
        else:
            X[k, _RelModelStructure(rel, x)] = 1
    return X


def _RelModelStructure(rel, x):
    """Compute symbolic model structure."""
    if IsSymbolic(rel):
        relVars = [str(e) for e in rel.atoms(sym.Symbol)]
    else:
        relVars = rel
    return [xi for xi in range(0, len(x)) if x[xi] in relVars]


def IsDifferentialConstraint(rel):
    """Return true if relation is differential constraint."""
    return isinstance(rel, list) and len(rel) == 3 and rel[2] is "diff"


def IsSymbolic(v):
    """Return true if variable is symbolic."""
    return isinstance(v, tuple(sym.core.all_classes))


def _ToEquations(rels):
    """Convert relation into symbolic equations."""
    def _ToEquation(rel):
        if IsSymbolic(rel) and not isinstance(rel, sym.Equality):
            return sym.Eq(rel, 0)
        else:
            return rel
    return list(map(lambda r: _ToEquation(r), rels))
    