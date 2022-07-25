import numpy as np
import sympy as sp
import time
from importlib.util import find_spec
import pathlib
import faultdiagnosistoolbox as fdt
from math import ceil, log10

if find_spec("black"):
    import black  # noqa: F811
    black_available = True
else:
    black_available = False


class ExtendedModel(fdt.DiagnosisModel):
    def __init__(self, modeldef, name=""):
        super().__init__(modeldef, name=name)

    def GetMSOModel(self, mso, residual_equation, name):
        if self.IsHighIndex(mso):
            raise RuntimeError("Only low-index MSOs are currently permitted")

        if self.IsHighIndex([ei for ei in mso if ei != residual_equation]):
            raise RuntimeError("""MSO is low-index but choice of residual equation gives high-index residual generator.
Choose more visely""")

        return MSOResGenModel(self, mso, residual_equation, name)


class MSOResGenModel(fdt.DiagnosisModel):
    def __init__(self, model, mso, residual_equation, name):
        mso_all_vars = set().union(*[ei.lhs.free_symbols for ei in model.syme[mso] if type(ei) is not list])
        mso_all_vars = [vi.name for vi in list(mso_all_vars)]

        mso_params = [pi for pi in model.parameters if pi in mso_all_vars]
        mso_vars = [xi for xi in model.x if xi in mso_all_vars] + ["r"]
        mso_known_vars = [zi for zi in model.z if zi in mso_all_vars]
        mso_diff_vars = [ei[1] for ei in model.syme[mso] if type(ei) is list and ei[2] == "diff"]
        mso_diff_vars = [vi for vi in model.x if vi in mso_diff_vars]

        init_model = model.copy()
        init_model.SubModel(mso, mso_all_vars)

        # Add redundancy variable r to the model
        red_idx = int(np.argwhere(mso == residual_equation))  # Index of redundant equation
        init_model.syme[red_idx] = sp.Eq(init_model.syme[red_idx].lhs - sp.symbols("r"), 0)
        init_model.x.append("r")
        init_model.X = np.column_stack((init_model.X, np.zeros(init_model.X.shape[0])))
        init_model.X[red_idx, -1] = 1

        self.x = init_model.x
        self.X = init_model.X
        self.f = init_model.f
        self.F = init_model.F
        self.z = init_model.z
        self.Z = init_model.Z
        self.e = init_model.e
        self.syme = init_model.syme
        self.modelType = init_model.modelType
        self.parameters = init_model.parameters
        self.name = name
        self.P = np.arange(0, len(self.x))
        self.Pfault = []

        self.z_orig = model.z  # Needed for correct indexing of original data
        self.mso_all_vars = mso_all_vars
        self.mso_params = mso_params
        self.mso_vars = mso_vars
        self.mso_known_vars = mso_known_vars
        self.mso_diff_vars = mso_diff_vars

        # Identify x1, x2, and x1d variables
        diff_constraints_idx = [eidx for eidx, _ in enumerate(self.e) if np.any(self.X[eidx] > 1)]

        self.x1_idx = [np.argwhere(ei == 2)[0, 0] for ei in self.X[diff_constraints_idx]]
        self.x1 = [self.x[xi] for xi in self.x1_idx]

        self.x1d_idx = [np.argwhere(ei == 3)[0, 0] for ei in self.X[diff_constraints_idx]]
        self.x1d = [self.x[xi] for xi in self.x1d_idx]

        self.x2_idx = [x for x, _ in enumerate(self.x) if x not in self.x1_idx + self.x1d_idx]
        self.x2 = [self.x[xi] for xi in self.x2_idx]

        self.hod_model = self.copy()
        self.hod_model.SubModel(diff_constraints_idx, self.x1_idx, remove=True)
        self.Gamma = self.hod_model.Matching(range(self.hod_model.X.shape[0]))
        self.hod_model.Gamma = self.Gamma

    def DAEResGen(self, diag_par, language="python", preamble=None, use_black=None):
        if use_black is None:
            use_black = black_available
        elif use_black and not black_available:
            raise RuntimeError("Install package black to activate automated source code formatting")

        self._diag_par = diag_par

        if language == "python":
            file_name = f"{self.name}.py"
            with open(file_name, "wt") as f:
                self.write_preamble(f, preamble)
                self.generate_python_dae_resgen(f)
                self.generate_utility_functions(f)
                self.generate_ic_code(f)

            if use_black:
                black.format_file_in_place(pathlib.Path(file_name), fast=True, write_back=black.WriteBack.YES, mode=black.FileMode())
            print(f"File {file_name} with class {self.name} generated.")
        else:
            raise NotImplementedError(f"Generation for language {language} not yetr implemented")

    def write_preamble(self, f, preamble):
        f.write('# File generated ' + time.strftime('%c') + '\n')

        f.write("from warnings import warn\n")
        f.write("import logging\n")
        f.write("import numpy as np\n")
        f.write("from scipy.interpolate import interp1d\n")
        f.write("from math import ceil, log10\n")
        f.write("from scipy.optimize import root\n")
        f.write("from scikits.odes import dae\n")

        if preamble:
            f.write(preamble)

        f.write("\n\n")

    def generate_python_dae_resgen(self, f):
        f.write(f"class {self.name}:\n")
        f.write("    def __init__(self, params, rtol=1e-3, atol=1e-3):\n")
        f.write("        self.vars = (")
        for vi in self.mso_vars:
            if vi != "r":
                f.write(f'"{vi}", ')
            else:
                f.write(f'"{vi}")\n')

        f.write("        self._params = (")
        for pi in self.mso_params:
            f.write(f'params["{pi}"]')
            if pi == self.mso_params[-1]:
                f.write(")\n")
            else:
                f.write(", ")
        self._alg_idx = np.sort([self.mso_vars.index(vi) for vi in self.mso_vars if vi not in self.mso_diff_vars])
        f.write("        self._alg_idx = np.array([")
        for v_idx in self._alg_idx:
            if v_idx != self._alg_idx[-1]:
                f.write(f"{v_idx}, ")
            else:
                f.write(f"{v_idx}])\n")
        f.write("\n")
        # Find index into w for x1
        n = len(self.mso_vars)
        f.write("        self._x1_idx = np.array([")
        for v_idx in self.x1_idx:
            if v_idx != self.x1_idx[-1]:
                f.write(f"{v_idx}, ")
            else:
                f.write(f"{v_idx}])\n")
        f.write("        self._x1d_idx = np.array([")
        for v_idx in self.x1d_idx:
            if v_idx != self.x1d_idx[-1]:
                f.write(f"{v_idx}, ")
            else:
                f.write(f"{v_idx}])\n")
        # Index for remaining variables (dx1, x2)
        nx1_idx = np.array([i for i in range(n) if i not in self.x1_idx])
        f.write("        self._nx1_idx = np.array([")
        for v_idx in nx1_idx:
            if v_idx != nx1_idx[-1]:
                f.write(f"{v_idx}, ")
            else:
                f.write(f"{v_idx}])\n")
        f.write("\n")
        f.write("        self.x1 = [")
        for v_idx in self.x1_idx:
            if v_idx != self.x1_idx[-1]:
                f.write(f"\"{self.mso_vars[v_idx]}\", ")
            else:
                f.write(f"\"{self.mso_vars[v_idx]}\"]\n")

        f.write("        self._Gamma_var = [")
        for gi in self.Gamma.matching:
            f.write("[")
            for xi in gi.col:
                f.write(f"\"{self.hod_model.x[xi]}\"")
                if xi != gi.col[-1]:
                    f.write(", ")
                else:
                    f.write("]")
            if gi != self.Gamma.matching[-1]:
                f.write(", ")
            else:
                f.write("]\n")

        f.write("        self._atol = atol\n")
        f.write("        self._rtol = rtol\n")
        f.write("\n")

        f.write('    def init(self, t, y, kind="linear"):\n')
        f.write('        self._sol = None\n')
        f.write('        self.res = None\n')
        f.write('        self._y_inter = interp1d(t, y, fill_value="extrapolate", kind=kind)\n')

        f.write('        self._integrator = dae("ida", self._sim_residual, old_api=False,\n')
        f.write('                               algebraic_vars_idx=self._alg_idx, compute_initcond="yp0",\n')
        f.write('                               rtol=self._rtol, atol=self._atol)\n\n')

        f.write("    def _sim_residual(self, t, w, wp, result):\n")
        f.write("        (")
        for pi in self.mso_vars:
            if pi != "r":
                f.write(f"{pi}, ")
            else:
                f.write("r) = w\n")

        for vi in self.mso_diff_vars:
            f.write(f"        {vi}_dot = wp[{self.mso_vars.index(vi)}]\n")

        f.write("\n")
        f.write("        z = self._y_inter(t)\n")
        for zi in self.mso_known_vars:
            f.write(f"        {zi} = z[{self.z_orig.index(zi)}]\n")

        f.write("\n")
        f.write("        (")
        for pi in self.mso_params:
            f.write(pi)
            if pi == self.mso_params[-1]:
                f.write(")")
            else:
                f.write(", ")
        f.write(" = self._params\n")

        f.write("\n")
        no_fault_subs = {sp.symbols(fi): 0 for fi in self.f}
        for i, ei in enumerate(self.syme):
            if type(ei) is sp.Equality:
                code = f"result[{i}] = {str(ei.lhs.subs(no_fault_subs))}"
            elif type(ei) is list and len(ei) == 3 and ei[2] == 'diff':
                code = f"result[{i}] = {ei[1]}_dot - {ei[0]}"
            f.write(f"        {code}\n")

    def generate_utility_functions(self, f):
        f.write("\n")
        f.write("    def state(self) -> list[str]:\n")
        f.write("        return self.x1\n")

        f.write("\n")
        f.write("    def solve(self, t, w0, wp0=None):\n")
        f.write("        if wp0 is not None:\n")
        f.write("            wp0 = 0 * w0\n")
        f.write("        self._sol = self._integrator.solve(t, w0, 0 * w0)\n")
        f.write("        self.flag = self._sol.flag\n")
        f.write("        self.message = self._sol.message\n")
        f.write("        if self.flag != 0:\n")
        f.write('            warn("Simulation ended with warning", RuntimeWarning)\n')
        f.write(f"        self.res = self._sol.values.y[:, {self.mso_vars.index('r')}]\n")

        f.write("\n")
        f.write("    def consistent_initial_conditions(self,\n")
        f.write("                                      initial_state_dict: dict,\n")
        f.write('                                      solver="hybr",\n')
        f.write("                                      solver_options={}) -> tuple[np.ndarray, np.ndarray]:\n")
        f.write("        # dx1 = f(x1, x2)\n")
        f.write("        # 0 &= g(x1, x2)\n")
        f.write("        # d/dt x1 = dx1\n")
        f.write("\n")
        f.write("        # w = res_gen.vars = (x1, dx1, x2)\n")
        f.write("        n = len(self.vars)\n")
        f.write("        n1 = len(self._x1_idx)\n")
        f.write("        n2 = n - n1\n")
        f.write("\n")
        f.write("        # Initialize random\n")
        f.write("        theta = np.random.default_rng().uniform(0, 1, n)\n")
        f.write("\n")
        f.write("        # Get initial values for x1(0)\n")
        f.write("        x10 = [initial_state_dict[self.vars[i]] for i in self._x1_idx]\n")
        f.write("        w = np.zeros(n)\n")
        f.write("        w[self._x1_idx] = x10\n")
        f.write("        wp = np.zeros(n)\n")
        f.write("\n")
        f.write("        # Residual equations, solve for (d/dt x1(0), dx1(0), x2(0))\n")
        f.write("        def initial_res(theta):\n")
        f.write("            w[self._nx1_idx] = theta[0:n2]\n")
        f.write("            wp[self._x1_idx] = theta[n2:n]\n")
        f.write("            r = np.zeros(n)\n")
        f.write("            self._sim_residual(0, w, wp, r)\n")
        f.write("            return r\n")
        f.write("\n")
        f.write("        sol = root(initial_res, theta, method=solver, options=solver_options)\n")
        f.write("\n")
        f.write("        if not sol.success:\n")
        f.write('            warn(f"Consistent initial condition computations exited in error: {sol.message}", RuntimeWarning)\n')
        f.write("\n")
        f.write("        w[self._nx1_idx] = sol.x[0:n2]\n")
        f.write("        wp[self._x1_idx] = sol.x[n2:n]\n")
        f.write("\n")
        f.write("        return w, wp\n")

    def get_hod_model(self):
        return self.hod_model

    def get_hod_matching(self):
        return self.Gamma

    def generate_ic_code(self, f):
        f_sub = {sp.symbols(fi): 0 for fi in self.hod_model.f}  # Substitute all fault variables
        f.write("\n")
        f.write("    def compute_initial_conditions(self, x1_0, z_0, diag_par, verbose=False):\n")
        block_num_size = ceil(log10(len(self.Gamma.matching))) + 1

        for k, g_i in enumerate(self.Gamma.matching):
            block_all_vars = set().union(*[ei.lhs.free_symbols for ei in self.hod_model.syme[g_i.row]])
            block_all_vars = block_all_vars - set(f_sub.keys())
            any_x1 = True if [x1i for x1i in self.x1 if sp.symbols(x1i) in block_all_vars] else False
            any_z = True if [zi for zi in self.hod_model.z if sp.symbols(zi) in block_all_vars] else False
            any_params = True if [pi for pi in self.hod_model.parameters if sp.symbols(pi) in block_all_vars] else False

            solve_vars = [self.hod_model.x[xidx] for xidx in g_i.col]
            pre_solve = [str(xi) for xi in block_all_vars if (str(xi) not in solve_vars and
                                                              str(xi) not in self.hod_model.parameters and
                                                              str(xi) not in self.hod_model.z and
                                                              str(xi) not in self.x1)]

            ind = " " * 8

            f.write(ind + f"def block_{k:0{block_num_size}}(x, subs):\n")
            f.write(ind + "    # Solve for: ")
            for xidx in g_i.col:
                f.write(f"{self.hod_model.x[xidx]}\n")
            if any_x1 or any_z or any_params or pre_solve:
                f.write(ind + "    ")
                if any_z:
                    f.write("z, ")
                else:
                    f.write("_, ")
                if any_x1:
                    f.write("x1, ")
                else:
                    f.write("_, ")
                if any_params:
                    f.write("params, ")
                else:
                    f.write("_, ")
                if pre_solve:
                    f.write("pre_solve, ")
                else:
                    f.write("_, ")
                f.write(" = subs\n")
            f.write(ind + "    ")
            for xv in solve_vars:
                f.write(f"{xv}, ")
            f.write(" = x\n")
            for zi in self.hod_model.z:
                if sp.symbols(zi) in block_all_vars:
                    f.write(ind + f"    {zi} = z[\"{zi}\"]\n")
            for x1i in self.x1:
                if sp.symbols(x1i) in block_all_vars:
                    f.write(ind + f"    {x1i} = x1[\"{x1i}\"]\n")
            for pi in self._diag_par:
                if sp.symbols(pi) in block_all_vars:
                    f.write(ind + f"    {pi} = params[\"{pi}\"]\n")
            for xi in pre_solve:
                if sp.symbols(xi) in block_all_vars:
                    f.write(ind + f"    {xi} = pre_solve[\"{xi}\"]\n")
            f.write("\n")
            f.write(ind + f"    result = np.zeros({len(g_i.row)})\n")
            for i, e_idx in enumerate(g_i.row):
                ei = self.hod_model.syme[e_idx].subs(f_sub)
                if type(ei) is sp.Equality:
                    code = f"result[{i}] = {str(ei.lhs)}"
                f.write(ind + f"    {code}\n")

            f.write(ind + "    return result\n")
            if g_i is not self.Gamma.matching[-1]:
                f.write("\n")

        f.write("\n")
        f.write(ind + "block_num_size = ceil(log10(len(self._Gamma_var))) + 1\n")
        f.write(ind + "pre_solve = {}\n")
        f.write("\n")
        f.write(ind + "for k, gi_var in enumerate(self._Gamma_var):\n")
        f.write(ind + '    root_fun = locals()[f"block_{k:0{block_num_size}}"]\n')
        f.write(ind + "    sol = root(lambda x: root_fun(x, (z_0, x1_0, diag_par, pre_solve)), 0)\n")
        f.write(ind + "    if not sol.success:\n")
        f.write(ind + '        logging.warning(f"Problems solving {gi_var} => {sol.x}: {sol.message}")\n')
        f.write(ind + "    elif verbose:\n")
        f.write(ind + '        print(f"Solving for: {gi_var} => {sol.x}")\n')
        f.write(ind + "    pre_solve.update({vi: xi_sol for xi_sol, vi in zip(sol.x, gi_var)})\n")
        f.write("\n")
        f.write(ind + "pre_solve.update(x1_0)\n")
        f.write("\n")
        f.write(ind + "w0 = np.array([pre_solve[vi] for vi in self.vars])\n")
        f.write(ind + "wp0 = np.zeros(len(w0))\n")
        f.write(ind + "wp0[self._x1_idx] = w0[self._x1d_idx]\n")
        f.write("\n")
        f.write(ind + "return w0, wp0\n")
