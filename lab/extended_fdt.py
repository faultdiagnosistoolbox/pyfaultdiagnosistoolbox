import numpy as np
import sympy as sp
import time
import faultdiagnosistoolbox as fdt


class ExtendedModel(fdt.DiagnosisModel):
    def __init__(self, modeldef, name=''):
        super().__init__(modeldef, name="")

    def DAEResGen(self, mso, residual_equation, name, language="python", preamble=None):
        if self.IsHighIndex(mso):
            raise RuntimeError("Only low-index MSOs are permitted")

        if self.IsHighIndex([ei for ei in mso if ei != residual_equation]):
            raise RuntimeError("""MSO is low-index but choice of residual equation gives high-index residual generator.
Choose more visely""")

        if language == "python":
            self.generate_python_dae_resgen(mso, residual_equation, name, preamble)
        else:
            raise NotImplementedError(f"Generation for language {language} not yetr implemented")

    def generate_python_dae_resgen(self, mso, residual_equation, name, preamble):
        file_name = f"{name}.py"

        with open(file_name, "wt") as f:
            mso_all_vars = set().union(*[ei.lhs.free_symbols for ei in self.syme[mso] if type(ei) is not list])
            mso_all_vars = [vi.name for vi in list(mso_all_vars)]

            mso_params = [pi for pi in self.parameters if pi in mso_all_vars]
            mso_vars = [xi for xi in self.x if xi in mso_all_vars] + ["r"]
            mso_known_vars = [zi for zi in self.z if zi in mso_all_vars]
            mso_diff_vars = [ei[1] for ei in self.syme[mso] if type(ei) is list and ei[2] == "diff"]
            mso_diff_vars = [vi for vi in self.x if vi in mso_diff_vars]

            f.write('# File generated ' + time.strftime('%c') + '\n')

            f.write("from warnings import warn\n")
            f.write("import numpy as np\n")
            f.write("from scipy.interpolate import interp1d\n")
            f.write("from scipy.optimize import root\n")
            f.write("from scikits.odes import dae\n")

            if preamble:
                f.write(preamble)

            f.write("\n\n")
            f.write(f"class {name}:\n")

            f.write("    def __init__(self, params, rtol=1e-3, atol=1e-3):\n")

            f.write("        self.vars = (")
            for vi in mso_vars:
                if vi != "r":
                    f.write(f'"{vi}", ')
                else:
                    f.write(f'"{vi}")\n')

            f.write("        self._params = (")
            for pi in mso_params:
                f.write(f'params["{pi}"]')
                if pi == mso_params[-1]:
                    f.write(")\n")
                else:
                    f.write(", ")
            self._alg_idx = np.sort([mso_vars.index(vi) for vi in mso_vars if vi not in mso_diff_vars])
            f.write("        self._alg_idx = np.array([")
            for v_idx in self._alg_idx:
                if v_idx != self._alg_idx[-1]:
                    f.write(f"{v_idx}, ")
                else:
                    f.write(f"{v_idx}])\n")
            f.write("\n")
            # Find index into w for x1
            n = len(mso_vars)
            x1_idx = np.array([mso_vars.index(self.x[i]) for
                               i in np.argwhere(np.any(self.X[mso] == 2, axis=0)).reshape(-1)])
            f.write("        self._x1_idx = np.array([")
            for v_idx in x1_idx:
                if v_idx != x1_idx[-1]:
                    f.write(f"{v_idx}, ")
                else:
                    f.write(f"{v_idx}])\n")
            # Index for remaining variables (dx1, x2)
            nx1_idx = np.array([i for i in range(n) if i not in x1_idx])
            f.write("        self._nx1_idx = np.array([")
            for v_idx in nx1_idx:
                if v_idx != nx1_idx[-1]:
                    f.write(f"{v_idx}, ")
                else:
                    f.write(f"{v_idx}])\n")
            f.write("\n")
            f.write("        self.x1 = [")
            for v_idx in x1_idx:
                if v_idx != x1_idx[-1]:
                    f.write(f"\"{mso_vars[v_idx]}\", ")
                else:
                    f.write(f"\"{mso_vars[v_idx]}\"]\n")

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
            for pi in mso_vars:
                if pi != "r":
                    f.write(f"{pi}, ")
                else:
                    f.write("r) = w\n")

            for vi in mso_diff_vars:
                f.write(f"        {vi}_dot = wp[{mso_vars.index(vi)}]\n")

            f.write("\n")
            f.write("        z = self._y_inter(t)\n")
            for zi in mso_known_vars:
                f.write(f"        {zi} = z[{self.z.index(zi)}]\n")

            f.write("\n")
            f.write("        (")
            for pi in mso_params:
                f.write(pi)
                if pi == mso_params[-1]:
                    f.write(")")
                else:
                    f.write(", ")
            f.write(" = self._params\n")

            f.write("\n")
            no_fault_subs = {sp.symbols(fi): 0 for fi in self.f}
            for i, e_idx in enumerate(mso):
                ei = self.syme[e_idx]
                res_append = " - r" if e_idx == residual_equation else ""
                if type(ei) is sp.Equality:
                    code = f"result[{i}] = {str(ei.lhs.subs(no_fault_subs))}{res_append}"
                elif type(ei) is list and len(ei) == 3 and ei[2] == 'diff':
                    code = f"result[{i}] = {ei[1]}_dot - {ei[0]}{res_append}"
                f.write(f"        {code}\n")

            f.write("\n")
            f.write("    def solve(self, t, w0, wp0=None):\n")
            f.write("        if wp0 is not None:\n")
            f.write("            wp0 = 0 * w0\n")
            f.write("        self._sol = self._integrator.solve(t, w0, 0 * w0)\n")
            f.write("        self.flag = self._sol.flag\n")
            f.write("        self.message = self._sol.message\n")
            f.write("        if self.flag != 0:\n")
            f.write('            warn("Simulation ended with warning", RuntimeWarning)\n')
            f.write(f"        self.res = self._sol.values.y[:, {mso_vars.index('r')}]\n")

            f.write("\n")
            f.write("    def consistent_initial_conditions(self, model,\n")
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

        print(f"File {file_name} with class {name} generated.")
