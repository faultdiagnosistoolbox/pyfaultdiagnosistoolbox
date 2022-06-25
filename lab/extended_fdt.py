import numpy as np
import sympy as sp
import time
import faultdiagnosistoolbox as fdt


class ExtendedModel(fdt.DiagnosisModel):
    def __init__(self, modeldef, name=''):
        super().__init__(modeldef, name="")

    def DAEResGen(self, mso, residual_equation, name, language="python"):
        if self.IsHighIndex(mso):
            raise RuntimeError("Only low-index MSOs are permitted")

        if self.IsHighIndex([ei for ei in mso if ei != residual_equation]):
            raise RuntimeError("""MSO is low-index but choice of residual equation gives high-index residual generator.
Choose more visely""")

        if language == "python":
            self.generate_python_dae_resgen(mso, residual_equation, name)
        else:
            raise NotImplementedError(f"Generation for language {language} not yetr implemented")

    def generate_python_dae_resgen(self, mso, residual_equation, name):
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
            f.write("from scikits.odes import dae\n")

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
            f.write("        self._atol = atol\n")
            f.write("        self._rtol = rtol\n")
            f.write("\n")

            f.write('    def init(self, t, y, kind="linear"):\n')
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
            f.write("    def solve(self, t, w0):\n")
            f.write("        r_sol = self._integrator.solve(t, w0, 0 * w0)\n")
            f.write("        self.flag = r_sol.flag\n")
            f.write("        self.message = r_sol.message\n")
            f.write("        if self.flag != 0:\n")
            f.write('            warn("Simulation ended with warning", RuntimeWarning)\n')
            f.write(f"        self.res = r_sol.values.y[:, {mso_vars.index('r')}]\n")

        print(f"File {file_name} with class {name} generated.")
