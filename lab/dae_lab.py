# %% Initial imports
from numpy.linalg import inv
import matplotlib.pyplot as plt
from controlpy.synthesis import controller_lqr
import faultdiagnosistoolbox as fdt
from seaborn import despine
import os
from faultdiagnosistoolbox.examples.ThreeTankSimulation import TankSimulation
import numpy as np
import sympy as sp
from scikits.odes import dae
from scipy.interpolate import interp1d
import faultdiagnosistoolbox.CodeGeneration as cg

# %% Utility function
def Jac(F, x):
    return sp.Matrix([[sp.diff(Fi, xj) for xj in x] for Fi in F])


# %% Define model
model = fdt.models.three_tank
model_params = {'Rv1': 1, 'Rv2': 1, 'Rv3': 1, 'CT1': 1, 'CT2': 1, 'CT3': 1}
model.Lint()

# Compute MSO amd MTES sets
msos = model.MSO()
mtes = model.MTES()
print(f"Found {len(msos)} MSO sets")
print(f"Found {len(mtes)} MTES sets")

# MSO causality sweep of two MSO sets
mso1 = msos[0]
mso2 = msos[1]
print(model.MSOCausalitySweep(mso1))
print(model.MSOCausalitySweep(mso2))

# %% Simulate system and residual generators
G = TankSimulation.ThreeTank_ss(model_params)
Lx, _, _ = controller_lqr(G['A'], G['B'], np.eye(3, 3), np.eye(1) * 0.5)
Lr = (1 / np.dot([1, 0, 0], -inv(G['A']-G['B'].dot(Lx)).dot(G['B'])))[0]


def ref(t):
    return 0.2 * np.sin(2 * np.pi * 1 / 10 * t) + 1


def controller(t, x):
    return np.transpose((-Lx.dot(np.transpose(x)) + Lr * ref(t)).clip(0))


# Simulate scenarios {NF, Rv1, Rv2, Rv3, CT1, CT2, CT3}
noise = 0
t = np.linspace(0, 20, 201)
x0 = [0, 0, 0]
sim = [TankSimulation.SimScenario(0, lambda t: 0*t, controller, model_params, t, x0),
       TankSimulation.SimScenario(1, lambda t: 0.3*TankSimulation.ramp(t, 6, 10), controller, model_params, t, x0),
       TankSimulation.SimScenario(2, lambda t: 0.3*TankSimulation.ramp(t, 6, 10), controller, model_params, t, x0),
       TankSimulation.SimScenario(3, lambda t: 0.3*TankSimulation.ramp(t, 6, 10), controller, model_params, t, x0),
       TankSimulation.SimScenario(4, lambda t: 0.3*TankSimulation.ramp(t, 6, 10), controller, model_params, t, x0),
       TankSimulation.SimScenario(5, lambda t: 0.3*TankSimulation.ramp(t, 6, 10), controller, model_params, t, x0),
       TankSimulation.SimScenario(6, lambda t: 0.3*TankSimulation.ramp(t, 6, 10), controller, model_params, t, x0)]

# Define measurements and add noise
nstd = np.array([0.01, 0, 0, 0, 0, 0.01, 0]) * noise
meas = [0, 5, 3]  # p1, q2, q0

M = np.zeros((sim[0].z0.shape[1], len(meas)))
M[meas, np.arange(0, len(meas))] = 1

for k in range(len(sim)):
    sim[k].z = sim[k].z0.dot(M) + np.random.normal(0, 1, (sim[k].z0.shape[0], len(meas))).dot(np.diag(nstd[meas]))


# %% Plot fault free and Rv1 scenarios
_, ax = plt.subplots(3, 1, num=40, clear=True)
ax[0].plot(sim[0].t, sim[0].z0[:, 0], 'b', label='p1')
ax[0].plot(sim[0].t, ref(sim[0].t), 'b--', label='p1ref')
ax[0].plot(sim[0].t, sim[0].z0[:, 1], 'r', label='p2')
ax[0].plot(sim[0].t, sim[0].z0[:, 2], 'g', label='p3')
ax[0].set_ylabel('Tank water pressure')
ax[0].set_title('Fault free simulation (no noise)')
ax[0].legend()
despine()

ax[1].plot(sim[0].t, sim[0].z0[:, 3])
ax[1].set_ylabel('q0')
despine()

ax[2].plot(sim[0].t, sim[0].f)
ax[2].set_ylabel('Fault signal')
ax[2].set_xlabel('t [s]')
despine()

_, ax = plt.subplots(3, 1, num=41, clear=True)
ax[0].plot(sim[1].t, sim[1].z0[:, 0], 'b', label='p1')
ax[0].plot(sim[1].t, ref(sim[1].t), 'b--', label='p1ref')
ax[0].plot(sim[1].t, sim[1].z0[:, 1], 'r', label='p2')
ax[0].plot(sim[1].t, sim[1].z0[:, 2], 'g', label='p3')
ax[0].set_ylabel('Tank water pressure')
ax[0].legend()
ax[0].set_title('Fault scenario Rv1  (without noise)')
despine()

ax[1].plot(sim[1].t, sim[1].z0[:, 3])
ax[1].set_ylabel('q0')
despine()

ax[2].plot(sim[1].t, sim[1].f)
ax[2].set_ylabel('Fault signal')
ax[2].set_xlabel('t [s]')
despine()


# %% Sequential residual generator in integral causality (no loops) (Fig. 3 in paper)
eqr2 = mso2
redeq2 = 6
m02 = [e for e in eqr2 if e != redeq2]
Gamma2 = model.Matching(m02)
model.SeqResGen(Gamma2, redeq2, 'ResGen2', batch=True, language='C')
print("")

# Compile generated source code
for resName in ['ResGen2']:
    print(f"Compiling residual generator: {resName} ... ", end='')
    compile_cmd = f"python {resName}_setup.py build_ext --inplace > /dev/null"
    if os.system(compile_cmd) == 0:
        print('Success!')
    else:
        print('Failure!')

# Simulate residual generators for all test cases and plot results
import ResGen2

m = len(sim)
N = len(sim[0].t)
Ts = sim[0].t[1] - sim[0].t[0]

r2 = np.zeros((m, N))

for fi in range(m):
    state2 = {'p1': sim[fi].z[0, 0], 'p2': sim[fi].z[0, 1]}
    r2[fi, :] = ResGen2.ResGen2(sim[fi].z, state2, model_params, Ts)

# Plot some results
_, ax = plt.subplots(num=10, clear=True)
for ri in range(m):
    ax.plot(sim[0].t, r2[ri, :])
despine()
ax.set_title('r2 (seq/integral/discrete)')

# %% Test that MSO-DAE is low-index
CT1, CT2, Rv1 = sp.symbols(["CT1", "CT2", "Rv1"])
p1, p2, dp1, dp2, q0, q1, q2, r, y1, y2, y3 = sp.symbols(["p1", "p2", "dp1", "dp2", "q0", "q1", "q2", "r", "y1", "y2", "y3"])
p1_dot, p2_dot = sp.symbols([r"\\dot{p}_1", r"\\dot{p}_2"])

w = sp.Matrix([p1, p2, dp1, dp2, q0, q1, q2, r])
whod = sp.Matrix([p1_dot, p2_dot, dp1, dp2, q0, q1, q2, r])

F = sp.Matrix([
    q2 - y2,
    p1 - y1 - r,
    -dp2 + (q1 - q2) / CT2,
    -p2_dot + dp2,
    -dp1 + (q0 - q1) / CT1,
    q0 - y3,
    -p1_dot + dp1,
    q1 - (p1 - p2) / Rv1])

Fhod = Jac(F, whod)
Fhod.det()

# %% Setup simulation code
sim_params = [model_params[pi] for pi in model.parameters]
def Resgen2_dae(_, w, wp, result, y, params):
    p1, p2, dp1, dp2, q0, q1, q2, r  = w
    p1_dot, p2_dot = wp[0:2]
    y1, y2, y3 = y
    Rv1, _, _, CT1, CT2, _ = params
    result[0] = q2 - y2
    result[1] = p1 - y1 - r
    result[2] = -dp2 + (q1 - q2) / CT2
    result[3] = -p2_dot + dp2
    result[4] = -dp1 + (q0 - q1) / CT1
    result[5] = q0 - y3
    result[6] = -p1_dot + dp1
    result[7] = q1 - (p1 - p2) / Rv1


# NF simulation
t = sim[0].t
y = sim[0].z
y_inter = interp1d(t, y.T, fill_value="extrapolate")
r2_sim = dae("ida", lambda t, w, wp, res: Resgen2_dae(t, w, wp, res, y_inter(t), sim_params),
             old_api=False, algebraic_vars_idx=np.arange(2, 8), compute_initcond='yp0',
             rtol=1e-3, atol=1e-3)
r2_init = np.hstack((y[0, 0:2], np.zeros(6)))
r2_sol = r2_sim.solve(t, r2_init, 0 * r2_init)

_, ax = plt.subplots(num=11, clear=True)
ax.plot(r2_sol.values.t, r2_sol.values.y[:, 7], label='DAE')
ax.plot(sim[0].t, r2[0, :], label='Sequential')
ax.legend()
despine()

# %% Simulate all modes
r2_dae = np.zeros((len(sim), len(t)))
for fi in range(m):
    y_inter = interp1d(t, sim[fi].z.T, fill_value="extrapolate", kind="linear")
    r2_sim = dae("ida", lambda t, w, wp, res: Resgen2_dae(t, w, wp, res, y_inter(t), sim_params),
                 old_api=False, algebraic_vars_idx=np.arange(2, 8), compute_initcond='yp0',
                 rtol=1e-3, atol=1e-3)
    r2_init = np.hstack((sim[fi].z[0, 0:2], np.zeros(6)))
    r2_sol = r2_sim.solve(t, r2_init, 0 * r2_init)
    r2_dae[fi, :] = r2_sol.values.y[:, 7]

_, ax = plt.subplots(1, 2, num=20, clear=True)
for ri in range(m):
    ax[0].plot(sim[0].t, r2[ri, :])
ax[0].set_xlabel("t [s]")
ax[0].set_title('r2 (seq/integral/discrete)')
despine()

for ri in range(m):
    ax[1].plot(sim[0].t, r2_dae[ri, :])
ax[1].set_xlabel("t [s]")
ax[1].set_title('r2 (dae)')
despine()


# %% Residual generator class
class ResGen2DAE:
    def __init__(self, params, t, y):
        self._y_inter = interp1d(t, y, fill_value="extrapolate")
        self._params = params
        alg_idx = np.arange(2, 8)
        self._integrator = dae("ida", lambda t, w, wp, res: self._sim_residual(t, w, wp, res),
                               old_api=False, algebraic_vars_idx=alg_idx, compute_initcond='yp0',
                               rtol=1e-3, atol=1e-3)

    def _sim_residual(self, t, w, wp, result):
        p1, p2, dp1, dp2, q0, q1, q2, r = w
        p1_dot, p2_dot = wp[0:2]
        y1, y2, y3 = self._y_inter(t)
        Rv1, _, _, CT1, CT2, _ = self._params
        result[0] = q2 - y2
        result[1] = p1 - y1 - r
        result[2] = -dp2 + (q1 - q2) / CT2
        result[3] = -p2_dot + dp2
        result[4] = -dp1 + (q0 - q1) / CT1
        result[5] = q0 - y3
        result[6] = -p1_dot + dp1
        result[7] = q1 - (p1 - p2) / Rv1

    def solve(self, t, w0):
        r_sol = self._integrator.solve(t, w0, 0 * w0)
        return r_sol.values.y[:, 7]

resgen_2 = ResGen2DAE(sim_params, t, y.T)
r2_class = resgen_2.solve(t, r2_init)

# %% Generate code
exprs = model.syme[mso2]
for i, e_idx in enumerate(mso2):
    ei = model.syme[e_idx]
    res_append = " - r" if e_idx==redeq2 else ""
    if type(ei) is sp.Equality:
        code = f"result[{i}] = {str(ei.lhs)}{res_append}"
    elif type(ei) is list and len(ei) == 3 and ei[2] == 'diff':
        code = f"result[{i}] = {ei[1]}_dot - {ei[0]}{res_append}"
    print(code)
