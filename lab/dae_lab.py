# %%
# %matplotlib
# %load_ext autoreload
# %autoreload 2
# %% Initial imports
from warnings import warn
from numpy.linalg import inv
import matplotlib.pyplot as plt
from controlpy.synthesis import controller_lqr
from seaborn import despine
import os
from faultdiagnosistoolbox.examples.ThreeTankSimulation import TankSimulation
import numpy as np
import sympy as sp
from scipy.optimize import fsolve, root
from extended_fdt import ExtendedModel
from faultdiagnosistoolbox.models._three_tank_model import modelDef as tt_model_def


# %% Utility function
def Jac(F, x):
    return sp.Matrix([[sp.diff(Fi, xj) for xj in x] for Fi in F])


# %% Define model
model = ExtendedModel(tt_model_def)
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

mso = mso2
redeq = mso[1]

# %% Verify that MSO-DAE is low-index
CT1, CT2, Rv1 = sp.symbols(["CT1", "CT2", "Rv1"])
p1, p2, dp1, dp2, q0, q1, q2, r, y1, y2, y3 = sp.symbols(["p1", "p2", "dp1", "dp2", "q0", "q1", "q2",
                                                          "r", "y1", "y2", "y3"])
p1_dot, p2_dot = sp.symbols([r"\\dot{p}_1", r"\\dot{p}_2"])
f = sp.symbols(['fV1', 'fV2', 'fV3', 'fT1', 'fT2', 'fT3'])

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

# %% Simulate system and residual generators
G = TankSimulation.ThreeTank_ss(model_params)
Lx, _, _ = controller_lqr(G.A, G.B, np.eye(3, 3), np.eye(1) * 0.5)
Lr = (1 / np.dot([1, 0, 0], -inv(G.A - G.B.dot(Lx)).dot(G.B)))[0]


def ref(t):
    return 0.2 * np.sin(2 * np.pi * 1 / 10 * t) + 1


def controller(t, x):
    return np.transpose((-Lx.dot(np.transpose(x)) + Lr * ref(t)).clip(0))


# Simulate scenarios {NF, Rv1, Rv2, Rv3, CT1, CT2, CT3}
noise = 0
t = np.linspace(0, 20, 201)
x0 = [0, 0, 0]
sim = [TankSimulation.SimScenario(0, lambda t: 0 * t, controller, model_params, t, x0),
       TankSimulation.SimScenario(1, lambda t: 0.3 * TankSimulation.ramp(t, 6, 10), controller, model_params, t, x0),
       TankSimulation.SimScenario(2, lambda t: 0.3 * TankSimulation.ramp(t, 6, 10), controller, model_params, t, x0),
       TankSimulation.SimScenario(3, lambda t: 0.3 * TankSimulation.ramp(t, 6, 10), controller, model_params, t, x0),
       TankSimulation.SimScenario(4, lambda t: 0.3 * TankSimulation.ramp(t, 6, 10), controller, model_params, t, x0),
       TankSimulation.SimScenario(5, lambda t: 0.3 * TankSimulation.ramp(t, 6, 10), controller, model_params, t, x0),
       TankSimulation.SimScenario(6, lambda t: 0.3 * TankSimulation.ramp(t, 6, 10), controller, model_params, t, x0)]

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


# %% Generate code for DAE and sequential residual generator in integral causality (no loops) (Fig. 3 in paper)
m02 = [e for e in mso if e != redeq]

# DAE residual generator
model.DAEResGen(mso, redeq, "ResGen2_dae")

# Sequential residual generator
Gamma2 = model.Matching(m02)
model.SeqResGen(Gamma2, redeq, 'ResGen2', batch=True, language='C')

# Compile generated source code
for resName in ['ResGen2']:
    print(f"Compiling residual generator: {resName} ... ", end='')
    compile_cmd = f"python {resName}_setup.py build_ext --inplace > /dev/null"
    if os.system(compile_cmd) == 0:
        print('Success!')
    else:
        print('Failure!')

# %% Simulate residual generators for all test cases
from ResGen2_dae import ResGen2_dae
from ResGen2 import ResGen2

m = len(sim)
N = len(sim[0].t)
Ts = sim[0].t[1] - sim[0].t[0]

r2_dae = ResGen2_dae(model_params)

r2_seq_sol = np.zeros((m, N))
r2_dae_sol = np.zeros((len(sim), len(t)))
for fi in range(m):
    state2 = {'p1': sim[fi].z[0, 0], 'p2': sim[fi].z[0, 1]}
    r2_seq_sol[fi, :] = ResGen2(sim[fi].z, state2, model_params, Ts)

    y = sim[fi].z.T
    r2_dae.init(t, y)
    r2_init = np.hstack((sim[fi].z[0, 0:2], np.zeros(6)))
    r2_dae.solve(t, r2_init)
    r2_dae_sol[fi] = r2_dae.res

# %% Plot some results
_, ax = plt.subplots(1, 2, num=20, clear=True)
for ri in range(m):
    ax[0].plot(sim[0].t, r2_seq_sol[ri])
ax[0].set_xlabel("t [s]")
ax[0].set_title('r2 (seq/integral/discrete)')
despine()

for ri in range(m):
    ax[1].plot(sim[0].t, r2_dae_sol[ri])
ax[1].set_xlabel("t [s]")
ax[1].set_title('r2 (dae)')
despine()


# %%
def consistent_initial_conditions(res_gen, model: ExtendedModel,
                                  initial_state_dict: dict,
                                  solver="hybr",
                                  solver_options={}) -> tuple[np.ndarray, np.ndarray]:
    # dx1 = f(x1, x2)
    # 0 &= g(x1, x2)
    # d/dt x1 = dx1

    # w = res_gen.vars = (x1, dx1, x2)

    # Find index into w for x1
    n = len(res_gen.vars)
    x1_idx = np.array([res_gen.vars.index(model.x[i]) for
                       i in np.argwhere(np.any(model.X[mso] == 2, axis=0)).reshape(-1)])
    n1 = len(x1_idx)

    # Index for remaining variables (dx1, x2)
    nx1_idx = np.array([i for i in range(n) if i not in x1_idx])
    n2 = n - n1

    # Initialize random
    theta = np.random.default_rng().uniform(0, 1, n)

    # Get initial values for x1(0)
    x10 = [initial_state_dict[res_gen.vars[i]] for i in x1_idx]
    w = np.zeros(n)
    w[x1_idx] = x10
    wp = np.zeros(n)

    # Residual equations, solve for (d/dt x1(0), dx1(0), x2(0))
    def initial_res(theta):
        w[nx1_idx] = theta[0:n2]
        wp[x1_idx] = theta[n2:n]
        r = np.zeros(n)
        res_gen._sim_residual(0, w, wp, r)
        return r

    sol = root(initial_res, theta, method=solver, options=solver_options)

    if not sol.success:
        warn(f"Consistent initial condition computations exited in error: {sol.message}", RuntimeWarning)

    w[nx1_idx] = sol.x[0:n2]
    wp[x1_idx] = sol.x[n2:n]

    return w, wp


t = sim[0].t
y = sim[0].z.T
init_dict = {'p1': y[0, 0], 'p2': y[1, 0]}

r2_dae = ResGen2_dae(model_params)
r2_dae.init(t, y)
w, wp = r2_dae.consistent_initial_conditions(model, init_dict)
r = np.zeros(8)
r2_dae._sim_residual(0, w, wp, r)
print(r)

r2_dae.solve(t, w, wp)



# %%
