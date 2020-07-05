# Residual generator design/simulation for Three-tank example from paper
#  Diagnosability Analysis Considering Causal Interpretations
#  for Differential Constraints, Erik Frisk, Anibal Bregon,
#  Jan Aaslund, Mattias Krysander, Belarmino Pulido, and Gautam Biswas
#  IEEE Transactions on Systems, Man and Cybernetics,
#  Part A: Systems and Humans (2012), Vol. 42, No. 5, 1216-1229.

from numpy.linalg import inv
import matplotlib.pyplot as plt
from controlpy.synthesis import controller_lqr
import faultdiagnosistoolbox as fdt
import seaborn as sns
import os
from TankSimulation import *

# %% Define model
model = fdt.models.three_tank
model_params = {'Rv1': 1, 'Rv2': 1, 'Rv3': 1, 'CT1': 1, 'CT2': 1, 'CT3': 1}

model.Lint()

# Plot model
plt.figure(10)
model.PlotModel()

# Isolability analysis
plt.figure(20)
model.IsolabilityAnalysis(plot=True)

plt.figure(21)
model.IsolabilityAnalysis(plot=True, causality='int')

plt.figure(22)
model.IsolabilityAnalysis(plot=True, causality='der')

plt.figure(23)
model.PlotDM(fault=True, eqclass=True)

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

# Simulate system and residual generators
G = ThreeTank_ss(model_params)
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
sim = [SimScenario(0, lambda t: 0*t, controller, model_params, t, x0),
       SimScenario(1, lambda t: 0.3*ramp(t, 6, 10), controller, model_params, t, x0),
       SimScenario(2, lambda t: 0.3*ramp(t, 6, 10), controller, model_params, t, x0),
       SimScenario(3, lambda t: 0.3*ramp(t, 6, 10), controller, model_params, t, x0),
       SimScenario(4, lambda t: 0.3*ramp(t, 6, 10), controller, model_params, t, x0),
       SimScenario(5, lambda t: 0.3*ramp(t, 6, 10), controller, model_params, t, x0),
       SimScenario(6, lambda t: 0.3*ramp(t, 6, 10), controller, model_params, t, x0)]

# Define measurements and add noise
nstd = np.array([0.01, 0, 0, 0, 0, 0.01, 0]) * noise
meas = [0, 5, 3]  # p1, q2, q0

M = np.zeros((sim[0].z0.shape[1], len(meas)))
M[meas, np.arange(0, len(meas))] = 1

for k in range(len(sim)):
    sim[k].z = sim[k].z0.dot(M) + np.random.normal(0, 1, (sim[k].z0.shape[0], len(meas))).dot(np.diag(nstd[meas]))


# Plot fault free and Rv1 scenarios
plt.figure(40, clear=True)
plt.subplot(3, 1, 1)
plt.plot(sim[0].t, sim[0].z0[:, 0], 'b', label='p1')
plt.plot(sim[0].t, ref(sim[0].t), 'b--', label='p1ref')
plt.plot(sim[0].t, sim[0].z0[:, 1], 'r', label='p2')
plt.plot(sim[0].t, sim[0].z0[:, 2], 'g', label='p3')
plt.ylabel('Tank water pressure')
plt.title('Fault free simulation (no noise)')
plt.legend()
sns.despine()

plt.subplot(3, 1, 2)
plt.plot(sim[0].t, sim[0].z0[:, 3])
plt.ylabel('q0')
sns.despine()

plt.subplot(3, 1, 3)
plt.plot(sim[0].t, sim[0].f)
plt.ylabel('Fault signal')
plt.xlabel('t [s]')
sns.despine()

plt.figure(41, clear=True)
plt.subplot(3, 1, 1)
plt.plot(sim[1].t, sim[1].z0[:, 0], 'b', label='p1')
plt.plot(sim[1].t, ref(sim[1].t), 'b--', label='p1ref')
plt.plot(sim[1].t, sim[1].z0[:, 1], 'r', label='p2')
plt.plot(sim[1].t, sim[1].z0[:, 2], 'g', label='p3')
plt.ylabel('Tank water pressure')
plt.legend()
plt.title('Fault scenario Rv1  (without noise)')
sns.despine()

plt.subplot(3, 1, 2)
plt.plot(sim[1].t, sim[1].z0[:, 3])
plt.ylabel('q0')
sns.despine()

plt.subplot(3, 1, 3)
plt.plot(sim[1].t, sim[1].f)
plt.ylabel('Fault signal')
plt.xlabel('t [s]')
sns.despine()


# Residual generator design
# Residual generator 1 in derivative causality (no loops), Fig. 2 in paper
eqr1 = mso2
redeq1 = 4
m01 = [e for e in eqr1 if e != redeq1]
Gamma1 = model.Matching(m01)
model.SeqResGen(Gamma1, redeq1, 'ResGen1', batch=True, language='C')
print("")

# Sequential residual generator in integral causality (no loops) (Fig. 3 in paper)
eqr2 = mso2
redeq2 = 6
m02 = [e for e in eqr2 if e != redeq2]
Gamma2 = model.Matching(m02)
model.SeqResGen(Gamma2, redeq2, 'ResGen2', batch=True, language='C')
print("")

# Sequential residual generator in mixed causality (Fig. 12 in paper)
eqr3 = mso1
redeq3 = 1
m03 = [e for e in eqr3 if e != redeq3]
Gamma3 = model.Matching(m03)
model.SeqResGen(Gamma3, redeq3, 'ResGen3', batch=True, language='C')

# Compile generated source code
for resName in ['ResGen1', 'ResGen2', 'ResGen3']:
    print(f"Compiling residual generator: {resName} ... ", end='')
    compile_cmd = f"python {resName}_setup.py build_ext --inplace > /dev/null"
    if os.system(compile_cmd) == 0:
        print('Success!')
    else:
        print('Failure!')

# Plot diagnosis system properties
plt.figure(50, clear=True)
ax = plt.subplot(1, 2, 1)
FSM = model.FSM([mso2, mso2, mso1], plot=True)
ax.set_yticklabels([f"r{k + 1}" for k in range(3)])
plt.title('Fault Signature Matrix')

plt.subplot(1, 2, 2)
model.IsolabilityAnalysisFSM(FSM, plot=True)
plt.title('Isolability matrix')

# %% Simulate residual generators for all test cases and plot results
import ResGen1
import ResGen2
import ResGen3
print("All residual generators imported/reloaded")

m = len(sim)
N = len(sim[0].t)
Ts = sim[0].t[1] - sim[0].t[0]

r1 = np.zeros((m, N))
r2 = np.zeros((m, N))
r3 = np.zeros((m, N))

for fi in range(m):
    state1 = {'p1': sim[fi].z[0, 0], 'p2': sim[fi].z[0, 1]}
    state2 = {'p1': sim[fi].z[0, 0], 'p2': sim[fi].z[0, 1]}
    state3 = {'p1': sim[fi].z[0, 0], 'p2': sim[fi].z[0, 1], 'p3': sim[fi].z[0, 2]}

    r1[fi, :] = ResGen1.ResGen1(sim[fi].z, state1, model_params, Ts)
    r2[fi, :] = ResGen2.ResGen2(sim[fi].z, state2, model_params, Ts)
    r3[fi, :] = ResGen3.ResGen3(sim[fi].z, state3, model_params, Ts)

    r1[fi, 0:2] = 0.0
    r3[fi, 0] = 0.0

# Plot some results
plt.figure(60)
plt.subplot(3, 1, 1)
for ri in range(m):
    plt.plot(sim[0].t, r1[ri, :])
sns.despine()
plt.title('r1 (seq/derivative/discrete)')

plt.subplot(3, 1, 2)
for ri in range(m):
    plt.plot(sim[0].t, r2[ri, :])
sns.despine()
plt.title('r2 (seq/integral/discrete)')

plt.subplot(3, 1, 3)
for ri in range(m):
    plt.plot(sim[0].t, r3[ri, :])
sns.despine()
plt.title('r3 (seq/mixed/discrete)')
plt.tight_layout()
