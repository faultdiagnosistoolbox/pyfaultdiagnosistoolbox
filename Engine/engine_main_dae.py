# %%
# %load_ext autoreload
# %autoreload 2

import matplotlib.pyplot as plt
import os
import sys
import pathlib
import numpy as np
import seaborn as sns
from extended_fdt import ExtendedModel
from faultdiagnosistoolbox import srank
from scipy.optimize import root

import sympy as sp

# %% Set up paths
if not('Engine' in sys.path):
    sys.path.append(os.getcwd() + '/Engine')

if not os.getcwd().endswith("Engine"):
    os.chdir('Engine')

from diagutil import loadmat
import GetMeasurementData as gm

sns.set(style='white', rc={'lines.linewidth': 0.75, 'axes.linewidth': 0.5})

# %% Model
from VEP4Engine import modeldef
from VEP_par import diag_par

model = ExtendedModel(modeldef, "VEP4")
model.Lint()

msos = model.MSO()

ts = [1649, 4011, 4016, 4017, 4066, 4074, 4477]
re = [73, 75, 75, 3, 76, 3, 2]
FSM = model.FSM([msos[ti] for ti in ts])

# %% Code generation
preamble = """
from math import sqrt, exp, pi
from extmodelfuns import a_eta_t_fun, b_eta_t_fun, eta_c_fun, W_ic_fun
from extmodelfuns import W_af_fun, W_es_fun, PSI_th_fun, PI_wg_fun, PSIli_wg_fun, Tflow_wg_fun
from extmodelfuns import W_c_fun, PHI_model_fun, PI_c_fun, W_t_fun, Tq_t_fun
from extmodelfuns import PI_t_fun

max_fun = max
min_fun = min
"""

for test, redIdx in zip(ts, re):
    mso = msos[test]
    red = mso[redIdx]
    resName = f"ResGen_{test}_{red}_dae"
    model.DAEResGen(mso, red, resName, preamble=preamble)

# %% Import generated residual generator modules
# from ResGen_1649_1_dae import ResGen_1649_1_dae
# from ResGen_4011_1_dae import ResGen_4011_1_dae
# from ResGen_4016_1_dae import ResGen_4016_1_dae
# from ResGen_4017_84_dae import ResGen_4017_84_dae
# from ResGen_4066_1_dae import ResGen_4066_1_dae
# from ResGen_4074_84_dae import ResGen_4074_84_dae
# from ResGen_4477_86_dae import ResGen_4477_86_dae

# res_gen_classes = [ResGen_1649_1_dae, ResGen_4011_1_dae, ResGen_4016_1_dae, ResGen_4017_84_dae, ResGen_4066_1_dae,
#                    ResGen_4074_84_dae, ResGen_4477_86_dae]

# res_gens = [ri(diag_par) for ri in res_gen_classes]

# %% Import measurement data
dataDir = f'{pathlib.Path.home()}/Research/Diagnosis/Data/EngineData/2016_EngineData/'
dataSets = {'NF': 'driving_cycle_FTP75_highway_no_fault_dataset1_16-01-20.mat',
            'fyp_im': 'driving_cycle_FTP75Highway_fault_y_pim_dataset_7_26-01-2016.mat',
            'fyw_af': 'driving_cycle_FTP75Highway_fault_y_waf_dataset_5_26-01-2016',
            'fyp_ic': 'driving_cycle_FTP75Highway_fault_y_pic_dataset_8_26-01-2016',
            'fyT_ic': 'driving_cycle_FTP75Highway_fault_y_Tic_dataset_6_26-01-2016'}
data = {}
print('Loading data ', end='')
for dd in dataSets:
    print('.', end='')
    data[dd] = gm.GetMeasurementData(dataDir + dataSets[dd])
print(f' Finished loading {len(data.keys())} datasets')
ftp75_hwfet = loadmat(dataDir + 'ftp75_hwfet.mat')
ftp75_hwfet = {'t': ftp75_hwfet['ftp75_hwfet']['t'], 'v': ftp75_hwfet['ftp75_hwfet']['v']}

# %% Run residual generators on measurement data
from ResGen_1649_1_dae import ResGen_1649_1_dae

ds = "NF"
state = data[ds]['state_init']
Ts = data[ds]['Ts']
z = data[ds]['z'].T
N = len(data[ds]['time'])
t = np.arange(0, N * Ts, Ts)

res_gen = ResGen_1649_1_dae(diag_par)
res_gen.init(t, z)
w, wp = res_gen.consistent_initial_conditions(model, data["NF"]["state_init"], solver="krylov", solver_options={"maxiter": 1000})

r = np.zeros(len(w))
res_gen._sim_residual(0, w, wp, r)
print(np.linalg.norm(r))

# w0 = np.zeros(len(res_gen.vars))
# for vi in data["NF"]["state_init"]:
#     w0[res_gen.vars.index(vi)] = data["NF"]["state_init"][vi]

res_gen.solve(t, w, wp)

# %% Use a matching to solve for consistent initial conditions
# MSO equations
# w = (x1, x1', x2), h.o.d = (dx1/dt, x1', x2)
# a) 0 = g(x1, x1', x2)
# b) 0 = h(x1, x2)
# c) dx1/dt = x1'

from ResGen_1649_1_dae import ResGen_1649_1_dae

res_gen = ResGen_1649_1_dae(diag_par)
res_gen.init(t, z)

mso = msos[ts[0]]
red = mso[re[0]]
m0 = [e for e in mso if e != red]

# Get highest-order-derivative form of exactly determined part of MSO
mso_hod_model = model.copy()
mso_hod_model.SubModel(m0, res_gen.vars)  # Get exactly determined part of MSO

# Identify x1, x2, and x1d variables
diff_constraints_idx = [eidx for eidx, _ in enumerate(mso_hod_model.e) if np.any(mso_hod_model.X[eidx] > 1)]
x1_idx = [xidx for xidx, _ in enumerate(mso_hod_model.x) if np.any(mso_hod_model.X[:, xidx] == 2)]
x1 = [mso_hod_model.x[xi] for xi in x1_idx]
x1d_idx = [xidx for xidx, _ in enumerate(mso_hod_model.x) if np.any(mso_hod_model.X[:, xidx] == 3)]
x1d = [mso_hod_model.x[xi] for xi in x1d_idx]
x2_idx = [x for x, _ in enumerate(mso_hod_model.x) if x not in x1_idx + x1d_idx]
x2 = [mso_hod_model.x[xi] for xi in x2_idx]

# Remove diff-constraints and variables x1
mso_hod_model.SubModel(diff_constraints_idx, x1_idx, remove=True)

# Verify full column rank
srank(mso_hod_model.X) == mso_hod_model.X.shape[0]

# Compute matching for hod
Gamma = mso_hod_model.Matching(range(mso_hod_model.X.shape[0]))

fig, ax = plt.subplots(num=10, clear=True)
mso_hod_model.PlotMatching(Gamma, verbose=True)


subs = {x1i: data["NF"]["state_init"][x1i] for x1i in x1}
subs.update({zi: yzi for zi, yzi in zip(model.z, data["NF"]["z"][0])})
subs.update({fi: 0 for fi in mso_hod_model.f})
subs.update(diag_par)

mso_hod_model.syme[Gamma.matching[1].row][0].subs(subs)

for gammai in Gamma.matching:
    print(str(mso_hod_model.syme[gammai.row][0].lhs.subs(subs)), end="")
    print(" => ", mso_hod_model.x[gammai.col[0]])
    # print(" => ", [mso_hod_model.x[xi] for xi in gammai.col])


def foo_res(theta):
    W_fc = theta
    return W_fc - 2.73945950655852e-5


root(foo_res, 0)

# %% Generate code for residual functions
# Lägg till lösta variabler eftersom i substitutionskedjan
with open("ic_blocks.py", "wt") as f:
    f.write("import numpy as np\n\n\n")
    for k, g_i in enumerate(Gamma.matching):
        block_all_vars = set().union(*[ei.lhs.free_symbols for ei in mso_hod_model.syme[g_i.row]])
        any_x1 = True if [x1i for x1i in x1 if sp.symbols(x1i) in block_all_vars] else False
        any_z = True if [zi for zi in mso_hod_model.z if sp.symbols(zi) in block_all_vars] else False
        any_params = True if [pi for zi in mso_hod_model.parameters if sp.symbols(zi) in block_all_vars] else False

        f.write(f"def block_{k}(x, subs):\n")
        f.write("    # Solve for: ")
        for xidx in g_i.col:
            f.write(f"{mso_hod_model.x[xidx]} \n")
        if any_x1 or any_z or any_params:
            f.write("    ")
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
            f.write(" = subs\n")
        f.write("    ")
        for xi in g_i.col:
            xv = mso_hod_model.x[xi]
            f.write(f"{xv}, ")
        f.write(" = x\n")
        for zi in mso_hod_model.z:
            if sp.symbols(zi) in block_all_vars:
                f.write(f"    {zi} = z[\"{zi}\"]\n")
        for x1i in x1:
            if sp.symbols(x1i) in block_all_vars:
                f.write(f"    {x1i} = x1[\"{x1i}\"]\n")
        for pi in diag_par:
            if sp.symbols(pi) in block_all_vars:
                f.write(f"    {pi} = params[\"{pi}\"]\n")
        f.write("\n")
        f.write(f"    result = np.zeros({len(g_i.row)})\n")
        for i, e_idx in enumerate(g_i.row):
            ei = mso_hod_model.syme[e_idx]
            if type(ei) is sp.Equality:
                code = f"result[{i}] = {str(ei.lhs)}"
            f.write(f"    {code}\n")

        f.write("    return result\n\n\n")
        # break
# %%
plt.show()

# %%
block_all_vars = set().union(*[ei.lhs.free_symbols for ei in mso_hod_model.syme[g_i.row]])
# %%
