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
import time
from math import log10, ceil
import sympy as sp
import logging
import black

# %% Set up paths
if not('Engine' in sys.path):
    sys.path.append(os.getcwd() + '/Engine')

if not os.getcwd().endswith("Engine"):
    os.chdir('Engine')

from diagutil import loadmat
import GetMeasurementData as gm

sns.set(style='white', rc={'lines.linewidth': 0.75, 'axes.linewidth': 0.5})

# %% Model
from VEP4Engine import modeldef  # noqa: E402
from VEP_par import diag_par  # noqa: E402

model = ExtendedModel(modeldef, "VEP4")
model.Lint()

msos = model.MSO()

ts = [1649, 4011, 4016, 4017, 4066, 4074, 4477]
re = [73, 75, 75, 3, 76, 3, 2]
FSM = model.FSM([msos[ti] for ti in ts])

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

ds = "NF"
state = data[ds]['state_init']
Ts = data[ds]['Ts']
z = data[ds]['z'].T
N = len(data[ds]['time'])
t = np.arange(0, N * Ts, Ts)

# %% Refactored code generation
preamble = """from math import sqrt, exp, pi
from extmodelfuns import a_eta_t_fun, b_eta_t_fun, eta_c_fun, W_ic_fun
from extmodelfuns import W_af_fun, W_es_fun, PSI_th_fun, PI_wg_fun, PSIli_wg_fun, Tflow_wg_fun
from extmodelfuns import W_c_fun, PHI_model_fun, PI_c_fun, W_t_fun, Tq_t_fun
from extmodelfuns import PI_t_fun

max_fun = max
min_fun = min
"""

mso = msos[ts[0]]
red = mso[re[0]]

mso_model = model.GetMSOModel(mso, red, "ResGen_1649_1_dae")

# Have a look at h.o.d-model and the corresponding matching
mso_hod_model = mso_model.get_hod_model()
Gamma = mso_hod_model.get_hod_matching()
srank(mso_hod_model.X) == mso_hod_model.X.shape[0]
fig, ax = plt.subplots(num=11, clear=True)
mso_hod_model.PlotMatching(Gamma, verbose=True)

# Print solving sequence
f_sub = {sp.symbols(fi): 0 for fi in mso_hod_model.f}  # Substitute all fault variables
for gammai in Gamma.matching:
    print(str(mso_hod_model.syme[gammai.row][0].lhs.subs(f_sub)), end="")
    print(" => ", mso_hod_model.x[gammai.col[0]])

mso_model.DAEResGen(diag_par, preamble=preamble, use_black=False)

from ResGen_1649_1_dae import ResGen_1649_1_dae  # noqa: E402
res_gen = ResGen_1649_1_dae(diag_par)
res_gen.init(t, z)

z_0 = {z: zv for zv, z in zip(data["NF"]["z"][0], model.z)}
x1_0 = {x1i: data["NF"]["state_init"][x1i] for x1i in res_gen.x1}
w0, wp0 = res_gen.compute_initial_conditions(x1_0, z_0, diag_par, verbose=True)

res = np.zeros(len(w0))
res_gen._sim_residual(0, w0, wp0, res)
np.sqrt(res.dot(res))

# # %% Code generation
# preamble = """from math import sqrt, exp, pi
# from extmodelfuns import a_eta_t_fun, b_eta_t_fun, eta_c_fun, W_ic_fun
# from extmodelfuns import W_af_fun, W_es_fun, PSI_th_fun, PI_wg_fun, PSIli_wg_fun, Tflow_wg_fun
# from extmodelfuns import W_c_fun, PHI_model_fun, PI_c_fun, W_t_fun, Tq_t_fun
# from extmodelfuns import PI_t_fun

# max_fun = max
# min_fun = min
# """

# mso = msos[ts[0]]
# red = mso[re[0]]
# model.DAEResGen(mso, red, "ResGen_1649_1_dae", diag_par, preamble=preamble, use_black=False)

# from ResGen_1649_1_dae import ResGen_1649_1_dae  # noqa: E402
# res_gen = ResGen_1649_1_dae(diag_par)

# # %% Use a matching to solve for consistent initial conditions
# # MSO equations
# # w = (x1, x1', x2), h.o.d = (dx1/dt, x1', x2)
# # a) 0 = g(x1, x1', x2)
# # b) 0 = h(x1, x2)
# # c) dx1/dt = x1'

# from ResGen_1649_1_dae import ResGen_1649_1_dae  # noqa: E402

# ds = "NF"
# state = data[ds]['state_init']
# Ts = data[ds]['Ts']
# z = data[ds]['z'].T
# N = len(data[ds]['time'])
# t = np.arange(0, N * Ts, Ts)

# res_gen = ResGen_1649_1_dae(diag_par)
# res_gen.init(t, z)

# mso = msos[ts[0]]
# red = mso[re[0]]
# m0 = [e for e in mso if e != red]

# # Get highest-order-derivative form of exactly determined part of MSO
# mso_model = model.copy()
# mso_model.SubModel(mso, res_gen.vars)

# red_idx = int(np.argwhere(mso == red))  # Index of redundant equation
# mso_model.syme[red_idx] = sp.Eq(mso_model.syme[red_idx].lhs - sp.symbols("r"), 0)
# mso_model.x.append("r")
# mso_model.X = np.column_stack((mso_model.X, np.zeros(mso_model.X.shape[0])))
# mso_model.X[red_idx, -1] = 1

# # Identify x1, x2, and x1d variables
# diff_constraints_idx = [eidx for eidx, _ in enumerate(mso_model.e) if np.any(mso_model.X[eidx] > 1)]

# x1_idx = [np.argwhere(ei == 2)[0, 0] for ei in mso_model.X[diff_constraints_idx]]
# x1 = [mso_model.x[xi] for xi in x1_idx]

# x1d_idx = [np.argwhere(ei == 3)[0, 0] for ei in mso_model.X[diff_constraints_idx]]
# x1d = [mso_model.x[xi] for xi in x1d_idx]

# x2_idx = [x for x, _ in enumerate(mso_model.x) if x not in x1_idx + x1d_idx]
# x2 = [mso_model.x[xi] for xi in x2_idx]

# # Remove diff-constraints and variables x1
# mso_hod_model = mso_model.copy()
# mso_hod_model.SubModel(diff_constraints_idx, x1_idx, remove=True)

# # Verify full column rank
# srank(mso_hod_model.X) == mso_hod_model.X.shape[0]

# # Compute matching for hod
# Gamma = mso_hod_model.Matching(range(mso_hod_model.X.shape[0]))

# fig, ax = plt.subplots(num=11, clear=True)
# mso_hod_model.PlotMatching(Gamma, verbose=True)

# f_sub = {sp.symbols(fi): 0 for fi in mso_hod_model.f}  # Substitute all fault variables
# for gammai in Gamma.matching:
#     print(str(mso_hod_model.syme[gammai.row][0].lhs.subs(f_sub)), end="")
#     print(" => ", mso_hod_model.x[gammai.col[0]])

# # %% Generate code for each block to solve for
# file_name = "ic_blocks.py"
# use_black = True  # Use document formatter
# block_num_size = ceil(log10(len(Gamma.matching))) + 1
# with open(file_name, "wt") as f:
#     f.write('# File generated ' + time.strftime('%c') + '\n\n')
#     f.write("import numpy as np\n")
#     f.write(preamble)
#     f.write("\n\n")

#     for k, g_i in enumerate(Gamma.matching):
#         block_all_vars = set().union(*[ei.lhs.free_symbols for ei in mso_hod_model.syme[g_i.row]])
#         block_all_vars = block_all_vars - set(f_sub.keys())
#         any_x1 = True if [x1i for x1i in x1 if sp.symbols(x1i) in block_all_vars] else False
#         any_z = True if [zi for zi in mso_hod_model.z if sp.symbols(zi) in block_all_vars] else False
#         any_params = True if [pi for pi in mso_hod_model.parameters if sp.symbols(pi) in block_all_vars] else False

#         solve_vars = [mso_hod_model.x[xidx] for xidx in g_i.col]
#         pre_solve = [str(xi) for xi in block_all_vars if (str(xi) not in solve_vars and
#                                                           str(xi) not in mso_hod_model.parameters and
#                                                           str(xi) not in mso_hod_model.z and
#                                                           str(xi) not in x1)]

#         print(f"Block {k}: Solve for ", [mso_hod_model.x[xi] for xi in g_i.col])

#         f.write(f"def block_{k:0{block_num_size}}(x, subs):\n")
#         f.write("    # Solve for: ")
#         for xidx in g_i.col:
#             f.write(f"{mso_hod_model.x[xidx]}\n")
#         if any_x1 or any_z or any_params or pre_solve:
#             f.write("    ")
#             if any_z:
#                 f.write("z, ")
#             else:
#                 f.write("_, ")
#             if any_x1:
#                 f.write("x1, ")
#             else:
#                 f.write("_, ")
#             if any_params:
#                 f.write("params, ")
#             else:
#                 f.write("_, ")
#             if pre_solve:
#                 f.write("pre_solve, ")
#             else:
#                 f.write("_, ")
#             f.write(" = subs\n")
#         f.write("    ")
#         for xv in solve_vars:
#             f.write(f"{xv}, ")
#         f.write(" = x\n")
#         for zi in mso_hod_model.z:
#             if sp.symbols(zi) in block_all_vars:
#                 f.write(f"    {zi} = z[\"{zi}\"]\n")
#         for x1i in x1:
#             if sp.symbols(x1i) in block_all_vars:
#                 f.write(f"    {x1i} = x1[\"{x1i}\"]\n")
#         for pi in diag_par:
#             if sp.symbols(pi) in block_all_vars:
#                 f.write(f"    {pi} = params[\"{pi}\"]\n")
#         for xi in pre_solve:
#             if sp.symbols(xi) in block_all_vars:
#                 f.write(f"    {xi} = pre_solve[\"{xi}\"]\n")
#         f.write("\n")
#         f.write(f"    result = np.zeros({len(g_i.row)})\n")
#         for i, e_idx in enumerate(g_i.row):
#             ei = mso_hod_model.syme[e_idx].subs(f_sub)
#             if type(ei) is sp.Equality:
#                 code = f"result[{i}] = {str(ei.lhs)}"
#                 # code = f"result[{i}] = {black.format_str(str(ei.lhs), mode=black.FileMode(line_length=180))}"
#             f.write(f"    {code}\n")
#             # f.write(f"    {black.format_str(code, mode=black.FileMode(line_length=180))}\n")

#         f.write("    return result\n")
#         if g_i is not Gamma.matching[-1]:
#             f.write("\n\n")

# if use_black:
#     black.format_file_in_place(pathlib.Path(file_name), fast=True, write_back=black.WriteBack.YES, mode=black.FileMode())

# # %% Solve for IC
# import ic_blocks  # noqa: E402

# z_0 = {z: zv for zv, z in zip(data["NF"]["z"][0], model.z)}
# # x1_0 = data["NF"]["state_init"]
# x1_0 = {x1i: data["NF"]["state_init"][x1i] for x1i in x1}
# pre_solve = {}
# block_num_size = ceil(log10(len(Gamma.matching))) + 1
# for k in range(len(Gamma.matching)):
#     root_fun = getattr(ic_blocks, f"block_{k:0{block_num_size}}")
#     sol = root(lambda x: root_fun(x, (z_0, x1_0, diag_par, pre_solve)), 0)
#     if not sol.success:
#         logging.warning(f"Problems solving {[mso_hod_model.x[xi] for xi in Gamma.matching[k].col]} => {sol.x}: {sol.message}")
#     else:
#         print(f"Solving for: {[mso_hod_model.x[xi] for xi in Gamma.matching[k].col]} => {sol.x}")
#     pre_solve.update({mso_hod_model.x[xi_idx]: xi_sol for xi_sol, xi_idx in zip(sol.x, Gamma.matching[k].col)})

# pre_solve.update(x1_0)

# w0 = np.array([pre_solve[vi] for vi in res_gen.vars])
# wp0 = np.zeros(len(w0))
# wp0[x1_idx] = w0[x1d_idx]

# res = np.zeros(len(w0))
# res_gen._sim_residual(0, w0, wp0, res)
# np.sqrt(res.dot(res))

# # %%
# plt.show()
