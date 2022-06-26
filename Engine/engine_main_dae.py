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

w0 = np.zeros(len(res_gen.vars))
for vi in data["NF"]["state_init"]:
    w0[res_gen.vars.index(vi)] = data["NF"]["state_init"][vi]

res_gen.solve(t, w0)

# %%
plt.show()

# %%
