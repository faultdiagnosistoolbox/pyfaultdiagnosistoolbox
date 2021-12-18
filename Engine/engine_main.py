# %load_ext autoreload
# %autoreload 2

from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
import sys
import pathlib
import numpy as np
from scipy.stats import gaussian_kde
import seaborn as sns
import platform

if not('Engine' in sys.path):
    sys.path.append(os.getcwd() + '/Engine')
os.chdir('Engine')

import batchfaultisolation as bfi
from diagutil import loadmat, RunResgenOnDataSets, Timer, PlotConfusionMatrix
import GetMeasurementData as gm

sns.set(style='white', rc={'lines.linewidth': 0.75, 'axes.linewidth': 0.5})

# # Model
import VEP4Engine
model = VEP4Engine.model
model.Lint()

# Plot the structural model
plt.figure(10, clear=True)
model.PlotModel(verbose=False)

# Plot isolability properties and the extended Dulmage-Mendelsohn
# decomposition with equivalence classes for the over-determined part
plt.figure(20, clear=True)
_ = model.IsolabilityAnalysis(plot=True)

plt.figure(21, clear=True)
_ = model.IsolabilityAnalysis(plot=True, causality='int')

plt.figure(22, clear=True)
_ = model.IsolabilityAnalysis(plot=True, causality='der')

plt.figure(23)
_ = model.PlotDM(fault=True, eqclass=True)

# MSO sets and test selection
print("Searching for MSO sets ...")
msos = model.MSO()
mtes = model.MTES()
print(f"  Found {len(msos)} MSO sets")
print(f"  Found {len(mtes)} MTES sets")

print(f"Checking index and observability properties of {len(msos)} MSO sets ...")
li = [model.IsLowIndex(m) for m in msos]
oi = [model.IsObservable(m) for m in msos]
print(f'  {sum(oi)} observable, {sum(li)} low (structural) differential index')


# Use tests selected from simple Kullback-Leibler selection strategy
ts = [1649, 4011, 4016, 4017, 4066, 4074, 4477]
re = [73, 75, 75, 3, 76, 3, 2]
for msoIdx, redIdx in zip(ts, re):
    mso = msos[msoIdx]
    red = mso[redIdx]
    m0 = [e for e in mso if e != red]
    Gamma = model.Matching(m0)
    print(f"MSO {msoIdx} with redundant equation {red}, causality: {Gamma.matchType}")
FSM = model.FSM([msos[ti] for ti in ts])

# Plot fault signature matrix and fault isolation matrix for selected set of tests
fIdx = [model.f.index(fi) for fi in ['fyw_af', 'fyp_im', 'fyp_ic', 'fyT_ic']]

fig, ax = plt.subplots(1, 2, num=30, clear=True)
ax[0].spy(FSM[:, fIdx], markersize=6, marker="o", color="b")
ax[0].set_xticks(range(len(fIdx)))
ax[0].set_xticklabels(np.array(model.f)[fIdx])
ax[0].set_yticks(range(len(ts)))
ax[0].set_yticklabels([f"MSO {k}" for k in ts])
ax[0].get_xaxis().set_ticks_position('bottom')
ax[0].set_title('Fault Signature Matrix')

plt.sca(ax[1])
IM = model.IsolabilityAnalysisArrs([msos[ti] for ti in ts])
ax[1].spy(IM[fIdx, :][:, fIdx], markersize=6, marker="o", color="b")
ax[1].set_xticks(range(len(fIdx)))
ax[1].set_xticklabels(np.array(model.f)[fIdx])
ax[1].set_yticks(range(len(fIdx)))
ax[1].set_yticklabels(np.array(model.f)[fIdx])
ax[1].get_xaxis().set_ticks_position('bottom')
ax[1].set_title('Fault isolation matrix')
fig.tight_layout()

# Code generation
for test, redIdx in zip(ts, re):
    mso = msos[test]
    red = mso[redIdx]
    m0 = [e for e in mso if e != red]
    resName = f"ResGen_{test}_{red}"
    Gamma = model.Matching(m0)
    model.SeqResGen(Gamma, red, resName, language='C', batch=True, user_functions=VEP4Engine.ext_funs,
                    external_headers=['extmodelfuns.h'], external_src=['extmodelfuns.cc'])
    print("")

# Compile sources
for test, redIdx in zip(ts, re):
    red = msos[test][redIdx]
    resName = f"ResGen_{test}_{red}"
    print(f"Compiling residual generator: {resName} ... ", end='')
    compile_cmd = f"python {resName}_setup.py build_ext --inplace"
    if platform.system == "Windows":
        compile_cmd = compile_cmd + " > nul"
    else:
        compile_cmd = compile_cmd + " > /dev/null"
    if os.system(compile_cmd) == 0:
        print('Success!')
    else:
        print('Failure!')

# Import generated residual generator modules
import ResGen_1649_1
import ResGen_4011_1
import ResGen_4016_1
import ResGen_4017_84
import ResGen_4066_1
import ResGen_4074_84
import ResGen_4477_86

res_gens = [ResGen_1649_1.ResGen_1649_1, ResGen_4011_1.ResGen_4011_1, ResGen_4016_1.ResGen_4016_1,
            ResGen_4017_84.ResGen_4017_84, ResGen_4066_1.ResGen_4066_1, ResGen_4074_84.ResGen_4074_84,
            ResGen_4477_86.ResGen_4477_86]

# # Import measurement data
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

_, ax = plt.subplots(num=90, clear=True)
ax.plot(ftp75_hwfet['t'] / 60.0, ftp75_hwfet['v'], linewidth=2)
sns.despine()
ax.set_xlabel('t [min]')
ax.set_ylabel('Velocity [km/h]')
_ = ax.set_title('EPA Highway Fuel Economy Test Cycle (HWFET)')

DS = 500  # Down sampling rate in plots
t_lim = [data['NF']['time'][0], data['NF']['time'][-1] / 60]

_, ax = plt.subplots(3, 3, num=80, clear=True)
yIdx = model.z.index('y_omega_e')
ax[0, 0].plot(data['NF']['time'][::DS] / 60.0, data['NF']['z'][::DS, yIdx] / (2 * np.pi) * 60)
ax[0, 0].set_ylabel('rpm')
ax[0, 0].set_title('Engine speed')
ax[0, 0].set_xlim(t_lim[0], t_lim[1])
ax[0, 0].get_yaxis().set_major_locator(MaxNLocator(4))
sns.despine()

yIdx = model.z.index('y_p_im')
ax[0, 1].plot(data['NF']['time'][::DS] / 60.0, data['NF']['z'][::DS, yIdx] / 1e3)
ax[0, 1].set_ylabel('kPa')
ax[0, 1].set_title('Intake manifold pressure')
ax[0, 1].set_xlim(t_lim[0], t_lim[1])
ax[0, 1].get_yaxis().set_major_locator(MaxNLocator(4))
sns.despine()

yIdx = model.z.index('y_W_af')
ax[0, 2].plot(data['NF']['time'][::DS] / 60.0, data['NF']['z'][::DS, yIdx])
ax[0, 2].set_ylabel('kg/s')
ax[0, 2].set_title('Air mass flow')
ax[0, 2].set_xlim(t_lim[0], t_lim[1])
sns.despine()
ax[0, 2].get_yaxis().set_major_locator(MaxNLocator(4))

yIdx = model.z.index('y_alpha_th')
ax[1, 0].plot(data['NF']['time'][::DS] / 60.0, data['NF']['z'][::DS, yIdx])
ax[1, 0].set_ylabel('%')
ax[1, 0].set_title('Throttle position')
ax[1, 0].set_xlim(t_lim[0], t_lim[1])
sns.despine()
ax[1, 0].get_yaxis().set_major_locator(MaxNLocator(4))

yIdx = model.z.index('y_wfc')
ax[1, 1].plot(data['NF']['time'][::DS] / 60.0, data['NF']['z'][::DS, yIdx] * 1.0e6)
ax[1, 1].set_ylabel('mg/s')
ax[1, 1].set_title('Injected fuel')
ax[1, 1].set_xlim(t_lim[0], t_lim[1])
sns.despine()
ax[1, 1].get_yaxis().set_major_locator(MaxNLocator(4))

yIdx = model.z.index('y_u_wg')
ax[1, 2].plot(data['NF']['time'][::DS] / 60.0, data['NF']['z'][::DS, yIdx])
ax[1, 2].set_ylabel('%')
ax[1, 2].set_title('Wastegate')
ax[1, 2].set_xlim(t_lim[0], t_lim[1])
sns.despine()
ax[1, 2].get_yaxis().set_major_locator(MaxNLocator(4))

yIdx = model.z.index('y_p_ic')
ax[2, 0].plot(data['NF']['time'][::DS] / 60.0, data['NF']['z'][::DS, yIdx] / 1e3)
ax[2, 0].set_xlabel('t [min]')
ax[2, 0].set_ylabel('kPa')
ax[2, 0].set_title('Intercooler pressure')
ax[2, 0].set_xlim(t_lim[0], t_lim[1])
sns.despine()
ax[2, 0].get_yaxis().set_major_locator(MaxNLocator(4))

yIdx = model.z.index('y_T_ic')
ax[2, 1].plot(data['NF']['time'][::DS] / 60.0, data['NF']['z'][::DS, yIdx])
ax[2, 1].set_xlabel('t [min]')
ax[2, 1].set_ylabel('K')
ax[2, 1].set_title('Intercooler temperature')
ax[2, 1].set_xlim(t_lim[0], t_lim[1])
sns.despine()
ax[2, 1].get_yaxis().set_major_locator(MaxNLocator(4))

yIdx = model.z.index('y_p_amb')
ax[2, 2].plot(data['NF']['time'][::DS] / 60.0, data['NF']['z'][::DS, yIdx] / 1.0e3)
ax[2, 2].set_xlabel('t [min]')
ax[2, 2].set_ylabel('kPa')
ax[2, 2].set_title('Ambient pressure')
ax[2, 2].set_xlim(t_lim[0], t_lim[1])
sns.despine()
ax[2, 2].get_yaxis().set_major_locator(MaxNLocator(4))
ax[2, 2].set_ylim(99, 101)
fig.tight_layout()
fig.subplots_adjust(top=0.9)
_ = fig.suptitle('Measurement data, no-fault dataset', fontsize=14, weight='bold')

# Run residual generators on measurement data
r = []
for k, ri in enumerate(res_gens):
    print(f"r{k + 1}: MSO {ri.__name__.split('_')[1]}")
    r.append(RunResgenOnDataSets(ri, data, VEP4Engine.diag_par))
    print('')

for ri in r:
    if 'NF' in ri:
        bias = np.mean(ri['NF']) * 0.0
        normFact = np.max(np.abs(ri['NF'] - bias)) * 1.1
        for FM in ri:
            ri[FM] = 1 / normFact * (ri[FM] - bias)

#  Simple threshold selection strategy
alpha = 1e-3
J = np.zeros(len(r))
for rIdx, ri in enumerate(r):
    N = len(ri['NF'])
    foo = np.sort(np.abs(ri['NF']))
    J[rIdx] = foo[np.ceil(N * (1 - alpha)).astype(np.int64)]

# Plot residuals
ds = 500
# dc = ['NF','fyp_im','fyw_af','fyp_ic','fyT_ic']
dc = ['NF', 'fyw_af']
for idx, fm in enumerate(dc):
    if fm in r[0]:  # Data set exists in first residual
        fig, ax = plt.subplots(3, 3, num=50 + idx, clear=True)
        for ridx, ri in enumerate(r):
            # ax: plt.Axes = fig.add_subplot(3, 3, ridx + 1)
            # ax_list.append(ax)
            ax_idx = (ridx // 3, ridx % 3)
            ax[ax_idx].plot(data[fm]['time'][::ds] / 60.0, ri[fm][::ds])
            ax[ax_idx].get_xaxis().set_major_locator(MaxNLocator(4))
            ax[ax_idx].get_yaxis().set_major_locator(MaxNLocator(4))
            if fm in model.f and FSM[ridx, model.f.index(fm)] == 1:
                for t1Idx, t2Idx in zip(data[fm]['fault_idx'][0:-1:2], data[fm]['fault_idx'][1::2]):
                    y1, y2 = ax[ax_idx].get_ylim()
                    t1 = data[fm]['time'][t1Idx] / 60.0
                    t2 = data[fm]['time'][t2Idx] / 60.0
                    ax[ax_idx].add_patch(
                        mpatches.Rectangle((t1, y1), t2 - t1, y2 - y1, facecolor='0.9', edgecolor='none'))
            ax[ax_idx].plot(data[fm]['time'][::ds] / 60.0, data[fm]['time'][::ds] * 0 + J[ridx], 'k--')
            ax[ax_idx].plot(data[fm]['time'][::ds] / 60.0, data[fm]['time'][::ds] * 0 - J[ridx], 'k--')

            if fm in model.f and FSM[ridx, model.f.index(fm)] == 1:
                ax[ax_idx].set_title(f'r{ridx + 1}: MSO {ts[ridx]} (*)', fontsize=10, weight='bold')
            else:
                ax[ax_idx].set_title(f'r{ridx + 1}: MSO {ts[ridx]}', fontsize=10, weight='bold')

            ax[ax_idx].set_xlim(0, np.max(data[fm]['time']) / 60)
            sns.despine()
        ax[2, 1].axis('off')
        ax[2, 2].axis('off')
        ax[1, 1].set_xlabel('t [min]')
        ax[1, 2].set_xlabel('t [min]')
        ax[2, 0].set_xlabel('t [min]')
        fig.tight_layout()
        fig.subplots_adjust(top=0.9)
        _ = fig.suptitle('Dataset: ' + fm, fontsize=12, weight='bold')

# Plot residual distributions
ds = 200
M = 50  # Number of data points in KDE plots
# dc = ['NF','fyp_im','fyw_af','fyp_ic','fyT_ic']
dc = ['NF', 'fyw_af']
for idx, fm in enumerate(dc):
    if fm in r[0]:  # Data set exists in first residual
        fig, ax = plt.subplots(3, 3, num=60 + idx, clear=True)
        for ridx, ri in enumerate(r):
            ax_idx = (ridx //3, ridx % 3)
            pNF = gaussian_kde(ri[fm][data[fm]['fault_vector'] == 0][::ds])

            rmin = np.min([np.min(ri['NF']), np.min(ri[fm]), -1])
            rmax = np.max([np.max(ri['NF']), np.max(ri[fm]), 1])
            r_range = np.arange(0, M) * (rmax - rmin) * 1.1 / (M - 1) + rmin - (rmax - rmin) * 0.05
            ax[ax_idx].plot(r_range, pNF(r_range), 'b', linewidth=2)
            if not fm == 'NF':
                pF = gaussian_kde(ri[fm][data[fm]['fault_vector'] == 1][::ds])
                ax[ax_idx].plot(r_range, pF(r_range), 'r', linewidth=2)
            if fm in model.f and FSM[ridx, model.f.index(fm)] == 1:
                ax[ax_idx].set_title(f'r{ridx+1}: MSO {ts[ridx]} (*)', fontsize=10, weight='bold')
            else:
                ax[ax_idx].set_title(f'r{ridx+1}: MSO {ts[ridx]}', fontsize=10, weight='bold')

            ax[ax_idx].get_xaxis().set_major_locator(MaxNLocator(3))
            ax[ax_idx].get_yaxis().set_major_locator(MaxNLocator(4))
            yl = ax[ax_idx].get_ylim()
            ax[ax_idx].plot([-1, -1], [0, 0.3 * yl[1]], 'k--')
            ax[ax_idx].plot([1, 1], [0, 0.3 * yl[1]], 'k--')
            sns.despine()
        ax[2, 1].axis('off')
        ax[2, 2].axis('off')
        fig.tight_layout()
        fig.subplots_adjust(top=0.9)
        fig.suptitle('Dataset: ' + fm, fontsize=12, weight='bold')


dc = ['NF', 'fyp_im', 'fyw_af', 'fyp_ic', 'fyT_ic']
# Compile with: python setup.py build_ext --inplace

t = Timer()
print(f'Running fault isolation strategy for {len(r)} residuals and {len(dc)} test cases ... ')
t.tic()
dx = bfi.BatchFaultIsolation(r, J, FSM[:, fIdx], dc, 1)
print(f'Finished in {t.toc():.2f} seconds')

# % C(i, j) = P(fi diagnos|fj injected fault)
nf = len(fIdx)
C = np.zeros((nf, nf))
for fiIdx, fi in enumerate(np.array(model.f)[fIdx]):
    for fjIdx, fj in enumerate(np.array(model.f)[fIdx]):
        tcIdx = dc.index(fj)
        C[fiIdx, fjIdx] = np.sum(dx[tcIdx][fiIdx, :][data[fj]['fault_vector'] > 0]) / np.sum(data[fj]['fault_vector'] > 0)

ds = 1000
for idx, fm in enumerate(dc):
    fig, ax = plt.subplots(2, 2, num=80 + idx, clear=True)

    for fiIdx in np.arange(0, nf):
        ax_idx = (fiIdx // 2, fiIdx % 2)
        ax[ax_idx].plot(data[fm]['time'][::ds] / 60, dx[idx][fiIdx, :][::ds], 'b')
        ax[ax_idx].get_yaxis().set_major_locator(MaxNLocator(4))
        ax[ax_idx].get_xaxis().set_major_locator(MaxNLocator(4))

        if fm in model.f and fm == model.f[fIdx[fiIdx]]:
            for t1Idx, t2Idx in zip(data[fm]['fault_idx'][0:-1:2], data[fm]['fault_idx'][1::2]):
                y1, y2 = ax[ax_idx].get_ylim()
                t1 = data[fm]['time'][t1Idx] / 60.0
                t2 = data[fm]['time'][t2Idx] / 60.0
                ax[ax_idx].add_patch(mpatches.Rectangle((t1, 0.05), t2 - t1, 0.9, facecolor='0.9', edgecolor='none'))
            perf = C[list(np.array(model.f)[fIdx]).index(fm), fiIdx]
        else:
            perf = np.sum(dx[idx][fiIdx, :] == 0) / len(data[fm]['fault_vector'])

        ax[ax_idx].set_title(f'{model.f[fIdx[fiIdx]]} (err={(1 - perf) * 100:.1f}%)', fontsize=10)
        sns.despine()
    ax[1, 0].set_xlabel('t [min]')
    ax[1, 1].set_xlabel('t [min]')
    fig.tight_layout()
    fig.subplots_adjust(top=0.9)
    fig.suptitle('Fault isolation: ' + fm, fontweight='bold', fontsize=12)

_, ax = plt.subplots(num=100, clear=True)
PlotConfusionMatrix(C)
ax.set_title('P(fi diag|fj injected)')
ax.set_xticks(np.arange(0, nf))
ax.set_xticklabels(np.array(model.f)[fIdx])
ax.set_yticks(np.arange(0, nf))
ax.set_yticklabels(np.array(model.f)[fIdx])
ax.set_xlabel('Injected fault')
ax.set_ylabel('Diagnosed fault')
_ = ax.set_title('Fault Isolation Performance Matrix')

plt.show()
