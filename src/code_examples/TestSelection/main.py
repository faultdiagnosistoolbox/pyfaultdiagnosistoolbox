# Residual Selection for Consistency Based Diagnosis Using Machine Learning Models

# Code corresponds to the paper "Residual Selection for Consistency Based Diagnosis Using Machine Learning Models"
# published at IFAC Safeprocess 2018 in Warszaw, Poland.
# 
# Note that the plots are not identical to the results in the paper where a Matlab implementation of the machine
# learning algorithms were used. However, the methodology is the same and the results are similar.

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from faultdiagnosistoolbox import RandomForestTestSelection, IsolabilityMatrix
from faultdiagnosistoolbox import DiagnosesAndConfusionMatrix, PlotConfusionMatrix

from example_util_functions import loadmat


# Load the data

# The data is loaded into a dictionary with 4 fields
# * modes - an array with names of no-fault and fault modes
# * res - An array with the 42 residuals
# * mode - a vector indicating which fault is active at each sample
# * fsm - A fault signature matrix based on model structure

data = loadmat('data.mat')['data']
nf = len(data['modes'])
nr = data['res'].shape[1]

# Pre process the data in two steps
# 1. Take absolute values of residuals
# 2. Threshold data (thdata)
# 
# The data is normalized so that a threshold at 1 corresponds to probability of false alarm of approximately 1%.

data['res'] = np.abs(data['res'])
thdata = data.copy()
thdata['res'] = thdata['res'] >= 1


# Plot the 7 first residuals in all fault modes. The residuals plotted in red are supposed to alarm for the fault
# according to the fault sensitivity matrix. (Fig. 4 in the paper)

fig, ax = plt.subplots(7, 8, num=10, clear=True, figsize=(10, 10))
for ri in range(7):
    for fm in range(nf):
        # plt.subplot(7, 8, ri*nf + fm + 1)
        if data['fsm'][ri, fm]==0:
            ax[ri, fm].plot(data['res'][data['mode'] == fm, ri], 'b', lw=0.3)
        else:
            ax[ri, fm].plot(data['res'][data['mode'] == fm, ri], 'r', lw=0.3)
        ax[ri, fm].tick_params(labelsize=6)
        ax[ri, fm].set_ylim(0, 3)
        sns.despine()
        if fm == 0:
            ax[ri, fm].set_ylabel(f'res-{ri + 1}', fontsize=8)
        if ri == 0:
            ax[ri, fm].set_title(data['modes'][fm], fontsize=8)
fig.tight_layout(w_pad=-0.75, h_pad=0)

# Basic analysis - performance of all 42 residuals
fsm = data['fsm']
fm = data['modes'][1:]

im = IsolabilityMatrix(fsm)
_, ax = plt.subplots(num=20, clear=True, figsize=(6, 6))
ax.spy(im[1:, 1:], marker='o', color='b')
ax.set_xticks(np.arange(len(fm)), fm)
ax.set_yticks(np.arange(len(fm)), fm)
ax.set_title('Isolability matrix')
ax.xaxis.tick_bottom()


_, C = DiagnosesAndConfusionMatrix(thdata)
fig, ax = plt.subplots(num=30, clear=True, figsize=(6, 6))
PlotConfusionMatrix(C, ax=ax)
ax.set_xticks(np.arange(nf), data['modes'])
ax.set_yticks(np.arange(nf), data['modes'])
ax.set_title('Fault Isolation Performance matrix, all 42 residuals')
fig.tight_layout()


# Test selection using Random Forest Classifiers

# First, build a random forest classifier based on the thresholded data.
# Here, 300 trees are trained in the tree ensemble.

res, C, rf, Crf = RandomForestTestSelection(thdata, n_estimators=300)

fig, ax = plt.subplots(num=31, clear=True, figsize=(6, 6))
PlotConfusionMatrix(Crf, ax=ax)
ax.set_xticks(np.arange(nf), data['modes'])
ax.set_yticks(np.arange(nf), data['modes'])
ax.set_title('Fault Isolation Performance matrix')
fig.tight_layout()


# Plot the variable importance, sorted, to get a ranking of predictor/residual usefullness in the classifier.
# Note that this classifier is not meant to be used in the diagnosis system. (Fig. 10 in the paper)

_, ax = plt.subplots(num=40, clear=True, figsize=(9, 6))
ax.plot(res['residualimportance'][res['sortidx']])
ax.set_yticks(fontsize=8)
ax.set_xticks(range(nr), res['sortidx'] + 1, fontsize=8, rotation=90)
ax.set_xlabel('Predictors')
ax.set_ylabel('Importance')
ax.set_title('Predictor importance')
sns.despine()


# Compute performance measures on false-alarm (FA), missed detection (MD), aggregated fault isolation (FI)
# and the probability of maximum isolability performance (FI-max) when selecting residuals according to the
# ranking computed above. Plot the three aggregated performance measures agains the number of selected residuals.
# (Fig. 11 in the paper)

num_res = [10, 12, 26, 27]
_, ax = plt.subplots(num=50, clear=True, figsize=(9, 7))
ax.plot(range(1, nr), res['pfa'], 'r', label='False alarm probability')
ax.plot(range(1, nr), res['pmd'], 'b', label='Missed detection probability')
ax.plot(range(1, nr), res['pfi'], 'y', label='False isolation probability')
for ni in num_res:
    ax.plot(ni+1, res['pfi'][ni], 'kx')

ax.legend()
ax.set_xlabel('Number of selected residuals')
ax.set_ylabel('Probability')
sns.despine()


# Plot the probability of maximum fault isolation performance for each fault. (Fig. 12 in the paper)

_, ax = plt.subplots(num=51, figsize=(10, 10))
for k in range(nf):
    ax.plot(res['pmfi'][:, k], label=thdata['modes'][k])
ax.legend(loc='upper right')
sns.despine()


# Compute and display confusion matrices corresponding to selecting 10, 12, 26, and 27 residuals.
# The results should be compared to the confusion matrix above where all 42 residuals were used.
# (Fig. 13 in the paper)

# Compare performance of 12 with 42 residuals.

ntests = 12
_, C = DiagnosesAndConfusionMatrix(thdata, residx=res['sortidx'][0:ntests])

fig, ax = plt.subplots(1, 2, num=90, clear=True, figsize=(12, 12))
PlotConfusionMatrix(C, ax=ax[0])
ax[0].set_title('No of tests: %d' % ntests)
ax[0].set_xticks(np.arange(nf))
ax[0].set_xticklabels(data['modes'])
ax[0].set_yticks(np.arange(nf))
ax[0].set_yticklabels(data['modes'])

_, C = DiagnosesAndConfusionMatrix(thdata)
PlotConfusionMatrix(C, ax[1])
ax[1].set_title('No of tests: 42')
ax[1].set_xticks(np.arange(nf))
ax[1].set_xticklabels(data['modes'])
ax[1].set_yticks(np.arange(nf))
ax[1].set_yticklabels(data['modes'])

fig.tight_layout()




