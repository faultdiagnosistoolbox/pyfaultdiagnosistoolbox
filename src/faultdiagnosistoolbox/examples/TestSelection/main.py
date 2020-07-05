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

plt.figure(10, clear=True, figsize=(10, 10))
for ri in range(7):
    for fm in range(nf):
        plt.subplot(7, 8, ri*nf + fm + 1)
        if data['fsm'][ri, fm]==0:
            plt.plot(data['res'][data['mode']==fm, ri], 'b', lw=0.3)
        else:
            plt.plot(data['res'][data['mode']==fm, ri], 'r', lw=0.3)
        plt.gca().tick_params(labelsize=6)
        plt.ylim(0, 3)
        sns.despine()
        if fm == 0:
            plt.ylabel(f'res-{ri + 1}', fontsize=8)
        if ri == 0:
            plt.title(data['modes'][fm], fontsize=8)
plt.tight_layout(w_pad=-0.75, h_pad=0)


# Basic analysis - performance of all 42 residuals
fsm = data['fsm']
fm = data['modes'][1:]

im = IsolabilityMatrix(fsm)
plt.figure(20, clear=True, figsize=(6, 6))
plt.spy(im[1:, 1:], marker='o', color='b')
plt.xticks(np.arange(len(fm)), fm)
plt.yticks(np.arange(len(fm)), fm)
plt.title('Isolability matrix')
plt.gca().xaxis.tick_bottom()


_, C = DiagnosesAndConfusionMatrix(thdata)
plt.figure(30, clear=True, figsize=(6, 6))
PlotConfusionMatrix(C)
plt.xticks(np.arange(nf), data['modes'])
plt.yticks(np.arange(nf), data['modes'])
plt.title('Fault Isolation Performance matrix, all 42 residuals')
plt.tight_layout()


# Test selection using Random Forest Classifiers

# First, build a random forest classifier based on the thresholded data.
# Here, 300 trees are trained in the tree ensemble.

res, C, rf, Crf = RandomForestTestSelection(thdata, n_estimators=100)

plt.figure(31, clear=True, figsize=(6, 6))
PlotConfusionMatrix(Crf)
plt.xticks(np.arange(nf), data['modes'])
plt.yticks(np.arange(nf), data['modes'])
plt.title('Fault Isolation Performance matrix')
plt.tight_layout()


# Plot the variable importance, sorted, to get a ranking of predictor/residual usefullness in the classifier.
# Note that this classifier is not meant to be used in the diagnosis system. (Fig. 10 in the paper)

plt.figure(40, clear=True, figsize=(9, 6))
plt.plot(res['residualimportance'][res['sortidx']])
plt.yticks(fontsize=8)
plt.xticks(range(nr), res['sortidx'] + 1, fontsize=8, rotation=90)
plt.xlabel('Predictors')
plt.ylabel('Importance')
plt.title('Predictor importance')
sns.despine()


# Compute performance measures on false-alarm (FA), missed detection (MD), aggregated fault isolation (FI)
# and the probability of maximum isolability performance (FI-max) when selecting residuals according to the
# ranking computed above. Plot the three aggregated performance measures agains the number of selected residuals.
# (Fig. 11 in the paper)

num_res = [10, 12, 26, 27]
plt.figure(50, clear=True, figsize=(9, 7))
plt.plot(range(1, nr), res['pfa'], 'r', label='False alarm probability')
plt.plot(range(1, nr), res['pmd'], 'b', label='Missed detection probability')
plt.plot(range(1, nr), res['pfi'], 'y', label='False isolation probability')
for ni in num_res:
    plt.plot(ni+1, res['pfi'][ni], 'kx')

plt.legend()
plt.xlabel('Number of selected residuals')
plt.ylabel('Probability')
sns.despine()


# Plot the probability of maximum fault isolation performance for each fault. (Fig. 12 in the paper)

plt.figure(figsize=(10, 10))
for k in range(nf):
    plt.plot(res['pmfi'][:, k], label=thdata['modes'][k])
plt.legend(loc='upper right')
sns.despine()


# Compute and display confusion matrices corresponding to selecting 10, 12, 26, and 27 residuals.
# The results should be compared to the confusion matrix above where all 42 residuals were used.
# (Fig. 13 in the paper)

# Compare performance of 12 with 42 residuals.

ntests = 12
plt.figure(90, clear=True, figsize=(12, 12))

plt.subplot(1, 2, 1)
_, C = DiagnosesAndConfusionMatrix(thdata, residx=res['sortidx'][0:ntests])
PlotConfusionMatrix(C)
plt.title('No of tests: %d' % ntests)
plt.xticks(np.arange(nf), data['modes'])
plt.yticks(np.arange(nf), data['modes'])

plt.subplot(1, 2, 2)
_, C = DiagnosesAndConfusionMatrix(thdata)
PlotConfusionMatrix(C)
plt.title('No of tests: 42')
plt.xticks(np.arange(nf), data['modes'])
plt.yticks(np.arange(nf), data['modes'])

plt.tight_layout()    




