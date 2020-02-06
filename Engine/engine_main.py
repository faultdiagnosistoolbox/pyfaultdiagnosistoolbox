import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
import sys
from importlib import reload
new_paths = ['../Misc/', '../src']
[sys.path.append(d) for d in new_paths if not d in sys.path];

import numpy as np
import scipy
from scipy.stats import gaussian_kde
from misc import BoxOff, loadmat, RunResgenOnDataSets
import GetMeasurementData as gm


# # Model

import VEP4Engine
reload(VEP4Engine)
model = VEP4Engine.model
model.Lint()


# Plot the structural model
plt.figure(10)
model.PlotModel(verbose=False)


# Plot isolability properties and the extended Dulmage-Mendelsohn decomposition with equivalence classes for the over-determined part
plt.figure(20)
model.IsolabilityAnalysis(plot=True)

plt.figure(21)
model.IsolabilityAnalysis(plot=True, causality='int')

plt.figure(22)
model.IsolabilityAnalysis(plot=True, causality='der')

plt.figure(23)
model.PlotDM(fault=True,eqclass=True)


# # MSO sets and test selection
print("Searchin for MSO sets ...")
msos = model.MSO()
print("  Found " + str(len(msos)) + " MSO sets")

print("Checking index and observability properties of " + str(len(msos)) + " MSO sets ...")
li = list(map(lambda m:model.IsLowIndex(m), msos))
oi = list(map(lambda m:model.IsObservable(m), msos))
print('  ' + str(sum(oi)) + ' observable, ' + str(sum(li)) + ' low (structural) differential index')


# Use tests selected from simple Kullback-Leibler selection strategy
ts = [ 1649, 4011, 4016, 4017, 4066, 4074, 4477]
re = [ 73  , 75  , 75  , 3   , 76  , 3   , 2]
for msoIdx,redIdx in zip(ts,re):
    mso = msos[msoIdx]
    red = mso[redIdx]
    m0 = [e for e in mso if e != red]
    Gamma = model.Matching(m0)
    print("MSO " + str(msoIdx) + " with redundant equation " + str(red) + ", causality: " + Gamma.matchType)
FSM = model.FSM(msos[ts])


# Plot fault signature matrix and fault isolation matrix for selected set of tests
fIdx = [model.f.index(fi) for fi in ['fyw_af','fyp_im','fyp_ic','fyT_ic']]
plt.figure(30)
plt.subplot(1,2,1)
plt.spy(FSM[:,fIdx],markersize=6,marker="o", color="b")
plt.xticks(np.arange(0,len(fIdx)), np.array(model.f)[fIdx])
plt.yticks(np.arange(0,len(ts)), ["MSO " + str(k) for k in ts])
plt.gca().xaxis.tick_bottom()
plt.title('Fault Signature Matrix');

plt.subplot(1,2,2)
IM=model.IsolabilityAnalysisArrs(msos[ts])
plt.spy(IM[fIdx,:][:,fIdx],markersize=6,marker="o", color="b")
plt.xticks(np.arange(0,len(fIdx)), np.array(model.f)[fIdx])
plt.yticks(np.arange(0,len(fIdx)), np.array(model.f)[fIdx])
plt.gca().xaxis.tick_bottom()
plt.title('Fault isolation matrix');
plt.tight_layout()


# # Code generation
for test,redIdx in zip(ts,re):
    mso = msos[test]
    red = mso[redIdx]
    m0 = [e for e in mso if e != red]
    resName = "ResGen_" + str(test) + "_" + str(red)
    Gamma = model.Matching(m0)
    model.SeqResGen(Gamma, red, resName, language='C', batch=True, user_functions=VEP4Engine.ext_funs,
               external_headers=['extmodelfuns.h'], external_src=['extmodelfuns.cc'])
    print("")


# Compile sources
for test,redIdx in zip(ts,re):
    red = msos[test][redIdx]
    resName = "ResGen_" + str(test) + "_" + str(red)
    sys.stdout.write("Compiling residual generator: " + resName + " ... ")
    if os.system("python " + resName + "_setup.py build_ext --inplace")==0:
        sys.stdout.write('Success!\n')
    else:
        sys.stdout.write('Failure!\n')


# Import generated residual generator modules
import ResGen_1649_1
import ResGen_4011_1
import ResGen_4016_1
import ResGen_4017_84
import ResGen_4066_1
import ResGen_4074_84
import ResGen_4477_86


# # Import measurement data
dataDir = '/Users/frisk/Diagnosis/Work/EngineDiagnosis/Work/Data/'
dataSets = {'NF' : 'driving_cycle_FTP75_highway_no_fault_dataset1_16-01-20.mat',
            'fyp_im' : 'driving_cycle_FTP75Highway_fault_y_pim_dataset_7_26-01-2016.mat',
            'fyw_af' : 'driving_cycle_FTP75Highway_fault_y_waf_dataset_5_26-01-2016',
            'fyp_ic' : 'driving_cycle_FTP75Highway_fault_y_pic_dataset_8_26-01-2016',
            'fyT_ic' : 'driving_cycle_FTP75Highway_fault_y_Tic_dataset_6_26-01-2016'}
data = {}
sys.stdout.write('Loading data ')
for dd in dataSets:
    sys.stdout.write('.')
    data[dd] = gm.GetMeasurementData(dataDir + dataSets[dd])
sys.stdout.write(' Finished loading ' + str(len(data.keys())) + ' datasets\n')
ftp75_hwfet = loadmat(dataDir + 'ftp75_hwfet.mat')
ftp75_hwfet = {'t': ftp75_hwfet['ftp75_hwfet']['t'], 'v':ftp75_hwfet['ftp75_hwfet']['v']}


plt.figure(90)
plt.plot(ftp75_hwfet['t']/60.0, ftp75_hwfet['v'], linewidth=2)
BoxOff()
plt.xlabel('t [min]')
plt.ylabel('Velocity [km/h]')
plt.title('EPA Highway Fuel Economy Test Cycle (HWFET)');


DS = 500 # Down sampling rate in plots
plt.figure(80)
plt.subplot( 3,3,1 )
yIdx = model.z.index('y_omega_e')
plt.plot( data['NF']['time'][::DS]/60.0, data['NF']['z'][::DS,yIdx]/(2*np.pi)*60 )
plt.ylabel('rpm')
plt.title('Engine speed')
plt.xlim([data['NF']['time'][0], data['NF']['time'][-1]/60])
BoxOff()
plt.gca().get_yaxis().set_major_locator(matplotlib.ticker.MaxNLocator(4))

plt.subplot( 3,3,2 )
yIdx = model.z.index('y_p_im')
plt.plot( data['NF']['time'][::DS]/60.0, data['NF']['z'][::DS,yIdx]/1e3 )
plt.ylabel('kPa')
plt.title('Intake manifold pressure')
plt.xlim([data['NF']['time'][0], data['NF']['time'][-1]/60])
BoxOff()
plt.gca().get_yaxis().set_major_locator(matplotlib.ticker.MaxNLocator(4))

plt.subplot( 3,3,3 )
yIdx = model.z.index('y_W_af')
plt.plot( data['NF']['time'][::DS]/60.0, data['NF']['z'][::DS,yIdx] )
plt.ylabel('kg/s')
plt.title('Air mass flow')
plt.xlim([data['NF']['time'][0], data['NF']['time'][-1]/60])
BoxOff()
plt.gca().get_yaxis().set_major_locator(matplotlib.ticker.MaxNLocator(4))

plt.subplot( 3,3,4 )
yIdx = model.z.index('y_alpha_th')
plt.plot( data['NF']['time'][::DS]/60.0, data['NF']['z'][::DS,yIdx] )
plt.ylabel('%')
plt.title('Throttle position')
plt.xlim([data['NF']['time'][0], data['NF']['time'][-1]/60])
BoxOff()
plt.gca().get_yaxis().set_major_locator(matplotlib.ticker.MaxNLocator(4))

plt.subplot( 3,3,5 )
yIdx = model.z.index('y_wfc')
plt.plot( data['NF']['time'][::DS]/60.0, data['NF']['z'][::DS,yIdx]*1.0e6 )
plt.ylabel('mg/s')
plt.title('Injected fuel')
plt.xlim([data['NF']['time'][0], data['NF']['time'][-1]/60])
BoxOff()
plt.gca().get_yaxis().set_major_locator(matplotlib.ticker.MaxNLocator(4))

plt.subplot( 3,3,6 )
yIdx = model.z.index('y_u_wg')
plt.plot( data['NF']['time'][::DS]/60.0, data['NF']['z'][::DS,yIdx] )
plt.ylabel('%')
plt.title('Wastegate')
plt.xlim([data['NF']['time'][0], data['NF']['time'][-1]/60])
BoxOff()
plt.gca().get_yaxis().set_major_locator(matplotlib.ticker.MaxNLocator(4))

plt.subplot( 3,3,7 )
yIdx = model.z.index('y_p_ic')
plt.plot( data['NF']['time'][::DS]/60.0, data['NF']['z'][::DS,yIdx]/1e3 )
plt.xlabel('t [min]')
plt.ylabel('kPa')
plt.title('Intercooler pressure')
plt.xlim([data['NF']['time'][0], data['NF']['time'][-1]/60])
BoxOff()
plt.gca().get_yaxis().set_major_locator(matplotlib.ticker.MaxNLocator(4))

plt.subplot( 3,3,8 )
yIdx = model.z.index('y_T_ic')
plt.plot( data['NF']['time'][::DS]/60.0, data['NF']['z'][::DS,yIdx] )
plt.xlabel('t [min]')
plt.ylabel('K')
plt.title('Intercooler temperature')
plt.xlim([data['NF']['time'][0], data['NF']['time'][-1]/60])
BoxOff()
plt.gca().get_yaxis().set_major_locator(matplotlib.ticker.MaxNLocator(4))

plt.subplot( 3,3,9 )
yIdx = model.z.index('y_p_amb')
plt.plot( data['NF']['time'][::DS]/60.0, data['NF']['z'][::DS,yIdx]/1.0e3 )
plt.xlabel('t [min]')
plt.ylabel('kPa')
plt.title('Ambient pressure')
plt.xlim([data['NF']['time'][0], data['NF']['time'][-1]/60])
BoxOff()
plt.gca().get_yaxis().set_major_locator(matplotlib.ticker.MaxNLocator(4))
plt.ylim(99,101)
plt.tight_layout()
plt.subplots_adjust(top=0.9)
plt.suptitle('Measurement data, no-fault dataset',fontsize=14, weight='bold');


# # Run residual generators on measurement data
print("r1: MSO 1649")
r1 = RunResgenOnDataSets( ResGen_1649_1.ResGen_1649_1, data, VEP4Engine.diag_par )
print("\nr2: MSO 4011")
r2 = RunResgenOnDataSets( ResGen_4011_1.ResGen_4011_1, data, VEP4Engine.diag_par )
print("\nr3: MSO 4016")
r3 = RunResgenOnDataSets( ResGen_4016_1.ResGen_4016_1, data, VEP4Engine.diag_par )
print("\nr4: MSO 4017")
r4 = RunResgenOnDataSets( ResGen_4017_84.ResGen_4017_84, data, VEP4Engine.diag_par )
print("\nr5: MSO 4066")
r5 = RunResgenOnDataSets( ResGen_4066_1.ResGen_4066_1, data, VEP4Engine.diag_par )
print("\nr6: MSO 4074")
r6 = RunResgenOnDataSets( ResGen_4074_84.ResGen_4074_84, data, VEP4Engine.diag_par )
print("\nr7: MSO 4477")
r7 = RunResgenOnDataSets( ResGen_4477_86.ResGen_4477_86, data, VEP4Engine.diag_par )
r = [r1, r2, r3, r4, r5, r6, r7]
for ri in r:
    if 'NF' in ri:
        bias = np.mean(ri['NF'])*0.0
        normFact = np.max(np.abs(ri['NF']-bias))*1.1
        for FM in ri:
                ri[FM] = 1/normFact*(ri[FM]-bias)


#  Simple threshold selection strategy

alpha = 1e-3
J = np.zeros(len(r))
for rIdx, ri in enumerate(r):
    N = len(ri['NF'])
    foo = np.sort(np.abs(ri['NF']))
    J[rIdx] = foo[np.ceil(N*(1-alpha)).astype(np.int64)]


# Plot residuals

ds = 500
# dc = ['NF','fyp_im','fyw_af','fyp_ic','fyT_ic']
dc = ['NF', 'fyw_af']
for idx,fm in enumerate(dc):
    if fm in r[0]: # Data set exists in first residual
        plt.figure(50+idx)
        for ridx, ri in enumerate(r):
            plt.subplot( 3,3, ridx+1)
            plt.plot( data[fm]['time'][::ds]/60.0, ri[fm][::ds])
            plt.gca().get_yaxis().set_major_locator(matplotlib.ticker.MaxNLocator(4))
            plt.gca().get_xaxis().set_major_locator(matplotlib.ticker.MaxNLocator(4))
            if fm in model.f and FSM[ridx,model.f.index(fm)]==1:
                for t1Idx,t2Idx in zip(data[fm]['fault_idx'][0:-1:2],data[fm]['fault_idx'][1::2]):
                    y1, y2 = plt.ylim()
                    t1 = data[fm]['time'][t1Idx]/60.0
                    t2 = data[fm]['time'][t2Idx]/60.0
                    plt.gca().add_patch(mpatches.Rectangle((t1,y1), t2-t1, y2-y1, facecolor='0.9', edgecolor='none'))
            plt.plot( data[fm]['time'][::ds]/60.0, data[fm]['time'][::ds]*0+J[ridx], 'k--')
            plt.plot( data[fm]['time'][::ds]/60.0, data[fm]['time'][::ds]*0-J[ridx], 'k--')
            
            if fm in model.f and FSM[ridx,model.f.index(fm)]==1:
                plt.title('r%d: MSO %d (*)' % (ridx+1,ts[ridx]), fontsize=10, weight='bold')
            else:
                plt.title('r%d: MSO %d' % (ridx+1,ts[ridx]), fontsize=10, weight='bold')
            plt.xlim((0,np.max(data[fm]['time'])/60))
            BoxOff()
        plt.subplot(3,3,7)
        plt.xlabel('t [min]');
        plt.subplot(3,3,5)
        plt.xlabel('t [min]');
        plt.subplot(3,3,6)
        plt.xlabel('t [min]');
        plt.tight_layout()
        plt.subplots_adjust(top=0.9)
        plt.suptitle('Dataset: ' + fm,fontsize=12, weight='bold');


# Plot residual distributions
ds = 200
M = 50 # Number of data points in KDE plots
# dc = ['NF','fyp_im','fyw_af','fyp_ic','fyT_ic']
dc = ['NF', 'fyw_af']
for idx,fm in enumerate(dc):
    if fm in r[0]: # Data set exists in first residual
        plt.figure(60+idx)
        for ridx, ri in enumerate(r):
            plt.subplot( 3,3, ridx+1)
            pNF = gaussian_kde(ri[fm][data[fm]['fault_vector']==0][::ds])
            if not fm is 'NF':
                pF = gaussian_kde(ri[fm][data[fm]['fault_vector']==1][::ds])
            rmin = np.min([np.min(ri['NF']), np.min(ri[fm]), -1])
            rmax = np.max([np.max(ri['NF']), np.max(ri[fm]), 1])
            r_range = np.arange(0,M)*(rmax-rmin)*1.1/(M-1) + rmin - (rmax-rmin)*0.05
            plt.plot( r_range, pNF(r_range), 'b', linewidth=2)
            if not fm is 'NF':
                plt.plot( r_range, pF(r_range), 'r', linewidth=2)
            if fm in model.f and FSM[ridx,model.f.index(fm)]==1:
                plt.title('r%d: MSO %d (*)' % (ridx+1,ts[ridx]), fontsize=10, weight='bold')
            else:
                plt.title('r%d: MSO %d' % (ridx+1,ts[ridx]), fontsize=10, weight='bold')
            plt.gca().get_xaxis().set_major_locator(matplotlib.ticker.MaxNLocator(3))
            plt.gca().get_yaxis().set_major_locator(matplotlib.ticker.MaxNLocator(4))
            yl = plt.ylim()
            plt.plot([-1,-1],[0,0.3*yl[1]], 'k--')
            plt.plot([1,1],[0,0.3*yl[1]], 'k--')
            BoxOff()
        plt.tight_layout()
        plt.subplots_adjust(top=0.9)
        plt.suptitle('Dataset: ' + fm,fontsize=12, weight='bold')





ds = 200
M = 50 # Number of data points in KDE plots
# dc = ['NF','fyp_im','fyw_af','fyp_ic','fyT_ic']
dc = ['NF', 'fyw_af']
for idx,fm in enumerate(dc):
    if fm in r[0]: # Data set exists in first residual
        _ = plt.figure(60+idx)
        for ridx, ri in enumerate(r):
            _ = plt.subplot( 3,3, ridx+1)
            pNF = gaussian_kde(ri[fm][data[fm]['fault_vector']==0][::ds])
            if not fm is 'NF':
                pF = gaussian_kde(ri[fm][data[fm]['fault_vector']==1][::ds])
            rmin = np.min([np.min(ri['NF']), np.min(ri[fm]), -1])
            rmax = np.max([np.max(ri['NF']), np.max(ri[fm]), 1])
            r_range = np.arange(0,M)*(rmax-rmin)*1.1/(M-1) + rmin - (rmax-rmin)*0.05
            _ = plt.plot( r_range, pNF(r_range), 'b', linewidth=2)
            if not fm is 'NF':
                _ = plt.plot( r_range, pF(r_range), 'r', linewidth=2)
            if fm in model.f and FSM[ridx,model.f.index(fm)]==1:
                _ = plt.title('r{}: MSO {} (*)'.format(ridx+1,ts[ridx]),
                              fontsize=10, weight='bold')
            else:
                plt.title('r%d: MSO %d' % (ridx+1,ts[ridx]), fontsize=10, weight='bold')
            plt.gca().get_xaxis().set_major_locator(matplotlib.ticker.MaxNLocator(3))
            plt.gca().get_yaxis().set_major_locator(matplotlib.ticker.MaxNLocator(4))
            yl = plt.ylim()
            plt.plot([-1,-1],[0,0.3*yl[1]], 'k--')
            plt.plot([1,1],[0,0.3*yl[1]], 'k--')
            BoxOff()
        plt.tight_layout()
        plt.subplots_adjust(top=0.9)
        plt.suptitle('Dataset: ' + fm,fontsize=12, weight='bold')

plt.show()
