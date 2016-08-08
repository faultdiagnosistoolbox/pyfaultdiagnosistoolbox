from VEP_par import TC_ENGINE
from misc import loadmat
from scipy.interpolate import interp1d
import numpy as np

def GetMeasurementData(fileName):
    """ GetMeasurementData(fileName)

    Input:
      fileName - filename of .mat-file with data

    Output:
      dictionary with the keys
      
        time - Global time vector
        z    - Known signal matrix, each signal as a column in the order
               y_p_ic     [Pa]
               y_p_im     [Pa]
               y_T_ic     [K]
               y_W_af     [kg/s]
               y_omega_e  [rad/s]
               y_alpha_th [%]
               y_u_wg     [%]
               y_wfc      [kg/s]
               y_T_amb    [K]
               y_p_amb    [Pa]
        fault_idx - Indices where faults are activated/deactivated
        fault_vector - 0/1-vecotr indicating active fault
        Ts - Sampling time
        state_init - dictionary with all states for a template state initialization
    """
    data = loadmat(fileName)
    if 'fault_idx' in data:
        fault_idx = data['fault_idx']
    else:
        fault_idx = []

    time3 = data['time_vectors']['time3']
    time4 = data['time_vectors']['time4']

    Ts = 0.001
    t0 = np.min(time3)
    time = t0 + np.arange(0,len(time3))*Ts

    W_fc_cyl1 = data['W_fc_cyl1']
    W_fc_cyl2 = data['W_fc_cyl2']
    W_fc_cyl3 = data['W_fc_cyl3']
    W_fc_cyl4 = data['W_fc_cyl4']

    y_T_amb = interp1d(time4, data['y_T_amb'],fill_value='extrapolate')(time) + 273.0 # [C] -> [K]

    y_p_ic = data['y_p_ic']*1.0e3 # [kPa] -> [Pa]
    y_p_im = data['y_p_im']*1.0e3 # [kPa] -> [Pa]
    y_T_ic = data['y_T_ic']+273.0 # [C] -> [K]
    y_W_af = data['y_W_af']/1.0e3 # [g/s] -> [kg/s]
    y_omega_e = data['y_omega_e']/(60.0/(2.0*np.pi)) # [rpm] -> [rad/s]
    y_alpha_th= data['y_alpha_th'] # [?]
    y_u_wg = np.zeros(len(y_alpha_th)) # [?]

    y_wfc = (W_fc_cyl1 + W_fc_cyl2 + W_fc_cyl3 + W_fc_cyl4)  
    y_wfc = interp1d(time4,y_wfc,fill_value='extrapolate')(time)/120.0*y_omega_e/1.0e6 # [kg/s]

    y_p_amb = data['y_p_amb']*1.0e3 # [KPa] -> [Pa]

    measurement = np.transpose(np.vstack(
        (y_p_ic,y_p_im, y_T_ic, y_W_af, y_omega_e, y_alpha_th, y_u_wg, y_wfc, y_T_amb, y_p_amb)))

    fault_vector = np.zeros(len(time))
    for t1,t2 in zip(fault_idx[0:-1:2],fault_idx[1::2]):
        fault_vector[t1:t2] = 1.0

    m_af_init = TC_ENGINE.controlVolumes.airFilter.p_init*TC_ENGINE.controlVolumes.airFilter.V/\
        TC_ENGINE.gasProp.air.R/TC_ENGINE.controlVolumes.airFilter.T_init

    m_c_init = TC_ENGINE.controlVolumes.Compressor.p_init*TC_ENGINE.controlVolumes.Compressor.V/\
        TC_ENGINE.gasProp.air.R/TC_ENGINE.controlVolumes.Compressor.T_init

    m_ic_init = TC_ENGINE.controlVolumes.interCooler.p_init*TC_ENGINE.controlVolumes.interCooler.V/\
        TC_ENGINE.gasProp.air.R/TC_ENGINE.controlVolumes.interCooler.T_init

    m_t_init = TC_ENGINE.controlVolumes.Turbine.p_init*TC_ENGINE.controlVolumes.Turbine.V/\
        TC_ENGINE.gasProp.exh.R/TC_ENGINE.controlVolumes.Turbine.T_init

    m_em_init = TC_ENGINE.controlVolumes.exhaustManifold.p_init*\
        TC_ENGINE.controlVolumes.exhaustManifold.V/TC_ENGINE.gasProp.exh.R/\
        TC_ENGINE.controlVolumes.exhaustManifold.T_init

    m_im_init = TC_ENGINE.controlVolumes.intakeManifold.p_init*\
        TC_ENGINE.controlVolumes.intakeManifold.V/\
        TC_ENGINE.gasProp.air.R/TC_ENGINE.controlVolumes.intakeManifold.T_init

    state_init = {
        'alpha_th' : 0,
        'm_af'     : m_af_init,
        'T_af'     : TC_ENGINE.controlVolumes.airFilter.T_init,
        'm_c'      : m_c_init,
        'T_c'      : TC_ENGINE.controlVolumes.Compressor.T_init,
        'm_ic'     : m_ic_init,
        'T_ic'     : TC_ENGINE.controlVolumes.interCooler.T_init,
        'wg_pos'   : 0,
        'm_t'      : m_t_init,
        'T_t'      : TC_ENGINE.controlVolumes.Turbine.T_init,
        'm_em'     : m_em_init,
        'T_em'     : TC_ENGINE.controlVolumes.exhaustManifold.T_init,
        'omega_tc' : TC_ENGINE.Turbo.omega_init,
        'm_im'     : m_im_init,
        'T_im'     : TC_ENGINE.controlVolumes.intakeManifold.T_init
    }

    res = {}
    res['time'] = time
    res['z'] = measurement
    res['fault_idx'] = fault_idx
    res['fault_vector'] = fault_vector
    res['Ts'] = Ts
    res['state_init'] = state_init
    return res
