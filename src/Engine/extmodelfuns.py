import numpy as np
from numpy import *
from scipy.interpolate import interp1d, interp2d

#syms igngrid(omega_e,p_im)
#syms a_eta_t_fun(W_t)
#syms b_eta_t_fun(W_t)

#syms eta_c_fun(eta_cmax,eta_cmin,W_ccorr,W_ccorrmax,PI_c,PI_cmax,Q_c11,Q_c12,Q_c22)
#syms cyl_temp_out(W_e)

#syms W_ic_fun(p_c,p_ic,plin_intercooler,H_intercooler,T_c)
#syms W_af_fun(p_amb,p_af,plin_af,H_af,T_amb)
#syms W_es_fun(p_amb,p_t,plin_exh,H_exhaust,T_t)
#syms PSI_th_fun(p_im,p_ic,gamma_air,PIli_th)
#syms PI_wg_fun(p_t,p_em,gamma_air)
#syms PSIli_wg_fun(PI_wg,PIli_wg,gamma_exh)
#syms Tflow_wg_fun(p_t,p_em,T_em,T_t)
#syms W_ig_fun(W_ac,AF_s,W_fc,q_HV,eta_otto,eta_ig_ch,omega_e)
#syms W_c_fun(PI_cnolim,p_af,U_c,D_c,T_af,R_air,PHI_model,fix_gain)
#syms PHI_model_fun(K1,K2,PSI_c)
#syms PI_c_fun(p_c,p_af)
#syms X_c_fun(W_ccorr,W_ccorrmax,PI_c,PI_cmax)
#syms W_t_fun(p_em,k1,PI_t,k2,T_em)
#syms Tq_t_fun(gamma_eg,cp_eg,eta_t,W_t,T_em,PI_t,omega_tc)
#syms PI_t_fun(p_t,p_em)
#syms reg_temp_fun(T_c,T_amb,W_c,Cic_1,Cic_2,Cic_3)

# Helper functions
#syms max_fun(choice_a,choice_b)
#syms min_fun(choice_a,choice_b)


def a_eta_t_fun( omega_tc ):
    a_omega = np.array([0, 9974928.83487578, 16575792.5370014, 32109852.6168100, 56054167.2758555, 76741191.7053602, 106900665.763039])
    omega_breakpoints = np.array([0, 6498.90800272609, 9788.46967029996, 13028.6036134574, 15822.8408393152, 18113.4807628027, 19969.5337025435])

    return interp1d(omega_breakpoints,a_omega)(omega_tc)

def b_eta_t_fun( omega_tc ):
    b_omega = np.array([0, -43392.0188412510, -103018.975128231, -262413.698439652, -533469.217487666, -777194.621957822, -1144332.45600859])
    omega_breakpoints = np.array([0, 6498.90800272609, 9788.46967029996, 13028.6036134574, 15822.8408393152, 18113.4807628027, 19969.5337025435])

    return interp1d(omega_breakpoints,b_omega)(omega_tc)


def eta_c_fun(eta_cmax,eta_cmin,W_ccorr,W_ccorrmax,PI_c,PI_cmax,Q_c11,Q_c12,Q_c22):
    X_c1 = W_ccorr-W_ccorrmax

    if PI_c>1: 
        X_c2 = np.sqrt(PI_c-1)+1-PI_cmax
    else:
        X_c2 = 1-PI_cmax
  
    eta_c_comparison = eta_cmax - (X_c1*X_c1*Q_c11+2*X_c1*X_c2*Q_c12+X_c2*X_c2*Q_c22)
  
    if eta_c_comparison >= 1:
        eta_cOut = 1 - 1e-9;
    elif eta_cmin < eta_c_comparison and eta_c_comparison < 1:
        eta_cOut = eta_c_comparison
    else:
        eta_cOut = eta_cmin
    return eta_cOut 


def W_ic_fun(p_c,p_ic,plin_intercooler,H_intercooler,T_c):
    if (p_c-p_ic-plin_intercooler)>0:
        return np.sqrt(np.abs(p_c/(H_intercooler*T_c)))*np.sqrt(p_c-p_ic)
    else:
        return np.sqrt(np.abs(p_c/(H_intercooler*T_c)))*(p_c-p_ic)/np.sqrt(plin_intercooler)

def W_af_fun(p_amb,p_af,plin_af,H_af,T_amb):
    if (p_amb-p_af-plin_af)>0:
        return np.sqrt(p_amb/(H_af*T_amb))*np.sqrt(p_amb-p_af)
    else:
        return np.sqrt(p_amb/(H_af*T_amb))*(p_amb-p_af)/np.sqrt(plin_af);

def W_es_fun(p_amb,p_t,plin_exh,H_exhaust,T_t):
    if plin_exh> (p_t-p_amb):
        return np.sqrt(p_t/(H_exhaust*T_t))*(p_t-p_amb)/np.sqrt(plin_exh)
    else:
        return np.sqrt(p_t/(H_exhaust*T_t))*np.sqrt(p_t-p_amb);

def PSI_th_fun(p_im,p_ic,gamma_air,PIli_th):
    PI = p_im/p_ic;

    if PI < PIli_th:
        PI_th = np.max([PI,(2/(gamma_air+1))**(gamma_air/(gamma_air-1))])
        return np.sqrt((2*gamma_air/(gamma_air-1))*(PI_th**(2/gamma_air)-PI_th**((gamma_air+1)/gamma_air)))
    else:
        PI_th = np.max([PIli_th,(2/(gamma_air+1))**(gamma_air/(gamma_air-1))])
        return np.sqrt((2*gamma_air/(gamma_air-1))*(PI_th**(2/gamma_air)-PI_th**((gamma_air+1)/gamma_air)))/(1-PIli_th)*(1 - PI)

def PI_wg_fun(p_t,p_em,gamma_air):
    if (-p_t/p_em+(2/(gamma_air+1))**(gamma_air/(gamma_air-1)))>0:
        return (2/(gamma_air+1))**(gamma_air/(gamma_air-1))
    else:
        return p_t/p_em

def PSI( PI, gamma):
    return np.sqrt( (2*gamma/(gamma-1))*(PI**(2/gamma) - PI**((gamma+1)/gamma)) )

def PSIli_wg_fun(PI_wg,PIli_wg,gamma_exh):
    if (-PI_wg+PIli_wg)>0:
        PI = np.max( [PI_wg, ( 2/(gamma_exh+1) )**( gamma_exh/(gamma_exh-1) )] )
        return PSI(PI,gamma_exh)
    else:
        PI = np.max( [PIli_wg, ( 2/(gamma_exh+1) )**( gamma_exh/(gamma_exh-1) ) ])
        return PSI(PI,gamma_exh)*(1-PI_wg)/(1-PIli_wg)

def Tflow_wg_fun(p_t,p_em,T_em,T_t):
    if p_em < p_t:
        return T_t
    else:
        return T_em

def W_ig_fun(W_ac,AF_s,W_fc,q_HV,eta_otto,eta_ig_ch,omega_e):
    return np.min([W_ac/(AF_s*W_fc),1])*W_fc*q_HV*eta_otto*eta_ig_ch/omega_e*(4*np.pi)

def W_c_fun(PI_cnolim,p_af,U_c,D_c,T_af,R_air,PHI_model,fix_gain):
    if (PI_cnolim - 1)>=0:
        return (p_af*U_c*pi*D_c**2/(T_af*2*R_air))*PHI_model
    else:
        return (p_af*U_c*np.pi*D_c**2/(T_af*2*R_air))*PHI_model*(1-(2-2*PI_cnolim)) + fix_gain*(2-2*PI_cnolim)

def PHI_model_fun(K1,K2,PSI_c):
    if ((1-K1*PSI_c**2)/K2)>0:
        return np.sqrt((1-K1*PSI_c**2)/K2)
    else:
        return 0.0

def PI_c_fun(p_c,p_af):
    if (p_c/p_af)>1:
        return p_c/p_af
    else:
        return 1.0 - 1e-9

def W_t_fun(p_em,k1,PI_t,k2,T_em):
    return (p_em*1e-3*k1*np.sqrt(np.max([1.0e-9,1-PI_t**k2])))/np.sqrt(T_em)

def Tq_t_fun(gamma_eg,cp_eg,eta_t,W_t,T_em,PI_t,omega_tc):
    return cp_eg*eta_t*W_t*T_em*np.max([0, 1-PI_t**((gamma_eg-1)/(gamma_eg))])/omega_tc

def PI_t_fun(p_t,p_em):
    if p_t/p_em > 1:
        return 1 - 1.0e-9
    else:
        return p_t/p_em

def reg_temp_fun(T_c,T_amb,W_c,Cic_1,Cic_2,Cic_3):
    if ((T_c-T_amb)*(Cic_1+(T_amb+T_c)*Cic_2/2+W_c*Cic_3)) > 273:
        return T_c - 273
    elif ((T_c-T_amb)*(Cic_1+(T_amb+T_c)*Cic_2/2+W_c*Cic_3)) < 0:
        return T_c
    else:
        return T_c - ((T_c-T_amb)*(Cic_1+(T_amb+T_c)*Cic_2/2+W_c*Cic_3))

def max_fun(choice_a,choice_b):
    return np.max([choice_a, choice_b])

def min_fun(choice_a,choice_b):
    return np.min([choice_a, choice_b])

def igngrid(imega_e,p_im):
    # Interpolate ignition angle as function of
    # engine speed, intake manifold pressure and ignition 2D look-up table
    #   Detailed explanation goes here
    dThGrid = np.array([
        [0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -5.68908e-02, -5.88293e-01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01],
        [0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -1.36223e-02, -8.07441e-01, -4.22760e+00, -1.19878e+01, -1.15178e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01],
        [0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -3.25697e-01, -2.23772e+00, -4.94829e+00, -8.77392e+00, -1.20997e+01, -1.39057e+01, -1.55571e+01, -1.70952e+01, -1.56320e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01],
        [0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -4.43190e-01, -3.46369e+00, -5.41470e+00, -7.72814e+00, -1.02632e+01, -1.21249e+01, -1.39619e+01, -1.56550e+01, -1.73184e+01, -1.88709e+01, -1.97153e+01, -2.02336e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01],
        [0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -8.62000e-02, -1.66926e+00, -3.61399e+00, -5.67795e+00, -8.55561e+00, -1.08785e+01, -1.25976e+01, -1.45814e+01, -1.59004e+01, -1.73517e+01, -1.83820e+01, -1.92232e+01, -1.97073e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01],
        [0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,  0.00000e+00, -5.77147e-01, -1.31053e+00, -2.70281e+00, -6.27204e+00, -9.24663e+00, -1.10989e+01, -1.30218e+01, -1.44008e+01, -1.58917e+01, -1.69460e+01, -1.78110e+01, -1.85908e+01, -1.97769e+01, -2.02998e+01, -2.02998e+01],
        [0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -4.01395e-02, -3.40622e-01, -1.02429e+00, -3.93728e+00, -6.99003e+00, -9.04656e+00, -1.09619e+01, -1.26314e+01, -1.38775e+01, -1.53915e+01, -1.61935e+01, -1.80491e+01, -1.92540e+01, -2.02998e+01, -2.02998e+01],
        [0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -5.78947e-02, -1.00000e-01, -5.88576e-01, -2.30318e+00, -5.23570e+00, -7.32383e+00, -9.28923e+00, -1.12614e+01, -1.25946e+01, -1.38110e+01, -1.48958e+01, -1.75274e+01, -1.87311e+01, -2.02998e+01, -2.02998e+01],
        [0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -6.49270e-02, -3.16327e-01, -1.54684e+00, -4.00289e+00, -6.57448e+00, -8.72839e+00, -1.04598e+01, -1.17786e+01, -1.31321e+01, -1.44884e+01, -1.69247e+01, -1.82082e+01, -2.02998e+01, -2.02998e+01],
        [0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -2.20683e-02, -5.60272e-03, -2.02359e-02, -4.28842e-01, -1.62681e+00, -3.93839e+00, -6.24253e+00, -8.42553e+00, -1.01078e+01, -1.16122e+01, -1.29991e+01, -1.43545e+01, -1.62633e+01, -1.76852e+01, -2.02998e+01, -2.02998e+01],
        [0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -4.69048e-02, 0.00000e+00, -7.04294e-02, -6.17368e-01, -2.06572e+00, -4.99934e+00, -6.87321e+00, -8.67248e+00, -1.02249e+01, -1.17765e+01, -1.32987e+01, -1.44689e+01, -1.57474e+01, -1.71623e+01, -1.83985e+01, -2.02998e+01],
        [0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -1.48693e-02, -6.65072e-02, -4.48565e-03, -5.98984e-01, -1.92382e+00, -4.29071e+00, -6.69579e+00, -8.83888e+00, -1.05148e+01, -1.21332e+01, -1.36353e+01, -1.45995e+01, -1.53986e+01, -1.66394e+01, -1.78724e+01, -2.02998e+01],
        [0.00000e+00, 0.00000e+00, -2.63158e-02, -1.20012e-02, 0.00000e+00, -2.63158e-02, -2.63158e-02, -4.02081e-01, -1.57529e+00, -4.18889e+00, -6.35175e+00, -8.26603e+00, -9.99547e+00, -1.14291e+01, -1.31123e+01, -1.40674e+01, -1.52962e+01, -1.61165e+01, -1.73494e+01, -2.02998e+01],
        [0.00000e+00, 0.00000e+00, -5.37205e-02, -4.72238e-02, -6.48339e-02, -2.20907e-01, -1.12724e-01, -5.75257e-01, -2.06047e+00, -4.28772e+00, -6.13409e+00, -7.66580e+00, -9.29463e+00, -1.06895e+01, -1.20211e+01, -1.35413e+01, -1.52724e+01, -1.55646e+01, -1.68265e+01, -2.02998e+01],
        [0.00000e+00, 0.00000e+00, 0.00000e+00, -7.26257e-02, -2.63158e-01, -1.19607e-01, -1.38158e-01, -1.11423e+00, -2.35216e+00, -4.35439e+00, -5.87412e+00, -7.24702e+00, -8.50691e+00, -9.65950e+00, -1.10283e+01, -1.28531e+01, -1.52554e+01, -1.51297e+01, -1.61995e+01, -2.02998e+01],
        [0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -3.23182e-02, -1.27059e-01, -2.97511e-01, -1.04808e+00, -2.28456e+00, -3.86254e+00, -5.20684e+00, -6.47333e+00, -7.76912e+00, -8.95529e+00, -1.00356e+01, -1.23924e+01, -1.46763e+01, -1.50338e+01, -1.55102e+01, -2.02998e+01],
        [0.00000e+00, 0.00000e+00, 0.00000e+00, -1.37788e-02, -3.15789e-02, -1.28444e-01, -2.47742e-01, -1.36697e+00, -3.01412e+00, -4.11608e+00, -5.27333e+00, -6.30427e+00, -7.85186e+00, -8.85404e+00, -1.00083e+01, -1.18490e+01, -1.40971e+01, -1.50068e+01, -2.02998e+01, -2.02998e+01],
        [0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -5.07458e-02, -2.29208e-01, -8.88953e-01, -2.22389e+00, -3.43865e+00, -4.62398e+00, -5.79535e+00, -7.11673e+00, -8.10389e+00, -9.60700e+00, -1.05154e+01, -1.14129e+01, -1.35179e+01, -1.50087e+01, -2.02998e+01, -2.02998e+01],
        [0.00000e+00, 0.00000e+00, -3.64035e-02, -1.01195e-01, -5.14225e-01, -1.15920e+00, -2.17880e+00, -3.20975e+00, -3.84711e+00, -4.66579e+00, -5.77850e+00, -6.98499e+00, -7.90927e+00, -1.04817e+01, -1.11104e+01, -1.13558e+01, -1.29467e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01],
        [0.00000e+00, 0.00000e+00, -4.21987e-02, -2.77493e-01, -1.96612e+00, -2.82002e+00, -2.64107e+00, -3.37338e+00, -4.32166e+00, -5.01108e+00, -5.65042e+00, -6.22493e+00, -7.43493e+00, -8.79125e+00, -1.00875e+01, -1.13837e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01]])
    pimGrid = 1.0e3*np.array([2.16000e+01, 3.25158e+01, 4.34316e+01, 5.43474e+01, 6.52632e+01, 7.61789e+01, 8.70947e+01, 9.80105e+01, 1.08926e+02, 1.19842e+02, 1.30758e+02, 1.41674e+02, 1.52589e+02, 1.63505e+02, 1.74421e+02, 1.85337e+02, 1.96253e+02, 2.07168e+02, 2.18084e+02, 2.29000e+02])
    NGrid = np.array([ 1.25000e+01, 1.71053e+01, 2.17105e+01, 2.63158e+01, 3.09211e+01, 3.55263e+01, 4.01316e+01, 4.47368e+01, 4.93421e+01, 5.39474e+01, 5.85526e+01, 6.31579e+01, 6.77632e+01, 7.23684e+01, 7.69737e+01, 8.15789e+01, 8.61842e+01, 9.07895e+01, 9.53947e+01, 1.00000e+02])

    return interp2(pimGrid,NGrid,dThGrid, kind='linear')(p_im,omega_e/(2*pi), fill_value=0)

