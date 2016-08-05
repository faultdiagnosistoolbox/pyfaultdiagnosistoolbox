import faultdiagnosistoolbox as fdt
import sympy as sym
from VEP_par import TC_ENGINE, diag_par

modeldef = {}
modeldef['type'] = 'Symbolic'
modeldef['x'] = [
  'W_af',                         # Restriktion Luftfilter
  'p_t','W_es',                   # Restriktion Avgassystem
  'p_ic','T_ic','W_ic',           # Restriktion Intercooler (T_c och p_c finns i kompressorn)
  'W_th','Aeff_th',               # Restriktion Throttle (anv?nder ?ven T_im, T_ic, p_im, p_ic)
  'W_wg','Aeff_wg',               # Restriktion Wastegate (anv?nder ?ven p_em,p_t,T_em,T_t)
  'T_af','p_af',                  # Kontrollvolym Luftfilter
  'W_c','T_cout','T_c','p_c',     # Kontrollvolym kompressor (Ska c_v vara ett tillst?nd?)
  'T_imcr','p_im','T_im',         # Kontrollvolym insugsr?r (Ska c_v vara ett tillst?nd?)
  'W_e','T_ti','T_em','p_em',     # Kontrollvolym avgasr?r (ska c_v vara ett tillst?nd?)
  'dh_is', 
  'W_twg','T_turb','T_t',   # Kontrollvolym turbin (ska c_v vara ett tillst?nd?)
  'alpha_th',                     # Aktuator Throttle
  'omega_e',                      # Massfl?de cylinder m. fyllnadsgrad
  'W_ig','W_fr',    
  'Tq_c','eta_c','omega_tc',      # Kompressormodell (Ska PHI_model ("dimensionless flow coeff.") vara med h?r?)
  'PSI_c','PI_c', 
  'W_t','Tq_t','eta_t','PI_t',    # Turbin- och wastegatemodell
  'u_wg','wg_pos',                # Aktuator Wastegate
  'DELTA_theta',                  # T?ndningseffektivitet
  'T_e', 
  'T_amb','p_amb' , 
  'PSI_th','PI_wg',   
  'PSIli_wg','m_af','m_c','m_ic','T_fwd_flow_ic', 
  'm_im','m_em','m_t', 
  'PI_cnolim','U_c','PHI_model','W_ccorr', 
  'Tq_e_cs','Tq_e_cb','Tq_e_f','Tq_e_p','eta_ign', 
  'W_i_p','FMEP','S_p','BMEP', 
  'C_eta_vol','T_in','eta_vol','W_ac','W_fc','T_tout', 
  'Tflow_wg', 
  'dmdt_af','dTdt_af','dmdt_c','dTdt_c',              
  'dmdt_ic','dTdt_ic','dmdt_im','dTdt_im', 
  'dmdt_em','dTdt_em','dmdt_t','dTdt_t', 
  'domegadt_tc','dwgdt_pos']

modeldef['f'] = [
    'fp_af', 
    'fw_af',     #L?ckage, variabel diameter, efter MAF, fore kompressor
    'fw_th',     #L?ckage, variabel diameter, efter gasspj?ll, i insugsror
    'fw_c',      #L?ckage, variabel diameter, efter kompressor, fore intercooler
     #,'fw_ic',     #L?ckage, variabel diameter, efter intercooler, f?re gasspj?ll
     'fc_vol',    #Insugsventilens stalldon har fastnat vid godtycklig position
     'fw_t',      #Okad friktion i turbinen.
     'fx_th',     #Matfel pa positionssensor f?r gasspjallslaget.
     'fyw_af',    #Matfel pa luftmassflodessensorn (MAF)
     'fyp_im',        #Matfel pa trycksensorn i insugsroret
     'fyp_ic',        #Matfel pa trycksensorn monterad i turboaggregatet (laddtrycket)
     'fyT_ic']       #Matfel pa temperatursensorn monterad i turboaggregatet.

modeldef['z'] = [
    #'y_p_af',
    'y_p_ic',
    'y_p_im','y_T_ic',
    'y_W_af','y_omega_e',
    'y_alpha_th',#'y_lambda',
    'y_u_wg','y_wfc', 
    'y_T_amb','y_p_amb']

modeldef['parameters'] = [
    'H_af','plin_af','H_exhaust','plin_exh', 
    'H_intercooler','plin_intercooler','PIli_th','gamma_air', 
    'PIli_wg','gamma_exh','R_air','cp_air','V_af', 
    'V_c','V_ic','V_im','R_exh','cp_exh', 
    'lambda_af','n_r','r_c','V_D','CEva_cool','AF_s','q_HV', 
    'eta_ig_ch','s','PI_bl','aux_dev_fric','B', 
    'D_c','K1','K2','fix_gain','R_a','cp_a', 
    'eta_cmax','eta_cmin','PI_c_max','T_std', 
    'k1','k2','cp_eg','k3','k4','D_t','gamma_eg','tau_wg', 
    'h_tot', 'Amax', 
    'A_0','A_1','A_2', 
    'J_tc','xi_fric_tc','V_em','V_t','c_2','c_3','p_std', 
    'eta_otto','C_tq1','C_tq2','xi_aux', 
    'a1','a2','a3','a4','a5','Q_c11','Q_c12','Q_c22','eta_tmin','eta_tmax','Cd','PI_cmax', 
    'cv_exh','cv_air','tau_th', 
    'W_ccorrmax','A_em', 
    'K_t','T0','h_ext','c0_em','c1_em','c2_em','D','L','Pr','my_exh','lambda_exh', 
    'Cic_1','Cic_2','Cic_3', 'TOL']

sym.var(modeldef['x']);
sym.var(modeldef['f']);
sym.var(modeldef['z']);
sym.var(modeldef['parameters']);

igngrid = sym.Function('igngrid')
a_eta_t_fun = sym.Function('a_eta_t_fun')
b_eta_t_fun = sym.Function('b_eta_t_fun')
eta_c_fun = sym.Function('eta_c_fun')
cyl_temp_out = sym.Function('cyl_temp_out')
W_ic_fun = sym.Function('W_ic_fun')
W_af_fun = sym.Function('W_af_fun')
W_es_fun = sym.Function('W_es_fun')
PSI_th_fun = sym.Function('PSI_th_fun')
PI_wg_fun = sym.Function('PI_wg_fun')
PSIli_wg_fun = sym.Function('PSIli_wg_fun')
Tflow_wg_fun = sym.Function('Tflow_wg_fun')
W_ig_fun = sym.Function('W_ig_fun')
W_c_fun = sym.Function('W_c_fun')
PHI_model_fun = sym.Function('PHI_model_fun')
PI_c_fun = sym.Function('PI_c_fun')
X_c_fun = sym.Function('X_c_fun')
W_t_fun = sym.Function('W_t_fun')
Tq_t_fun = sym.Function('Tq_t_fun')
PI_t_fun = sym.Function('PI_t_fun')
reg_temp_fun = sym.Function('reg_temp_fun')

max_fun = sym.Function('max_fun')
min_fun = sym.Function('min_fun')


ext_funs = {
    'syms igngrid'  : 'syms igngrid',
    'a_eta_t_fun'   : 'a_eta_t_fun',
    'b_eta_t_fun'   : 'b_eta_t_fun',
    'eta_c_fun'     : 'eta_c_fun',
    'cyl_temp_out'  : 'cyl_temp_out',
    'W_ic_fun'      : 'W_ic_fun',
    'W_af_fun'      : 'W_af_fun',
    'W_es_fun'      : 'W_es_fun',
    'PSI_th_fun'    : 'PSI_th_fun',
    'PI_wg_fun'     : 'PI_wg_fun',
    'PSIli_wg_fun'  : 'PSIli_wg_fun',
    'Tflow_wg_fun'  : 'Tflow_wg_fun',
    'W_ig_fun'      : 'W_ig_fun',
    'W_c_fun'       : 'W_c_fun',
    'PHI_model_fun' : 'PHI_model_fun',
    'PI_c_fun'      : 'PI_c_fun',
    'X_c_fun'       : 'X_c_fun',
    'W_t_fun'       : 'W_t_fun',
    'Tq_t_fun'      : 'Tq_t_fun',
    'PI_t_fun'      : 'PI_t_fun',
    'reg_temp_fun'  : 'reg_temp_fun',
    'max_fun'       : 'max_fun',
    'min_fun'       : 'min_fun' }

modeldef['rels'] = [
    # Restriction intercooler
    -W_ic + W_ic_fun(p_c,p_ic,plin_intercooler,H_intercooler,T_c),    
    # Restriction air-filter
    -W_af + W_af_fun(p_amb,p_af,plin_af,H_af,T_amb)+fp_af,
    
    # Restriction exhaust
    -W_es + W_es_fun(p_amb,p_t,plin_exh,H_exhaust,T_t),
    
    # Restriction throttle
    -PSI_th + PSI_th_fun(p_im,p_ic,gamma_air,PIli_th),
    -W_th + (Aeff_th*p_ic)/(sym.sqrt(R_air*T_ic))*PSI_th,
    -T_imcr + T_ic,
    
    # Restriction wastegate:
    -PI_wg + PI_wg_fun(p_t,p_em,gamma_air),
    -PSIli_wg + PSIli_wg_fun(PI_wg,PIli_wg,gamma_exh),
    -W_wg + (Aeff_wg*p_em)/(sym.sqrt(R_exh*Tflow_wg))*PSIli_wg,
    -Tflow_wg + Tflow_wg_fun(p_t,p_em,T_em,T_t),
    
    # Control volume air filter
    -p_af + m_af*R_air*T_af/V_af,
    -dmdt_af + W_af - W_c + fw_af,
    -dTdt_af + (W_af*cv_air*(T_amb-T_af)+R_air*(T_amb*W_af-T_af*W_c))/(m_af*cv_air),
    fdt.DiffConstraint('dmdt_af','m_af'),
    fdt.DiffConstraint('dTdt_af','T_af'),
    
    # Control volume, compressor
    -p_c + m_c*R_air*T_c/V_c,
    -dmdt_c + W_c - W_ic + fw_c,
    -dTdt_c + (W_c*cv_air*(T_cout-T_c)+R_air*(T_cout*W_c-T_c*W_ic))/(m_c*cv_air),
    fdt.DiffConstraint('dmdt_c','m_c'),
    fdt.DiffConstraint('dTdt_c','T_c'),
    
    # Control volume air cooler
    -p_ic + m_ic*R_air*T_ic/V_ic,
    -dmdt_ic + W_ic - W_th,
    -dTdt_ic + (W_ic*cv_air*(T_fwd_flow_ic-T_ic)+R_air*(T_fwd_flow_ic*W_ic-T_ic*W_th))/(m_ic*cv_air),
    fdt.DiffConstraint('dmdt_ic','m_ic'),
    fdt.DiffConstraint('dTdt_ic','T_ic'),
    
    # Control volume, air intake
    -p_im + m_im*R_air*T_im/V_im,
    -dmdt_im + (W_th - W_ac) + fw_th,
    -dTdt_im + (W_th*cv_air*(T_imcr-T_im) +R_air*(T_imcr*W_th-T_im*W_ac))/(m_im*cv_air),
    fdt.DiffConstraint('dmdt_im','m_im'),
    fdt.DiffConstraint('dTdt_im','T_im'),
    
    # Control volume, exhaust
    -p_em + m_em*R_exh*T_em/V_em,
    -dmdt_em + W_e - W_twg,
    -dTdt_em + (W_e*cv_exh*(T_ti-T_em)+R_exh*(T_ti*W_e-T_em*W_twg))/(m_em*cv_exh),
    fdt.DiffConstraint('dmdt_em','m_em'),
    fdt.DiffConstraint('dTdt_em','T_em'),
    
    # Control volume, turbine
    -p_t + m_t*R_exh*T_t/V_t,
    -dmdt_t + W_twg - W_es,
    -dTdt_t + (W_twg*cv_exh*(T_turb-T_t)+R_exh*(T_turb*W_twg-T_t*W_es))/(m_t*cv_exh),
    fdt.DiffConstraint('dmdt_t','m_t'),
    fdt.DiffConstraint('dTdt_t','T_t'),
    
    # Intercooler, temperature model
    -T_fwd_flow_ic + max_fun(T_amb, T_c - max_fun(TOL,(T_c-T_amb)*(Cic_1 + (T_amb+T_c)*Cic_2/2+W_ic*Cic_3))),
    
    # Throttle actuator
    -Aeff_th + A_0 + A_1*(alpha_th*86/100+4)*sym.pi/180 + A_2*((alpha_th*86/100+4)*sym.pi/180)**2  ,
#     ,dalphadt_th == 1/tau_th*(Beta_ped*100-alpha_th)
#     ,DiffConstraint(dalphadt_th,alpha_th)
    
    # Mass flow with volumetric efficiency and fuel flow into the cylinder
    -W_ac + eta_vol*V_D*1/(2*sym.pi)*omega_e*p_im/(n_r*R_air*T_im),
    -T_in + T_im - CEva_cool*(1/lambda_af - 1),
    #,eta_vol == (1+fc_vol)*C_eta_vol*((r_c-(p_em/p_im)**(1/gamma_air))*T_im)/((r_c-1)*T_in) 
    #((r_c-(p_em/p_im)**(1/gamma_air))*T_im)/((r_c-1)*T_in) 
    -eta_vol + (1+fc_vol)*C_eta_vol*((r_c-(p_em/p_im)**(1/gamma_exh))*T_im)/((r_c-1)*T_in),
    -C_eta_vol + a1*(max_fun(p_im,2.5e4)*1e-3)**4 + a2*(max_fun(p_im,2.5e4)*1e-3)**3 + a3*(max_fun(p_im,2.5e4)*1e-3)**2 + a4*(max_fun(p_im,2.5e4)*1e-3) + a5,
    
    # Bransleflode in till cylindern
    #,W_fc == W_ac/(lambda_af*AF_s) ## Anv?nder insignal ist?llet
    
    # Engine torque
    -Tq_e_cs + Tq_e_cb-Tq_e_f-Tq_e_p,
    -Tq_e_cb + eta_ign*W_ig/(4*sym.pi),
    -Tq_e_f + W_fr/(4*sym.pi),
    -Tq_e_p + W_i_p/(4*sym.pi),
    -W_ig + W_ig_fun(W_ac,AF_s,W_fc,q_HV,eta_otto,eta_ig_ch,omega_e),
    -W_fr + V_D*FMEP,
    -FMEP + xi_aux*((0.464+0.0072*S_p**(1.8))*PI_bl*10**5+0.0215*BMEP)*(0.075/B)**0.5,
    -S_p + omega_e*0.5/(2*sym.pi*s),
    -BMEP + (C_tq1+C_tq2*p_im),
    -W_i_p + V_D*(p_em-p_im),
    -W_e + W_fc + W_ac,
    
    # Engine-out temperature
    -T_e + T0+W_e*K_t,
    
    # Exhaust temperature drop
    -T_ti + T_amb+(T_e-T_amb)*sym.exp(-h_tot*A_em/(W_e/4*cp_exh)),
    
    # Compressor
    -Tq_c + cp_air*W_c*(T_cout-T_af)/omega_tc,
    -W_c + W_c_fun(PI_cnolim,p_af,U_c,D_c,T_af,R_air,PHI_model,fix_gain),
    -PHI_model + PHI_model_fun(K1,K2,PSI_c),
    -PSI_c + (2*cp_air*((PI_c**((gamma_air-1)/gamma_air))-1)*T_af)/(U_c**2),
    -PI_c + PI_c_fun(p_c,p_af),
    -PI_cnolim + p_c/p_af,
    -U_c + omega_tc*D_c/2,
    -T_cout + T_af*(1+(PI_c**((gamma_air-1)/gamma_air)-1)/eta_c),
    -eta_c + eta_c_fun(eta_cmax,eta_cmin,W_ccorr,W_ccorrmax,PI_c,PI_cmax,Q_c11,Q_c12,Q_c22),
    -W_ccorr + W_c*sym.sqrt(T_af/T_std)/(p_af/p_std),
    
    # Turbine
    -W_t + W_t_fun(p_em,k1,PI_t,k2,T_em),
    -Tq_t + Tq_t_fun(gamma_eg,cp_eg,eta_t,W_t,T_em,PI_t,omega_tc),
    #,eta_t == (a_eta_t_fun(W_t)*k1.*sym.sqrt(1-PI_t**k2)+b_eta_t_fun(W_t))/dh_is #
    -eta_t + max_fun(TOL, min_fun(1-TOL, (a_eta_t_fun(omega_tc)*sym.sqrt(T_em)*W_t/(p_em*1e-3)+b_eta_t_fun(omega_tc))/dh_is)),
    -dh_is + cp_exh*T_em*(1-PI_t**((gamma_exh-1)/gamma_exh)),
# Efficiency:
    -T_tout + T_em*(1-eta_t*(1-PI_t**((gamma_eg-1)/gamma_eg))),
    -PI_t + PI_t_fun(p_t,p_em),
    
    # Adiabatic mixer
    -T_turb + (W_t*T_tout + W_wg*Tflow_wg)/(W_t + W_wg),
    -W_twg + W_t + W_wg,
    
    # Wastegate actuator
    -Aeff_wg + Amax*Cd*wg_pos,
    -dwgdt_pos + 1/tau_wg*(u_wg-wg_pos),
    fdt.DiffConstraint('dwgdt_pos','wg_pos'),
    
    # Ignition efficiency
    -eta_ign + 1-(DELTA_theta/100)**2*c_2-(DELTA_theta/100)**3*c_3,
    
    # Ignition angle modifier
    -DELTA_theta + igngrid(omega_e,p_im),
    
    # Turbine compressor, dynamics
    -domegadt_tc + (1-fw_t)*1/J_tc*(Tq_t-Tq_c-omega_tc*xi_fric_tc),
    #,DiffConstraint(dWdt_t,W_t)
    #,DiffConstraint(dWdt_c,W_c)
    fdt.DiffConstraint('domegadt_tc','omega_tc'),
    
    # Measurements
    #,y_p_af == p_af               # No measurement, computed
    -y_p_ic + p_ic + fyp_ic,        # Pressure before throttle
    -y_p_im + p_im + fyp_im,        # Pressure in the intake manifold
    -y_T_ic + T_ic + fyT_ic,        # Temperature before throttle
    -y_W_af + W_af + fyw_af,        # Air mass-flow after air filter
    -y_omega_e + omega_e,           # Engine speed
    -y_alpha_th + alpha_th + fx_th, # Throttle position
    -y_u_wg + u_wg,                 # Wastegate reference
    -y_wfc + W_fc,                  # Computed total fuel mass flow
    -y_T_amb + T_amb,               # Ambient temperature
    -y_p_amb + p_amb]               # Ambient pressure

model = fdt.DiagnosisModel(modeldef, name='VEP4');
