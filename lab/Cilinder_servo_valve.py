import faultdiagnosistoolbox as fdt
import sympy as sym

modelDef={}
modelDef['type'] = 'Symbolic'

modelDef['x'] = ['q_A1', 'q_A2', 'p_A1', 'p_A2', 'V_A1', 'V_A2', 'q_leak_A',
                 'vp_A', 'dvp_A', 'xp_A', 'dxp_A', 'ap_A', 'v_A', 'x_A', 'F_A', 'dV_A1', 'dV_A2',
                 'dp_A1', 'dp_A2', 'c_leak_f1_A', 'Bp_f2_A', 'beta_f3_A',
                'K_f4_A', 'q_AT', 'q_AP', 'p_AT', 'p_AP', 'y', 'y_f5_A']
modelDef['z'] = ['s1','s2','s3','s4','s5','s6','s7']
modelDef['f'] = ['f1_A','f2_A','f3_A','f4_A','f5_A']
modelDef['parameters'] = ['A_A1','A_A2','V0_A1','V0_A2','Bp_A','beta_A',
                       'c_leak_A','K_A','m_A','k_f1_A','k_f2_A','k_f3_A',
                       'k_f4_A', 'k_f5_A', 'Cq_A', 'rho_A']

sym.var(modelDef['x'])
sym.var(modelDef['f'])
sym.var(modelDef['z'])
sym.var(modelDef['parameters'])

modelDef['rels'] = [
    -V_A1 + V0_A1 + A_A1*(-xp_A),
    -V_A2 + V0_A2 + A_A2*(xp_A),
    -q_leak_A + c_leak_f1_A*(p_A1 - p_A2),
    -q_A1 + q_leak_A + dV_A1 + dp_A1*V_A1/beta_f3_A,
    -q_A2 - q_leak_A + dV_A2 + dp_A2*V_A2/beta_f3_A, 
    -p_A2*A_A2 + p_A1*A_A1 + m_A*ap_A + Bp_f2_A*vp_A + K_f4_A*(xp_A - x_A),
    -F_A + K_f4_A*(xp_A - x_A),
    -vp_A + dxp_A,
    -ap_A + dvp_A,
    -q_AP + Cq_A*y_f5_A*sym.sqrt(2.0/rho_A) * sym.sqrt(sym.functions.Abs(p_AP - p_A1)) * sym.sign(p_AP - p_A1), 
    -q_A1 + Cq_A*y_f5_A*sym.sqrt(2.0/rho_A) * sym.sqrt(sym.functions.Abs(p_AP - p_A1)) * sym.sign(p_AP - p_A1), 
    -q_AT + -Cq_A*y_f5_A*sym.sqrt(2.0/rho_A) * sym.sqrt(sym.functions.Abs(p_A2 - p_AT)) * sym.sign(p_A2 - p_AT),
    -q_A2 + -Cq_A*y_f5_A*sym.sqrt(2.0/rho_A) * sym.sqrt(sym.functions.Abs(p_A2 - p_AT)) * sym.sign(p_A2 - p_AT), 
    -c_leak_f1_A + c_leak_A*(1 + k_f1_A*f1_A),
    -Bp_f2_A + Bp_A *(1 + k_f2_A*f2_A), 
    -beta_f3_A + beta_A*(1 + k_f3_A*f3_A),
    -K_f4_A + K_A*(1 + k_f4_A*f4_A), 
    -y_f5_A + y*(1 + k_f5_A*f5_A), 
    -s1 + q_AP,
    -s2 + p_AP,
    -s3 + q_AT,
    -s4 + p_AT,
    -s5 + x_A,
    -s6 + y,
    -s7 + F_A,
    fdt.DiffConstraint('dV_A1','V_A1'),
    fdt.DiffConstraint('dV_A2','V_A2'),
    fdt.DiffConstraint('dp_A1','p_A1'),
    fdt.DiffConstraint('dp_A2','p_A2'),
    fdt.DiffConstraint('dxp_A','xp_A'),
    fdt.DiffConstraint('dvp_A','vp_A'),
    fdt.DiffConstraint('v_A','x_A')]

model = fdt.DiagnosisModel(modelDef, name='Actuator')


# Emilias original model
modelDef={}
modelDef['type'] = 'Symbolic'

modelDef['x'] = ['q_A1', 'q_A2', 'p_A1', 'p_A2', 'V_A1', 'V_A2', 'q_leak_A',
              'vp_A', 'xp_A', 'ap_A', 'v_A', 'x_A', 'F_A', 'dV_A1', 'dV_A2',
              'dp_A1', 'dp_A2', 'c_leak_f1_A', 'Bp_f2_A', 'beta_f3_A',
              'K_f4_A', 'q_AT', 'q_AP', 'p_AT', 'p_AP', 'y', 'y_f5_A']
modelDef['z'] = ['s1','s2','s3','s4','s5','s6','s7']
modelDef['f'] = ['f1_A','f2_A','f3_A','f4_A','f5_A']
modelDef['parameters'] = ['A_A1','A_A2','V0_A1','V0_A2','Bp_A','beta_A',
                       'c_leak_A','K_A','m_A','k_f1_A','k_f2_A','k_f3_A',
                       'k_f4_A', 'k_f5_A', 'Cq_A', 'rho_A']

sym.var(modelDef['x'])
sym.var(modelDef['f'])
sym.var(modelDef['z'])
sym.var(modelDef['parameters'])

modelDef['rels'] = [
    -V_A1 + V0_A1 + A_A1*(-xp_A),
    -V_A2 + V0_A2 + A_A2*(xp_A),
    -q_leak_A + c_leak_f1_A*(p_A1 - p_A2),
    -q_A1 + q_leak_A + dV_A1 + dp_A1*V_A1/beta_f3_A,
    -q_A2 - q_leak_A + dV_A2 + dp_A2*V_A2/beta_f3_A, 
    -p_A2*A_A2 + p_A1*A_A1 + m_A*ap_A + Bp_f2_A*vp_A + K_f4_A*(xp_A - x_A),
    -F_A + K_f4_A*(xp_A - x_A),
    -q_AP + Cq_A*y_f5_A*sym.sqrt(2.0/rho_A) * sym.sqrt(sym.functions.Abs(p_AP - p_A1)) * sym.sign(p_AP - p_A1), 
    -q_A1 + Cq_A*y_f5_A*sym.sqrt(2.0/rho_A) * sym.sqrt(sym.functions.Abs(p_AP - p_A1)) * sym.sign(p_AP - p_A1), 
    -q_AT + -Cq_A*y_f5_A*sym.sqrt(2.0/rho_A) * sym.sqrt(sym.functions.Abs(p_A2 - p_AT)) * sym.sign(p_A2 - p_AT),
    -q_A2 + -Cq_A*y_f5_A*sym.sqrt(2.0/rho_A) * sym.sqrt(sym.functions.Abs(p_A2 - p_AT)) * sym.sign(p_A2 - p_AT), 
    -c_leak_f1_A + c_leak_A*(1 + k_f1_A*f1_A),
    -Bp_f2_A + Bp_A *(1 + k_f2_A*f2_A), 
    -beta_f3_A + beta_A*(1 + k_f3_A*f3_A),
    -K_f4_A + K_A*(1 + k_f4_A*f4_A), 
    -y_f5_A + y*(1 + k_f5_A*f5_A), 
    -s1 + q_AP,
    -s2 + p_AP,
    -s3 + q_AT,
    -s4 + p_AT,
    -s5 + x_A,
    -s6 + y,
    -s7 + F_A,
    fdt.DiffConstraint('dV_A1','V_A1'),
    fdt.DiffConstraint('dV_A2','V_A2'),
    fdt.DiffConstraint('dp_A1','p_A1'),
    fdt.DiffConstraint('dp_A2','p_A2'),
    fdt.DiffConstraint('vp_A','xp_A'),
    fdt.DiffConstraint('ap_A','vp_A'),
    fdt.DiffConstraint('v_A','x_A')]

model_emilia = fdt.DiagnosisModel(modelDef, name='Actuator (original)')