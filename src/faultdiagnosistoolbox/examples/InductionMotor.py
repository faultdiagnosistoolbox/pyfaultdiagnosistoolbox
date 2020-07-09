import matplotlib.pyplot as plt
import faultdiagnosistoolbox as fdt
import sympy as sym

# Small induction motor example
# The model equations are taken from
#
# Aguilera, F., et al. "Current-sensor fault detection and isolation
# for induction-motor drives using a geometric approach."
# Control Engineering Practice 53 (2016): 35-46.

# %% Define model
model_def = {'type': 'Symbolic',
             'x': ['i_a', 'i_b', 'lambda_a', 'lambda_b', 'w',
                   'di_a', 'di_b', 'dlambda_a', 'dlambda_b', 'dw', 'q_a', 'q_b', 'Tl'],
             'f': ['f_a', 'f_b'], 'z': ['u_a', 'u_b', 'y1', 'y2', 'y3'],
             'parameters': ['a', 'b', 'c', 'd', 'L_M', 'k', 'c_f', 'c_t']}

sym.var(model_def['x'])
sym.var(model_def['f'])
sym.var(model_def['z'])
sym.var(model_def['parameters'])

model_def['rels'] = [
  -q_a + w*lambda_a,
  -q_b + w*lambda_b, 
  -di_a + -a*i_a + b*c*lambda_a + b*q_b + d*u_a,
  -di_b + -a*i_b + b*c*lambda_b + b*q_a + d*u_b,
  -dlambda_a + L_M*c*i_a - c*lambda_a-q_b, 
  -dlambda_b + L_M*c*i_b - c*lambda_b-q_a, 
  -dw + -k*c_f*w + k*c_t*(i_a*lambda_b - i_b*lambda_a) - k*Tl,
  fdt.DiffConstraint('di_a','i_a'),
  fdt.DiffConstraint('di_b','i_b'),
  fdt.DiffConstraint('dlambda_a','lambda_a'),
  fdt.DiffConstraint('dlambda_b','lambda_b'),
  fdt.DiffConstraint('dw','w'),
  -y1 + i_a + f_a,
  -y2 + i_b + f_b,
  -y3 + w]

model = fdt.DiagnosisModel(model_def, name ='Induction motor')

model.Lint()

# Plot model
plt.figure(10)
model.PlotModel()
ax = plt.gca()
ax.set_xticklabels(latex_var['x'] + latex_var['f'] + latex_var['z'], rotation=0)

# Isolability analysis
plt.figure(20)
model.IsolabilityAnalysis(plot=True, causality='der')
plt.title('Isolability with derivative causality')

plt.figure(21)
model.IsolabilityAnalysis(plot=True, causality='int')
plt.title('Isolability with integral causality')

plt.figure(22)
model.IsolabilityAnalysis(plot=True)
plt.title('Isolability with mixed causality')

plt.figure(23)
model.PlotDM(fault=True, eqclass=True)

# Find set of MSOS and MTES
msos = model.MSO()
mtes = model.MTES()
print(f"Found {len(msos)} MSO sets and {len(mtes)} MTES sets.")

# Check observability and low index for MSO sets
oi_mso = [model.IsObservable(m_i) for m_i in msos]
li_mso = [model.IsLowIndex(m_i) for m_i in msos]
print(f'Out of {len(msos)} MSO sets, {sum(oi_mso)} observable, {sum(li_mso)} low (structural) differential index')

# Check observability and low index for MTES sets
oi_mtes = [model.IsObservable(m_i) for m_i in mtes]
li_mtes = [model.IsLowIndex(m_i) for m_i in mtes]
print(f'Out of {len(mtes)} MTES sets, {sum(oi_mtes)} observable, {sum(li_mtes)} low (structural) differential index')

# Isolability analysis and FSM of set of MSO and MTES sets
FSM = model.FSM(msos)
plt.figure(30, clear=True)
plt.subplot(1, 2, 1)
_ = model.IsolabilityAnalysisArrs(msos, plot=True)
plt.title('Isolability matrix')

ax = plt.subplot(1, 2, 2)
plt.spy(FSM, markersize=10, marker='o')
ax.xaxis.tick_bottom()
ax.set_xticks(range(model.nf()))
ax.set_xticklabels(model.f)
ax.set_yticks(range(len(msos)))
ax.set_yticklabels([f"MSO {i + 1}" for i in range(len(msos))])
ax.set_xlabel('Fault')
plt.title('Fault Signature Matrix (MSO sets)')
plt.tight_layout()

FSM = model.FSM(mtes)
plt.figure(31, clear=True)
plt.subplot(1, 2, 1)
_ = model.IsolabilityAnalysisArrs(mtes, plot=True)
plt.title('Isolability matrix')

ax = plt.subplot(122)
plt.spy(FSM, markersize=10, marker='o')
ax.xaxis.tick_bottom()
ax.set_xticks(range(model.nf()))
ax.set_xticklabels(model.f)
ax.set_yticks(range(len(mtes)))
ax.set_yticklabels([f"MTES {i + 1}" for i in range(len(mtes))])
ax.set_xlabel('Fault')
ax.set_title('Fault Signature Matrix (MTES sets)')
plt.tight_layout()

# Generate C-code for residual generator based on first MTES
model.MSOCausalitySweep(mtes[0])
red_eq = mtes[0][10]
model.syme[red_eq]
M0 = [e for e in mtes[0] if e != red_eq]

Gamma = model.Matching(M0)
model.SeqResGen(Gamma, red_eq, 'ResGen', batch=True, language='C')

plt.show()

