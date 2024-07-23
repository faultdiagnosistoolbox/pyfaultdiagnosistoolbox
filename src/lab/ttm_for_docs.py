# %%
import faultdiagnosistoolbox as fdt
import sympy as sym
import matplotlib.pyplot as plt
from seaborn import despine

# %%
modelDef = {}
modelDef["type"] = "Symbolic"
modelDef["x"] = ["p1", "p2", "p3", "q0", "q1", "q2", "q3", "dp1", "dp2", "dp3"]
modelDef["f"] = ["fV1", "fV2", "fV3", "fT1", "fT2", "fT3"]
modelDef["z"] = ["y1", "y2", "y3"]
modelDef["parameters"] = ["Rv1", "Rv2", "Rv3", "CT1", "CT2", "CT3"]

p1, p2, p3, q0, q1, q2, q3, dp1, dp2, dp3 = sym.symbols(modelDef["x"])
fV1, fV2, fV3, fT1, fT2, fT3 = sym.symbols(modelDef["f"])
y1, y2, y3 = sym.symbols(modelDef["z"])
Rv1, Rv2, Rv3, CT1, CT2, CT3 = sym.symbols(modelDef["parameters"])

modelDef["rels"] = [
    -q1 + 1 / Rv1 * (p1 - p2) + fV1,
    -q2 + 1 / Rv2 * (p2 - p3) + fV2,
    -q3 + 1 / Rv3 * p3 + fV3,
    -dp1 + 1 / CT1 * (q0 - q1) + fT1,
    -dp2 + 1 / CT2 * (q1 - q2) + fT2,
    -dp3 + 1 / CT3 * (q2 - q3) + fT3,
    -y1 + p1,
    -y2 + q2,
    -y3 + q0,
    fdt.DiffConstraint("dp1", "p1"),
    fdt.DiffConstraint("dp2", "p2"),
    fdt.DiffConstraint("dp3", "p3"),
]

model = fdt.DiagnosisModel(modelDef, name="Three Tank System")

# %%
model.PlotModel()

# %% Basic usage
df, ndf = model.DetectabilityAnalysis()

# %% Isolability analysis
fig_ia, ax = plt.subplots(num=20, clear=True)
_ = model.IsolabilityAnalysis(ax=ax)
ax.set_xlabel("Diagnosed fault")
ax.set_ylabel("Injected fault")
despine(ax=ax)
# fig_ia.savefig("ia.png")

fig_ia_int, ax = plt.subplots(num=21, clear=True)
_ = model.IsolabilityAnalysis(ax=ax, causality="int")
ax.set_xlabel("Diagnosed fault")
ax.set_ylabel("Injected fault")
despine(ax=ax)
# fig_ia_int.savefig("ia_int.png")

fig_ia_der, ax = plt.subplots(num=22, clear=True)
_ = model.IsolabilityAnalysis(ax=ax, causality="der")
ax.set_xlabel("Diagnosed fault")
ax.set_ylabel("Injected fault")
despine(ax=ax)
# fig_ia_der.savefig("ia_der.png")


fig_ia, ax = plt.subplots(1, 3, num=20, clear=True, layout="compressed")
_ = model.IsolabilityAnalysis(ax=ax[0])
ax[0].set_xlabel("Diagnosed fault")
ax[0].set_ylabel("Injected fault")
ax[0].set_title("Causality: all")
despine(ax=ax[0])

_ = model.IsolabilityAnalysis(ax=ax[1], causality="int")
ax[1].set_xlabel("Diagnosed fault")
ax[1].set_ylabel("Injected fault")
ax[1].set_title("Causality: int")
despine(ax=ax[1])

_ = model.IsolabilityAnalysis(ax=ax[2], causality="der")
ax[2].set_xlabel("Diagnosed fault")
ax[2].set_ylabel("Injected fault")
ax[2].set_title("Causality: der")
despine(ax=ax[2])
# fig_ia.savefig("ia.png", bbox_inches='tight')


# %% Dulmage-Mendelsohn decomposition
fig_dmperm, ax = plt.subplots(num=23, clear=True)
_ = model.PlotDM(ax=ax, fault=True, eqclass=True)
despine(ax=ax)
# fig_dmperm.savefig("dmperm.png")


# %% Sensor placement

sens, _ = model.SensorPlacementIsolability()
model2 = model.copy()
model2.AddSensors(sens[0])

fig_sp, ax = plt.subplots(num=30, clear=True, layout="constrained")
_ = model2.IsolabilityAnalysis(ax=ax)
ax.set_xlabel("Diagnosed fault")
ax.set_ylabel("Injected fault")
ax.set_title("Causality: all")
despine(ax=ax)

# %% MSO:s and MTES:s
msos = model.MSO()
mtes = model.MTES()

li = [m for m in msos if model.IsLowIndex(m)]
oi = [m for m in msos if model.IsObservable(m)]

model.MSOCausalitySweep(msos[0])
model.MSOCausalitySweep(msos[1])

ts = model.TestSelection(li)
ts_msos = [li[ts_i] for ts_i in ts]
FSM = model.FSM(model.FSM(ts_msos))
model.IsolabilityAnalysisFSM(FSM)

# %% Causality of MSOs
mso = ts_msos[0]
print(mso)
model.MSOCausalitySweep(mso)

res = mso[1]  # y1 = p1
M0 = [ei for ei in mso if ei != res]
Gamma = model.Matching(M0)

model.SeqResGen(Gamma, res, "residual")
