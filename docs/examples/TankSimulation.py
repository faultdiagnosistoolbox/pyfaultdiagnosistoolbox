from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp


def ramp(t, t0, t1):
    return (t > t0) * (t <= t1) * (t - t0) / (t1 - t0) + (t > t1)


def fThreeTank(x, u, f, params):
    # Define parameter variables
    Rv1 = params["Rv1"]
    Rv2 = params["Rv2"]
    Rv3 = params["Rv3"]
    CT1 = params["CT1"]
    CT2 = params["CT2"]
    CT3 = params["CT3"]

    # State variables
    p1 = x[0]
    p2 = x[1]
    p3 = x[2]

    # Control input
    q0 = u

    # Fault variables
    fRv1 = f[0]
    fRv2 = f[1]
    fRv3 = f[2]
    fCT1 = f[3]
    fCT2 = f[4]
    fCT3 = f[5]

    # Model equations
    q1 = 1 / Rv1 * (p1 - p2) + fRv1
    q2 = 1 / Rv2 * (p2 - p3) + fRv2
    q3 = 1 / Rv3 * p3 + fRv3

    p1d = 1 / CT1 * (q0 - q1) + fCT1
    p2d = 1 / CT2 * (q1 - q2) + fCT2
    p3d = 1 / CT3 * (q2 - q3) + fCT3

    return [p1d, p2d, p3d]


def ThreeTank_ss(params):
    # Define parameter variables
    Rv1 = params["Rv1"]
    Rv2 = params["Rv2"]
    Rv3 = params["Rv3"]
    CT1 = params["CT1"]
    CT2 = params["CT2"]
    CT3 = params["CT3"]

    # State-space matrices
    A = np.array(
        [
            [-1 / CT1 / Rv1, 1 / CT1 / Rv1, 0],
            [1 / CT2 / Rv1, -1 / CT2 * (1 / Rv1 + 1 / Rv2), 1 / CT2 / Rv2],
            [0, 1 / CT3 / Rv2, -1 / CT3 * (1 / Rv2 + 1 / Rv3)],
        ]
    )
    B = np.reshape(np.array([1 / CT1, 0, 0], dtype=np.float64), (3, 1))
    C = np.eye(3)
    D = np.zeros((3, 1))

    return {"A": A, "B": B, "C": C, "D": D}


@dataclass
class SimResult:
    t: np.ndarray  # time
    f: np.ndarray  # fault signal
    w: np.ndarray  # all variables (p1, p2, p3, q0, q1, q2, q3)
    z: np.ndarray  # measurements (p1, q2, q0)


def SimScenario(fi, fault, controller, params, t, x0):
    Fi = np.zeros(6)
    if fi > 0:
        Fi[fi - 1] = 1

    res = solve_ivp(lambda t, x: fThreeTank(x, controller(t, x), fault(t) * Fi, params), [t[0], t[-1]], x0, t_eval=t)
    f = fault(res.t)[:, None] * Fi

    q0 = controller(res.t, res.y)
    q1 = 1 / params["Rv1"] * (res.y[0] - res.y[1]) + f[:, 0]
    q2 = 1 / params["Rv2"] * (res.y[1] - res.y[2]) + f[:, 1]
    q3 = 1 / params["Rv3"] * res.y[2] + f[:, 2]

    y = np.column_stack((res.y[0], q2, q0))

    return SimResult(res.t, fault(res.t), np.row_stack((res.y, q0, q1, q2, q3)), y)
