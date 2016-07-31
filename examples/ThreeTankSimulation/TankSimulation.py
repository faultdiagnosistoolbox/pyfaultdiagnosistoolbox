import numpy as np
from scipy.integrate import ode

class MyStruct:
    pass

def ramp(t,t0,t1):
    return (t>t0)*(t<=t1)*(t-t0)/(t1-t0)+(t>t1)

def fThreeTank(x,u,f,params):
    #Define parameter variables
    Rv1 = params['Rv1']
    Rv2 = params['Rv2']
    Rv3 = params['Rv3']
    CT1 = params['CT1']
    CT2 = params['CT2']
    CT3 = params['CT3']

    # State variables
    p1 = x[0]
    p2 = x[1]
    p3 = x[2]
  
    # Control input
    q0 = u;
  
    # Fault variables
    fRv1 = f[0]
    fRv2 = f[1]
    fRv3 = f[2]
    fCT1 = f[3]
    fCT2 = f[4]
    fCT3 = f[5]
  
    # Model equations  
    q1 = 1/Rv1*(p1-p2)  + fRv1
    q2 = 1/Rv2*(p2-p3)  + fRv2
    q3 = 1/Rv3*p3       + fRv3
  
    p1d = 1/CT1*(q0-q1) + fCT1
    p2d = 1/CT2*(q1-q2) + fCT2
    p3d = 1/CT3*(q2-q3) + fCT3
  
    return [p1d,p2d,p3d]

def ThreeTank_ss(params):
    # Define parameter variables
    Rv1 = params['Rv1']
    Rv2 = params['Rv2']
    Rv3 = params['Rv3']
    CT1 = params['CT1']
    CT2 = params['CT2']
    CT3 = params['CT3']

    # State-space matrices
    A = np.array([[-1/CT1/Rv1, 1/CT1/Rv1, 0],
         [1/CT2/Rv1, -1/CT2*(1/Rv1+1/Rv2), 1/CT2/Rv2],
         [0, 1/CT3/Rv2, -1/CT3*(1/Rv2+1/Rv3)]], dtype=np.float64)
    B = np.reshape(np.array([1/CT1,0,0], dtype=np.float64),(3,1))
    C = np.eye(3, dtype=np.float64)
    D = np.zeros((3,1), dtype=np.float64)
    
    return {'A':A, 'B':B, 'C':C, 'D':D}

def SimScenario(fi,fault, controlLaw, params, t, x0):
    Fi = np.zeros(6,dtype=np.float64)
    if fi>0:
        Fi[fi-1]=1
    
    Tend = np.max(t)
    N = len(t)
    dt = t[1]-t[0] # Ugly
    
    tsim = np.zeros(N)
    x = np.zeros((N,3))
    x[0,:] = x0

    f = lambda t,x: fThreeTank(x,controlLaw(t,x),fault(t)*Fi,params)
    r = ode(f).set_integrator('vode', method='bdf')
    r.set_initial_value(x0, 0)

    k = 0
    while r.successful() and k < N-1:
        tsim[k+1] = r.t + dt
        x[k+1,:] = r.integrate(r.t+dt)
        k = k+1

    sim = MyStruct();
    sim.t = tsim
    if fi==0:
        sim.f = 0*sim.t
    else:
        sim.f = fault(sim.t)
    sim.Fi = fi
    
    f = fault(sim.t).reshape((N,1))*([Fi])

    #q0 = np.transpose(-Lx.dot(np.transpose(x)) + Lr*ref(sim.t))
    q0 = controlLaw(sim.t, x)
    q1 = (1/params['Rv1']*(x[:,0]-x[:,1]) + f[:,0]).reshape((N,1))
    q2 = (1/params['Rv2']*(x[:,1]-x[:,2]) + f[:,1]).reshape((N,1))
    q3 = (1/params['Rv3']*x[:,2] + f[:,2]).reshape((N,1))

    sim.z0 = np.concatenate((x, q0, q1, q2, q3), axis=1)
    return sim

