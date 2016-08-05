import numpy as np
class VEPPar:
    pass

# Parameter file for VEP engine based on TSRT10 project FUDGE 2014 and
# summer project on Vehicular Systems 2015. 

# This file

TC_ENGINE = VEPPar()

TC_ENGINE.name = 'VEP_par'

## Ambient conditions
TC_ENGINE.ambient = VEPPar()
TC_ENGINE.ambient.p = 1.01300e+05 
TC_ENGINE.ambient.T = 2.93000e+02 

## Fuel properties
TC_ENGINE.fuel = VEPPar()
TC_ENGINE.fuel.q_lhv = 4.40000e+07 
TC_ENGINE.fuel.AFs = 1.40000e+01 

## Gas properties
TC_ENGINE.gasProp = VEPPar()
TC_ENGINE.gasProp.air = VEPPar()
TC_ENGINE.gasProp.throttle = VEPPar()
TC_ENGINE.gasProp.exh = VEPPar()
TC_ENGINE.gasProp.air.cp = 1.00520e+03 
TC_ENGINE.gasProp.air.gamma = 1.40000e+00 
TC_ENGINE.gasProp.air.R = 2.87200e+02 
TC_ENGINE.gasProp.throttle.cp = 1.00520e+03 
TC_ENGINE.gasProp.throttle.gamma = 2.00000e+00 
TC_ENGINE.gasProp.throttle.R = 2.87200e+02 
TC_ENGINE.gasProp.exh.cp = 1.25667e+03 
TC_ENGINE.gasProp.exh.gamma = 1.30000e+00 
TC_ENGINE.gasProp.exh.R = 2.90000e+02 

## Air filter
TC_ENGINE.airFilter = VEPPar()
TC_ENGINE.airFilter.H = 1.59598e+08 
TC_ENGINE.airFilter.p_lin = 1.60000e+03 

## Compressor
TC_ENGINE.Compressor = VEPPar()
TC_ENGINE.Compressor.massFlow = VEPPar()
TC_ENGINE.Compressor.massFlow.fit = VEPPar()
TC_ENGINE.Compressor.map = VEPPar()
TC_ENGINE.Compressor.efficiency = VEPPar()
TC_ENGINE.Compressor.efficiency.fit = VEPPar()
TC_ENGINE.Compressor.ExtendedEllipse = VEPPar()

TC_ENGINE.Compressor.D = 3.94000e-02 
TC_ENGINE.Compressor.massFlow.fit.k1 = 0.31972 #Old model: 5.67000e-01 
TC_ENGINE.Compressor.massFlow.fit.k2 = 18.1518 #Old model: 5.72000e-02 
TC_ENGINE.Compressor.map.pnom = 1.01300e+05 
TC_ENGINE.Compressor.map.Tnom = 2.93150e+02 
TC_ENGINE.Compressor.efficiency.fit.maxEfficiencyWcompCorr = 1.34730e-01 
TC_ENGINE.Compressor.efficiency.fit.maxEfficiencyPi_c = 2.09997e+00 
TC_ENGINE.Compressor.efficiency.fit.maxEfficiency = 7.36335e-01 
TC_ENGINE.Compressor.efficiency.fit.minEfficiency = 3.00000e-01 
TC_ENGINE.Compressor.efficiency.fit.Q_11 = 3.50283e+01 
TC_ENGINE.Compressor.efficiency.fit.Q_12 = -3.31011e+00 
TC_ENGINE.Compressor.efficiency.fit.Q_22 = 7.30446e-01 
TC_ENGINE.Compressor.MassFlowCoeff = [0.349075268491121, -0.186697794872490, 0.026600750345879] 
TC_ENGINE.Compressor.EnthalpyCoeff = 1e3*np.array([-0.000123211728717, -0.641987364359664, -0.000000024068063, 0.000277267255442, -0.355191900649545, -0.000002158117386, 0.027753777417142, 0.000001340414244, -0.002658906021361, 1.575702400427787])
TC_ENGINE.Compressor.fix_gain = 2.00000e-01 
TC_ENGINE.Compressor.ExtendedEllipse.Param=[0.477914879022786, 0.755152812979778, 1.90269061831234, 0.230576385816686, 0.802079743307562, 0.105421683936242, 0.442052718736814, 3.40214390269208, 0.722748681312551, 1.92853286250695, 0.741490025247638, 2.28900674245182, 1.82973314559249, 1.46991129576191, 0.885296164620165, 0.0668431406594969, 1.34044235115134, 2.91656174433059, 0.555838491751263] 
TC_ENGINE.Compressor.ExtendedEllipse.PI_max_map = 3.588530347611503 
TC_ENGINE.Compressor.ExtendedEllipse.Wc_max_map = 0.208537901446384 
TC_ENGINE.Compressor.ExtendedEllipse.Nc_max_map = 2.037183271576263e+05 
TC_ENGINE.Compressor.ExtendedEllipse.Tref = 298 
TC_ENGINE.Compressor.ExtendedEllipse.Pref =101300 

## Intercooler
TC_ENGINE.interCooler = VEPPar()
TC_ENGINE.interCooler.temperature = VEPPar()

TC_ENGINE.interCooler.epsilon = 8.00000e-01 
TC_ENGINE.interCooler.H = 3.46780e+08 
TC_ENGINE.interCooler.p_lin = 8.00000e+02 
TC_ENGINE.interCooler.temperature.fit = [ -5.64164e-01, 4.65003e-03, -2.99972e+00]

## Throttle
TC_ENGINE.throttle = VEPPar()
TC_ENGINE.throttle.gamma = 2.00000e+00 
TC_ENGINE.throttle.tau_th = 1.00000e-01 
TC_ENGINE.throttle.a_0 = 5.18849e-06 #7.9041e-05 #5.18849e-06 
TC_ENGINE.throttle.a_1 = -6.38593e-05 #-5.9427e-04 #-6.38593e-05 
TC_ENGINE.throttle.a_2 = 1.23068e-03 #0.0020 #1.23068e-03 
TC_ENGINE.throttle.PI_lin = 9.90000e-01 

## Engine geometry
TC_ENGINE.geometry = VEPPar()
TC_ENGINE.geometry.bore = 8.17000e-02 
TC_ENGINE.geometry.stroke = 9.32000e-02 
TC_ENGINE.geometry.connectingRodLength = 1.43800e-01 
TC_ENGINE.geometry.V_d = 1.95300e-03 
TC_ENGINE.geometry.r_c = 1.08000e+01 
TC_ENGINE.geometry.n_cyl = 4.00000e+00 
TC_ENGINE.geometry.n_r = 2.00000e+00 
TC_ENGINE.geometry.VehicleArea = 2.00000e+00 

## Ignition map and ignition efficiency
TC_ENGINE.Ignition = VEPPar()
TC_ENGINE.Ignition.eta_ign = VEPPar()

TC_ENGINE.Ignition.dThGrid = [
  0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -5.68908e-02, -5.88293e-01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01 
, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -1.36223e-02, -8.07441e-01, -4.22760e+00, -1.19878e+01, -1.15178e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01 
, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -3.25697e-01, -2.23772e+00, -4.94829e+00, -8.77392e+00, -1.20997e+01, -1.39057e+01, -1.55571e+01, -1.70952e+01, -1.56320e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01 
, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -4.43190e-01, -3.46369e+00, -5.41470e+00, -7.72814e+00, -1.02632e+01, -1.21249e+01, -1.39619e+01, -1.56550e+01, -1.73184e+01, -1.88709e+01, -1.97153e+01, -2.02336e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01 
, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -8.62000e-02, -1.66926e+00, -3.61399e+00, -5.67795e+00, -8.55561e+00, -1.08785e+01, -1.25976e+01, -1.45814e+01, -1.59004e+01, -1.73517e+01, -1.83820e+01, -1.92232e+01, -1.97073e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01 
, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -5.77147e-01, -1.31053e+00, -2.70281e+00, -6.27204e+00, -9.24663e+00, -1.10989e+01, -1.30218e+01, -1.44008e+01, -1.58917e+01, -1.69460e+01, -1.78110e+01, -1.85908e+01, -1.97769e+01, -2.02998e+01, -2.02998e+01 
, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -4.01395e-02, -3.40622e-01, -1.02429e+00, -3.93728e+00, -6.99003e+00, -9.04656e+00, -1.09619e+01, -1.26314e+01, -1.38775e+01, -1.53915e+01, -1.61935e+01, -1.80491e+01, -1.92540e+01, -2.02998e+01, -2.02998e+01 
, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -5.78947e-02, -1.00000e-01, -5.88576e-01, -2.30318e+00, -5.23570e+00, -7.32383e+00, -9.28923e+00, -1.12614e+01, -1.25946e+01, -1.38110e+01, -1.48958e+01, -1.75274e+01, -1.87311e+01, -2.02998e+01, -2.02998e+01 
, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -6.49270e-02, -3.16327e-01, -1.54684e+00, -4.00289e+00, -6.57448e+00, -8.72839e+00, -1.04598e+01, -1.17786e+01, -1.31321e+01, -1.44884e+01, -1.69247e+01, -1.82082e+01, -2.02998e+01, -2.02998e+01 
, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -2.20683e-02, -5.60272e-03, -2.02359e-02, -4.28842e-01, -1.62681e+00, -3.93839e+00, -6.24253e+00, -8.42553e+00, -1.01078e+01, -1.16122e+01, -1.29991e+01, -1.43545e+01, -1.62633e+01, -1.76852e+01, -2.02998e+01, -2.02998e+01 
, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -4.69048e-02, 0.00000e+00, -7.04294e-02, -6.17368e-01, -2.06572e+00, -4.99934e+00, -6.87321e+00, -8.67248e+00, -1.02249e+01, -1.17765e+01, -1.32987e+01, -1.44689e+01, -1.57474e+01, -1.71623e+01, -1.83985e+01, -2.02998e+01 
, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -1.48693e-02, -6.65072e-02, -4.48565e-03, -5.98984e-01, -1.92382e+00, -4.29071e+00, -6.69579e+00, -8.83888e+00, -1.05148e+01, -1.21332e+01, -1.36353e+01, -1.45995e+01, -1.53986e+01, -1.66394e+01, -1.78724e+01, -2.02998e+01 
, 0.00000e+00, 0.00000e+00, -2.63158e-02, -1.20012e-02, 0.00000e+00, -2.63158e-02, -2.63158e-02, -4.02081e-01, -1.57529e+00, -4.18889e+00, -6.35175e+00, -8.26603e+00, -9.99547e+00, -1.14291e+01, -1.31123e+01, -1.40674e+01, -1.52962e+01, -1.61165e+01, -1.73494e+01, -2.02998e+01 
, 0.00000e+00, 0.00000e+00, -5.37205e-02, -4.72238e-02, -6.48339e-02, -2.20907e-01, -1.12724e-01, -5.75257e-01, -2.06047e+00, -4.28772e+00, -6.13409e+00, -7.66580e+00, -9.29463e+00, -1.06895e+01, -1.20211e+01, -1.35413e+01, -1.52724e+01, -1.55646e+01, -1.68265e+01, -2.02998e+01 
, 0.00000e+00, 0.00000e+00, 0.00000e+00, -7.26257e-02, -2.63158e-01, -1.19607e-01, -1.38158e-01, -1.11423e+00, -2.35216e+00, -4.35439e+00, -5.87412e+00, -7.24702e+00, -8.50691e+00, -9.65950e+00, -1.10283e+01, -1.28531e+01, -1.52554e+01, -1.51297e+01, -1.61995e+01, -2.02998e+01 
, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -3.23182e-02, -1.27059e-01, -2.97511e-01, -1.04808e+00, -2.28456e+00, -3.86254e+00, -5.20684e+00, -6.47333e+00, -7.76912e+00, -8.95529e+00, -1.00356e+01, -1.23924e+01, -1.46763e+01, -1.50338e+01, -1.55102e+01, -2.02998e+01 
, 0.00000e+00, 0.00000e+00, 0.00000e+00, -1.37788e-02, -3.15789e-02, -1.28444e-01, -2.47742e-01, -1.36697e+00, -3.01412e+00, -4.11608e+00, -5.27333e+00, -6.30427e+00, -7.85186e+00, -8.85404e+00, -1.00083e+01, -1.18490e+01, -1.40971e+01, -1.50068e+01, -2.02998e+01, -2.02998e+01 
, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -5.07458e-02, -2.29208e-01, -8.88953e-01, -2.22389e+00, -3.43865e+00, -4.62398e+00, -5.79535e+00, -7.11673e+00, -8.10389e+00, -9.60700e+00, -1.05154e+01, -1.14129e+01, -1.35179e+01, -1.50087e+01, -2.02998e+01, -2.02998e+01 
, 0.00000e+00, 0.00000e+00, -3.64035e-02, -1.01195e-01, -5.14225e-01, -1.15920e+00, -2.17880e+00, -3.20975e+00, -3.84711e+00, -4.66579e+00, -5.77850e+00, -6.98499e+00, -7.90927e+00, -1.04817e+01, -1.11104e+01, -1.13558e+01, -1.29467e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01 
, 0.00000e+00, 0.00000e+00, -4.21987e-02, -2.77493e-01, -1.96612e+00, -2.82002e+00, -2.64107e+00, -3.37338e+00, -4.32166e+00, -5.01108e+00, -5.65042e+00, -6.22493e+00, -7.43493e+00, -8.79125e+00, -1.00875e+01, -1.13837e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01, -2.02998e+01] 

TC_ENGINE.Ignition.pimGrid = [ 2.16000e+01, 3.25158e+01, 4.34316e+01, 5.43474e+01, 6.52632e+01, 7.61789e+01, 8.70947e+01, 9.80105e+01, 1.08926e+02, 1.19842e+02, 1.30758e+02, 1.41674e+02, 1.52589e+02, 1.63505e+02, 1.74421e+02, 1.85337e+02, 1.96253e+02, 2.07168e+02, 2.18084e+02, 2.29000e+02] 
TC_ENGINE.Ignition.NGrid = [ 1.25000e+01, 1.71053e+01, 2.17105e+01, 2.63158e+01, 3.09211e+01, 3.55263e+01, 4.01316e+01, 4.47368e+01, 4.93421e+01, 5.39474e+01, 5.85526e+01, 6.31579e+01, 6.77632e+01, 7.23684e+01, 7.69737e+01, 8.15789e+01, 8.61842e+01, 9.07895e+01, 9.53947e+01, 1.00000e+02] 
TC_ENGINE.Ignition.eta_ign.c2 = 3.21700e+00 
TC_ENGINE.Ignition.eta_ign.c3 = 4.72200e+00 

## Volymetric efficiency
TC_ENGINE.air2cylinder = VEPPar()
TC_ENGINE.air2cylinder.fuelEvapPar = 1.10000e+02 
TC_ENGINE.air2cylinder.eta_vol_curve = [ 1.73380e-09, -8.39040e-07, 1.36250e-04, -7.40090e-03, 9.35810e-01]

## Torque model parameters
TC_ENGINE.torque = VEPPar()
TC_ENGINE.torque.bmep = VEPPar()
TC_ENGINE.torque.bmep.fit = VEPPar()

TC_ENGINE.torque.bmep.fit.C_tq1 = -2.76952e+05 
TC_ENGINE.torque.bmep.fit.C_tq2 = 1.25822e+01 
TC_ENGINE.torque.etaOtto = 5.65188e-01 
TC_ENGINE.torque.eta_ig_ch = 6.76945e-01 
TC_ENGINE.torque.Pi_bl = 1.85000e+00 
TC_ENGINE.torque.aux_dev_fric = 4.87111e-01 

## Control volumes
TC_ENGINE.controlVolumes = VEPPar()
TC_ENGINE.controlVolumes.airFilter = VEPPar()
TC_ENGINE.controlVolumes.Compressor = VEPPar()
TC_ENGINE.controlVolumes.interCooler = VEPPar()
TC_ENGINE.controlVolumes.intakeManifold = VEPPar()
TC_ENGINE.controlVolumes.exhaustManifold = VEPPar()
TC_ENGINE.controlVolumes.Turbine = VEPPar()
TC_ENGINE.controlVolumes.EGRCooler = VEPPar()

TC_ENGINE.controlVolumes.airFilter.V = 1.50850e-02 
TC_ENGINE.controlVolumes.airFilter.T_init = 2.93000e+02 
TC_ENGINE.controlVolumes.airFilter.p_init = 1.01000e+05 
TC_ENGINE.controlVolumes.Compressor.V = 6.00000e-03 
TC_ENGINE.controlVolumes.Compressor.T_init = 2.93000e+02 
TC_ENGINE.controlVolumes.Compressor.p_init = 1.01000e+05 
TC_ENGINE.controlVolumes.interCooler.V = 4.70000e-03 
TC_ENGINE.controlVolumes.interCooler.T_init = 2.93000e+02 
TC_ENGINE.controlVolumes.interCooler.p_init = 1.01000e+05 
TC_ENGINE.controlVolumes.intakeManifold.V = 1.80000e-03 
TC_ENGINE.controlVolumes.intakeManifold.T_init = 2.93000e+02 
TC_ENGINE.controlVolumes.intakeManifold.p_init = 1.01000e+05 
TC_ENGINE.controlVolumes.exhaustManifold.V = 1.59000e-03 
TC_ENGINE.controlVolumes.exhaustManifold.p_init = 1.05000e+05 
TC_ENGINE.controlVolumes.exhaustManifold.T_init = 9.00000e+02 
TC_ENGINE.controlVolumes.Turbine.V = 2.00000e-02 
TC_ENGINE.controlVolumes.Turbine.p_init = 1.02000e+05 
TC_ENGINE.controlVolumes.Turbine.T_init = 3.00000e+02 
TC_ENGINE.controlVolumes.EGRCooler.V = 1.00000e-03 
TC_ENGINE.controlVolumes.EGRCooler.T_init = 2.93000e+02 
TC_ENGINE.controlVolumes.EGRCooler.p_init = 1.01000e+05 

## Exhaust manifold
TC_ENGINE.exhMan = VEPPar()
TC_ENGINE.exhTemp = VEPPar()

TC_ENGINE.exhMan.eManDiam = 4.00000e-02 
TC_ENGINE.exhMan.eManLen = 1.70000e-01 
TC_ENGINE.exhMan.nrOfParallelPipes = 4.00000e+00 
TC_ENGINE.exhMan.h_ext = 9.64300e+01  #Used by old model
TC_ENGINE.exhMan.h_tot = 0.03646e+3 
TC_ENGINE.exhMan.exh_mu = 4.23228e-05 
TC_ENGINE.exhMan.exh_lambda = 6.84994e-02 
TC_ENGINE.exhTemp.zeroFlowTemp = 1.14683e+3 
TC_ENGINE.exhTemp.tempChangeWithFlow = 0.84336e+3 

## Turbo shaft
TC_ENGINE.Turbo = VEPPar()
TC_ENGINE.Turbo.damp_tc = 9.71870e-07 
TC_ENGINE.Turbo.inertia_tc = 1.25200e-04 
TC_ENGINE.Turbo.omega_init = 1.50000e+03 
TC_ENGINE.Turbo.omega_max = 2.28965e+04 
TC_ENGINE.Turbo.omega_min = 0 #1.00000e+03 

## Turbine
TC_ENGINE.Turbine = VEPPar()
TC_ENGINE.Turbine.massFlow = VEPPar()
TC_ENGINE.Turbine.massFlow.fit = VEPPar()
TC_ENGINE.Turbine.efficiency = VEPPar()
TC_ENGINE.Turbine.efficiency.fit = VEPPar()

TC_ENGINE.Turbine.D = 3.83000e-02 
TC_ENGINE.Turbine.massFlow.fit.k1 = 2.92000e-02 #0.015874410938786 #0.0126 #2.92000e-02 
TC_ENGINE.Turbine.massFlow.fit.k2 = 1.931710380902679 #2.3411 #2.34213e+00 
TC_ENGINE.Turbine.massFlow.fit.kvec = [-6.70830242304781e-06, 0.575543912119626, 0.000112194039663471, -0.447500246017316]  # for physical model with speeds
TC_ENGINE.Turbine.efficiency.fit.maxEfficiency = 6.37445e-01 
TC_ENGINE.Turbine.efficiency.fit.maxEfficiencyBSR = 4.12740e-01 
TC_ENGINE.Turbine.efficiency.fit.minEfficiency = 3.00000e-01 
TC_ENGINE.Turbine.efficiency.fit.etaT_curve = [ 5.70200e+00, -4.77600e+01, 1.57900e+02, -2.57400e+02, 2.06500e+02, -6.45500e+01]
TC_ENGINE.Turbine.efficiency.fit.PiT_limits = [ 1.00000e+00, 2.10000e+00] 
TC_ENGINE.Turbine.efficiency.fit.a_omega = [0, 9974928.83487578, 16575792.5370014, 32109852.6168100, 56054167.2758555, 76741191.7053602, 106900665.763039] 
TC_ENGINE.Turbine.efficiency.fit.b_omega = [0, -43392.0188412510, -103018.975128231, -262413.698439652, -533469.217487666, -777194.621957822, -1144332.45600859]
TC_ENGINE.Turbine.efficiency.fit.omega_breakpoints = [0, 6498.90800272609, 9788.46967029996, 13028.6036134574, 15822.8408393152, 18113.4807628027, 19969.5337025435]

## Exhaust system
TC_ENGINE.exhaustSystem = VEPPar()
TC_ENGINE.exhaustSystem.H = 2.1027e+08 #Old model, potentially wrong: 5.26092e+07 
TC_ENGINE.exhaustSystem.p_lin = 1.00000e+02 

## Bypass velve
TC_ENGINE.Bypass = VEPPar()
TC_ENGINE.Bypass.Amax = 2.80000e-03 
TC_ENGINE.Bypass.Cd = 9.00000e-01 
TC_ENGINE.Bypass.tau_bp = 1.00000e-01 
TC_ENGINE.Bypass.PI_lin = 9.90000e-01 

## Wastegate
TC_ENGINE.Wastegate = VEPPar()
TC_ENGINE.Wastegate.Amax = 4.90870e-04 #5.78e-04 #6.7126e-04 #4.90870e-04 
TC_ENGINE.Wastegate.Cd = 1.00000e+00 
TC_ENGINE.Wastegate.tau_wg = 1.00000e-01 
TC_ENGINE.Wastegate.PI_lin = 9.90000e-01 



lambda_cyl = 1
c_em = [0.26, 0.6, 0, 0]
Prandtl = 0.7

diag_par = {
    'H_af' : TC_ENGINE.airFilter.H, 
    'plin_af' : TC_ENGINE.airFilter.p_lin, 
    'H_exhaust' : TC_ENGINE.exhaustSystem.H, 
    'plin_exh' : TC_ENGINE.exhaustSystem.p_lin, 
    'H_intercooler' : TC_ENGINE.interCooler.H, 
    'plin_intercooler' : TC_ENGINE.interCooler.p_lin, 
    'PIli_th' : TC_ENGINE.throttle.PI_lin, 
    'gamma_air' : TC_ENGINE.gasProp.air.gamma, 
    'PIli_wg' : TC_ENGINE.Wastegate.PI_lin, 
    'gamma_exh' : TC_ENGINE.gasProp.exh.gamma, 
    'R_air' : TC_ENGINE.gasProp.air.R, 
    'cp_air' : TC_ENGINE.gasProp.air.cp, 
    'V_af' : TC_ENGINE.controlVolumes.airFilter.V, 
    'V_c' : TC_ENGINE.controlVolumes.Compressor.V, 
    'V_ic' : TC_ENGINE.controlVolumes.interCooler.V, 
    'V_im' : TC_ENGINE.controlVolumes.intakeManifold.V, 
    'R_exh' : TC_ENGINE.gasProp.exh.R, 
    'cp_exh' : TC_ENGINE.gasProp.exh.cp, 
    'A_0' : TC_ENGINE.throttle.a_0, 
    'A_1' : TC_ENGINE.throttle.a_1, 
    'A_2' : TC_ENGINE.throttle.a_2, 
    'lambda' : lambda_cyl, 
    'n_r' : TC_ENGINE.geometry.n_r, 
    'r_c' : TC_ENGINE.geometry.r_c, 
    'V_D' : TC_ENGINE.geometry.V_d, 
    'CEva_cool' : TC_ENGINE.air2cylinder.fuelEvapPar, 
    'AF_s' : TC_ENGINE.fuel.AFs, 
    'q_HV' : TC_ENGINE.fuel.q_lhv, 
    'C_tq1' : TC_ENGINE.torque.bmep.fit.C_tq1, 
    'C_tq2' : TC_ENGINE.torque.bmep.fit.C_tq2, 
    's' : TC_ENGINE.geometry.stroke, 
    'PI_bl' : TC_ENGINE.torque.Pi_bl, 
    'aux_dev_fric' : TC_ENGINE.torque.aux_dev_fric, 
    'B' : TC_ENGINE.geometry.bore, 
    'D_c' : TC_ENGINE.Compressor.D, 
    'K1' : TC_ENGINE.Compressor.massFlow.fit.k1, 
    'K2' : TC_ENGINE.Compressor.massFlow.fit.k2, 
    'fix_gain' : TC_ENGINE.Compressor.fix_gain, 
    'R_a' : TC_ENGINE.gasProp.air.R, 
    'cp_a' : TC_ENGINE.gasProp.throttle.cp, 
    'Q_c11' : TC_ENGINE.Compressor.efficiency.fit.Q_11, 
    'Q_c12' : TC_ENGINE.Compressor.efficiency.fit.Q_12, 
    'Q_c22' : TC_ENGINE.Compressor.efficiency.fit.Q_22, 
    'eta_cmax' : TC_ENGINE.Compressor.efficiency.fit.maxEfficiency, 
    'eta_cmin' : TC_ENGINE.Compressor.efficiency.fit.minEfficiency, 
    'W_ccorrmax' : TC_ENGINE.Compressor.efficiency.fit.maxEfficiencyWcompCorr, 
    'PI_c_max' : TC_ENGINE.Compressor.efficiency.fit.maxEfficiencyPi_c, 
    'T_std' : TC_ENGINE.Compressor.map.Tnom, 
    'k1' : TC_ENGINE.Turbine.massFlow.fit.k1, 
    'k2' : TC_ENGINE.Turbine.massFlow.fit.k2, 
    'cp_eg' : TC_ENGINE.gasProp.exh.cp, 
    'k3' : TC_ENGINE.Turbine.efficiency.fit.maxEfficiency, 
    'k4' : TC_ENGINE.Turbine.efficiency.fit.maxEfficiencyBSR, 
    'D_t' : TC_ENGINE.Turbine.D, 
    'gamma_eg' : TC_ENGINE.gasProp.exh.gamma, 
    'tau_wg' : TC_ENGINE.Wastegate.tau_wg, 
    'Amax' : TC_ENGINE.Wastegate.Amax, 
    'c_2' : TC_ENGINE.Ignition.eta_ign.c2, 
    'c_3' : TC_ENGINE.Ignition.eta_ign.c3, 
    'xi_fric_tc' : TC_ENGINE.Turbo.damp_tc, 
    'J_tc' : TC_ENGINE.Turbo.inertia_tc, 
    'V_em' :  TC_ENGINE.controlVolumes.exhaustManifold.V,  
    'V_t' :  TC_ENGINE.controlVolumes.Turbine.V, 
    'p_std' : TC_ENGINE.Compressor.map.pnom, 
    'eta_otto' : TC_ENGINE.torque.etaOtto, 
    'eta_ig_ch' : TC_ENGINE.torque.eta_ig_ch, 
    'xi_aux' : TC_ENGINE.torque.aux_dev_fric, 
    'a1' : TC_ENGINE.air2cylinder.eta_vol_curve[0], 
    'a2' : TC_ENGINE.air2cylinder.eta_vol_curve[1], 
    'a3' : TC_ENGINE.air2cylinder.eta_vol_curve[2], 
    'a4' : TC_ENGINE.air2cylinder.eta_vol_curve[3], 
    'a5' : TC_ENGINE.air2cylinder.eta_vol_curve[4], 
    'eta_tmin' : TC_ENGINE.Turbine.efficiency.fit.minEfficiency, 
    'eta_tmax' : TC_ENGINE.Turbine.efficiency.fit.maxEfficiency, 
    'Cd' : TC_ENGINE.Wastegate.Cd, 
    'PI_cmax' : TC_ENGINE.Compressor.efficiency.fit.maxEfficiencyPi_c, 
    'tau_th' : TC_ENGINE.throttle.tau_th,  
    'A_em' :  TC_ENGINE.exhMan.eManDiam*TC_ENGINE.exhMan.eManLen*np.pi,  # Area f?r exhaust manifold
    'K_t' : TC_ENGINE.exhTemp.tempChangeWithFlow, 
    'T0' : TC_ENGINE.exhTemp.zeroFlowTemp, 
    'h_ext' : TC_ENGINE.exhMan.h_ext, 
    'h_tot' : TC_ENGINE.exhMan.h_tot, 
    'c0_em' : c_em[0], 
    'c1_em' : c_em[1], 
    'c2_em' : c_em[2], 
    'D' : TC_ENGINE.exhMan.eManDiam, 
    'L' : TC_ENGINE.exhMan.eManLen, 
    'Pr' : Prandtl, 
    'my_exh' : TC_ENGINE.exhMan.exh_mu, 
    'lambda_exh' : TC_ENGINE.exhMan.exh_lambda,  
    'Cic_1' : TC_ENGINE.interCooler.temperature.fit[0], 
    'Cic_2' : TC_ENGINE.interCooler.temperature.fit[1], 
    'Cic_3' : TC_ENGINE.interCooler.temperature.fit[2], 
    'cv_exh' : TC_ENGINE.gasProp.exh.R/(TC_ENGINE.gasProp.exh.gamma-1), 
    'cv_air' : TC_ENGINE.gasProp.air.R/(TC_ENGINE.gasProp.air.gamma-1), 
    'omega_tc_max' : TC_ENGINE.Turbo.omega_max, 
    'omega_tc_min' : TC_ENGINE.Turbo.omega_min,
    'TOL' : 1e-16}
