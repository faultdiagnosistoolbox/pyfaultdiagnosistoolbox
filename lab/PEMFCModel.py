import faultdiagnosistoolbox as fdt
import sympy as sym

# %% Define model object
modelDef={}
modelDef['type'] = 'Symbolic'

modelDef['x'] = ['AHamb' , 'AHcmp' , 'cpair' , 'd_p' , 'dpcritic' , 'EAct' , 'EDif' , 
        'ENernst' , 'EOhm' , 'eta_cmp' , 'Idens' , 'Ilim' , 'lambdaan' , 'lambdaca' , 
        'lambdamem' , 'Mair' , 'Man' , 'Mca' , 'mH2an' , 'mH2sman' , 
        'mH2Oan' , 'mH2Oca' , 'mH2Osman' , 'mH2O_smandes' , 'mH2Osmca',
        'mH2O_smcades' , 'mN2ca' , 'mN2smca' , 'mO2ca' , 'mO2smca' , 'Msman',
        'Msmca' , 'pairamb' , 'pan' , 'pca' , 'Pcmp' , 'pH2an' , 'pH2sman',
        'pH2Oamb' , 'pH2Oan' , 'pH2Oca' , 'pH2Osman' ,'pH2Osmca' , 'pN2ca',
        'pN2smca' , 'pO2ca' , 'pO2smca' , 'psatamb' , 'psatfc' , 'psatsman' ,
        'psatsmca' , 'psman' , 'psmca' , 'RHan' , 'RHca' , 'sigmamem' ,
        'Tcain' , 'Tcmp' , 'Tfc' , 'Vfc' , 'Waircmp' , 'Wanin' , 'Wanout' ,
        'Wcain' , 'Wcaout' , 'Wcmp' , 'WH2anin' , 'WH2anout' , 'WH2smanout' ,
        'WH2Oanin' , 'WH2ancons' , 'WH2Oanout' , 'WH2Ocagen' , 'WH2Ocain' ,
        'WH2Ocaout' , 'WH2Ocmp' , 'WH2Omem' , 'WH2Osmaninj' , 'WH2Osmanout' ,
        'WH2O_smcain' , 'WH2Osmcainj' , 'WH2Osmcaout' , 'WN2cain' ,
        'WN2caout' , 'WN2smcain' , 'WN2smcaout' , 'WO2cacons' , 'WO2cain' ,
        'WO2caout' , 'WO2smcain' , 'WO2smcaout' , 'Wsmanin' , 'Wsmanout' ,
        'Wsmcaout' , 'xH2an' , 'xH2sman' , 'xH2Oan' , 'xH2Oca' , 'xH2Osman' ,
        'xH2Osmca' , 'xN2ca' , 'xN2smca' , 'xO2ca' , 'xO2smca' , 'I0' ,
        'rhoair', 'dmO2ca', 'dmN2ca', 'dmH2Oca', 'dmH2an', 'dmH2Oan', 'dTfc',
        'dmO2smca', 'dmN2smca', 'dmH2Osmca', 'CDcain', 'CDcaout', 'dmH2sman',
        'dmH2Osman']


modelDef['f'] = ['f_cmp', 'f_inj' , 'f_leak' , 'f_vin', 'f_ohm' ,
                 'f_ecsa', 'f_vout']
modelDef['z'] = [ 'betacmp' , 'ncmp', 'Tanin', 'RHsmandes' , 'RHsmcades',
                  'Load'  , 'CDcain0', 'CDanin' ,'CDcaout0' , 'CDanout' ,
                  'CDsmanin']

modelDef['parameters'] = ['Aanin' , 'Acain' , 'Aanout' , 'Acaout' , 'Afc' ,
        'Asmanin' , 'cpH2' , 'cpH2O', 'cpN2' , 'cpO2' , 'DO2N2' , 'epsca' ,
        'etaem' , 'F' , 'I_0Pt' , 'kk' , 'Kfc' , 'MH2' , 'MH2O' , 'MN2' ,
        'MO2' , 'omegatil' , 'pamb' , 'pH2tank' , 'Qdot' , 'R' , 'RHamb' ,
        'Tamb' , 'tauinj' , 'tca' , 'tmem' , 'Van' , 'Vca' , 'Vm' , 'Vsman' ,
        'Vsmca' , 'xO2amb' , 'xN2amb' , 'lambdamem_map' , 'lambdamem_map_x1',
        'lambdamem_map_x2', 'lambdamem_map_x3', 'lambdamem_map_x4' ,
        'WH2Omem_map' , 'WH2Omem_map_x1', 'WH2Omem_map_x2', 'WH2Omem_map_x3',
        'WH2Omem_map_x4' , 'Wcmp_map' , 'Wcmp_map_x1', 'Wcmp_map_x2',
        'eta_map' , 'eta_map_x1', 'eta_map_x2', 'ECSA' , 'Hf_H2O' ]

sym.var(modelDef['x'])
sym.var(modelDef['f'])
sym.var(modelDef['z'])
sym.var(modelDef['parameters'])

lambdamem_fun = sym.Function('lambdamem_fun')
WH2Omem_fun = sym.Function('WH2Omem_fun')
Wcmp_fun = sym.Function('Wcmp_fun')
eta_fun = sym.Function('eta_fun')
ifstatement_fun = sym.Function('ifstatement_fun')

ext_funs = {
        'lambdamem_fun' : 'lambdamem_fun',
        'WH2Omem_fun' : 'WH2Omem_fun',
        'Wcmp_fun' : 'Wcmp_fun',
        'eta_fun' : 'eta_fun',
        'ifstatement_fun' : 'ifstatement_fun'}

modelDef['rels'] = [
  -Vfc + ENernst-EAct-EOhm-EDif,
  -ENernst + 1.23-0.9e-3*(Tfc-298)+Tfc*R/1000/2/F*sym.log(pH2an*sym.sqrt(pO2ca)),
  -EAct + Tfc*R/1000/F*sym.log(Idens/I0),                    
  -EOhm + tmem*Idens/sigmamem,
  -EDif + omegatil*Idens*Tfc*sym.log(Ilim/(Ilim-Idens)),
  -sigmamem + ((0.005139*lambdamem-0.00326)*sym.exp(1268*(1/303-1/Tfc)))*(1-f_ohm),
  -Ilim + -2*F*DO2N2*epsca**1.5/Vm/tca*(Tfc/273)**0.823*sym.log(1-xO2ca),
  -xO2ca + pO2ca/pca,
  -lambdamem + lambdamem_fun(lambdamem_map,lambdamem_map_x1,lambdamem_map_x2,lambdamem_map_x3,lambdamem_map_x4,lambdaan,lambdaca,Idens,d_p),
  -WH2Omem + WH2Omem_fun(WH2Omem_map,WH2Omem_map_x1,WH2Omem_map_x2,WH2Omem_map_x3,WH2Omem_map_x4,lambdaan,lambdaca,Idens,d_p)*MH2O/1000,
  -d_p + pca-pan,
  -lambdaca + ifstatement_fun(-(RHca-1),0.3+10.8*RHca-16*RHca**2+14.1*RHca**3,-30.41+61.98*RHca-25.96*RHca**2+3.7*RHca**3),
  -lambdaan + ifstatement_fun(-(RHan-1),0.3+10.8*RHan-16*RHan**2+14.1*RHan**3,-30.41+61.98*RHan-25.96*RHan**2+3.7*RHan**3),
  -RHca + pH2Oca/psatfc,
  -RHan + pH2Oan/psatfc,
  -psatfc + 10**(-1.69e-10*Tfc**4+3.85e-7*Tfc**3-3.39e-4*Tfc**2+0.143*Tfc-20.92)*1e-2,
  -dmO2ca + WO2cain-WO2caout-WO2cacons,
  -dmN2ca + WN2cain-WN2caout,
  -dmH2Oca + WH2Ocain-WH2Ocaout+WH2Ocagen+WH2Omem,
  -WO2cacons + Idens*Afc/4/F*MO2/1000,
  -WH2Ocagen + Idens*Afc/2/F*MH2O/1000,
  -pO2ca + mO2ca*R*Tfc/MO2/Vca*1e-5,
  -pN2ca + mN2ca*R*Tfc/MN2/Vca*1e-5,
  -pH2Oca + mH2Oca*R*Tfc/MH2O/Vca*1e-5,
  -pca + pO2ca+pN2ca+pH2Oca,
  -xN2ca + pN2ca/pca,
  -xH2Oca + pH2Oca/pca,
  -Mca + MO2*xO2ca+MN2*xN2ca+MH2O*xH2Oca,
  -dmH2an + WH2anin-WH2anout-WH2ancons,
  -dmH2Oan + WH2Oanin-WH2Oanout-WH2Omem,
  -WH2ancons + Idens*Afc/2/F*MH2/1000,
  -pH2an + mH2an*R*Tfc/MH2/Van*1e-5,
  -pH2Oan + mH2Oan*R*Tfc/MH2O/Van*1e-5,
  -pan + pH2an+pH2Oan,
  -xH2an + pH2an/pan,
  -xH2Oan + pH2Oan/pan,
  -Man + MH2*xH2an+MH2O*xH2Oan,
  -dTfc + ((WH2Oanin+WH2Ocain-WH2Oanout-WH2Ocaout)*Hf_H2O+(WH2anin*cpH2+WH2Oanin*cpH2O)*(Tanin-298.15)+(WO2cain*cpO2+WH2Ocain*cpH2O+WN2cain*cpN2)*(Tcain-298.15)-(WH2anout*cpH2+WH2Oanout*cpH2O+WO2caout*cpO2+WH2Ocaout*cpH2O+WN2caout*cpN2)*(Tfc-298.15)-Vfc*Load-Qdot)/Kfc,
  -dpcritic + (2/(kk+1))**(kk/(kk-1)),
  -Wsmcaout + ifstatement_fun(pca/psmca-dpcritic,
    CDcain*Acain/sym.sqrt(R/Msmca*Tcain)*psmca*1e5*(pca/psmca)**(1/kk)*sym.sqrt(2*kk/(kk-1)*(1-(pca/psmca)**((kk-1)/kk))),  
    CDcain*Acain/sym.sqrt(R/Msmca*Tcain)*psmca*1e5*sym.sqrt(kk)*(2/(kk+1))**((kk+1)/2/(kk-1))), 
  -WO2cain + Wcain*xO2smca*MO2/Msmca,
  -WN2cain + Wcain*xN2smca*MN2/Msmca,
  -WH2Ocain + Wcain*xH2Osmca*MH2O/Msmca,
  -Wsmanout + ifstatement_fun(pan/psman-dpcritic,
    CDanin*Aanin/sym.sqrt(R/Msman*Tanin)*psman*1e5*(pan/psman)**(1/kk)*sym.sqrt(2*kk/(kk-1)*(1-(pan/psman)**((kk-1)/kk))), 
    CDanin*Aanin/sym.sqrt(R/Msman*Tanin)*psman*1e5*sym.sqrt(kk)*(2/(kk+1))**((kk+1)/2/(kk-1))), 
  -WH2anin + Wanin*xH2sman*MH2/Msman,
  -WH2Oanin + Wanin*xH2Osman*MH2O/Msman,
  -Wcaout + ifstatement_fun(pamb/pca-dpcritic,
    CDcaout*Acaout/sym.sqrt(R/Mca*Tfc)*pca*1e5*(pamb/pca)**(1/kk)*sym.sqrt(2*kk/(kk-1)*(1-(pamb/pca)**((kk-1)/kk))),  
    CDcaout*Acaout/sym.sqrt(R/Mca*Tfc)*pca*1e5*sym.sqrt(kk)*(2/(kk+1))**((kk+1)/2/(kk-1))), 
  -WO2caout + Wcaout*xO2ca*MO2/Mca,
  -WN2caout + Wcaout*xN2ca*MN2/Mca,
  -WH2Ocaout + Wcaout*xH2Oca*MH2O/Mca,
  -Wanout + ifstatement_fun(pamb/pan-dpcritic,
    CDanout*Aanout/sym.sqrt(R/Man*Tfc)*pan*1e5*(pamb/pan)**(1/kk)*sym.sqrt(2*kk/(kk-1)*(1-(pamb/pan)**((kk-1)/kk))), 
    CDanout*Aanout/sym.sqrt(R/Man*Tfc)*pan*1e5*sym.sqrt(kk)*(2/(kk+1))**((kk+1)/2/(kk-1))), 
  -WH2anout + Wanout*xH2an*MH2/Man,
  -WH2Oanout + Wanout*xH2Oan*MH2O/Man,
  -WO2smcain + Waircmp*xO2amb*MO2/Mair,
  -WN2smcain + Waircmp*xN2amb*MN2/Mair,
  -WH2O_smcain + WH2Ocmp,
  -dmO2smca + WO2smcain-WO2smcaout,
  -dmN2smca + WN2smcain-WN2smcaout,
  -dmH2Osmca + WH2O_smcain-WH2Osmcaout+WH2Osmcainj,
  -pO2smca + mO2smca*R*Tcmp/MO2/Vsmca*1e-5,
  -pN2smca + mN2smca*R*Tcmp/MN2/Vsmca*1e-5,
  -pH2Osmca + mH2Osmca*R*Tcmp/MH2O/Vsmca*1e-5,
  -psmca + pO2smca+pN2smca+pH2Osmca,
  -xO2smca + pO2smca/psmca,
  -xN2smca + pN2smca/psmca,
  -xH2Osmca + pH2Osmca/psmca,
  -Msmca + xO2smca*MO2+xN2smca*MN2+xH2Osmca*MH2O,
  -Wcain + Wsmcaout*(1-f_leak),
  -WO2smcaout + Wsmcaout*xO2smca*MO2/Msmca,
  -WN2smcaout + Wsmcaout*xN2smca*MN2/Msmca,
  -WH2Osmcaout + Wsmcaout*xH2Osmca*MH2O/Msmca,
  -psatsmca + 10**(-1.69e-10*Tcmp**4+3.85e-7*Tcmp**3-3.39e-4*Tcmp**2+0.143*Tcmp-20.92)*1e-2,
  -mH2O_smcades + psatsmca*1e5*(RHsmcades-pH2Osmca/psatsmca)*Vsmca*MH2O/R/Tcmp,
  -WH2Osmcainj + mH2O_smcades/tauinj,
  -Wsmanin + ifstatement_fun(psman/pH2tank-dpcritic,
    CDsmanin*Asmanin/sym.sqrt(R/Msman*Tanin)*pH2tank*1e5*(psman/pH2tank)**(1/kk)*sym.sqrt(2*kk/(kk-1)*(1-(psman/pH2tank)**((kk-1)/kk))), 
    CDsmanin*Asmanin/sym.sqrt(R/Msman*Tanin)*pH2tank*1e5*sym.sqrt(kk)*(2/(kk+1))**((kk+1)/2/(kk-1))), 
  -dmH2sman + Wsmanin-WH2smanout,
  -dmH2Osman + -WH2Osmanout+WH2Osmaninj,
  -pH2sman + mH2sman*R*Tanin/MH2/Vsman*1e-5,
  -pH2Osman + mH2Osman*R*Tanin/MH2O/Vsman*1e-5,
  -psman + pH2sman+pH2Osman,
  -xH2sman + pH2sman/psman,
  -xH2Osman + pH2Osman/psman,
  -Msman + xH2sman*MH2+xH2Osman*MH2O,
  -Wanin + Wsmanout,
  -WH2smanout + Wsmanout*xH2sman*MH2/Msman,
  -WH2Osmanout + Wsmanout*xH2Osman*MH2O/Msman,
  -psatsman + 10**(-1.69e-10*Tanin**4+3.85e-7*Tanin**3-3.39e-4*Tanin**2+0.143*Tanin-20.92)*1e-2,
  -mH2O_smandes + psatsman*1e5*(RHsmandes-pH2Osman/psatsman)*Vsman*MH2O/R/Tanin,
  -WH2Osmaninj + mH2O_smandes/tauinj*(1-f_inj),
  -psatamb + 10**(-1.69e-10*Tamb**4+3.85e-7*Tamb**3-3.39e-4*Tamb**2+0.143*Tamb-20.92)*1e-2,
  -pH2Oamb + RHamb*psatamb,
  -pairamb + pamb-pH2Oamb,
  -Mair + xO2amb*MO2+xN2amb*MN2,
  -cpair + xO2amb*cpO2+xN2amb*cpN2,
  -AHamb + pH2Oamb/pairamb,
  -Wcmp + Wcmp_fun(Wcmp_map,Wcmp_map_x1,Wcmp_map_x2,betacmp,ncmp)*rhoair/60/1000,
  -eta_cmp + eta_fun(eta_map,eta_map_x1,eta_map_x2,betacmp,ncmp)*(1-f_cmp),
  -Pcmp + Wcmp*Tamb*cpair/eta_cmp/etaem*(betacmp**((kk-1)/kk)-1),
  -Tcmp + Tamb*(1+1/eta_cmp*(betacmp**((kk-1)/kk)-1)),
  -Tcain + Tcmp,
  -AHcmp + AHamb,
  -Waircmp + Wcmp/(1+AHcmp),
  -WH2Ocmp + Wcmp-Waircmp,
  -I0 + I_0Pt*ECSA*(1-f_ecsa),
  -Idens + Load/Afc,
  -CDcain + CDcain0*(1-f_vin),
  -CDcaout + CDcaout0*(1-f_vout),
  -rhoair + pairamb*betacmp*1e5/Tcmp/R*Mair,
  fdt.DiffConstraint('dmO2ca','mO2ca'),
  fdt.DiffConstraint('dmN2ca','mN2ca'),
  fdt.DiffConstraint('dmH2Oca','mH2Oca'),
  fdt.DiffConstraint('dmH2an','mH2an'),
  fdt.DiffConstraint('dmH2Oan','mH2Oan'),
  fdt.DiffConstraint('dTfc','Tfc'),
  fdt.DiffConstraint('dmO2smca','mO2smca'),
  fdt.DiffConstraint('dmN2smca','mN2smca'),
  fdt.DiffConstraint('dmH2Osmca','mH2Osmca'),
  fdt.DiffConstraint('dmH2sman','mH2sman'),
  fdt.DiffConstraint('dmH2Osman','mH2Osman')]

# %%
model = fdt.DiagnosisModel( modelDef, name='PEM Fuel Cell' )

