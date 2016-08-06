#include <math.h>
#include "extmodelfuns.h"

#define max(a,b) a>b ? a : b;
#define min(a,b) a>b ? b : a;

#define TOL 1e-16

double robsqrt(double x) {return sqrt(fabs(x));};

double max_fun(double choice_a, double choice_b) {return max(choice_a, choice_b);};
double min_fun(double choice_a, double choice_b) {return min(choice_a, choice_b);};

double PSI( double PI, double gamma)
{
  return robsqrt( ( 2*gamma/(gamma-1) )*( pow(PI,2/gamma) - pow(PI,(gamma+1)/gamma) ) );
}

double
W_af_fun( double p_amb, double p_af, double plin_af, double H_af, double T_amb )
{
  double W_af;
  if( p_amb - p_af - plin_af > 0 ) {
    W_af = robsqrt( p_amb/(H_af*T_amb) )*robsqrt(p_amb-p_af);
  } else {
    W_af = robsqrt( p_amb/(H_af*T_amb) )*(p_amb-p_af)/robsqrt(plin_af);
  }
  return W_af;
}

double PSI_th_fun(double p_im, double p_ic, double gamma_air, double PIli_th)
{
  double PI = p_im/p_ic;
  double PI_th;
  double PSI_th;

  if( PI < PIli_th ) {
    PI_th = max( PI, pow( 2/(gamma_air+1),gamma_air/(gamma_air-1) ) );
    PSI_th = robsqrt( ( 2*gamma_air/(gamma_air-1) )*( pow(PI_th,2/gamma_air)-pow(PI_th,(gamma_air+1)/gamma_air) ) );
  } else {
    PI_th = max( PIli_th, pow( 2/(gamma_air+1),gamma_air/(gamma_air-1) ) );
    PSI_th = robsqrt( ( 2*gamma_air/(gamma_air-1) )*( pow( PI_th,2/gamma_air )-pow( PI_th,(gamma_air+1)/gamma_air ) ) )/(1-PIli_th)*(1 - PI);
  }
  return PSI_th;
}

double PI_c_fun(double p_c, double p_af)
{
  
  if( p_c/p_af>1 ) {
    return p_c/p_af;
  } else {
    return 1 - TOL;
  }
}
    
double PHI_model_fun(double K1, double K2, double PSI_c)
{
  double PHI_model;
  if ( (1-K1*PSI_c*PSI_c)/K2 > 0 ) {
    PHI_model = robsqrt((1-K1*PSI_c*PSI_c)/K2);
  } else {
    PHI_model = 0;
  }
return PHI_model;
}

double W_c_fun(double PI_cnolim, double p_af, double U_c, double D_c, double T_af, double R_air, double PHI_model, double fix_gain)
{
  double Wc;
  if( PI_cnolim - 1 >= 0 ) {
    Wc =  ( p_af*U_c*M_PI*D_c*D_c/(T_af*2*R_air) )*PHI_model;
  } else {
    Wc = ( p_af*U_c*M_PI*D_c*D_c/(T_af*2*R_air) )*PHI_model*(1-(2-2*PI_cnolim)) + fix_gain*(2-2*PI_cnolim);
  }
  return Wc;
}


double eta_c_fun(double eta_cmax, double eta_cmin, double W_ccorr, double W_ccorrmax, double PI_c, double PI_cmax, double Q_c11, double Q_c12, double Q_c22) 
{
  double X_c1, X_c2;
  double eta_c_comparison;

  X_c1 = W_ccorr-W_ccorrmax; 

  if( PI_c > 1 ) {
    X_c2 = robsqrt(PI_c-1)+1-PI_cmax;
  } else {
    X_c2 = 1-PI_cmax;
  }
  
  eta_c_comparison = eta_cmax - ( X_c1*X_c1*Q_c11 + 2*X_c1*X_c2*Q_c12 + X_c2*X_c2*Q_c22 );
  
  if( eta_c_comparison >= 1 ) {
    return 1 - TOL;
  } else if( eta_cmin < eta_c_comparison && eta_c_comparison < 1 ) {
    return eta_c_comparison;
  } else {
    return eta_cmin;
  } 
}

double W_ic_fun(double p_c, double p_ic, double plin_intercooler, double H_intercooler, double T_c)
{
  if( p_c-p_ic-plin_intercooler > 0 ) {
    return robsqrt(p_c/(H_intercooler*T_c))*robsqrt(p_c-p_ic);
  } else {
    return robsqrt(p_c/(H_intercooler*T_c))*(p_c-p_ic)/robsqrt(plin_intercooler);
  }
}

double Tflow_wg_fun(double p_t, double p_em, double T_em, double T_t)
{
  if( p_em < p_t ) {
    return T_t;
  } else {
    return T_em;
  }
}

double W_es_fun(double p_amb, double p_t, double plin_exh, double H_exhaust, double T_t)
{
  if( plin_exh > p_t-p_amb ) {
    return robsqrt( p_t/(H_exhaust*T_t) )*(p_t-p_amb)/robsqrt(plin_exh);
  } else {
    return robsqrt( p_t/(H_exhaust*T_t) )*robsqrt(p_t-p_amb);
  }
}

double PI_wg_fun( double p_t, double p_em, double gamma_air )
{
  double x = pow( 2/(gamma_air+1) , gamma_air/(gamma_air-1) );
  if( -p_t/p_em + x > 0 ) {
    return x;
  } else {
    return p_t/p_em;
  }
} 

double PSIli_wg_fun( double PI_wg, double PIli_wg, double gamma_exh)
{
  double PI;
  double x = pow( 2/(gamma_exh+1), gamma_exh/(gamma_exh - 1) );
  if( PIli_wg - PI_wg > 0 ) {
    PI = max(PI_wg, x);	
    return PSI(PI,gamma_exh);
  } else {
    PI = max(PIli_wg, x);	
    return PSI(PI,gamma_exh)*(1-PI_wg)/(1-PIli_wg);
  }
}

double PI_t_fun( double p_t, double p_em )
{
  double x = p_t/p_em;
  if( x > 1 ) {
    return 1 - TOL;
  } else {
    return x;
  }
} 

double W_t_fun(double p_em, double k1, double PI_t, double k2, double T_em)
{
  double x = max( TOL , 1-pow(PI_t, k2) );
  return ( p_em*0.001*k1*robsqrt(x) )/robsqrt(T_em);
}

double eta_t_fun( double k3, double BSR, double k4, double eta_min )
{
  double eta_t = k3*(1-((BSR - k4)/k4)*((BSR - k4)/k4));
  
  if (eta_t > k3) {
    return k3;
  } else if (eta_t < eta_min) {
    return eta_min;
  } else { 
    return eta_t;
  } 
}

double Tq_t_fun(double gamma_eg, double cp_eg, double eta_t, double W_t, double T_em, double PI_t, double omega_tc)
{
  double x = max(0, 1-pow(PI_t, (gamma_eg-1)/gamma_eg) );
  return cp_eg*eta_t*W_t*T_em*x/omega_tc;
}


double igngrid(double omega_e, double p_im) 
{
  double dThGrid[20][20] = {
    {   0.00000E00,   0.00000E00,  0.00000E00,   0.00000E00, -5.68908E-02, 
	-5.88293E-01,  -2.02998E01, -2.02998E01,  -2.02998E01,  -2.02998E01, 
	-2.02998E01,  -2.02998E01, -2.02998E01,  -2.02998E01,  -2.02998E01, 
	-2.02998E01,  -2.02998E01, -2.02998E01,  -2.02998E01,  -2.02998E01 }, 
    {   0.00000E00,   0.00000E00,  0.00000E00,   0.00000E00, -1.36223E-02, 
	-8.07441E-01,  -4.22760E00, -1.19878E01,  -1.15178E01,  -2.02998E01,
	-2.02998E01,  -2.02998E01, -2.02998E01,  -2.02998E01,  -2.02998E01, 
	-2.02998E01,  -2.02998E01, -2.02998E01,  -2.02998E01,  -2.02998E01 }, 
    {   0.00000E00,   0.00000E00,  0.00000E00,   0.00000E00, -3.25697E-01,
	-2.23772E00,  -4.94829E00, -8.77392E00,  -1.20997E01,  -1.39057E01, 
	-1.55571E01,  -1.70952E01, -1.56320E01,  -2.02998E01,  -2.02998E01, 
	-2.02998E01,  -2.02998E01, -2.02998E01,  -2.02998E01,  -2.02998E01 },
    {   0.00000E00,   0.00000E00,  0.00000E00,   0.00000E00, -4.43190E-01, 
	-3.46369E00,  -5.41470E00, -7.72814E00,  -1.02632E01,  -1.21249E01,
	-1.39619E01,  -1.56550E01, -1.73184E01,  -1.88709E01,  -1.97153E01, 
	-2.02336E01,  -2.02998E01, -2.02998E01,  -2.02998E01,  -2.02998E01 }, 
    {   0.00000E00,   0.00000E00,  0.00000E00,   0.00000E00, -8.62000E-02, 
	-1.66926E00,  -3.61399E00, -5.67795E00,  -8.55561E00,  -1.08785E01, 
	-1.25976E01,  -1.45814E01, -1.59004E01,  -1.73517E01,  -1.83820E01,
	-1.92232E01,  -1.97073E01, -2.02998E01,  -2.02998E01,  -2.02998E01 },
    {   0.00000E00,   0.00000E00,  0.00000E00,   0.00000E00,   0.00000E00, 
	-5.77147E-01,  -1.31053E00, -2.70281E00,  -6.27204E00,  -9.24663E00, 
	-1.10989E01,  -1.30218E01, -1.44008E01,  -1.58917E01,  -1.69460E01, 
	-1.78110E01,  -1.85908E01, -1.97769E01,  -2.02998E01,  -2.02998E01 },
    {   0.00000E00,   0.00000E00,  0.00000E00,   0.00000E00,   0.00000E00, 
	-4.01395E-02, -3.40622E-01, -1.02429E00,  -3.93728E00,  -6.99003E00, 
	-9.04656E00,  -1.09619E01, -1.26314E01,  -1.38775E01,  -1.53915E01, 
	-1.61935E01,  -1.80491E01, -1.92540E01,  -2.02998E01,  -2.02998E01 },
    {   0.00000E00,   0.00000E00,  0.00000E00,   0.00000E00,   0.00000E00, 
	-5.78947E-02, -1.00000E-01,-5.88576E-01,  -2.30318E00,  -5.23570E00, 
	-7.32383E00,  -9.28923E00, -1.12614E01,  -1.25946E01,  -1.38110E01, 
	-1.48958E01,  -1.75274E01, -1.87311E01,  -2.02998E01,  -2.02998E01 },
    {   0.00000E00,   0.00000E00,  0.00000E00,   0.00000E00,   0.00000E00,
	0.00000E00, -6.49270E-02,-3.16327E-01,  -1.54684E00,  -4.00289E00, 
	-6.57448E00,  -8.72839E00, -1.04598E01,  -1.17786E01,  -1.31321E01, 
	-1.44884E01,  -1.69247E01, -1.82082E01,  -2.02998E01,  -2.02998E01 },
    {   0.00000E00,   0.00000E00,  0.00000E00,   0.00000E00, -2.20683E-02, 
	-5.60272E-03, -2.02359E-02,-4.28842E-01,  -1.62681E00,  -3.93839E00, 
	-6.24253E00,  -8.42553E00, -1.01078E01,  -1.16122E01,  -1.29991E01, 
	-1.43545E01,  -1.62633E01, -1.76852E01,  -2.02998E01,  -2.02998E01 },
    {   0.00000E00,   0.00000E00,  0.00000E00,   0.00000E00, -4.69048E-02,  
	0.00000E00, -7.04294E-02,-6.17368E-01,  -2.06572E00,  -4.99934E00, 
	-6.87321E00,  -8.67248E00, -1.02249E01,  -1.17765E01,  -1.32987E01, 
	-1.44689E01,  -1.57474E01, -1.71623E01,  -1.83985E01,  -2.02998E01 },
    {   0.00000E00,   0.00000E00,  0.00000E00,   0.00000E00, -1.48693E-02, 
	-6.65072E-02, -4.48565E-03,-5.98984E-01,  -1.92382E00,  -4.29071E00, 
	-6.69579E00,  -8.83888E00, -1.05148E01,  -1.21332E01,  -1.36353E01, 
	-1.45995E01,  -1.53986E01, -1.66394E01,  -1.78724E01,  -2.02998E01 },
    {   0.00000E00,   0.00000E00,-2.63158E-02, -1.20012E-02,   0.00000E00, 
	-2.63158E-02, -2.63158E-02,-4.02081E-01,  -1.57529E00,  -4.18889E00, 
	-6.35175E00,  -8.26603E00, -9.99547E00,  -1.14291E01,  -1.31123E01, 
	-1.40674E01,  -1.52962E01, -1.61165E01,  -1.73494E01,  -2.02998E01 },
    {   0.00000E00,   0.00000E00,-5.37205E-02, -4.72238E-02, -6.48339E-02, 
	-2.20907E-01, -1.12724E-01,-5.75257E-01,  -2.06047E00,  -4.28772E00, 
	-6.13409E00,  -7.66580E00, -9.29463E00,  -1.06895E01,  -1.20211E01, 
	-1.35413E01,  -1.52724E01, -1.55646E01,  -1.68265E01,  -2.02998E01 },
    {   0.00000E00,   0.00000E00,  0.00000E00, -7.26257E-02, -2.63158E-01, 
	-1.19607E-01, -1.38158E-01, -1.11423E00,  -2.35216E00,  -4.35439E00, 
	-5.87412E00,  -7.24702E00, -8.50691E00,  -9.65950E00,  -1.10283E01, 
	-1.28531E01,  -1.52554E01, -1.51297E01,  -1.61995E01,  -2.02998E01 },
    {   0.00000E00,   0.00000E00,  0.00000E00,   0.00000E00, -3.23182E-02, 
	-1.27059E-01, -2.97511E-01, -1.04808E00,  -2.28456E00,  -3.86254E00, 
	-5.20684E00,  -6.47333E00, -7.76912E00,  -8.95529E00,  -1.00356E01, 
	-1.23924E01,  -1.46763E01, -1.50338E01,  -1.55102E01,  -2.02998E01 },
    {   0.00000E00,   0.00000E00,  0.00000E00, -1.37788E-02, -3.15789E-02, 
	-1.28444E-01, -2.47742E-01, -1.36697E00,  -3.01412E00,  -4.11608E00, 
	-5.27333E00,  -6.30427E00, -7.85186E00,  -8.85404E00,  -1.00083E01, 
	-1.18490E01,  -1.40971E01, -1.50068E01,  -2.02998E01,  -2.02998E01 },
    {   0.00000E00,   0.00000E00,  0.00000E00,   0.00000E00, -5.07458E-02, 
	-2.29208E-01, -8.88953E-01, -2.22389E00,  -3.43865E00,  -4.62398E00, 
	-5.79535E00,  -7.11673E00, -8.10389E00,  -9.60700E00,  -1.05154E01, 
	-1.14129E01,  -1.35179E01, -1.50087E01,  -2.02998E01,  -2.02998E01 }, 
    {   0.00000E00,   0.00000E00,-3.64035E-02, -1.01195E-01, -5.14225E-01, 
	-1.15920E00,  -2.17880E00, -3.20975E00,  -3.84711E00,  -4.66579E00, 
	-5.77850E00,  -6.98499E00, -7.90927E00,  -1.04817E01,  -1.11104E01,  
	-1.13558E01,  -1.29467E01, -2.02998E01,  -2.02998E01, -2.02998E01 },
    {   0.00000E00,   0.00000E00,-4.21987E-02, -2.77493E-01, -1.96612E00, 
	-2.82002E00,  -2.64107E00, -3.37338E00,  -4.32166E00, -5.01108E00, 
	-5.65042E00,  -6.22493E00, -7.43493E00,  -8.79125E00, -1.00875E01, 
	-1.13837E01,  -2.02998E01, -2.02998E01,  -2.02998E01, -2.02998E01} };

  double pimGrid[20] =     { 2.16000E04, 3.25158E04, 4.34316E04, 5.43474E04, 6.52632E04,
			     7.61789E04, 8.70947E04, 9.80105E04, 1.08926E05, 1.19842E05,
			     1.30758E05, 1.41674E05, 1.52589E05, 1.63505E05, 1.74421E05,
			     1.85337E05, 1.96253E05, 2.07168E05, 2.18084E05, 2.29000E05 };
  
  double nGrid[20]   =     { 1.25000E01, 1.71053E01, 2.17105E01, 2.63158E01, 3.09211E01,
			     3.55263E01, 4.01316E01, 4.47368E01, 4.93421E01, 5.39474E01, 
			     5.85526E01, 6.31579E01, 6.77632E01, 7.23684E01, 7.69737E01,
			     8.15789E01, 8.61842E01, 9.07895E01, 9.53947E01, 1.00000E02 };
  
  int row, col;
  double rowf, colf;

  //  double row_min = pimGrid[0], row_max = pimGrid[19], col_min = nGrid[0], col_max = nGrid[19];
  
  for(row = 1;row < 20; row++) {
    if(omega_e < pimGrid[row]) {
      break;	
    }
  }

  for(col = 1;col < 20; col++) {
    if(p_im < nGrid[col]) {
      break;	
    }
  }
  
  rowf = (omega_e - pimGrid[row-1])/(pimGrid[row] - pimGrid[row-1]);
  colf = (p_im - nGrid[row-1])/(nGrid[row] - nGrid[row-1]);
  
  return( (1-rowf)*(1-colf)*dThGrid[row-1][col-1] + rowf*(1-colf)*dThGrid[row][col-1] + 
	  (1-rowf)*colf*dThGrid[row-1][col] + rowf*colf*dThGrid[row][col] );
}


double PI_th_fun( double p_im, double p_ic, double gamma_air )
{
  double x = pow( 2/(gamma_air+1) , gamma_air/(gamma_air-1) );
  if( p_im/p_ic - x> 0 ) {
    return p_im/p_ic;
  } else {
    return x;
  }
} 

double PSIli_th_fun(double PI_th, double PIli_th, double gamma_air)
{
  double psi;
  if( PIli_th - PI_th > 0 ) {
    psi = PSI(PI_th,gamma_air);
    return psi;
  } else {
    psi = PSI(PIli_th,gamma_air);
    return psi*(1-PI_th)/(1-PIli_th);
  }
}

double reg_temp_fun( double T_c, double T_amb, double W_c, double Cic_1, double Cic_2, double Cic_3)
{
  double T_ic;
  
  if( ((T_c-T_amb)*(Cic_1+(T_amb+T_c)*Cic_2/2+W_c*Cic_3)) > 273.15 ) {
    T_ic = T_c - 273.15;
  } else if( ((T_c-T_amb)*(Cic_1+(T_amb+T_c)*Cic_2/2+W_c*Cic_3)) < 0 ) {
    T_ic = T_c;
  } else {
    T_ic =  T_c - (T_c-T_amb)*(Cic_1+(T_amb+T_c)*Cic_2/2+W_c*Cic_3);
  }
  return T_ic;
}


double W_ig_fun(double W_ac, double AF_s, double W_fc, double q_HV, double eta_otto, double eta_ig_ch, double omega_e)
{
  double x = min(W_ac/(AF_s*W_fc),1.0);
  return x*W_fc*q_HV*eta_otto*eta_ig_ch/omega_e*(4*M_PI);
}

double a_eta_t_fun(double omega_tc)
{
  double a_omega[7] = {0,9974928.83487578,16575792.5370014,32109852.6168100,56054167.2758555,76741191.7053602,106900665.763039};
  double omega_breakpoints[7] = {0,6498.90800272609,9788.46967029996,13028.6036134574,15822.8408393152,18113.4807628027,19969.5337025435};

  if( omega_tc < 0) return a_omega[0];

  for( int ii = 1;ii <= 6; ii++ )
  {
    if( omega_tc < omega_breakpoints[ii])
    {
        return a_omega[ii-1] + (a_omega[ii]-a_omega[ii-1])/(omega_breakpoints[ii] - omega_breakpoints[ii-1])*(omega_tc-omega_breakpoints[ii-1]);
    }
  }

  return a_omega[6];
}

double b_eta_t_fun(double omega_tc)
{
  double b_omega[7] = {0,-43392.0188412510,-103018.975128231,-262413.698439652,-533469.217487666,-777194.621957822,-1144332.45600859};
  double omega_breakpoints[7] = {0,6498.90800272609,9788.46967029996,13028.6036134574,15822.8408393152,18113.4807628027,19969.5337025435};

  if( omega_tc < 0 ) return b_omega[0];

  for( int ii = 1;ii <= 6; ii++ )
  {
    if( omega_tc < omega_breakpoints[ii])
    {
        return b_omega[ii-1] + (b_omega[ii]-b_omega[ii-1])/(omega_breakpoints[ii] - omega_breakpoints[ii-1])*(omega_tc-omega_breakpoints[ii-1]);
    }
  }

  return b_omega[6];
}
