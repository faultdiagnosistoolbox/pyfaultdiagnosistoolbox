#ifndef _EXTMODELFUNS_H_
#define _EXTMODELFUNS_H_

#ifdef _WIN32
#define _USE_MATH_DEFINES
#include <cmath>
#endif


double eta_c_fun(double eta_cmax, double eta_cmin, double W_ccorr, double W_ccorrmax, double PI_c, double PI_cmax, double Q_c11, double Q_c12, double Q_c22);
double eta_t_fun( double k3, double BSR, double k4, double eta_min );
double PHI_model_fun(double K1, double K2, double PSI_c);
double PI_c_fun(double p_c, double p_af);
double PI_t_fun( double p_t, double p_em );
double PI_wg_fun( double p_t, double p_em, double gamma_air );
double PSI( double PI, double gamma);
double PSI_th_fun(double p_im, double p_ic, double gamma_air, double PIli_th);
double PSIli_wg_fun( double PI_wg, double PIli_wg, double gamma_exh);
double Tflow_wg_fun(double p_t, double p_em, double T_em, double T_t);
double Tq_t_fun(double gamma_eg, double cp_eg, double eta_t, double W_t, double T_em, double PI_t, double omega_tc);
double W_af_fun( double p_amb,double p_af,double plin_af,double H_af,double T_amb );
double W_c_fun(double PI_cnolim, double p_af, double U_c, double D_c, double T_af, double R_air, double PHI_model, double fix_gain);
double W_ic_fun(double p_c, double p_ic, double plin_intercooler, double H_intercooler, double T_c);
double W_es_fun(double p_amb, double p_t, double plin_exh, double H_exhaust, double T_t);
double W_t_fun(double p_em, double k1, double PI_t, double k2, double T_em);

double igngrid(double omega_e, double p_im);
double PI_th_fun( double p_im, double p_ic, double gamma_air );
double PSIli_th_fun(double PI_th, double PIli_th, double gamma_air);
double reg_temp_fun( double T_c, double T_amb, double W_c, double Cic_1, double Cic2, double Cic3);
double W_i_g_fun( double W_ac, double AF_s, double W_fc, double q_HV, double eta_otto, double eta_ig_ch, double omega_e);
double W_ic_fun(double p_c, double p_ic, double plin_intercooler, double H_intercooler, double T_c);
double W_ig_fun(double W_ac, double AF_s, double W_fc, double q_HV, double eta_otto, double eta_ig_ch, double omega_e);

double a_eta_t_fun(double omega_tc);
double b_eta_t_fun(double omega_tc);

double max_fun(double choice_a, double choice_b);
double min_fun(double choice_a, double choice_b);

#endif
