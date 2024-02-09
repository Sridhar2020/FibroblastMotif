#include "TP06_OpSplit_2D.h"

double calculate_MacFib_current( double *U, double dt)
{
/* Indices for u array */
  int fib_volt = 21;
  int rkv = 22;
  int skv = 23;
   
  /*Fibroblast declarations*/
  double fib_alpha_K1;
  double fib_beta_K1;
  double fib_rec_IK1;
  double fib_IKv;
  double fib_IK1;
  double fib_INaK;
  double fib_IbNa;
  double fib_current;
  double fib_to_myo;
  double myo_to_fib;
  double fib_r_bar;
  double fib_s_bar;
  double tr1;
  double tr2;
  double fib_tau_r;
  double fib_tau_s;


 
  /*Parameters for fibroblasts*/
 double FIB_GKV = 0.155; //0.25; // units: ns/pF
 double FIB_GK1 = 0.4822; //units: ns/pF
 double FIB_I_NAK = 2.002; //units pA/pF
 double FIB_G_BNA = 0.0095; //units ns/pF
 double fib_Ek = -87.0; //units mV
 double fib_B = -200.0; //units mV
 double V_rev = -150.0; //units mV
 double FIB_K_mK = 1.0; //units mmol/L
 double FIB_K_mna = 11.0; //units mmol/L
 double FIB_Ko = 5.3581; //units mmol/L
 double FIB_Nai = 8.5547; //units mmol/L
 double FIB_CAP = 0.050; //nF //units pF 
 double FIB_Nao = 130.0;
 
       const double Rgas = 8314.472;      /* Universal Gas Constant (J/kmol*K) */
  	const double Frdy = 96485.3415;  /* Faraday's Constant (C/mol) */
  	const double Temp = 310.0;    /* Temperature (K) 37C */
  	double RTonF = (Rgas * Temp) / Frdy;

   /* Fibroblast sodium reversal potentials */
        double Ena_f = RTonF*log(FIB_Nao/FIB_Nai); //Nain
	

    /*inward rectifying potassium*/ 
    fib_alpha_K1 = 0.1/(1+exp(0.06*(U[fib_volt]-fib_Ek-200)));
    fib_beta_K1= 3*exp(0.0002*(U[fib_volt]-fib_Ek+100))+exp(0.1*(U[fib_volt]-fib_Ek-10))/(1+exp(-0.5*(U[fib_volt]-fib_Ek)));
    fib_rec_IK1= fib_alpha_K1/(fib_alpha_K1 + fib_beta_K1);
    fib_IK1 = FIB_GK1*fib_rec_IK1*(U[fib_volt]-fib_Ek);

    /*Na-K pump in the fibroblast*/
     fib_INaK = FIB_I_NAK*(FIB_Ko/(FIB_Ko+FIB_K_mK))*(pow(FIB_Nai,1.5)/(pow(FIB_Nai,1.5)+pow(FIB_K_mna,1.5)))*((U[fib_volt]-V_rev)/(U[fib_volt] - fib_B));  

   /*inward Kv current*/
     fib_r_bar = 1/(1 + exp(-(U[fib_volt] + 20.0 -30)/11));
     fib_s_bar = 1/(1 + exp((U[fib_volt] + 23.0 -30 )/7));
     tr1 = (U[fib_volt] + 20 -30)/25.9;
     tr2 = (U[fib_volt] + 23 -30)/22.7;  
     fib_tau_r = 20.3 + 138*exp(-(tr1*tr1));
     fib_tau_s = 1574 + 5268*exp(-(tr2*tr2));	    
     U[rkv] = fib_r_bar - (fib_r_bar-U[rkv])*exp(-dt/fib_tau_r);
     U[skv] = fib_s_bar - (fib_s_bar-U[skv])*exp(-dt/fib_tau_s); 
     fib_IKv = FIB_GKV*U[rkv]*U[skv]*(U[fib_volt]-fib_Ek);       
   
  /*Background sodium current*/
     fib_IbNa = FIB_G_BNA*(U[fib_volt] - Ena_f);    
   
    //Fibroblast current  
    
     fib_current = fib_IKv + fib_IK1 + fib_INaK + fib_IbNa;
   return (fib_current);
}
