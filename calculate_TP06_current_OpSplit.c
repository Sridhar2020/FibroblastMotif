
#include "TP06_OpSplit_2D.h"

double calculate_TP06_current_OpSplit( double *U, double dt, double **lookup, int celltype, double stimCurrent )
{

/* This function returns the total current flow for element n */
/* Includes parameters for four variants as described in TP06
 * paper, with different values for GKr, GKs, GpCa, GpK, and tau_f.
 * Parameter set 2 is equivalent to 'default' epicardial cells. */

	int Vmhi, Vmlo, Vmlo_diff;
	int gain = GAIN;
	int offset = VMOFFSET;

  /* Indices for u array */
  int V =      1;
  int M =      2;
  int H =      3;
  int J =      4;
  int R =      5;
  int S =      6;
  int D =      7;
  int F =      8;
  int F2 =     9;
  int FCass = 10;
  int Xr1 =   11;
  int Xr2 =   12;
  int Xs =    13;
  int RR =    14;
  int OO =    15;
  int CaSS =  16;
  int CaSR =  17;
  int Cai =   18;
  int Nai =   19;
  int Ki =    20;
  

  /* indices for lookup table */
	int na_h_exp       = 1;
	int na_h_inf       = 2;
	int na_j_exp       = 3;
	int na_j_inf       = 4;
	int na_m_exp       = 5;
	int na_m_inf       = 6;
	int ca_d_exp       = 7;
	int ca_d_inf       = 8;
	int ca_f_exp       = 9;
	int ca_f_inf       = 10;
	int ca_f2_exp      = 11;
	int ca_f2_inf      = 12;
	int k_xr1_exp      = 13;
	int k_xr1_inf      = 14;
	int k_xr2_exp      = 15;
	int k_xr2_inf      = 16;
	int k_xs_exp       = 17;
	int k_xs_inf       = 18;
	int to_r_epi_inf   = 19;
	int to_r_epi_exp   = 20;
	int to_s_epi_inf   = 21;
	int to_s_epi_exp   = 22;
	int to_r_endo_inf  = 23;
	int to_r_endo_exp  = 24;
	int to_s_endo_inf  = 25;
	int to_s_endo_exp  = 26;
	int to_r_M_inf     = 27;
	int to_r_M_exp     = 28;
	int to_s_M_inf     = 29;
	int to_s_M_exp     = 30;

  /* Terms for Solution of Conductance and Reversal Potential */
 	const double Rgas = 8314.472;      /* Universal Gas Constant (J/kmol*K) */
  	const double Frdy = 96485.3415;  /* Faraday's Constant (C/mol) */
  	const double Temp = 310.0;    /* Temperature (K) 37C */
  	double RTonF = (Rgas * Temp) / Frdy;
  	double FonRT = 1.0/RTonF;
  	double VmoRTonF = U[V]/RTonF;

//	Cellular capacitance
	const double CAPACITANCE=0.185;
	
//	External concentrations
	const double Ko = 5.4;
	const double KoNorm = 5.4;
	const double Cao=2.0;
	const double Nao=140.0;
	const double Nao3 = 2744000.0;

//	Intracellular volumes
	const double Vc= 0.016404;
	const double Vsr=0.001094;
	const double Vss=0.00005468;

//	Calcium buffering dynamics
	const double Bufc=0.2;
	const double Kbufc=0.001;
	const double Bufsr=10.0;
	const double Kbufsr=0.3;
	const double Bufss=0.4;
	const double Kbufss=0.00025;

//	Intracellular calcium flux dynamics
	const double Vmaxup=0.006375;
	const double Kup=0.00025;
	const double Vrel=0.102;//40.8;
	const double k1bar=0.15;
	const double k2bar=0.045;
	const double k3=0.060;
	const double k4=0.005;//0.000015;
	const double EC=1.5;
	const double maxsr=2.5;
	const double minsr=1.0;
	const double Vleak=0.00036;
	const double Vxfer=0.0038;

//	Parameters for currents
//	Parameters for IKr
	double Gkr;
	const double GkrPar1=0.134;
	const double GkrPar2=0.153;
	const double GkrPar3=0.172;
	const double GkrPar4=0.172;
	const double GkrR2 = 0.134 * 1.75;
	const double GkrR1 = 0.134 * 0.8;

//	Parameters for Iks
	const double pKNa=0.03;
	double Gks;
	const double GksEpi=0.392;
	const double GksEndo=0.392;
	const double GksMcell=0.098;
	double GksPar1=0.270;
	const double GksPar2=0.392;
	const double GksPar3=0.441;
	const double GksPar4=0.441;
	const double GksR2 = 0.270 * 1.75;
	const double GksR1 = 0.270 * 0.8;

//	Parameters for Ik1
	const double GK1=5.405;
//	Parameters for Ito
	double Gto;
	const double GtoEpi=0.294;
	const double GtoEndo=0.073;
	const double GtoMcell=0.294;
//	Parameters for INa
	const double GNa=14.838; // nS/PF
	const double INa_Vshift = 0.0; // -3.4 for acid conditions;
//	Parameters for IbNa
	const double GbNa=0.00029;
//	Parameters for INaK
	const double KmK=1.0;
	const double KmNa=40.0;
	const double knak=2.724;
//	Parameters for ICaL
	const double GCaL=0.00003980;
	const double GCaL_atp = 1.0; //0.87 for atpi 3.0 mM
	const double GCaL_pH = 1.0;//0.922 for pHi = pHo = 7.09
//	Parameters for IbCa
	const double GbCa=0.000592;
//	Parameters for INaCa
	const double naca_pH = 1.0; // 0.743 for acid
	const double knaca=1000;
	const double KmNai=87.5;
	const double KmCa=1.38;
	const double ksat=0.1;
	const double nn=0.35;
//	Parameters for IpCa
	double GpCa;
	const double GpCaPar1=0.0619;
	const double GpCaPar2=0.1238;
	const double GpCaPar3=0.3714;
	const double GpCaPar4=0.8666;
	const double KpCa=0.0005;
//	Parameters for IpK;
	double GpK;
	const double GpKPar1=0.0730;
	const double GpKPar2=0.0146;
	const double GpKPar3=0.0073;
	const double GpKPar4=0.00219;
	const double inverseVcF2=1.0/(2.0*Vc*Frdy);
	const double inverseVcF=1.0/(Vc*Frdy);
	const double inversevssF2=1.0/(2.0*Vss*Frdy);

//  Parameters for IKATP
	double ekatp;              /* K reversal potential (mV) */
	double gkbaratp;           /* Conductance of the ATP-sensitive K channel (nS/uF) */
	const double gkatp = 0;//3.9 ; // (nS/pF) or (mS/uF) = 0.039 nS/nF;/* Maximum conductance
	                              //of the ATP-sensitive K channel (nS/uF) */
	double patp;     		   /* Percentage availibility of open channels */
	const double natp = 0.24;  /* K dependence of ATP-sensitive K current */
	const double atpi = 6.8;   // 4.6 ischaemia or 6.8 normal  /* Intracellular ATP concentraion (mM) */
	const double hatp = 2.0;   /* Hill coefficient */
	const double katp = 0.042; // 0.25 ischaemia, 0.042 normal /* Half-maximal saturation point of
		                       // ATP-sensitive K current (mM) */

        // ATP    6.8   6.5   6.0   5.5   5.0   4.5   4.0
        // kATP   0.042 0.070 0.117 0.164 0.212 0.259 0.306

  	// currents
	double IKr;
	double IKs;
	double IK1;
	double Ito;
	double INa;
	double IbNa;
	double ICaL;
	double IbCa;
	double INaCa;
	double IpCa;
	double IpK;
	double INaK;
	double Irel;
	double Ileak;
	double Iup;
	double Ixfer;
	double k1;
	double k2;
	double kCaSR;
	double IKatp;

  /* activation and inactivation parameters */
	double tau_h, h_inf;
  	double tau_j, j_inf;
  	double tau_m, m_inf;
  	double d_inf, tau_d;
    double f_inf, tau_f;
	// adjust for different parameters
          const double tau_f_multiplier = 0.6; // par 1
	//const double tau_f_multiplier = 1.0; // par 2
	//const double tau_f_multiplier = 1.5; // par 3
	//const double tau_f_multiplier = 2.0; // par 4
	//const double tau_f_multiplierR2 = 1.0; // par R2
    //const double tau_f_multiplierR1 = 0.5; // par R1
    double f2_inf, tau_f2;
    double fCass_inf, tau_fCass;
	double xr1_inf, tau_xr1;
 	double xr2_inf, tau_xr2;
  	double xs_inf, tau_xs;
	double r_inf, tau_r;
	double s_inf, tau_s;

  /* Reversal potentials */
  	double Ena = RTonF*log(Nao/U[Nai]);
	double Ek = RTonF*(log((Ko/U[Ki])));
	double Eks = RTonF*(log((Ko+pKNa*Nao)/(U[Ki]+pKNa*U[Nai])));
	double Eca = 0.5*RTonF*(log((Cao/U[Cai])));

	double naca1, naca2, naca3;
	double Ak1;
	double Bk1;
	double rec_iK1;
	double rec_iNaK;
	double rec_ipK;

  /* ion concentrations and parameters for calculating them */
	double dRR;
	double CaCSQN;
	double dCaSR;
	double bjsr, cjsr;
	double CaSSBuf;
	double dCaSS;
	double bcss, ccss;
	double CaBuf;
	double dCai;
	double bc, cc;
	double dNai;
	double dKi;

	/* other variables */
	double alpha_h;  /* Na alpha-h rate constant (ms^-1) */
	double beta_h;   /* Na beta-h rate constant (ms^-1) */
	double alpha_j;  /* Na alpha-j rate constant (ms^-1) */
	double beta_j;   /* Na beta-j rate constant (ms^-1) */
	double alpha_m;  /* Na alpha-m rate constant (ms^-1) */
 	double beta_m;   /* Na beta-m rate constant (ms^-1) */
    double ad, bd, cd;
 	double af, bf, cf;
	double af2, bf2, cf2;
	double axr1, bxr1;
	double axr2, bxr2;
	double axs, bxs;

/* inhomogeneity parameter for 200 x 3 geometry */
	//int m;

	//if (n<=200) m = n;
	//else if (n<=400) m = n-200;
	//else m = n-400;

	//GksPar1 = GksPar1 + ((double) m)/200.0;

/* lookup table indices */

	// for INa shift by up to -3.4 mV
  	Vmhi = ceil(U[V]+INa_Vshift) * gain + offset;
  	Vmlo_diff = U[V]+INa_Vshift - floor(U[V]+INa_Vshift);
  	Vmlo = floor(U[V]+INa_Vshift) * gain + offset;

  /* Inward current iNa */
  /* lookup table code */
  	/*m_inf = lookup[na_m_inf][Vmlo] + Vmlo_diff * (lookup[na_m_inf][Vmhi] - lookup[na_m_inf][Vmlo]);
  	tau_m = lookup[na_m_exp][Vmlo] + Vmlo_diff * (lookup[na_m_exp][Vmhi] - lookup[na_m_exp][Vmlo]);

  	h_inf = lookup[na_h_inf][Vmlo] + Vmlo_diff * (lookup[na_h_inf][Vmhi] - lookup[na_h_inf][Vmlo]);
  	tau_h = lookup[na_h_exp][Vmlo] + Vmlo_diff * (lookup[na_h_exp][Vmhi] - lookup[na_h_exp][Vmlo]);

  	j_inf = lookup[na_j_inf][Vmlo] + Vmlo_diff * (lookup[na_j_inf][Vmhi] - lookup[na_j_inf][Vmlo]);
  	tau_j = lookup[na_j_exp][Vmlo] + Vmlo_diff * (lookup[na_j_exp][Vmhi] - lookup[na_j_exp][Vmlo]);
       */	
  	alpha_m = 1.0/(1.0+exp((-60.0-U[V])/5.0));
	beta_m = 0.1/(1.0+exp((U[V]+35.0)/5.0))+0.10/(1.0+exp((U[V]-50.0)/200.0));
	tau_m = alpha_m * beta_m;
	m_inf = 1.0/((1.0+exp((-56.86-U[V])/9.03))*(1.0+exp((-56.86-U[V])/9.03)));

	if (U[V] >= -40.0)
	    {
	    alpha_h = 0.0;
	    beta_h = (0.77/(0.13*(1.0 + exp(-(U[V] + 10.66)/11.1))));
	    }
	else
	    {
	    alpha_h = (0.057 * exp(-(U[V] + 80.0)/6.8));
	    beta_h = (2.7 * exp(0.079 * U[V])+(3.1e5) * exp(0.3485 * U[V]));
	    }
	tau_h = 1.0/(alpha_h + beta_h);
	h_inf = 1.0/((1.0 + exp((U[V] + 71.55)/7.43))*(1.0 + exp((U[V] + 71.55)/7.43)));

	if(U[V] >= -40.0)
	    {
	    alpha_j = 0.;
	    beta_j = (0.6 * exp((0.057) * U[V])/(1.0 + exp(-0.1 * (U[V] + 32.0))));
	    }
	else
	    {
		alpha_j = (((-2.5428e4)*exp(0.2444*U[V])-(6.948e-6)*exp(-0.04391*U[V]))*(U[V]+37.78)/(1.0+exp(0.311*(U[V]+79.23))));
		beta_j = (0.02424*exp(-0.01052*U[V])/(1.0+exp(-0.1378*(U[V]+40.14))));
	    }
	tau_j = 1.0/(alpha_j + beta_j);
	j_inf = h_inf;
       	 
 
  	U[M] = m_inf - ( m_inf - U[M] ) * exp( -dt / tau_m );
  	U[H] = h_inf - ( h_inf - U[H] ) * exp( -dt / tau_h );
	U[J] = j_inf - ( j_inf - U[J] ) * exp( -dt / tau_j );

  	INa = GNa*U[M]*U[M]*U[M]*U[H]*U[J]*(U[V]-Ena);

  	Vmhi = ceil(U[V]) * gain + offset;
  	Vmlo_diff = U[V] - floor(U[V]);
  	Vmlo = floor(U[V]) * gain + offset;

  /* Currents in Ca channels */
  /* lookup table code */
      
  	/*d_inf = lookup[ca_d_inf][Vmlo] + Vmlo_diff * (lookup[ca_d_inf][Vmhi] - lookup[ca_d_inf][Vmlo]);
  	tau_d = lookup[ca_d_exp][Vmlo] + Vmlo_diff * (lookup[ca_d_exp][Vmhi] - lookup[ca_d_exp][Vmlo]);

  	f_inf = lookup[ca_f_inf][Vmlo] + Vmlo_diff * (lookup[ca_f_inf][Vmhi] - lookup[ca_f_inf][Vmlo]);
  	tau_f = lookup[ca_f_exp][Vmlo] + Vmlo_diff * (lookup[ca_f_exp][Vmhi] - lookup[ca_f_exp][Vmlo]);

	if (U[V] >= 0) tau_f *= tau_f_multiplier; // Vm >= 0 added 20/12/2010

  	f2_inf = lookup[ca_f2_inf][Vmlo] + Vmlo_diff * (lookup[ca_f2_inf][Vmhi] - lookup[ca_f2_inf][Vmlo]);
  	tau_f2 = lookup[ca_f2_exp][Vmlo] + Vmlo_diff * (lookup[ca_f2_exp][Vmhi] - lookup[ca_f2_exp][Vmlo]);
       */
   	
    d_inf = 1.0/(1.0+exp((-8.0-U[V])/7.5));
    ad = 1.4/(1.0+exp((-35.0-U[V])/13.0))+0.25;
    bd=1.4/(1.0+exp((U[V]+5.0)/5.0));
	cd=1.0/(1.0+exp((50.0-U[V])/20.0));
    tau_d = ad * bd * cd;

    f_inf = 1.0/(1.0+exp((U[V]+20.0)/7.0));
    af=1102.5*exp(-(U[V]+27.0)*(U[V]+27.0)/225.0);
	bf=200.0/(1.0+exp((13.0-U[V])/10.0));
	cf=(180.0/(1.0+exp((U[V]+30.0)/10.0)))+20.0;
   if(U[V] < 0)
   {
    tau_f = af + bf + cf;
   }
  else if(U[V] >= 0 )
   {
    tau_f= tau_f_multiplier*(af + bf + cf);
   } 
  
    f2_inf = 0.67/(1.0+exp((U[V]+35.0)/7.0))+0.33;
    af2=562.0*exp(-(U[V]+27.0)*(U[V]+27.0)/240.0);
	bf2=31.0/(1.0+exp((25.0-U[V])/10.0));
	cf2=80.0/(1.0+exp((U[V]+30.0)/10.0));
	tau_f2 = af2 + af2 + af2;
  

	fCass_inf = 0.6/(1.0+(U[CaSS]/0.05)*(U[CaSS]/0.05))+0.4;
	tau_fCass = 80.0/(1.0+(U[CaSS]/0.05)*(U[CaSS]/0.05))+2.0;

     U[D] = d_inf - (d_inf - U[D]) * exp( -dt / tau_d );
  	U[F] = f_inf - (f_inf - U[F]) * exp( -dt / tau_f );
  	U[F2] = f2_inf - (f2_inf - U[F2]) * exp( -dt / tau_f2 );
  	U[FCass] = fCass_inf - (fCass_inf - U[FCass]) * exp( -dt / tau_fCass );

  	ICaL = GCaL_pH*GCaL_atp*GCaL*U[D]*U[F]*U[F2]*U[FCass]*4.0*(U[V]-15.0)*(Frdy/RTonF)*(0.25*exp(2.0*(U[V]-15.0)/RTonF)*U[CaSS]-Cao)/(exp(2.0*(U[V]-15.0)/RTonF)-1.0);

  /* Rapidly inactivating K current */
  /* lookup table code */
	/*xr1_inf = lookup[k_xr1_inf][Vmlo] + Vmlo_diff * (lookup[k_xr1_inf][Vmhi] - lookup[k_xr1_inf][Vmlo]);
	tau_xr1 = lookup[k_xr1_exp][Vmlo] + Vmlo_diff * (lookup[k_xr1_exp][Vmhi] - lookup[k_xr1_exp][Vmlo]);

	xr2_inf = lookup[k_xr2_inf][Vmlo] + Vmlo_diff * (lookup[k_xr2_inf][Vmhi] - lookup[k_xr2_inf][Vmlo]);
	tau_xr2 = lookup[k_xr2_exp][Vmlo] + Vmlo_diff * (lookup[k_xr2_exp][Vmhi] - lookup[k_xr2_exp][Vmlo]);
        */ 
  	
    xr1_inf = 1.0/(1.0+exp((-26.0-U[V])/7.0));
    axr1 = 450.0/(1.0+exp((-45.0-U[V])/10.0));
    bxr1 = 6.0/(1.0+exp((U[V]+30.0)/11.5));
    tau_xr1 = axr1 * bxr1;

    xr2_inf = 1.0/(1.0+exp((U[V]+88.0)/24.0));
    axr2 = 3.0/(1.0+exp((-60.0-U[V])/20.0));
    bxr2 = 1.12/(1.0+exp((U[V]-60.0)/20.0));
    tau_xr2 = axr2 * bxr2;
     
 
 	Gkr = GkrPar1; //GkrPar4;
 	U[Xr1] = xr1_inf - (xr1_inf - U[Xr1]) * exp( -dt / tau_xr1 );
 	U[Xr2] = xr2_inf - (xr2_inf - U[Xr2]) * exp( -dt / tau_xr2 );
	IKr = Gkr*sqrt(Ko/5.4)*U[Xr1]*U[Xr2]*(U[V]-Ek);

  /* Slowly inactivating K current */
  /* celltype is 0 (endo) 2 (epi) or 1 (M) */
  /* lookup table code */
         
  	/*xs_inf = lookup[k_xs_inf][Vmlo] + Vmlo_diff * (lookup[k_xs_inf][Vmhi] - lookup[k_xs_inf][Vmlo]);
	tau_xs = lookup[k_xs_exp][Vmlo] + Vmlo_diff * (lookup[k_xs_exp][Vmhi] - lookup[k_xs_exp][Vmlo]);
        */       
      	
 	xs_inf = 1.0/(1.0+exp((-5.0-U[V])/14.0));
	axs = (1400.0/(sqrt(1.0+exp((5.0-U[V])/6.0))));
	bxs = (1.0/(1.0+exp((U[V]-35.0)/15.0)));
	tau_xs = axs * bxs + 80.0;
     	

	U[Xs] = xs_inf - (xs_inf - U[Xs]) * exp( -dt / tau_xs );

	//if (celltype == 0) { Gks = GksEndo; }
	//else if (celltype == 1) { Gks = GksMcell; }
	//else  {Gks = GksEpi; }
	Gks = GksPar1; //GksPar4;
	IKs = Gks*U[Xs]*U[Xs]*(U[V]-Eks);

  /* Time independent K current */
  	Ak1 = 0.1/(1.0+exp(0.06*(U[V]-Ek-200.0)));
	Bk1 = (3.0*exp(0.0002*(U[V]-Ek+100.0))+exp(0.1*(U[V]-Ek-10.0)))/(1.0+exp(-0.5*(U[V]-Ek)));
	rec_iK1 = Ak1/(Ak1+Bk1);
	IK1 = GK1*rec_iK1*(U[V] - Ek);

  /* Plateau K current */
  	GpK=GpKPar1; //GpKPar4;
	rec_ipK = 1.0/(1.0+exp((25.0-U[V])/5.98));
	IpK=GpK*rec_ipK*(U[V]-Ek);

  /* transient outward current */
  /* celltype is 0 (endo) 2 (epi) or 1 (M) */

  /* EPI only code */

		/* lookup table code */
          /*      r_inf = lookup[to_r_epi_inf][Vmlo] + Vmlo_diff * (lookup[to_r_epi_inf][Vmhi] - lookup[to_r_epi_inf][Vmlo]);
		tau_r = lookup[to_r_epi_exp][Vmlo] + Vmlo_diff * (lookup[to_r_epi_exp][Vmhi] - lookup[to_r_epi_exp][Vmlo]);
		s_inf = lookup[to_s_epi_inf][Vmlo] + Vmlo_diff * (lookup[to_s_epi_inf][Vmhi] - lookup[to_s_epi_inf][Vmlo]);
		tau_s = lookup[to_s_epi_exp][Vmlo] + Vmlo_diff * (lookup[to_s_epi_exp][Vmhi] - lookup[to_s_epi_exp][Vmlo]);
	  */    	
	       
               r_inf = 1.0/(1.0+exp((20.0-U[V])/6.0));
		s_inf = 1.0/(1.0+exp((U[V]+20.0)/5.0));
		tau_r = 9.5*exp(-(U[V]+40.0)*(U[V]+40.0)/1800.0)+0.8;
		tau_s = 85.0*exp(-(U[V]+45.0)*(U[V]+45.0)/320.0)+5.0/(1.0+exp((U[V]-20.0)/5.0))+3.0;
	       	
		Gto = GtoEpi;//GtoEpi;

	U[S] = s_inf - (s_inf - U[S]) * exp(-dt / tau_s);
	U[R] = r_inf - (r_inf - U[R]) * exp(-dt / tau_r);
	Ito = Gto*U[R]*U[S]*(U[V]-Ek);


  /* ATP dependent K current */
	ekatp = (RTonF * log(Ko/Ki));
	patp = 1.0/(1.0+(pow((atpi/katp),hatp)));
	gkbaratp = gkatp*patp*(pow((Ko/KoNorm),natp));

	IKatp = gkbaratp*(U[V]-ekatp);
        IKatp = 0.0; 

  /* Na Ca exchanger */
	naca1 = knaca*(1.0/(KmNai*KmNai*KmNai+Nao3))*(1.0/(KmCa+Cao));
	naca2 = (1.0/(1.0+ksat*exp((nn-1.0)*VmoRTonF)));
	naca3 = (exp(nn*VmoRTonF)*U[Nai]*U[Nai]*U[Nai]*Cao-exp((nn-1.0)*VmoRTonF)*Nao3*U[Cai]*2.5);
	INaCa = naca_pH * naca1 * naca2 * naca3;

	//INaCa=knaca*(1.0/(KmNai*KmNai*KmNai+Nao3))*(1.0/(KmCa+Cao))*(1.0/(1.0+ksat*exp((nn-1.0)*VmoRTonF)))*(exp(nn*VmoRTonF)*U[Nai]*U[Nai]*U[Nai]*Cao-exp((nn-1.0)*VmoRTonF)*Nao3*Cai*2.5);

  /* Background Na current */
	IbNa=GbNa*(U[V]-Ena);

  /* iNaK */
	rec_iNaK = (1.0/(1.0+0.1245*exp(-0.1*VmoRTonF)+0.0353*exp(-VmoRTonF)));
	INaK=knak*(Ko/(Ko+KmK))*(U[Nai]/(U[Nai]+KmNa))*rec_iNaK;

  /* Plateau Ca current */
 	GpCa=GpCaPar1; //GpCaPar4;
  	IpCa=GpCa*U[Cai]/(KpCa+U[Cai]);

  /* Background Ca current */
	IbCa=GbCa*(U[V]-Eca);

  /* intracellular ion concentrations */

	kCaSR = maxsr-((maxsr-minsr)/(1.0+(EC/U[CaSR])*(EC/U[CaSR])));
	k1 = k1bar/kCaSR;
	k2 = k2bar*kCaSR;
	dRR = k4 * (1.0-U[RR]) - k2*U[CaSS]*U[RR];
	U[RR] += dt*dRR;
	U[OO] = k1*U[CaSS]*U[CaSS]*U[RR]/(k3+k1*U[CaSS]*U[CaSS]);
	Irel = Vrel*U[OO]*(U[CaSR]-U[CaSS]);
	Ileak = Vleak*(U[CaSR]-U[Cai]);
	Iup = Vmaxup/(1.0+((Kup*Kup)/(U[Cai]*U[Cai])));
	Ixfer = Vxfer*(U[CaSS] - U[Cai]);

	CaCSQN = Bufsr*U[CaSR]/(U[CaSR]+Kbufsr);
	dCaSR = dt*(Iup-Irel-Ileak);
	bjsr = Bufsr-CaCSQN-dCaSR-U[CaSR]+Kbufsr;
	cjsr = Kbufsr*(CaCSQN+dCaSR+U[CaSR]);
	U[CaSR] = (sqrt(bjsr*bjsr+4.0*cjsr)-bjsr)/2.0;

	CaSSBuf=Bufss*U[CaSS]/(U[CaSS]+Kbufss);
	dCaSS = dt*(-Ixfer*(Vc/Vss)+Irel*(Vsr/Vss)+(-ICaL*inversevssF2*CAPACITANCE));
	bcss = Bufss-CaSSBuf-dCaSS-U[CaSS]+Kbufss;
	ccss = Kbufss*(CaSSBuf+dCaSS+U[CaSS]);
	U[CaSS] = (sqrt(bcss*bcss+4.0*ccss)-bcss)/2.0;

	CaBuf = Bufc*U[Cai]/(U[Cai]+Kbufc);
	dCai = dt*((-(IbCa+IpCa-2.0*INaCa)*inverseVcF2*CAPACITANCE)-(Iup-Ileak)*(Vsr/Vc)+Ixfer);
	bc = Bufc-CaBuf-dCai-U[Cai]+Kbufc;
	cc = Kbufc*(CaBuf+dCai+U[Cai]);
	U[Cai] = (sqrt(bc*bc+4.0*cc)-bc)/2.0;

	dNai=-(INa+IbNa+3.0*INaK+3.0*INaCa)*inverseVcF*CAPACITANCE;
	U[Nai] += dt*dNai;

	dKi=-(stimCurrent+IK1+Ito+IKr+IKs-2.0*INaK+IpK)*inverseVcF*CAPACITANCE;
	U[Ki] += dt*dKi;

//	if (n==1) printf("%g %g %g %g %g %g\n",INa,Ito,ICaL,IKr,IKs,U[Cai]);
	return( IKr + IKs + IK1 + Ito + IKatp + INa + IbNa + ICaL + IbCa + INaK + INaCa + IpCa + IpK + stimCurrent );

}
