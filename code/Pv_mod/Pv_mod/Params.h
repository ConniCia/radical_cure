/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
///                                                                       ///
///  Individual-based Plasmodium vivax transmission model.                ///
///                                                                       ///
///  Dr Michael White                                                     ///
///  Institut Pasteur                                                     ///
///  michael.white@pasteur.fr                                             ///
///                                                                       ///
///  This contains the Params struct and parameter reading code.          ///
///                                                                       ///
///  This header is common to all model code.                             ///
///                                                                       ///
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

// Include guard
#ifndef PVIVAX_MODEL_PARAMS
#define PVIVAX_MODEL_PARAMS

#include <vector>

using namespace std;


///////////////////////////////////////////////////////////////////
//                                                               //
// 2.1.1. Define model constants                                 //
//                                                               //
///////////////////////////////////////////////////////////////////

#define t_step 1            // Time step for update in humans
#define mosq_steps 20       // Number of mosquito steps per human step

#define N_H_comp 6          // Number of human compartments (indexed by p)
#define N_M_comp 6          // Number of mossquito compartments (indexed by p)

#define N_gen 3             // Number of gender categories (indexed by g)
#define N_age 58            // Number of age categories for calculation of equilibrium set up (indexed by i)
#define N_het 11            // Number of heterogeneity categories for calculation of equilibrium set up (indexed by j)
#define K_max 10            // Maximum umber of hypnozoites (indexed by k)
#define N_int 11            // Number of interventions

#define N_mosq 2            // Number of mosquito species to be modelled (indexed by v)

#define N_eq_setup 5       // Number of times for initialising equilibrium (indexed by eq)

const double CONST_LOG_2 = 0.6931471805599453094172321214581766L;


struct SimTimes {
	double start, end, burnin;
};


///////////////////////////////////////////////////////////////////
//                                                               //
// 2.1.2. Define structure for parameters                        //
//                                                               //
///////////////////////////////////////////////////////////////////

struct Params
{
	//////////////////////////////////////////////////////////////////////////
	//  Class member functions
	//////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////
	// Read parameters from input files

	SimTimes read(const char *parameter_File, const char *mosquito_File[N_mosq]);


	//////////////////////////////////////////////////////////////////////////
	//  Data
	//////////////////////////////////////////////////////////////////////////

	int N_pop;                     // human population size


	/////////////////////////////////////
	// Equilibrium EIR (no seasonality)

	double EIR_dom_equil;          // EIR at equilibrium (domestic)
	double EIR_occ_equil;          // EIR at equilibrium (occupational)


	/////////////////////////////////
	// Human demography

	double mu_H;                   // human death rate
	double age_mean;               // mean age of human population
	double age_max;                // maximum age in human population
	double het_max;                // maximum heterogeneity

	double P_occup;                // Proportion of males with occupational exposure
	double risk_occup;             // Extra risk from occupational mosquitoes

	double occup_age_low;          // lowest age of occupational exposure
	double occup_age_high;         // highest age of occupational exposure

	double P_preg;                 // Proportion pregnant females within age range
	double dur_preg;               // Duration of pregnancy and breastfeeding
	double preg_rate;              // Rate at which age eligible women become pregnant

	double preg_age_low;           // lowest age of pregnancy
	double preg_age_high;          // highest age of pregnancy


	////////////////////////////////////////////////////////
	// Heterogeneity in exposure and age-dependent biting

	double age_0;                  // age-dependent biting parameter
	double rho_age;                // age-dependent biting parameter

	double sig_het;                // heterogeneity in exposure - standard deviation on a log scale


	/////////////////////////////////
	// Transmission probabilities

	double bb;                     // mosquito -> human transmission probability
	double c_PCR;                  // human -> mosquito transmission probability - PCR detectable
	double c_LM;                   // human -> mosquito transmission probability - LM detectable
	double c_D;                    // human -> mosquito transmission probability - clinical disease
	double c_T;                    // human -> mosquito transmission probability - treatment

	double d_latent;               // duration of latency in the liver
	int H_track;                   // Number of steps for tracking lagged lam_H (due to duration of latency)


	/////////////////////////////////
	// Human recovery parameters

	double r_LM;                   // rate of recovery from LM detectable infection
	double r_D;                    // rate of recovery from symptomatic disease
	double r_T;                    // rate of progression through treatment
	double r_P;                    // rate of loss of prophylaxis

	double d_PCR_min;              // minimum duration of PCR-detectable infection - full immunity
	double d_PCR_max;              // maximum duration of PCR-detectable infection - no immunity
	double A_PCR_50pc;             // scale parameter for effect of anti-parasite immunity on PCR-detectable infection
	double K_PCR;                  // shape parameter for effect of anti-parasite immunity on PCR-detectable infection

	double r_PCR;                  // recovery rate from PCR detectable infection - depends on level of immunity


	////////////////////////////////
	// Anti-parasite immunity parameters

	double u_par;                  // scale parameter for acquisition of anti-parasite immunity
	double r_par;                  // rate of decay of anti-parasite immunity
	double a_par_rel;              // following a relapse: increment of anti-parasite immunity

	double phi_LM_max;             // probability of LM-detectable infection with no immunity
	double phi_LM_min;             // probability of LM-detectable infection with full immunity
	double phi_LM_rel;             // following a relapse: factor of phi_LM
	double A_LM_50pc;              // anti-parasite immunity scale parameter
	double K_LM;                   // anti-parasite immunity shape parameter

	double phi_LM;                 // probability that PCR detectable infection becomes detectable by light microscopy (LM) - DYNAMICALLY UPDATED


	////////////////////////////////
	// Clinical immunity parameters

	double u_clin;                 // scale parameter for acquisition of clinical immunity
	double r_clin;                 // rate of decay of clinical immunity
	double a_clin_rel;             // following a relapse: increment of clinical immunity

	double phi_D_max;              // probability of clinical episode with no immunity
	double phi_D_min;              // probability of clinical episode with full immunity
	double phi_D_rel;              // following a relapse: factor of phi_D
	double A_D_50pc;               // clinical immunity scale parameter
	double K_D;                    // clinical immunity shape parameter

	double phi_D;                  // probability that LM detectable infection progresses to symptomatic disease - DYNAMICALLY UPDATED


	/////////////////////////////////////////////
	// maternal immunity

	double P_mat;                  // New-born immunity relative to mother's
	double d_mat;                  // Inverse of decay rate of maternal immunity


	/////////////////////////////////
	// Relapse parameters

	double ff;                     // relapse rate
	double gamma_L;                // liver clearance rate


	/////////////////////////////////
	// Baseline case management parameters

	int CM_regimen;             // treatment regimen: 0 = blood-stage (baseline default); 1 = primaquine; 2 = tafenoquine
	double CM_cover;            // proportion of symptomatic cases treated

	double CM_CQ_eff;	        // chloroquine efficacy (on its own)
	double CM_CQ_eff_wPQ;	    // chloroquine efficacy (co-administered with primaquine)
	double CM_CQ_proph;		    // duration of chloroquine prophylaxis

	double CM_PQ_eff;		    // primaquine efficacy
	double CM_PQ_proph;		    // duration of primaquine_prophylaxis
	double CM_PQ_adhere;		// primaquine adherence
	double CM_PQ_lowage;	    // minimum ag for primaquine
	int CM_PQ_G6PD_risk;	    // G6PD risk for primaquine
	int CM_PQ_CYP2D6_risk;	    // CYP2D6 risk for primaquine
	int	CM_PQ_preg_risk;	    // G6PD risk for pregnancy

	double CM_TQ_eff;		    // tafenoquine efficacy
	double CM_TQ_proph;		    // duration of tafenoquine prophylaxis
	double CM_TQ_adhere;		// tafenoquine adherence
	double CM_TQ_lowage;	    // minimum ag for tafenoquine
	int CM_TQ_G6PD_risk;	    // G6PD risk for tafenoquine
	int CM_TQ_CYP2D6_risk;	    // CYP2D6 risk for tafenoquine
	int	CM_TQ_preg_risk;	    // G6PD risk for pregnancy

	int	CM_G6PD_test; 	        // G6PD testing


	double CM_CQ_coveff;        // product of coverage and efficacy for chloroquine
	double CM_PQ_coveff;        // product of coverage and efficacy for primaquine
	double CM_TQ_coveff;        // product of coverage and efficacy for tafenoquine


	////////////////////////////////
	// G6PD and CYP2D6

	double G6PD_prev;              // Prevalence of G6PD deficiency
	double mu_G6PD_nor;            // Mean G6PD activity in normals
	double sig_G6PD_nor;           // Standard deviation in G6PD activity in normals
	double mu_G6PD_het;            // Mean G6PD activity in heterozygous deficients
	double sig_G6PD_het;           // Standard deviation in G6PD activity in heterozygous deficients
	double mu_G6PD_def;            // Mean G6PD activity in homozygous deficients
	double sig_G6PD_def;           // Standard deviation in G6PD activity in homozygous deficients

	double CYP2D6_prev;			   // Prevalence of CYP2D6 phenotype

	double low_G6PD_activity[N_gen];   // Proportion with G6PD activity less than 30%


	///////////////////////////////////////////////////////////////////
	// Correlation between interventions and intervention rounds

	double rho_round_LLIN;      // Between round correlation of LLINS
	double rho_round_IRS;       // Between round correlation of IRS
	double rho_LLIN_IRS;        // Correlation between LLIN and IRS coverage

	double rho_round_MDA;       // Between round correlation of MDA
	double rho_MDA_VC;          // Correlation between MDA and vector control (LLINs or IRS)

	double sig_round_LLIN;      // Derived parameter for correlation between rounds of LLINs
	double sig_round_IRS;       // Derived parameter for correlation between rounds of IRS
	double sig_round_MDA;       // Derived parameter for correlation between rounds of MDA

	float V_int[N_int][N_int];
	float V_int_dummy[N_int][N_int];


	/////////////////////////////////
	// Entomological parameters

	double mm_0[N_mosq];                   // number of mosquitoes per human (An. farauti)
	double aa[N_mosq];                     // mosquito biting rate (in the absence of vector control)
	double mu_M[N_mosq];                   // mosquito death rate
	double tau_M[N_mosq];                  // duration of sporogony

	double lam_M[N_mosq];                  // Force of infection on mosquites - updated dynamically

	int M_track;                           // Number of steps for tracking lagged lam_M*S_M (needed for lag due to duration of sporogony)
	vector<vector<double>> lam_S_M_track;  // Lagged force of infection on moquitoes


	////////////////////////////////
	// Seasonality parameters

	double dry_seas[N_mosq];        // Proportion of dry season transmission compared to mean
	double kappa_seas[N_mosq];      // Shape parameter for seasonality
	double t_peak_seas[N_mosq];     // Offset for seasonal transmission
	double denom_seas[N_mosq];      // Denominator for seasonality


	///////////////////////////////
	// Larval parameters

	double d_E_larvae;              // Development time of early larval instars
	double d_L_larvae;              // Development time of late larval instars
	double d_pupae;                 // Development time of pupae
	double mu_E0;                   // Mortality rate of early larval instars (low density)
	double mu_L0;                   // Mortality rate of late larval instars (low density)
	double mu_P;                    // Mortality rate of pupae
	double beta_larvae;             // Number of eggs laid per day per mosquito
	double gamma_larvae;            // Effect of density dependence on late instars relative to early instars

	double omega_larvae[N_mosq];    // Useful pre-calculated quantity

	double Karry[N_mosq];           // Larval carry capacity

	double eps_max[N_mosq];         // Number of eggs per day


	////////////////////////////////
	// Intervention 0 parameters
	// Treatment regimen = 0; blood-stage

	double CM0_CQ_eff;                 // chloroquine efficacy
	double CM0_CQ_proph;               // chloroquine prophylaxis


	////////////////////////////////
	// Intervention 1 parameters
	// Treatment regimen = 1; primaquine

	double CM1_CQ_eff;               // chloroquine efficacy (on its own)
	double CM1_CQ_eff_wPQ;           // chloroquine efficacy (co-administered with PQ)
	double CM1_CQ_proph;             // duration of prophylaxis of chloroquine
	double CM1_PQ_eff;               // primaquine efficacy
	double CM1_PQ_proph;             // duration of primaquine prophylaxis (number of days + 1)
	double CM1_PQ_adhere;            // adherence to full primaquine regimen
	double CM1_PQ_lowage;            // youngest age for primaquine treatment regimen
	int	   CM1_PQ_G6PD_risk;         // Risk in G6PD-deficient individuals (0 = no risk; 1 = risk)
	int	   CM1_PQ_CYP2D6_risk;       // Risk of not working in low CYP2D6 metabolizers  (0 = no risk; 1 = risk, i.e. doesn't work)
	int	   CM1_PQ_preg_risk;         // Risk in pregnant women  (0 = no risk; 1 = risk)
	int    CM1_G6PD_test;		  	 // Is G6PD testing implemented in case management (0 = no; 1 = yes)


	////////////////////////////////
	// Intervention 2 parameters
	// Treatment regimen = 2; tafenoquine

	double CM2_CQ_eff;               // chloroquine drug efficacy (on its own)
	double CM2_CQ_eff_wPQ;           // chloroquine drug efficacy (co-administered with PQ)
	double CM2_CQ_proph;             // duration of prophylaxis of chloroquine
	double CM2_PQ_eff;               // primaquine efficacy
	double CM2_PQ_proph;             // duration of primaquine prophylaxis (number of days + 1)
	double CM2_PQ_adhere;            // adherence to full primaquine regimen
	double CM2_PQ_lowage;            // youngest age for primaquine treatment regimen
	int	   CM2_PQ_G6PD_risk;         // Risk in G6PD-deficient individuals (0 = no risk; 1 = risk)
	int	   CM2_PQ_CYP2D6_risk;       // Risk of not working in low CYP2D6 metabolizers  (0 = no risk; 1 = risk, i.e. doesn't work)
	int	   CM2_PQ_preg_risk;         // Risk in pregnant women  (0 = no risk; 1 = risk)
	double CM2_TQ_eff;               // tafenoquine efficacy
	double CM2_TQ_proph;             // duration of tafenoquine prophylaxis (number of days + 1)
	double CM2_TQ_adhere;            // adherence to full tafenoquine regimen
	double CM2_TQ_lowage;            // youngest age for tafenoquine treatment regimen
	int	   CM2_TQ_G6PD_risk;         // Risk in G6PD-deficient individuals (0 = no risk; 1 = risk)
	int	   CM2_TQ_CYP2D6_risk;       // Risk of not working in low CYP2D6 metabolizers  (0 = no risk; 1 = risk, i.e. doesn't work)
	int	   CM2_TQ_preg_risk;         // Risk in pregnant women  (0 = no risk; 1 = risk)
	int    CM2_G6PD_test;		     // Is G6PD testing implemented in case management (0 = no; 1 = yes)


	////////////////////////////////
	// Intervention 5 parameters
	// MDA chloroquine parameters

	double MDA0_cover;            // Coverage of chloroquine
	double MDA0_CQ_eff;           // Efficacy of chloroquine
	double MDA0_CQ_proph;         // Duration of chloroquine prophylaxis


	////////////////////////////////
	// Intervention 6 parameters
	// MDA PQ parameters

	double MDA1_cover;            // Coverage
	double MDA1_CQ_eff;           // Efficacy of chloroquine (on its own)
	double MDA1_CQ_eff_wPQ;       // Efficacy of chloroquine (co-administered with PQ)
	double MDA1_CQ_proph;         // Duration of chloroquine prophylaxis
	double MDA1_PQ_eff;           // Efficacy of primaquine
	double MDA1_PQ_proph;         // Duration of primaquine prophylaxis
	double MDA1_PQ_adhere;        // Primaquine adherence
	double MDA1_PQ_lowage;        // Lower age limit for primaquine treatment (in days)
	int    MDA1_PQ_G6PD_risk;     // Risk in G6PD - deficient individuals
	int    MDA1_PQ_CYP2D6_risk;   // Risk of not working in low CYP2D6 metabolizers
	int    MDA1_PQ_preg_risk;     // Risk in pregnant women
	int    MDA1_G6PD_test;	      // G6PD testing for PQ eligible cases  (0 = no test; 1 = test)


	////////////////////////////////
	// Intervention 7 parameters
	// MDA TQ parameters

	double MDA2_cover;            // Coverage
	double MDA2_CQ_eff;           // Efficacy of chloroquine (on its own)
	double MDA2_CQ_eff_wPQ;       // Efficacy of chloroquine (co-administered with PQ)
	double MDA2_CQ_proph;         // Duration of chloroquine prophylaxis
	double MDA2_PQ_eff;           // Efficacy of primaquine
	double MDA2_PQ_proph;         // Duration of primaquine prophylaxis
	double MDA2_PQ_adhere;        // Primaquine adherence
	double MDA2_PQ_lowage;        // Lower age limit for primaquine treatment (in days)
	int    MDA2_PQ_G6PD_risk;     // Risk in G6PD - deficient individuals
	int    MDA2_PQ_CYP2D6_risk;   // Risk of not working in low CYP2D6 metabolizers
	int    MDA2_PQ_preg_risk;     // Risk in pregnant women
	double MDA2_TQ_eff;           // Efficacy of tafenoquine
	double MDA2_TQ_proph;         // Duration of tafenoquine prophylaxis
	double MDA2_TQ_adhere;        // Tafenoquine adherence
	double MDA2_TQ_lowage;        // Lower age limit for tafenoquine treatment (in days)
	int    MDA2_TQ_G6PD_risk;     // Risk in G6PD - deficient individuals
	int    MDA2_TQ_CYP2D6_risk;   // Risk of not working in low CYP2D6 metabolizers
	int    MDA2_TQ_preg_risk;     // Risk in pregnant women
	int    MDA2_G6PD_test;		  // G6PD testing for PQ eligible cases  (0 = no test; 1 = test)


	///////////////////////////////
	// Intervention 8 parameters
	// MSAT with chloroquine parameters

	double MSAT0_cover;           // Coverage
	double MSAT0_RDT_PCR;         // What diagnostic tool : 1 for RDT(= LM); 2 for PCR
	double MSAT0_sens;            // Sensitivity of diagnostic tool
	double MSAT0_CQ_eff;          // Efficacy of chloroquine
	double MSAT0_CQ_proph;        // Duration of chloroquine prophylaxis


	///////////////////////////////
	// Intervention 9 parameters
	// MSAT with PQ parameters

	double MSAT1_cover;            // Coverage
	double MSAT1_RDT_PCR;          // What diagnostic tool : 1 for RDT(= LM); 2 for PCR
	double MSAT1_sens;             // Sensitivity of diagnostic tool
	double MSAT1_CQ_eff;           // Efficacy of chloroquine (on its own)
	double MSAT1_CQ_eff_wPQ;       // Efficacy of chloroquine (co-administered with PQ)
	double MSAT1_CQ_proph;         // Duration of chloroquine prophylaxis
	double MSAT1_PQ_eff;           // Efficacy of primaquine
	double MSAT1_PQ_proph;         // Duration of primaquine prophylaxis
	double MSAT1_PQ_adhere;        // Primaquine adherence
	double MSAT1_PQ_lowage;        // Lower age limit for primaquine treatment (in days)
	int    MSAT1_PQ_G6PD_risk;     // Risk in G6PD - deficient individuals
	int    MSAT1_PQ_CYP2D6_risk;   // Risk of not working in low CYP2D6 metabolizers
	int    MSAT1_PQ_preg_risk;     // Risk in pregnant women
	int    MSAT1_G6PD_test;	       // G6PD testing for PQ eligible cases  (0 = no test; 1 = test)


	///////////////////////////////
	// Intervention 10 parameters
	// MSAT with PQ parameters

	double MSAT2_cover;            // Coverage
	double MSAT2_RDT_PCR;          // What diagnostic tool : 1 for RDT(= LM); 2 for PCR
	double MSAT2_sens;             // Sensitivity of diagnostic tool
	double MSAT2_CQ_eff;           // Efficacy of chloroquine (on its own)
	double MSAT2_CQ_eff_wPQ;       // Efficacy of chloroquine (co-administered with primaquine)
	double MSAT2_CQ_proph;         // Duration of chloroquine prophylaxis
	double MSAT2_PQ_eff;           // Efficacy of primaquine
	double MSAT2_PQ_proph;         // Duration of primaquine prophylaxis
	double MSAT2_PQ_adhere;        // Primaquine adherence
	double MSAT2_PQ_lowage;        // Lower age limit for primaquine treatment (in days)
	int    MSAT2_PQ_G6PD_risk;     // Risk in G6PD - deficient individuals
	int    MSAT2_PQ_CYP2D6_risk;   // Risk of not working in low CYP2D6 metabolizers
	int    MSAT2_PQ_preg_risk;     // Risk in pregnant women
	double MSAT2_TQ_eff;           // Efficacy of tafenoquine
	double MSAT2_TQ_proph;         // Duration of tafenoquine prophylaxis
	double MSAT2_TQ_adhere;        // Tafenoquine adherence
	double MSAT2_TQ_lowage;        // Lower age limit for tafenoquine treatment (in days)
	int    MSAT2_TQ_G6PD_risk;     // Risk in G6PD - deficient individuals
	int    MSAT2_TQ_CYP2D6_risk;   // Risk of not working in low CYP2D6 metabolizers
	int    MSAT2_TQ_preg_risk;     // Risk in pregnant women
	int    MSAT2_G6PD_test;	       // G6PD testing for PQ eligible cases  (0 = no test; 1 = test)


	///////////////////////////////
	// Intervention 11 parameters
	// STAT with PQ parameters

	double STAT1_cover;            // Coverage
	double STAT1_sens;             // Sensitivity of diagnostic tool
	double STAT1_spec;             // Sensitivity of diagnostic tool
	double STAT1_RDT_PCR;          // What diagnostic tool (extra parasitological testing)): 0 = nothing; 1 for RDT(= LM); 2 for PCR
	double STAT1_CQ_eff;           // Efficacy of chloroquine drugs (on its own)
	double STAT1_CQ_eff_wPQ;       // Efficacy of chloroquine drugs (co-administered with PQ)
	double STAT1_CQ_proph;         // Duration of chloroquine prophylaxis
	double STAT1_PQ_eff;           // Efficacy of primaquine
	double STAT1_PQ_proph;         // Duration of primaquine prophylaxis
	double STAT1_PQ_adhere;        // Primaquine adherence
	double STAT1_PQ_lowage;        // Lower age limit for primaquine treatment (in days)
	int    STAT1_PQ_G6PD_risk;     // Risk in G6PD - deficient individuals
	int    STAT1_PQ_CYP2D6_risk;   // Risk of not working in low CYP2D6 metabolizers
	int    STAT1_PQ_preg_risk;     // Risk in pregnant women
	int    STAT1_G6PD_test;	       // G6PD testing for PQ eligible cases  (0 = no test; 1 = test)


	///////////////////////////////
	// Intervention 12 parameters
	// STAT with TQ parameters

	double STAT2_cover;            // Coverage
	double STAT2_sens;             // Sensitivity of diagnostic tool
	double STAT2_spec;             // Sensitivity of diagnostic tool
	double STAT2_RDT_PCR;          // What diagnostic tool (extra parasitological testing)): 0 = nothing; 1 for RDT(= LM); 2 for PCR
	double STAT2_CQ_eff;           // Efficacy of chloroquine (on its own)
	double STAT2_CQ_eff_wPQ;       // Efficacy of chloroquine (co-administered with PQ)
	double STAT2_CQ_proph;         // Duration of chloroquine prophylaxis
	double STAT2_PQ_eff;           // Efficacy of primaquine
	double STAT2_PQ_proph;         // Duration of primaquine prophylaxis
	double STAT2_PQ_adhere;        // Primaquine adherence
	double STAT2_PQ_lowage;        // Lower age limit for primaquine treatment (in days)
	int    STAT2_PQ_G6PD_risk;     // Risk in G6PD - deficient individuals
	int    STAT2_PQ_CYP2D6_risk;   // Risk of not working in low CYP2D6 metabolizers
	int    STAT2_PQ_preg_risk;     // Risk in pregnant women
	double STAT2_TQ_eff;           // Efficacy of tafenoquine
	double STAT2_TQ_proph;         // Duration of tafenoquine prophylaxis
	double STAT2_TQ_adhere;        // Tafenoquine adherence
	double STAT2_TQ_lowage;        // Lower age limit for tafenoquine treatment (in days)
	int    STAT2_TQ_G6PD_risk;     // Risk in G6PD - deficient individuals
	int    STAT2_TQ_CYP2D6_risk;   // Risk of not working in low CYP2D6 metabolizers
	int    STAT2_TQ_preg_risk;     // Risk in pregnant women
	int    STAT2_G6PD_test;	       // G6PD testing for PQ eligible cases  (0 = no test; 1 = test)


	/////////////////////////////////////////////////////////
	// Interventions 3 & 4 parameters
	// Baseline entomological parameters required for LLINs

	double LLIN_half_life;      // Half-life of loss of LLINs
	double P_LLIN_loss;         // Daily probability of losing LLIN - pre-calculated for efficiency
	double PYR_half_life;       // Half-life of pyrethroid decay on LLINs
	double PYR_decay;           // Daily pyrethroid decay - pre-calculated for efficiency

	double r_LLIN_0[N_mosq];    // Probability mosquito repelled (with full insecticide activity)
	double r_LLIN_net[N_mosq];  // Probability mosquito repelled due to barrier effect of net (no insecticide)
	double s_LLIN_0[N_mosq];    // Probability mosquito feeds successfully
	double d_LLIN_0[N_mosq];    // Probability mosquito dies during feeding attempt

	double IRS_half_life;       // Half-life of IRS insecticide decay
	double IRS_decay;           // Daily IRS insecticide decay - pre-calculated for efficiency

	double r_IRS_0[N_mosq];     // Probability mosquito repelled by IRS (full insecticide activity)
	double d_IRS_0[N_mosq];     // Probability mosquito killed by IRS (full insecticide activity)
	double s_IRS_0[N_mosq];     // Probability mosquito survives feeding attempt with IRS (full insecticide activity)

	double Q_0[N_mosq];         // Human Blood Index (proportion of blood meals taken on humans)
	double CHI_endo[N_mosq];	// Endophily - proportion of mosquitoes resting indoors after feeding (no intervention)
	double PSI_indoors[N_mosq]; // Proportion of bites taken on humans indoors
	double PSI_bed[N_mosq];     // Proportion of bites taken on humans in bed

	double delta_1;	            // Time spent foraging for a blood meal
	double delta_2;	            // Time_spent_digesting_blood_meal
	double delta;	            // Duration of gonotrophic cycle

	double p_1[N_mosq];         // Daily? death probability when foraging for a blood meal
	double p_2[N_mosq];         // Daily? death probability when digesting blood meal


	/////////////////////////////////////////////////////////
	// Intervention 13 parameters
	// Ivermectin parameters

	double d_IVM_0;          // Probability mosquito dies during feeding attempt
	double IVM_half_life;    // Half-life of loss of LLINs

	double P_IVM_decay;      // Daily pyrethroid decay - pre-calculated for efficiency


	//////////////////////////////////////////////////////
	// Pre-multiplication of quantities for efficiency

	double A_par_decay;         // Decay factor for BS immunity
	double A_clin_decay;        // Decay factor for clinical immunity
	double mat_decay;           // Decay factor for maternal immunity

	double age_0_inv;           // Inverse of age-dependent biting parameter

	double A_PCR_50pc_inv;      // Immune scalar for PCR-detectable infection
	double A_LM_50pc_inv;       // Immune scalar for LM-detectable infection
	double A_D_50pc_inv;        // Immune scalar for clinical disease

	double P_dead;              // Probability of dying in each time step
	double Preg_daily;          // Probability of woman 18-40 years becoming pregnant per day

	double P_PYR_decay;         // Proportional decay of pyrethroid per time step (pre-calculated for efficiency)
	double P_IRS_decay;         // Proportional decay of IRS insecticide per time step (pre-calculated for efficiency)


	////////////////////////////////////////
	// Matrices for hypnozoite transitions

	double D_MAT[K_max + 1][K_max + 1];
	double OD_MAT[K_max + 1][K_max + 1];
	double K_MAT[K_max + 1][K_max + 1];
	double L_MAT[K_max + 1][K_max + 1];
	double H_MAT[K_max + 1][K_max + 1];
};

#endif
