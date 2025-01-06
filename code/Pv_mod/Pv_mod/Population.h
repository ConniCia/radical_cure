/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
///                                                                       ///
///  Individual-based Plasmodium vivax transmission model.                ///
///                                                                       ///
///  Dr Michael White                                                     ///
///  Institut Pasteur                                                     ///
///  michael.white@pasteur.fr                                             ///
///                                                                       ///
///  This file and accompanying .cpp contain:                             ///
///                                                                       ///
///  A class called Population stores all individuals.                    ///
///                                                                       ///
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

// Include guard
#ifndef PVIVAX_MODEL_POPULATION
#define PVIVAX_MODEL_POPULATION

#include "Individual.h"


/////////////////////////////////////////////////////////////////////////////////////////
//                                                                                     //
// 5.1.1. Define a structure for a population of individuals                           //
//                                                                                     //
//        Note that this object only stores information at a fixed point in time       //
//        and as such is memoryless                                                    //
//                                                                                     //
/////////////////////////////////////////////////////////////////////////////////////////

class Population
{
public:
	//////////////////////////////////////////////////////////////////////////
	//  5.1.1.1.  Class member functions                                    //
	//////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////
	// Initialise a population of individuals at equilibrium

	void pop_setup(Params& theta);

	void pop_at_equil(Params& theta);

	void ind_at_equil(Params& theta);


	////////////////////////////////////////////////////
	// Update human individuals

	void human_step(Params& theta);

	////////////////////////////////////////////////////
	// Summarise population outputs

	void summary();


private:
	void gauher(Params& theta);

public:
	//////////////////////////////////////////////////////////////////////////
	//  5.1.1.2.  Data                                                      //
	//////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////
	// Human population

	int N_pop;                      // Population size - we have balanced demography at the moment so this will be effectively constant

	vector<Individual> people;      // A vector of individuals

	vector<vector<double>> pi_n;    // Proportion of bites on humans that person n receives
	vector<vector<double>> lam_n;   // Biting rate of a single mosquito on person n


	///////////////////////////////////////////////
	//  Mosquito population
	//
	//  Depends dynamically on vector control

	double yM[N_mosq][N_M_comp];      // mosquito state

	double SUM_pi_w[N_mosq];
	double SUM_pi_z[N_mosq];


	double delta_1_VC[N_mosq];        // duration spent foraging
	double delta_VC[N_mosq];          // duration of gonotrophic cycle = duration between blood meals

	double Z_VC[N_mosq];              // average probability of mosquito repeating during a single attempt
	double W_VC[N_mosq];              // average probability of successfully feeding on a human during a single attempt
	double Q_VC[N_mosq];              // human blood index
	double p_1_VC[N_mosq];            // probability of surviving foraging
	double mu_M_VC[N_mosq];           // mosquito death rate
	double aa_VC[N_mosq];             // mosquito biting rate
	double exp_muM_tauM_VC[N_mosq];   // probability of surviving sporogony
	double beta_VC[N_mosq];           // egg oviposition rate


	////////////////////////////////////////////////////////////////
	//  Objects for storing summary output of the population

	int yH[N_H_comp];   // Human compartmental states

	int prev_all[15];   // Contains {N_pop, PvPR_PCR, PvPR_LM, Pv_clin, PvHR, PvHR_batches, new_PCR, new_LM, new_D, new_CQ, new_PQ, new_TQ, G6PD_test, PQ_effective, TQ_effective}
	int prev_U5[15];    // Contains {N_pop, PvPR_PCR, PvPR_LM, Pv_clin, PvHR, PvHR_batches, new_PCR, new_LM, new_D, new_CQ, new_PQ, new_TQ, G6PD_test, PQ_effective, TQ_effective}
	int prev_U10[15];   // Contains {N_pop, PvPR_PCR, PvPR_LM, Pv_clin, PvHR, PvHR_batches, new_PCR, new_LM, new_D, new_CQ, new_PQ, new_TQ, G6PD_test, PQ_effective, TQ_effective}

	double EIR_dom_t;        // EIR (domestic)
	double EIR_occ_t;        // EIR (occupational)
	int LLIN_cov_t;          // LLIN coverage
	int IRS_cov_t;           // IRS coverage
	int IVM_cov_t;           // IVM coverage
	int CQ_treat_t;          // Indicator for chloroquine (or other blood-stage drug) treatment
	int PQ_treat_t;          // Indicator for primaquine treatment
	int TQ_treat_t;          // Indicator for tafenoquine treatment
	int pregnant_t;          // Coverage with front-line treatment (primaquine or tafenoquine)


	int PvHR_t;
	int PvHR_batch_t;
	int I_PCR_new_t;
	int I_LM_new_t;
	int I_D_new_t;
	int T_new_t;
	int G6PD_tests_t;
	int CQ_relapse_t;
	int CQ_bite_t;
	int PQ_G6PD_def_t;
	int TQ_G6PD_def_t;


	int PQ_effective_t;     // Effective 8-aminoquinoline treatment (PQ)
	int PQ_overtreat_t;     // Over-treatment with 8-aminoquinolines (PQ) - defined as someone without hypnozoites being treated
	int PQ_overtreat_9m_t;  // Over-treatment with 8-aminoquinolines (PQ) - defined as someone without BS infection in last 9 mths being treated

	int TQ_effective_t;     // Effective 8-aminoquinoline treatment (TQ)
	int TQ_overtreat_t;     // Over-treatment with 8-aminoquinolines (TQ) - defined as someone without hypnozoites being treated
	int TQ_overtreat_9m_t;  // Over-treatment with 8-aminoquinolines (TQ) - defined as someone without BS infection in last 9 mths being treated

	int cases_M_O16_t;       // Detected cases in males over 16
	int cases_M_U16_t;       // Detected cases in males under 16
	int cases_F_O16_t;       // Detected cases in females over 16
	int cases_F_U16_t;       // Detected cases in females under 16
	int cases_preg_t;        // Detected cases in pregnant women

	double A_par_mean_t;     // Average anti-parasite immunity
	double A_clin_mean_t;    // Average anti-clinical immunity


	///////////////////////////////////////////////////////////////
	//  Equilibrium settings of population
	//  These are only needed for initialising the simulation

	////////////////////////////////////////
	// Age and heterogeneity demographics

	double age_bounds[N_age + 1];

	double age_demog[N_age];
	double age_bite[N_age];
	double age_mids[N_age];

	double x_het_bounds[N_het + 1];

	double x_het[N_het];
	double w_het[N_het];

	double x_age_het[N_age][N_het];
	double w_age_het[N_age][N_het];

	double w_gen_het_at_birth[N_gen][N_het];

	double x_gen_age_het[N_mosq][N_gen][N_age][N_het];
	double w_gen_age_het[N_mosq][N_gen][N_age][N_het];


	////////////////////////////////////////
	// Ageing rates

	double r_age[N_age];

	double P_age_bite;     // Proportional change for age-dependent biting


	////////////////////////////////////////
	// Equilibrium vectors

	vector<vector<vector<vector<vector<double>>>>> yH_eq;
	vector<vector<vector<double>>> lam_eq;
	vector<vector<vector<vector<double>>>> A_par_eq;
	vector<vector<vector<double>>> A_par_eq_mean;
	vector<vector<vector<vector<double>>>> A_clin_eq;
	vector<vector<vector<double>>> A_clin_eq_mean;
	vector<vector<vector<vector<double>>>> phi_LM_eq;
	vector<vector<vector<vector<double>>>> phi_D_eq;
	vector<vector<vector<vector<double>>>> r_PCR_eq;

	vector<vector<vector<vector<double>>>> yH_eq_cumsum;

	double denom_w_occup;


	////////////////////////////////////////
	// Equilibrium moves to zero HPZ due to PQ

	int equil_setup_count;

	vector<vector<vector<vector<vector<double>>>>> theta_HPZzero;

	vector<vector<vector<vector<double>>>> theta_HPZzero_weight;


	////////////////////////////////////////
	// Treatment effectiveness

	vector<vector<double>> treat_PQeff;


	////////////////////////////////////////
	// Indicator for age group of 20 year old woman

	int index_age_20;
	int index_preg_age_low;
	int index_preg_age_high;
	int index_occup_age_low;
	int index_occup_age_high;
};

#endif
