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
///  4. INDIVIDUAL-BASED MODEL                                            ///
///     Details of the stochastic individual-based model for each         ///
///     person. Transitions occur with a fixed time step according to     ///
///     compting hazards                                                  ///
///                                                                       ///
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

// Include guard
#ifndef PVIVAX_MODEL_INDIVIDUAL
#define PVIVAX_MODEL_INDIVIDUAL

#include "Params.h"


/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// 4.1.1.  Define a class for humans                                       //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

class Individual
{
public:
	//////////////////////////////////////////////////////////////////////////
	//  Class constructors and destructors
	//////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////
	// 4.1.1.1. Class constructor

	Individual(double a, double zeta)
	{
		age = a;
		zeta_het = zeta;
	}


	////////////////////////////////////////////////////
	// Copy and move constructors

	// Delete unwanted copy constructors
	Individual(Individual&) = delete;
	void operator= (Individual&) = delete;

	// Allow default move constructors
	Individual(Individual&&) = default;
	Individual& operator= (Individual&&) = default;


	//////////////////////////////////////////////////////////////////////////
	//  Class member functions
	//////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////
	// 4.1.1.2. Function declarations within the human class

	void state_mover(Params& theta, double lam_bite);
	void ager(Params& theta);
	void case_management(Params& theta);
	void vec_con_updater(Params& theta);

	//////////////////////////////////////////////////////////////////////////
	//  Data
	//////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////
	// 4.1.1.4. Person-specific age and exposure classifiers

	double age;                      // Person's age
	double zeta_het;                 // Heterogeneity in exposure to mosquitoes

	int gender;                      // 0 = male (domestic); 1 = male (occupational); 2 = female
	int G6PD_def;                    // Is the person G6PD deficient? 0 = normal; 1 = deficient (homozygous); 2 = deficient (heterozygous - female only)
	double G6PD_activity;            // G6PD activity - quantitative measurement (true value)
	double G6PD_read;                // G6PD activity - quantitative measurement (from dx)
	bool CYP2D6;                     // Does the person have CYP2D6 phenotype? 0 = normal; 1 = low metabolizer


	////////////////////////////////////////////////////
	//  4.1.1.4. Child-bearing age.
	//           Indicator for people of 18-22 years of age. New-born
	//           children acquire a proportion of the immunity of a
	//           women aged 20 years.

	bool preg_age;
	bool pregnant;
	double preg_timer;


	double lam_bite_lag;             // Lagged force of infection due to moquito bites
	vector<double> lam_bite_track;   // Tracking lagged force of infection due to moquito bites

	double lam_rel_lag;              // Lagged force of infection due to relapses
	vector<double> lam_rel_track;    // Tracking lagged force of infection due to relapses

	double lam_H_lag;                // Lagged total force of infection


	////////////////////////////////////////////////////
	// 4.1.1.5. Indicators for compartments

	bool S;
	bool I_PCR;
	bool I_LM;
	bool I_D;
	bool T;
	bool P;

	/////////////////////////////////////////////////////////////////
	//  4.1.1.6.  Number of batches of hypnozoites. Must be an integer.

	int Hyp;

	////////////////////////////////////////////////////
	// Indicator for competing hazards move

	int CH_move;


	////////////////////////////////////////////////////
	//  4.1.1.7.  Indicators for new events

	bool I_PCR_new;         // New PCR-detectable infection (I_LM, I_D & T included here)
	bool I_LM_new;          // New LM-detectable infection (I_D & T included here)
	bool I_D_new;           // New clinical episode (untreated)
	bool T_new;             // New clinical episode (treated)
	bool CQ_treat;          // Indicator for chloroquine (or other blood-stage drug) treatment
	bool PQ_treat;          // Indicator for primaquine treatment
	bool TQ_treat;          // Indicator for tafenoquine treatment

	bool G6PD_test;			// New G6PD test administered

	bool CQ_effective;     // Effective blood-stage treatment (CQ)

	bool PQ_effective;     // Effective 8-aminoquinoline treatment (PQ)
	bool PQ_overtreat;     // Over-treatment with 8-aminoquinolines (PQ) - defined as someone without hypnozoites being treated
	bool PQ_overtreat_9m;  // Over-treatment with 8-aminoquinolines (PQ) - defined as someone without BS infection in last 9 mths being treated

	bool TQ_effective;     // Effective 8-aminoquinoline treatment (TQ)
	bool TQ_overtreat;     // Over-treatment with 8-aminoquinolines (TQ) - defined as someone without hypnozoites being treated
	bool TQ_overtreat_9m;  // Over-treatment with 8-aminoquinolines (TQ) - defined as someone without BS infection in last 9 mths being treated


	bool CQ_relapse;
	bool CQ_bite;
	bool PQ_G6PD_def;
	bool TQ_G6PD_def;

	////////////////////////////////////////////////////
	//  4.1.1.8.  Person-specific levels of immunity and indicators
	//            for whether immunity is suppressed

	double A_par;
	double A_clin;

	double A_par_mat;
	double A_clin_mat;

	bool A_par_boost;
	bool A_clin_boost;

	double A_par_timer;
	double A_clin_timer;

	bool AQ8_proph;
	double AQ8_proph_timer;


	//////////////////////////////////////////////////////////
	// 4.1.1.9. Person-specific intervention access parameter
	//
	// Note that while there are 9 interventions in total, 2 of these are
	// are related to first-line treatment. Therefore there are N_int = 6
	// 'pulsed' interventions, i.e. those distributed by campaigns where
	// access will be an issue.
	//
	// Of course there will be differential access to first-line treatment,
	// but that's a story for another day.

	double zz_int[N_int];


	///////////////////////////////
	// LLINs

	bool LLIN;                 // Does the person have an LLIN?
	double LLIN_age;           // How old is the LLIN

	double r_LLIN[N_mosq];     // repellency
	double d_LLIN[N_mosq];     // killing effect
	double s_LLIN[N_mosq];     // survival


	///////////////////////////////
	// IRS

	bool IRS;                 // Is the person protected by IRS
	double IRS_age;           // How long ago was their house sprayed?

	double r_IRS[N_mosq];     // repellency
	double d_IRS[N_mosq];     // killing effect
	double s_IRS[N_mosq];     // survival


	///////////////////////////////
	// Ivermectin (IVM)

	bool IVM;         // Is the person given IVM
	double IVM_age;   // How long ago were they given IVM?

	double d_IVM;     // killing effect
	double s_IVM;     // survival


	///////////////////////////////
	// STAT

	double T_last_BS;         // Tracking of time since last PCR-detectable blood-stage infection

	///////////////////////////////////////////////////
	// Individual-level effect of vector control

	double z_VC[N_mosq];      // probability of mosquito being repelled from this individual during a single feeding attempt
	double y_VC[N_mosq];      // probability of mosquito feeding on this individual during a single attempt
	double w_VC[N_mosq];      // probability of mosquito feeding and surviving on this individual during a single feeding attempt
};

#endif
