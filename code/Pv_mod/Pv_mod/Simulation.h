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
///  7. MAIN - SIMULATION                                                 ///
///     Here we read in parameters from files (model parameters,          ///
///     mosquito parameters, intervention parameters).                    ///
///     A population of individuals is created at equilibrium and         ///
///     then simulated.                                                   ///
///                                                                       ///
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

// Include guard
#ifndef PVIVAX_MODEL_SIMULATION
#define PVIVAX_MODEL_SIMULATION

#include "Intervention.h"


/////////////////////////////////////////////////////////////////////////////////////////
//                                                                                     //
// 7.1.1. Define a class for storing the output of a simulation                        //
//                                                                                     //
/////////////////////////////////////////////////////////////////////////////////////////

class Simulation
{
public:
	//////////////////////////////////////////////////////////////////////////
	//  7.1.1.1.  Functions                                                 //
	//////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////
	// Set up simulation times and storage for outputs

	Simulation(SimTimes times);


	/////////////////////////////////////
	// Simulate the model and store the output

	void run(Params& theta, Population& POP, Intervention& INTVEN);


	/////////////////////////////////////
	// Write output file

	void write_output(const char *output_File);


	//////////////////////////////////////////////////////////////////////////
	//  7.1.1.2.  Data                                                      //
	//////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////
	// Vector of simulation times

	int N_time;
	SimTimes times;

	vector<double> t_vec;


	//////////////////////////////////////////
	// Tracking output

	vector<vector<int>> yH_t;
	vector<vector<vector<double>>> yM_t;


	vector<double> EIR_dom_t;
	vector<double> EIR_occ_t;


	vector<vector<int>> prev_all;   // Contains {N_pop, PvPR_PCR, PvPR_LM, Pv_clin, PvHR, PvHR_batches, new_PCR, new_LM, new_D, new_T}
	vector<vector<int>> prev_U5;    // Contains {N_pop, PvPR_PCR, PvPR_LM, Pv_clin, PvHR, PvHR_batches, new_PCR, new_LM, new_D, new_T}
	vector<vector<int>> prev_U10;   // Contains {N_pop, PvPR_PCR, PvPR_LM, Pv_clin, PvHR, PvHR_batches, new_PCR, new_LM, new_D, new_T}


	//////////////////////////////////////////
	// Tracking coverage over time

	vector<int> LLIN_cov_t;
	vector<int> IRS_cov_t;
	vector<int> IVM_cov_t;
	vector<int> CQ_treat_t;
	vector<int> PQ_treat_t;
	vector<int> TQ_treat_t;
	vector<int> pregnant_t;


	vector<int> PvHR_t;
	vector<int> PvHR_batch_t;
	vector<int> I_PCR_new_t;
	vector<int> I_LM_new_t;
	vector<int> I_D_new_t;
	vector<int> T_new_t;
	vector<int> G6PD_tests_t;
	vector<int> PQ_effective_t;
	vector<int> TQ_effective_t;
	vector<int> CQ_relapse_t;
	vector<int> CQ_bite_t;
	vector<int> PQ_G6PD_def_t;
	vector<int> TQ_G6PD_def_t;


	vector<int> PQ_overtreat_t;
	vector<int> PQ_overtreat_9m_t;

	vector<int> TQ_overtreat_t;
	vector<int> TQ_overtreat_9m_t;

	vector<int> cases_M_O16_t;
	vector<int> cases_M_U16_t;
	vector<int> cases_F_O16_t;
	vector<int> cases_F_U16_t;
	vector<int> cases_preg_t;


	//////////////////////////////////////////
	// Tracking immunity over time

	vector<double> A_par_mean_t;
	vector<double> A_clin_mean_t;
};

#endif
