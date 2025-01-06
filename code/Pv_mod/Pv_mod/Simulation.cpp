/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
///                                                                       ///
///  Individual-based Plasmodium vivax transmission model.                ///
///                                                                       ///
///  Dr Michael White                                                     ///
///  Institut Pasteur                                                     ///
///  michael.white@pasteur.fr                                             ///
///                                                                       ///
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

#include "Simulation.h"

#include <iostream>
#include <fstream>
#include <cmath>


////////////////////////////////////////////////////////////
//                                                        //
//  7.2.1.  Function declarations                         //
//                                                        //
////////////////////////////////////////////////////////////

// Defined in Mosquito.cpp
void mosquito_step(double t, Params& theta, Population& POP);


Simulation::Simulation(SimTimes times) :
	times(times)
{
	/////////////////////////////////////////////////////////////////////////
	// 7.2.1.1.  Vector of simulation times

	// Number of time steps for simulation:
	N_time = (1 / t_step)*(times.burnin + times.end - times.start) * 365;

	for (int i = 0; i < N_time; i++)
	{
		t_vec.push_back((double)(times.start * 365 - times.burnin * 365 + i * t_step));
	}


	/////////////////////////////////////////////////////////////////////////
	// 7.2.1.2. Create storage for output

	yH_t.resize(N_time);
	for (int i = 0; i < N_time; i++)
	{
		yH_t[i].resize(N_H_comp);
	}


	yM_t.resize(N_time);
	for (int i = 0; i < N_time; i++)
	{
		yM_t[i].resize(N_mosq);
		for (int v = 0; v < N_mosq; v++)
		{
			yM_t[i][v].resize(N_M_comp);
		}
	}


	// prev_all.resize(N_time);
	// for (int i = 0; i < N_time; i++)
	// {
	// 	prev_all[i].resize(15);
	// }
  //
	// prev_U5.resize(N_time);
	// for (int i = 0; i < N_time; i++)
	// {
	// 	prev_U5[i].resize(15);
	// }
	//
	// prev_U10.resize(N_time);
	// for (int i = 0; i < N_time; i++)
	// {
	// 	prev_U10[i].resize(15);
	// }

	EIR_dom_t.resize(N_time);
	EIR_occ_t.resize(N_time);

	LLIN_cov_t.resize(N_time);
	IRS_cov_t.resize(N_time);
	IVM_cov_t.resize(N_time);
	CQ_treat_t.resize(N_time);
	PQ_treat_t.resize(N_time);
	TQ_treat_t.resize(N_time);
	pregnant_t.resize(N_time);


	PvHR_t.resize(N_time);
	PvHR_batch_t.resize(N_time);
	I_PCR_new_t.resize(N_time);
	I_LM_new_t.resize(N_time);
	I_D_new_t.resize(N_time);
	T_new_t.resize(N_time);
	G6PD_tests_t.resize(N_time);
	CQ_relapse_t.resize(N_time);
	CQ_bite_t.resize(N_time);
	PQ_effective_t.resize(N_time);
	TQ_effective_t.resize(N_time);
	PQ_G6PD_def_t.resize(N_time);
	TQ_G6PD_def_t.resize(N_time);


	PQ_overtreat_t.resize(N_time);
	PQ_overtreat_9m_t.resize(N_time);

	TQ_overtreat_t.resize(N_time);
	TQ_overtreat_9m_t.resize(N_time);

	cases_M_O16_t.resize(N_time);
	cases_M_U16_t.resize(N_time);
	cases_F_O16_t.resize(N_time);
	cases_F_U16_t.resize(N_time);
	cases_preg_t.resize(N_time);

	A_par_mean_t.resize(N_time);
	A_clin_mean_t.resize(N_time);
}


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  7.2.2. Simulate the model and store the output in SIM                   //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

void Simulation::run(Params& theta, Population& POP, Intervention& INTVEN)
{

	for (int i = 0; i < N_time; i++)
	{
		if (t_vec[i] / 365.0 - floor(t_vec[i] / 365.0) < 0.5*t_step / 365.0)
		{
			cout << "time = " << t_vec[i] / 365.0 << "\t" << 100.0*(t_vec[i] - t_vec[0]) / (double(t_step*N_time)) << "% complete" << endl;
		}

		POP.human_step(theta);

		mosquito_step(t_vec[i], theta, POP);

		INTVEN.distribute(t_vec[i], theta, POP);

		POP.summary();

		//////////////////////////////////////
		// Fill out Simulation object

		for (int k = 0; k < N_H_comp; k++)
		{
			yH_t[i][k] = POP.yH[k];
		}

		for (int k = 0; k < N_M_comp; k++)
		{
			for (int v = 0; v < N_mosq; v++)
			{
				yM_t[i][v][k] = POP.yM[v][k];
			}
		}

		// for (int k = 0; k < 15; k++)
		// {
		// 	prev_all[i][k] = POP.prev_all[k];
		// 	prev_U5[i][k] = POP.prev_U5[k];
		// 	prev_U10[i][k] = POP.prev_U10[k];
		// }


		// LLIN_cov_t[i] = POP.LLIN_cov_t;
		// IRS_cov_t[i]  = POP.IRS_cov_t;
		// IVM_cov_t[i]  = POP.IVM_cov_t;

		PvHR_t[i]         = POP.PvHR_t;
		PvHR_batch_t[i]   = POP.PvHR_batch_t;
		I_PCR_new_t[i]    = POP.I_PCR_new_t;
		I_LM_new_t[i]     = POP.I_LM_new_t;
		I_D_new_t[i]      = POP.I_D_new_t;
		T_new_t[i]        = POP.T_new_t;
		CQ_relapse_t[i]   = POP.CQ_relapse_t;
		CQ_bite_t[i]      = POP.CQ_bite_t;
		G6PD_tests_t[i]   = POP.G6PD_tests_t;
		PQ_effective_t[i] = POP.PQ_effective_t;
		TQ_effective_t[i] = POP.TQ_effective_t;
		PQ_G6PD_def_t[i]  = POP.PQ_G6PD_def_t;
		TQ_G6PD_def_t[i]  = POP.TQ_G6PD_def_t;

		CQ_treat_t[i] = POP.CQ_treat_t;
		PQ_treat_t[i] = POP.PQ_treat_t;
		TQ_treat_t[i] = POP.TQ_treat_t;
		// pregnant_t[i] = POP.pregnant_t;
		//
		// PQ_overtreat_t[i]    = POP.PQ_overtreat_t;
		// PQ_overtreat_9m_t[i] = POP.PQ_overtreat_9m_t;
		//
		// TQ_overtreat_t[i]    = POP.TQ_overtreat_t;
		// TQ_overtreat_9m_t[i] = POP.TQ_overtreat_9m_t;
		//
		EIR_dom_t[i] = EIR_dom_t[i] + POP.aa_VC[0] * POP.yM[0][5];
		EIR_occ_t[i] = EIR_occ_t[i] + POP.aa_VC[1] * POP.yM[1][5];
		//
		// cases_M_O16_t[i] = POP.cases_M_O16_t;
		// cases_M_U16_t[i] = POP.cases_M_U16_t;
		// cases_F_O16_t[i] = POP.cases_F_O16_t;
		// cases_F_U16_t[i] = POP.cases_F_U16_t;
		// cases_preg_t[i]  = POP.cases_preg_t;

		A_par_mean_t[i] = POP.A_par_mean_t;
		A_clin_mean_t[i] = POP.A_clin_mean_t;

	}
}

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  7.2.3. Write output to file                                             //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

void Simulation::write_output(const char *output_File)
{
	cout << "Start writing output to file......" << endl;
	cout << endl;

	ofstream output_Stream(output_File);

	for (int i = (int)(1 / t_step)*(times.burnin) * 365; i < N_time; i++)
	{

	  // days
		output_Stream << t_vec[i] << "\t";


	  // mosquitoes
	  for (int v = 0; v < N_mosq; v++)
	  {
	    for (int k = 0; k < N_M_comp; k++)
	    {
	      output_Stream << yM_t[i][v][k] << "\t";
	    }
	  }

		output_Stream << EIR_dom_t[i] << "\t";
		output_Stream << EIR_occ_t[i] << "\t";


		// humans
		for (int k = 0; k < N_H_comp; k++)
		{
		  output_Stream << yH_t[i][k] << "\t";
		}

		// output_Stream << LLIN_cov_t[i] << "\t";
		// output_Stream << IRS_cov_t[i] << "\t";
		// output_Stream << IVM_cov_t[i] << "\t";
		output_Stream << PvHR_t[i]         << "\t";
		output_Stream << PvHR_batch_t[i]   << "\t";
		output_Stream << I_PCR_new_t[i]    << "\t";
		output_Stream << I_LM_new_t[i]     << "\t";
		output_Stream << I_D_new_t[i]      << "\t";
		output_Stream << T_new_t[i]        << "\t";
		output_Stream << CQ_treat_t[i]     << "\t";
		output_Stream << CQ_relapse_t[i]   << "\t";
		output_Stream << CQ_bite_t[i]      << "\t";
		output_Stream << G6PD_tests_t[i]   << "\t";
		output_Stream << PQ_treat_t[i]     << "\t";
		output_Stream << PQ_effective_t[i] << "\t";
		output_Stream << PQ_G6PD_def_t[i]  << "\t";
		output_Stream << TQ_treat_t[i]     << "\t";
		output_Stream << TQ_effective_t[i] << "\t";
		output_Stream << TQ_G6PD_def_t[i]  << "\t";


		// for (int k = 0; k < 15; k++)
		// {
		//   output_Stream << prev_all[i][k] << "\t";
		// }

		// Write output for age categories U5 and U10
		/*for (int k = 0; k<15; k++)
		 {
		 output_Stream << prev_U5[i][k] << "\t";
		 }
		 for (int k = 0; k<15; k++)
		 {
		 output_Stream << prev_U10[i][k] << "\t";
		 }*/

		// output_Stream << PQ_overtreat_t[i] << "\t";
		// output_Stream << PQ_overtreat_9m_t[i] << "\t";
		//
		// output_Stream << TQ_overtreat_t[i] << "\t";
		// output_Stream << TQ_overtreat_9m_t[i] << "\t";

		// // Write number of pregnant women
		// output_Stream << pregnant_t[i] << "\t";
		//
		// output_Stream << cases_M_O16_t[i] << "\t";
		// output_Stream << cases_M_U16_t[i] << "\t";
		// output_Stream << cases_F_O16_t[i] << "\t";
		// output_Stream << cases_F_U16_t[i] << "\t";
		// output_Stream << cases_preg_t[i] << "\t";

		// Write A_par_mean_t and A_clin_mean_t
		output_Stream << A_par_mean_t[i] << "\t";
		output_Stream << A_clin_mean_t[i];

		output_Stream << endl;
	}

	output_Stream.close();


	cout << "Output successfully written to file......" << endl;
	cout << endl;
}
