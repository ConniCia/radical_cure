/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
///                                                                       ///
///  Individual-based Plasmodium vivax transmission model.                ///
///                                                                       ///
///  Dr Michael White                                                     ///
///  Institut Pasteur                                                     ///
///  michael.white@pasteur.fr                                             ///
///                                                                       ///
///  EQUILIBRIUM SETUP                                                    ///
///  This set of functions calculates the equilibrium set up of the       ///
///  population. It is only called once while the population is           ///
///  being initialised.                                                   ///
///                                                                       ///
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

#include "Population.h"

#include <iostream>
#include <cmath>
#include "randlib.h"


////////////////////////////////////////////////////////////
//                                                        //
//  5.2.1.  Function declarations                         //
//                                                        //
////////////////////////////////////////////////////////////

void ludcmp(vector<vector<double>> &a, int n_dim, vector<int> &indx, double &d);
void lubksb(vector<vector<double>> &a, int n_dim, vector<int> &indx, vector<double> &b);
void matrix_inv(vector<vector<double>> &a, int n, vector<vector<double>> &a_inv);
void inv_MM_bb(vector<vector<double>> &MM, vector<double> &bb, vector<double> &xx, int n_dim);
void MM_ij( int g, int i, int j, Params& theta, double r_age[], vector<vector<double>> &MM,
	        vector<vector<vector<double>>> lam_eq, vector<vector<vector<vector<double>>>> phi_LM_eq,
	        vector<vector<vector<vector<double>>>> phi_D_eq, vector<vector<vector<vector<double>>>> r_PCR_eq,
	        vector<vector<double>> treat_PQeff);


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  5.2.2.  Update the vector of human classes                              //
//                                                                          //
//          THINK CAREFULLY ABOUT THE ORDERING OF EVENTS                    //
//////////////////////////////////////////////////////////////////////////////

void Population::human_step(Params& theta)
{

	//////////////////////////////////////////////////////////////////////////
	// 5.2.2.1. Temporary objects for setting up individuals' intervention
	//          access characteristics

	float GMN_parm[(N_int)*(N_int + 3) / 2 + 1];
	float GMN_work[N_int];
	float GMN_zero[N_int];
	float zz_GMN[N_int];

	for (int k = 0; k < N_int; k++)
	{
		GMN_zero[k] = 0.0;
	}


	///////////////////////////////////////////////
	// 5.2.2.2. Apply ageing

	for (int n = 0; n < N_pop; n++)
	{
		people[n].ager(theta);
	}


	///////////////////////////////////////////////
	// 5.2.2.3. Deaths
	//
	// Look again at how things are erased from vectors.

	int N_dead = 0;

	for (size_t n = 0; n < people.size(); n++)
	{
		/////////////////////////////////////////////
		// Everyone has an equal probability of dying

		if (theta.P_dead > genunf(0, 1))
		{
			people.erase(people.begin() + n);

			pi_n.erase(pi_n.begin() + n);
			lam_n.erase(lam_n.begin() + n);

			N_dead = N_dead + 1;
			n = n - 1;      // If we erase something, the next one moves into it's place so we don't want to step forward.
		}
		else {

			///////////////////////////////////////////
			// People die once they reach the maximum age

			if (people[n].age > theta.age_max)
			{
				people.erase(people.begin() + n);

				pi_n.erase(pi_n.begin() + n);
				lam_n.erase(lam_n.begin() + n);

				N_dead = N_dead + 1;
				n = n - 1;       // If we erase something, the next one moves into it's place so we don't want to step forward.
			}
		}
	}

	/////////////////////////////////////////////////////////////
	// 5.2.2.4. Births - set up to ensure balanced population.
	//          Can be adjusted to account for changing demography.

	double zeta_start, het_dif_track, q_rand;

	vector<double> zero_push(N_mosq);
	for (int v = 0; v < N_mosq; v++)
	{
		zero_push[v] = 0.0;
	}


	for (int n = 0; n < N_dead; n++)
	{
		zeta_start = exp(gennor(-0.5*theta.sig_het*theta.sig_het, theta.sig_het));

		while (zeta_start > theta.het_max)
		{
			zeta_start = exp(gennor(-0.5*theta.sig_het*theta.sig_het, theta.sig_het));
		}

		Individual HH(0.0, zeta_start);

		HH.S     = 1;
		HH.I_PCR = 0;
		HH.I_LM  = 0;
		HH.I_D   = 0;
		HH.T     = 0;
		HH.P     = 0;

		HH.A_par  = 0.0;
		HH.A_clin = 0.0;

		HH.A_par_boost  = 0;
		HH.A_clin_boost = 0;

		HH.A_par_timer  = -1.0;
		HH.A_clin_timer = -1.0;

		HH.AQ8_proph = 0;
		HH.AQ8_proph_timer = -1.0;

		HH.Hyp = 0;

		if (genunf(0.0, 1.0) < 0.5)
		{
			if (genunf(0.0, 1.0) > theta.P_occup)
			{
				// Domestic males

				HH.gender = 0;
			}
			else {
				// Occupational males

				HH.gender = 1;
			}
		}
		else {

			// females

			HH.gender = 2;
		}

		if ((HH.gender == 0) || (HH.gender == 1))
		{

			if (genunf(0.0, 1.0) < theta.G6PD_prev)
			{
				// males, deficient (hemizygous)

				HH.G6PD_def = 1;
				HH.G6PD_activity = gennor(theta.mu_G6PD_def, theta.sig_G6PD_def);
			}
			else {
				// males, normal (hemizygous)

				HH.G6PD_def = 0;
				HH.G6PD_activity = gennor(theta.mu_G6PD_nor, theta.sig_G6PD_nor);
			}
		}
		else {

			q_rand = genunf(0.0, 1.0);

			if (q_rand <= theta.G6PD_prev*theta.G6PD_prev)
			{
				// females, deficient (homozygous)

				HH.G6PD_def = 1;
				HH.G6PD_activity = gennor(theta.mu_G6PD_def, theta.sig_G6PD_def);
			}

			if ((q_rand > theta.G6PD_prev*theta.G6PD_prev) && (q_rand <= theta.G6PD_prev*theta.G6PD_prev + 2 * theta.G6PD_prev*(1.0 - theta.G6PD_prev)))
			{
				// females, deficient (heterozygous)

				HH.G6PD_def = 2;
				HH.G6PD_activity = gennor(theta.mu_G6PD_het, theta.sig_G6PD_het);
			}

			if (q_rand > (theta.G6PD_prev*theta.G6PD_prev + 2 * theta.G6PD_prev*(1.0 - theta.G6PD_prev)))
			{
				// females, normal (homozygous)

				HH.G6PD_def = 0;
				HH.G6PD_activity = gennor(theta.mu_G6PD_nor, theta.sig_G6PD_nor);
			}
		}

		if (genunf(0.0, 1.0) < theta.CYP2D6_prev)
		{
			HH.CYP2D6 = 1;
		}
		else {
			HH.CYP2D6 = 0;
		}

		HH.preg_age = 0;
		HH.pregnant = 0;
		HH.preg_timer = 0.0;


		/////////////////////////////////
		// Assign levels of maternally-acquired immunity
		// by finding women of child-bearing age with the
		// closest level of heterogeneity

		HH.A_par_mat = 0.0;
		HH.A_clin_mat = 0.0;

		het_dif_track = 1e10;

		for (size_t j = 0; j < people.size(); j++)
		{
			if (people[j].preg_age == 1)
			{
				if (abs(HH.zeta_het - people[j].zeta_het) < het_dif_track)
				{
					HH.A_par_mat = theta.P_mat*people[j].A_par_mat;
					HH.A_clin_mat = theta.P_mat*people[j].A_clin_mat;

					het_dif_track = (HH.zeta_het - people[j].zeta_het)*(HH.zeta_het - people[j].zeta_het);
				}
			}
		}


		///////////////////////////////////////////////////
		// Lagged exposure equals zero - they're not born yet!

		for (int k = 0; k < theta.H_track; k++)
		{
			HH.lam_bite_track.push_back(0.0);
			HH.lam_rel_track.push_back(0.0);
		}


		///////////////////////////////////////////////////
		// Assign intervention access scores

		for (int p = 0; p < N_int; p++)
		{
			for (int q = 0; q < N_int; q++)
			{
				theta.V_int_dummy[p][q] = theta.V_int[p][q];
			}
		}

		setgmn(GMN_zero, *theta.V_int_dummy, N_int, GMN_parm);

		genmn(GMN_parm, zz_GMN, GMN_work);

		for (int k = 0; k < N_int; k++)
		{
			HH.zz_int[k] = zz_GMN[k];
		}


		///////////////////////////////////////////////////
		// Born with no interventions

		HH.LLIN = 0;
		HH.IRS = 0;
		HH.IVM = 0;

		for (int v = 0; v < N_mosq; v++)
		{
			HH.w_VC[v] = 1.0;
			HH.y_VC[v] = 1.0;
			HH.z_VC[v] = 0.0;
		}


		/////////////////////////////////////////////////////////////////
		// 2.4.5. Push the created individual onto the vector of people

		people.push_back(move(HH));

		pi_n.push_back(zero_push);
		lam_n.push_back(zero_push);
	}



	///////////////////////////////////////////////////
	// Update individual-level vector control

	for (int n = 0; n < N_pop; n++)
	{
		people[n].vec_con_updater(theta);
	}


	///////////////////////////////////////////////////
	// 5.2.2.5. Update proportion of bites
	//
	//          Note the ordering of n and g loops. Need to
	//          check if this makes a difference for speed.
	//
	//          Should be able to make this quicker


	for (int n = 0; n < N_pop; n++)
	{
		pi_n[n][0] = people[n].zeta_het*(1.0 - theta.rho_age*exp(-people[n].age*theta.age_0_inv));

		pi_n[n][1] = 0.0;

		if (people[n].gender == 1)
		{
			if( (people[n].age > theta.occup_age_low) && ((people[n].age <= theta.occup_age_high)) )
			{
				pi_n[n][1] = people[n].zeta_het*(1.0 - theta.rho_age*exp(-people[n].age*theta.age_0_inv));
			}
		}
	}

	double SIGMA_PI[N_mosq];
	for (int v = 0; v < N_mosq; v++)
	{
		SIGMA_PI[v] = 0.0;
	}

	for (int n = 0; n < N_pop; n++)
	{
		for (int v = 0; v < N_mosq; v++)
		{
			SIGMA_PI[v] = SIGMA_PI[v] + pi_n[n][v];
		}
	}

	for (int v = 0; v < N_mosq; v++)
	{
		SIGMA_PI[v] = 1.0 / SIGMA_PI[v];
	}

	for (int n = 0; n < N_pop; n++)
	{
		for (int v = 0; v < N_mosq; v++)
		{
			pi_n[n][v] = pi_n[n][v] * SIGMA_PI[v];
		}
	}


	///////////////////////////////////////////////////
	// 5.2.2.6. Update population-level vector control quantities

	for (int v = 0; v < N_mosq; v++)
	{
		SUM_pi_w[v] = 0;
	}

	for (int n = 0; n < N_pop; n++)
	{
		for (int v = 0; v < N_mosq; v++)
		{
			SUM_pi_w[v] = SUM_pi_w[v] + pi_n[n][v] * people[n].w_VC[v];
		}
	}


	for (int v = 0; v < N_mosq; v++)
	{
		W_VC[v] = 1.0 - theta.Q_0[v] + theta.Q_0[v] * SUM_pi_w[v];
		Z_VC[v] = theta.Q_0[v] * SUM_pi_z[v];

		delta_1_VC[v] = theta.delta_1 / (1.0 - Z_VC[v]);
		delta_VC[v] = delta_1_VC[v] + theta.delta_2;

		p_1_VC[v] = theta.p_1[v] * W_VC[v] / (1.0 - Z_VC[v] * theta.p_1[v]);

		mu_M_VC[v] = -log(p_1_VC[v] * theta.p_2[v]) / delta_VC[v];

		Q_VC[v] = 1.0 - (1.0 - theta.Q_0[v]) / W_VC[v];

		aa_VC[v] = Q_VC[v] / delta_VC[v];

		exp_muM_tauM_VC[v] = exp(-mu_M_VC[v] * theta.tau_M[v]);
		beta_VC[v] = theta.eps_max[v] * mu_M_VC[v] / (exp(delta_VC[v] * mu_M_VC[v]) - 1.0);
	}


	///////////////////////////////////////////////////
	// 5.2.2.7. Update individual-level force of infection on humans

	for (int n = 0; n < N_pop; n++)
	{
		for (int v = 0; v < N_mosq; v++)
		{
			lam_n[n][v] = aa_VC[v] * pi_n[n][v] * people[n].w_VC[v] / SUM_pi_w[v];
		}
	}


	///////////////////////////////////////////////////
	// 5.2.2.8. Implement moves between compartments
	//
	// TO DO: Can take some multiplications out of the loop.

	double lam_bite_base[N_mosq];
	double lam_bite_n;     // better notation (this is lam_bite)

	for (int v = 0; v < N_mosq; v++)
	{
		lam_bite_base[v] = (double(N_pop))*theta.bb*yM[v][5];
	}

	for (int n = 0; n < N_pop; n++)
	{
		lam_bite_n = 0.0;

		for (int v = 0; v < N_mosq; v++)
		{
			lam_bite_n = lam_bite_n + lam_n[n][v] * lam_bite_base[v];
		}

		people[n].state_mover(theta, lam_bite_n);
	}

}


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  5.2.3.  Summarise the output from the population                        //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

void Population::summary()
{
	for (int k = 0; k < N_H_comp; k++)
	{
		yH[k] = 0.0;
	}

	// for (int k = 0; k < 15; k++)
	// {
	// 	prev_all[k] = 0.0;
	// 	prev_U5[k] = 0.0;
	// 	prev_U10[k] = 0.0;
	// }


	for (int n = 0; n < N_pop; n++)
	{
		////////////////////////////////////////
		// Numbers in each compartment

		yH[0] = yH[0] + people[n].S;
		yH[1] = yH[1] + people[n].I_PCR;
		yH[2] = yH[2] + people[n].I_LM;
		yH[3] = yH[3] + people[n].I_D;
		yH[4] = yH[4] + people[n].T;
		yH[5] = yH[5] + people[n].P;


		//////////////////////////////////////////////
		//////////////////////////////////////////////
		// Summary - full population

		// ////////////////////////////////////////
		// // Prevalence
		//
		// prev_all[0] = prev_all[0] + 1;                                                     // Numbers - denominator
		// prev_all[1] = prev_all[1] + people[n].I_PCR + people[n].I_LM +
		// 	                      + people[n].I_D + people[n].T;                           // PCR detectable infections
		// prev_all[2] = prev_all[2] + people[n].I_LM + people[n].I_D + people[n].T;          // LM detectable infections
		// prev_all[3] = prev_all[3] + people[n].I_D + people[n].T;                           // Clinical episodes
		//
		// if (people[n].Hyp > 0)
		// {
		// 	prev_all[4] = prev_all[4] + 1;                     // Hypnozoite positive
		//
		// 	prev_all[5] = prev_all[5] + people[n].Hyp;    // Number of batches of hypnozoites
		// }
		//
		// ////////////////////////////////////////
		// // Incidence
		//
		// prev_all[6]  = prev_all[6]  + people[n].I_PCR_new;
		// prev_all[7]  = prev_all[7]  + people[n].I_LM_new;
		// prev_all[8]  = prev_all[8]  + people[n].I_D_new;
		// prev_all[9]  = prev_all[9]  + people[n].CQ_treat;
		// prev_all[10] = prev_all[10] + people[n].PQ_treat;
		// prev_all[11] = prev_all[11] + people[n].TQ_treat;
		//
		// prev_all[12] = prev_all[12] + people[n].G6PD_test;
		// prev_all[13] = prev_all[13] + people[n].PQ_effective;
		// prev_all[14] = prev_all[14] + people[n].TQ_effective;

		//////////////////////////////////////////////
		//////////////////////////////////////////////
		// Summary - under 5's

		// if (people[n].age < 1825.0)
		// {
		// 	////////////////////////////////////////
		// 	// Prevalence
		//
		// 	prev_U5[0] = prev_U5[0] + 1;                                                   // Numbers - denominator
		// 	prev_U5[1] = prev_U5[1] + people[n].I_PCR + people[n].I_LM
		// 		                    + people[n].I_D + people[n].T;                         // PCR detectable infections
		// 	prev_U5[2] = prev_U5[2] + people[n].I_LM + people[n].I_D + people[n].T;        // LM detectable infections
		// 	prev_U5[3] = prev_U5[3] + people[n].I_D + people[n].T;                         // Clinical episodes
		//
		// 	if (people[n].Hyp > 0)
		// 	{
		// 		prev_U5[4] = prev_U5[4] + 1;                     // Hypnozoite positive
		//
		// 		prev_U5[5] = prev_U5[5] + people[n].Hyp;    // Number of batches of hypnozoites
		// 	}
		//
		//
		// 	////////////////////////////////////////
		// 	// Incidence
		//
		// 	prev_U5[6]  = prev_U5[6]  + people[n].I_PCR_new;
		// 	prev_U5[7]  = prev_U5[7]  + people[n].I_LM_new;
		// 	prev_U5[8]  = prev_U5[8]  + people[n].I_D_new;
		// 	prev_U5[9]  = prev_U5[9]  + people[n].CQ_treat;
		// 	prev_U5[10] = prev_U5[10] + people[n].PQ_treat;
		// 	prev_U5[11] = prev_U5[11] + people[n].TQ_treat;
		//
		// 	prev_U5[12] = prev_U5[12] + people[n].G6PD_test;
		// 	prev_U5[13] = prev_U5[13] + people[n].PQ_effective;
		// 	prev_U5[14] = prev_U5[14] + people[n].TQ_effective;
		// }

		//////////////////////////////////////////////
		//////////////////////////////////////////////
		// Summary - under 10's

		// if (people[n].age < 3650.0)
		// {
		// 	////////////////////////////////////////
		// 	// Prevalence
		//
		// 	prev_U10[0] = prev_U10[0] + 1;                                               // Numbers - denominator
		// 	prev_U10[1] = prev_U10[1] + people[n].I_PCR + people[n].I_LM
		// 		                      + people[n].I_D + people[n].T;                     // PCR detectable infections
		// 	prev_U10[2] = prev_U10[2] + people[n].I_LM + people[n].I_D + people[n].T;    // LM detectable infections
		// 	prev_U10[3] = prev_U10[3] + people[n].I_D + people[n].T;                     // Clinical episodes
		//
		// 	if (people[n].Hyp > 0)
		// 	{
		// 		prev_U10[4] = prev_U10[4] + 1;                     // Hypnozoite positive
		//
		// 		prev_U10[5] = prev_U10[5] + people[n].Hyp;    // Number of batches of hypnozoites
		// 	}
		//
		//
		// 	////////////////////////////////////////
		// 	// Incidence
		//
		// 	prev_U10[6]  = prev_U10[6]  + people[n].I_PCR_new;
		// 	prev_U10[7]  = prev_U10[7]  + people[n].I_LM_new;
		// 	prev_U10[8]  = prev_U10[8]  + people[n].I_D_new;
		// 	prev_U10[9]  = prev_U10[9]  + people[n].CQ_treat;
		// 	prev_U10[10] = prev_U10[10] + people[n].PQ_treat;
		// 	prev_U10[11] = prev_U10[11] + people[n].TQ_treat;
		//
		// 	prev_U10[12] = prev_U10[12] + people[n].G6PD_test;
		// 	prev_U10[13] = prev_U10[13] + people[n].PQ_effective;
		// 	prev_U10[14] = prev_U10[14] + people[n].TQ_effective;
		// }
	}


	//////////////////////////////
	// Intervention coverage

	// LLIN_cov_t = 0;
	// IRS_cov_t = 0;
	// IVM_cov_t  = 0;
	PvHR_t         = 0;
  PvHR_batch_t   = 0;
  I_PCR_new_t    = 0;
  I_LM_new_t     = 0;
  I_D_new_t      = 0;
  T_new_t        = 0;
	CQ_treat_t     = 0;
	CQ_relapse_t   = 0;
	CQ_bite_t      = 0;
	G6PD_tests_t   = 0;
	PQ_treat_t     = 0;
	TQ_treat_t     = 0;
	PQ_effective_t = 0;
	TQ_effective_t = 0;
	PQ_G6PD_def_t  = 0;
	TQ_G6PD_def_t  = 0;

	// pregnant_t = 0;
	//
	// PQ_overtreat_t    = 0;
	// PQ_overtreat_9m_t = 0;
	//
	// TQ_overtreat_t = 0;
	// TQ_overtreat_9m_t = 0;
	//
	for (int n = 0; n < N_pop; n++)
	{
	// 	LLIN_cov_t = LLIN_cov_t + people[n].LLIN;
	// 	IRS_cov_t  = IRS_cov_t + people[n].IRS;
	// 	IVM_cov_t  = IVM_cov_t  + people[n].IVM;
	if (people[n].Hyp > 0)
	{
	  PvHR_t       = PvHR_t       + 1;              // Hypnozoite positive
	  PvHR_batch_t = PvHR_batch_t + people[n].Hyp;  // Number of batches of hypnozoites
	}
	I_PCR_new_t    = I_PCR_new_t    + people[n].I_PCR_new;
	I_LM_new_t     = I_LM_new_t     + people[n].I_LM_new;
	I_D_new_t      = I_D_new_t      + people[n].I_D_new;
	T_new_t        = T_new_t        + people[n].T_new;
	if (people[n].CQ_treat > 0)
	{
	  CQ_treat_t     = CQ_treat_t     + people[n].CQ_treat;
	  CQ_relapse_t   = CQ_relapse_t   + people[n].CQ_relapse;
	  CQ_bite_t      = CQ_bite_t      + people[n].CQ_bite;
	}
	G6PD_tests_t   = G6PD_tests_t   + people[n].G6PD_test;
	if (people[n].PQ_treat > 0)
	{
	  PQ_treat_t     = PQ_treat_t     + people[n].PQ_treat;
	  PQ_effective_t = PQ_effective_t + people[n].PQ_effective;
	  PQ_G6PD_def_t  = PQ_G6PD_def_t  + people[n].PQ_G6PD_def;
	}
	if (people[n].TQ_treat > 0)
	{
	  TQ_treat_t     = TQ_treat_t     + people[n].TQ_treat;
	  TQ_effective_t = TQ_effective_t + people[n].TQ_effective;
	  TQ_G6PD_def_t  = TQ_G6PD_def_t  + people[n].TQ_G6PD_def;
	}

	// 	pregnant_t = pregnant_t + people[n].pregnant;
	//
	// 	PQ_overtreat_t    = PQ_overtreat_t    + people[n].PQ_overtreat;
	// 	PQ_overtreat_9m_t = PQ_overtreat_9m_t + people[n].PQ_overtreat_9m;
	//
	// 	TQ_overtreat_t    = TQ_overtreat_t    + people[n].TQ_overtreat;
	// 	TQ_overtreat_9m_t = TQ_overtreat_9m_t + people[n].TQ_overtreat_9m;
	}


	//////////////////////////////
	// Age and gender stratified cases

	// cases_M_O16_t = 0;       // Detected cases in males over 16
	// cases_M_U16_t = 0;       // Detected cases in males under 16
	// cases_F_O16_t = 0;       // Detected cases in females over 16
	// cases_F_U16_t = 0;       // Detected cases in females under 16
	// cases_preg_t  = 0;       // Detected cases in pregnant women
	//
	// for (int n = 0; n < N_pop; n++)
	// {
	// 	if( (people[n].gender == 0) || (people[n].gender == 1) )
	// 	{
	// 		if (people[n].age > 5840.0)
	// 		{
	// 			cases_M_O16_t = cases_M_O16_t + people[n].CQ_treat;
	// 		} else {
	// 			cases_M_U16_t = cases_M_U16_t + people[n].CQ_treat;
	// 		}
	// 	}
	//
	// 	if( people[n].gender == 2 )
	// 	{
	// 		if (people[n].age > 5840.0)
	// 		{
	// 			cases_F_O16_t = cases_F_O16_t + people[n].CQ_treat;
	// 		} else {
	// 			cases_F_U16_t = cases_F_U16_t + people[n].CQ_treat;
	// 		}
	//
	// 		if (people[n].pregnant == 1)
	// 		{
	// 			cases_preg_t = cases_preg_t + people[n].CQ_treat;
	// 		}
	// 	}
	// }


	//////////////////////////////
	// Immunity

	double A_par_mean = 0.0, A_clin_mean = 0.0;

	for (int n = 0; n < N_pop; n++)
	{
		A_par_mean  = A_par_mean  + people[n].A_par;
		A_clin_mean = A_clin_mean + people[n].A_clin;
	}

	A_par_mean_t  = A_par_mean / ((double)N_pop);
	A_clin_mean_t = A_clin_mean / ((double)N_pop);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
//          //                                                                                        //
//  5.2.4.  //  LU decomposition of a matrix                                                          //
//          //  Based on ludcmp.cpp from Numerical Recipes in C++                                     //
//          //                                                                                        //
////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                    //
//  Given a matrix a[1..n][1..n], this routine replaces it by the LU decomposition of a rowwise       //
//  permutation of itself. a and n are input. a is output, arranged as in equation (2.3.14) above;    //
//  indx[1..n] is an output vector that records the row permutation effected by the partial           //
//  pivoting; d is output as ±1 depending on whether the number of row interchanges was even          //
//  or odd, respectively. This routine is used in combination with lubksb to solve linear equations   //
//  or invert a matrix                                                                                //
//                                                                                                    //
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

void ludcmp(vector<vector<double>> &a, int n_dim, vector<int> &indx, double &d)
{
	const double TINY = 1.0e-20;
	int i, imax, j, k;
	double big, dum, sum, temp;


	vector<double> vv(n_dim);
	d = 1.0;
	for (i = 0; i < n_dim; i++) {
		big = 0.0;
		for (j = 0; j < n_dim; j++)
			if ((temp = fabs(a[i][j])) > big) big = temp;
		if (big == 0.0) throw("Singular matrix in routine ludcmp");
		vv[i] = 1.0 / big;
	}
	for (j = 0; j < n_dim; j++) {
		for (i = 0; i < j; i++) {
			sum = a[i][j];
			for (k = 0; k < i; k++) sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
		}
		big = 0.0;
		for (i = j; i < n_dim; i++) {
			sum = a[i][j];
			for (k = 0; k < j; k++) sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
			if ((dum = vv[i] * fabs(sum)) >= big) {
				big = dum;
				imax = i;
			}
		}
		if (j != imax) {
			for (k = 0; k < n_dim; k++) {
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
			}
			d = -d;
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (a[j][j] == 0.0) a[j][j] = TINY;
		if (j != n_dim - 1) {
			dum = 1.0 / (a[j][j]);
			for (i = j + 1; i < n_dim; i++) a[i][j] *= dum;
		}
	}
}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//          //                                                      //
//  5.2.5.  //  Matrix back substitution                            //
//          //  Based on lubksb.cpp from Numerical Recipes in C++   //
//          //                                                      //
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

void lubksb(vector<vector<double>> &a, int n_dim, vector<int> &indx, vector<double> &b)
{
	int i, ii = 0, ip, j;
	double sum;


	for (i = 0; i < n_dim; i++) {
		ip = indx[i];
		sum = b[ip];
		b[ip] = b[i];
		if (ii != 0)
			for (j = ii - 1; j < i; j++) sum -= a[i][j] * b[j];
		else if (sum != 0.0)
			ii = i + 1;
		b[i] = sum;
	}
	for (i = n_dim - 1; i >= 0; i--) {
		sum = b[i];
		for (j = i + 1; j < n_dim; j++) sum -= a[i][j] * b[j];
		b[i] = sum / a[i][i];
	}
}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//          //                                                      //
//  5.2.6.  //  Matrix inversion and calculation of determinant.    //
//          //  Based on ludcmp.cpp from Numerical Recipes in C++   //
//          //                                                      //
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

void matrix_inv(vector<vector<double>> &a, int n, vector<vector<double>> &a_inv)
{
	vector<int> a_index(n);
	vector<double> col(n);
	double d;

	ludcmp(a, n, a_index, d);

	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < n; i++)
		{
			col[i] = 0.0;
		}
		col[j] = 1.0;

		lubksb(a, n, a_index, col);

		for (int i = 0; i < n; i++)
		{
			a_inv[i][j] = col[i];
		}
	}

}


///////////////////////////////////////////////////
///////////////////////////////////////////////////
//          //                                   //
//  5.2.7.  //  Matrix multiplication            //
//          //  Calculates inv(MM*)xx            //
//          //                                   //
///////////////////////////////////////////////////
///////////////////////////////////////////////////

void inv_MM_bb(vector<vector<double>> &MM, vector<double> &bb, vector<double> &xx, int n_dim)
{
	///////////////////////////////////////////////
	// 5.2.7.1. calculate inv(MM)

	vector<vector<double>> MM_inv;
	MM_inv.resize(n_dim);
	for (int k = 0; k < n_dim; k++)
	{
		MM_inv[k].resize(n_dim);
	}

	matrix_inv(MM, n_dim, MM_inv);


	///////////////////////////////////////////////
	// 5.2.7.2. calculate xx = MM_inv*bb

	for (int i = 0; i < n_dim; i++)
	{
		xx[i] = 0.0;

		for (int j = 0; j < n_dim; j++)
		{
			xx[i] = xx[i] + MM_inv[i][j] * bb[j];
		}
	}

}


///////////////////////////////////////////////////
///////////////////////////////////////////////////
//          //                                   //
//  5.2.8.  //  Equilibrium matrix               //
//          //                                   //
///////////////////////////////////////////////////
///////////////////////////////////////////////////

void MM_ij( int g, int i, int j, Params& theta, double r_age[], vector<vector<double>> &MM,
	        vector<vector<vector<double>>> lam_eq, vector<vector<vector<vector<double>>>> phi_LM_eq,
	        vector<vector<vector<vector<double>>>> phi_D_eq, vector<vector<vector<vector<double>>>> r_PCR_eq,
	        vector<vector<double>> treat_PQeff )
{
	//////////////////////////////////////////////
	// 5.2.8.1. Initialise matrix with zeros

	for (int i1 = 0; i1 < (N_H_comp * (K_max + 1)); i1++)
	{
		for (int j1 = 0; j1 < (N_H_comp * (K_max + 1)); j1++)
		{
			MM[i1][j1] = 0.0;
		}
	}


	//////////////////////////////////////////////
	// 5.2.8.2. Fill out non-zero elements

	for (int k1 = 0; k1 < (K_max + 1); k1++)
	{
		for (int k2 = 0; k2 < (K_max + 1); k2++)
		{
			MM[0*(K_max + 1) + k1][0*(K_max + 1) + k2] = - lam_eq[g][i][j]*theta.D_MAT[k1][k2] - theta.ff*theta.K_MAT[k1][k2]
				                                         + theta.gamma_L*theta.L_MAT[k1][k2] - theta.mu_H*theta.D_MAT[k1][k2] - r_age[i] * theta.D_MAT[k1][k2];
			MM[0*(K_max + 1) + k1][1*(K_max + 1) + k2] = + r_PCR_eq[g][i][j][k2] * theta.D_MAT[k1][k2];
			MM[0*(K_max + 1) + k1][5*(K_max + 1) + k2] = + theta.r_P*theta.D_MAT[k1][k2];

			MM[1*(K_max + 1) + k1][0*(K_max + 1) + k2] = + lam_eq[g][i][j]*(1.0 - phi_LM_eq[g][i][j][k2])*theta.OD_MAT[k1][k2] + theta.ff*(1.0 - phi_LM_eq[g][i][j][k2])*theta.K_MAT[k1][k2];
			MM[1*(K_max + 1) + k1][1*(K_max + 1) + k2] = - lam_eq[g][i][j]*theta.D_MAT[k1][k2] - theta.ff*theta.K_MAT[k1][k2] - r_PCR_eq[g][i][j][k2]*theta.D_MAT[k1][k2]
				                                         + lam_eq[g][i][j]*(1.0 - phi_LM_eq[g][i][j][k2])*theta.OD_MAT[k1][k2] + theta.ff*(1.0 - phi_LM_eq[g][i][j][k2])*theta.K_MAT[k1][k2]
				                                         + theta.gamma_L*theta.L_MAT[k1][k2] - theta.mu_H*theta.D_MAT[k1][k2] - r_age[i] * theta.D_MAT[k1][k2];
			MM[1*(K_max + 1) + k1][2*(K_max + 1) + k2] = + theta.r_LM*theta.D_MAT[k1][k2];

			MM[2*(K_max + 1) + k1][0*(K_max + 1) + k2] = + lam_eq[g][i][j]*phi_LM_eq[g][i][j][k2]*(1.0 - phi_D_eq[g][i][j][k2])*theta.OD_MAT[k1][k2] + theta.ff*phi_LM_eq[g][i][j][k2]*(1.0 - phi_D_eq[g][i][j][k2])*theta.K_MAT[k1][k2];
			MM[2*(K_max + 1) + k1][1*(K_max + 1) + k2] = + lam_eq[g][i][j]*phi_LM_eq[g][i][j][k2]*(1.0 - phi_D_eq[g][i][j][k2])*theta.OD_MAT[k1][k2] + theta.ff*phi_LM_eq[g][i][j][k2]*(1.0 - phi_D_eq[g][i][j][k2])*theta.K_MAT[k1][k2];
			MM[2*(K_max + 1) + k1][2*(K_max + 1) + k2] = - lam_eq[g][i][j]*theta.D_MAT[k1][k2] - theta.ff*theta.K_MAT[k1][k2] - theta.r_LM*theta.D_MAT[k1][k2]
				                                         + lam_eq[g][i][j]*(1.0 - phi_D_eq[g][i][j][k2])*theta.OD_MAT[k1][k2] + theta.ff*(1.0 - phi_D_eq[g][i][j][k2])*theta.K_MAT[k1][k2]
				                                         + theta.gamma_L*theta.L_MAT[k1][k2] - theta.mu_H*theta.D_MAT[k1][k2] - r_age[i]*theta.D_MAT[k1][k2];
			MM[2*(K_max + 1) + k1][3*(K_max + 1) + k2] = + theta.r_D*theta.D_MAT[k1][k2];

			MM[3*(K_max + 1) + k1][0*(K_max + 1) + k2] = + lam_eq[g][i][j]*phi_LM_eq[g][i][j][k2]*phi_D_eq[g][i][j][k2]*(1.0 - theta.CM_CQ_coveff)*theta.OD_MAT[k1][k2]
				                                         + theta.ff*phi_LM_eq[g][i][j][k2]*phi_D_eq[g][i][j][k2]*(1.0 - theta.CM_CQ_coveff)*theta.K_MAT[k1][k2];
			MM[3*(K_max + 1) + k1][1*(K_max + 1) + k2] = + lam_eq[g][i][j]*phi_LM_eq[g][i][j][k2]*phi_D_eq[g][i][j][k2]*(1.0 - theta.CM_CQ_coveff)*theta.OD_MAT[k1][k2]
				                                         + theta.ff*phi_LM_eq[g][i][j][k2]*phi_D_eq[g][i][j][k2]*(1.0 - theta.CM_CQ_coveff)*theta.K_MAT[k1][k2];
			MM[3*(K_max + 1) + k1][2*(K_max + 1) + k2] = + lam_eq[g][i][j]*phi_D_eq[g][i][j][k2]*(1.0 - theta.CM_CQ_coveff)*theta.OD_MAT[k1][k2] + theta.ff*phi_D_eq[g][i][j][k2]*(1.0 - theta.CM_CQ_coveff)*theta.K_MAT[k1][k2];
			MM[3*(K_max + 1) + k1][3*(K_max + 1) + k2] = - lam_eq[g][i][j]*theta.D_MAT[k1][k2] - theta.r_D*theta.D_MAT[k1][k2] + lam_eq[g][i][j]*theta.OD_MAT[k1][k2]
				                                         + theta.gamma_L*theta.L_MAT[k1][k2] - theta.mu_H*theta.D_MAT[k1][k2] - r_age[i]*theta.D_MAT[k1][k2];

			MM[4*(K_max + 1) + k1][0*(K_max + 1) + k2] = + lam_eq[g][i][j]*phi_LM_eq[g][i][j][k2]*phi_D_eq[g][i][j][k2]*theta.CM_CQ_coveff*(1.0 - treat_PQeff[g][i])*theta.OD_MAT[k1][k2]
				                                         + theta.ff*phi_LM_eq[g][i][j][k2]*phi_D_eq[g][i][j][k2]*theta.CM_CQ_coveff*(1.0 - treat_PQeff[g][i])*theta.K_MAT[k1][k2];
			MM[4*(K_max + 1) + k1][1*(K_max + 1) + k2] = + lam_eq[g][i][j]*phi_LM_eq[g][i][j][k2]*phi_D_eq[g][i][j][k2]*theta.CM_CQ_coveff*(1.0 - treat_PQeff[g][i])*theta.OD_MAT[k1][k2]
				                                         + theta.ff*phi_LM_eq[g][i][j][k2]*phi_D_eq[g][i][j][k2]*theta.CM_CQ_coveff*(1.0 - treat_PQeff[g][i])*theta.K_MAT[k1][k2];
			MM[4*(K_max + 1) + k1][2*(K_max + 1) + k2] = + lam_eq[g][i][j]*phi_D_eq[g][i][j][k2]*theta.CM_CQ_coveff*(1.0 - treat_PQeff[g][i])*theta.OD_MAT[k1][k2]
				                                         + theta.ff*phi_D_eq[g][i][j][k2]*theta.CM_CQ_coveff*(1.0 - treat_PQeff[g][i])*theta.K_MAT[k1][k2];
			MM[4*(K_max + 1) + k1][4*(K_max + 1) + k2] = - lam_eq[g][i][j]*theta.D_MAT[k1][k2] - theta.r_T*theta.D_MAT[k1][k2] + lam_eq[g][i][j]*theta.OD_MAT[k1][k2]
				                                         + theta.gamma_L*theta.L_MAT[k1][k2] - theta.mu_H*theta.D_MAT[k1][k2] - r_age[i]*theta.D_MAT[k1][k2];

			MM[5*(K_max + 1) + k1][4*(K_max + 1) + k2] = + theta.r_T*theta.D_MAT[k1][k2];
			MM[5*(K_max + 1) + k1][5*(K_max + 1) + k2] = - lam_eq[g][i][j]*theta.D_MAT[k1][k2] - theta.r_P*theta.D_MAT[k1][k2] + lam_eq[g][i][j]*theta.OD_MAT[k1][k2]
				                                         + theta.gamma_L*theta.L_MAT[k1][k2] - theta.mu_H*theta.D_MAT[k1][k2] - r_age[i]*theta.D_MAT[k1][k2];
		}
	}


	/////////////////////////////////////////////
	// 5.2.8.3. Account for transitions into
	// the treatment state with zero hypnozoites

	for (int k1 = 0; k1 < (K_max + 1); k1++)
	{
		for (int k2 = 0; k2 < (K_max + 1); k2++)
		{
			MM[4*(K_max + 1) + 0][0*(K_max + 1) + k2] = MM[4*(K_max + 1) + 0][0*(K_max + 1) + k2]
				                                        + lam_eq[g][i][j]*phi_LM_eq[g][i][j][k2]*phi_D_eq[g][i][j][k2]*theta.CM_CQ_coveff*treat_PQeff[g][i]*theta.OD_MAT[k1][k2]
				                                        + theta.ff*phi_LM_eq[g][i][j][k2]*phi_D_eq[g][i][j][k2]*theta.CM_CQ_coveff*treat_PQeff[g][i]*theta.K_MAT[k1][k2];
			MM[4*(K_max + 1) + 0][1*(K_max + 1) + k2] = MM[4*(K_max + 1) + 0][1*(K_max + 1) + k2]
				                                        + lam_eq[g][i][j]*phi_LM_eq[g][i][j][k2]*phi_D_eq[g][i][j][k2]*theta.CM_CQ_coveff*treat_PQeff[g][i]*theta.OD_MAT[k1][k2]
				                                        + theta.ff*phi_LM_eq[g][i][j][k2]*phi_D_eq[g][i][j][k2]*theta.CM_CQ_coveff*treat_PQeff[g][i]*theta.K_MAT[k1][k2];
			MM[4*(K_max + 1) + 0][2*(K_max + 1) + k2] = MM[4*(K_max + 1) + 0][2*(K_max + 1) + k2]
				                                        + lam_eq[g][i][j]*phi_D_eq[g][i][j][k2]*theta.CM_CQ_coveff*treat_PQeff[g][i]*theta.OD_MAT[k1][k2]
				                                        + theta.ff*phi_D_eq[g][i][j][k2]*theta.CM_CQ_coveff*treat_PQeff[g][i]*theta.K_MAT[k1][k2];
		}
	}
}


/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//          //                                                 //
//  5.2.9.  //  Guass-Hermite weights for Gaussian quadrature  //
//          //  integration with Normal distribution.          //
//          //  Code is adapted from gauher.cpp from NR3       //
//          //                                                 //
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

void Population::gauher(Params& theta)
{
	double x[N_het];
	double w[N_het];

	/////////////////////////////
	// PART 1

	const double EPS = 1.0e-14, PIM4 = 0.7511255444649425;
	const int MAXIT = 10;
	int i, its, j, m;
	double p1, p2, p3, pp, z, z1;

	m = (N_het + 1) / 2;
	for (i = 0; i < m; i++) {
		if (i == 0) {
			z = sqrt(double(2 * N_het + 1)) - 1.85575*pow(double(2 * N_het + 1), -0.16667);
		}
		else if (i == 1) {
			z -= 1.14*pow(double(N_het), 0.426) / z;
		}
		else if (i == 2) {
			z = 1.86*z - 0.86*x[0];
		}
		else if (i == 3) {
			z = 1.91*z - 0.91*x[1];
		}
		else {
			z = 2.0*z - x[i - 2];
		}
		for (its = 0; its < MAXIT; its++) {
			p1 = PIM4;
			p2 = 0.0;
			for (j = 0; j < N_het; j++) {
				p3 = p2;
				p2 = p1;
				p1 = z * sqrt(2.0 / (j + 1))*p2 - sqrt(double(j) / (j + 1))*p3;
			}
			pp = sqrt(double(2 * N_het))*p2;
			z1 = z;
			z = z1 - p1 / pp;
			if (fabs(z - z1) <= EPS) break;
		}
		if (its >= MAXIT) throw("too many iterations in gauher");
		x[i] = z;
		x[N_het - 1 - i] = -z;
		w[i] = 2.0 / (pp*pp);
		w[N_het - 1 - i] = w[i];
	}


	/////////////////////////////
	// PART 2

	double w_sum = 0.0;

	for (int j = 0; j < N_het; j++)
	{
		w_sum = w_sum + w[j];
	}

	////////////////////////////////
	// Note that N_het-1-j here instead of j to ensure x_het increasing

	for (int j = 0; j < N_het; j++)
	{
		x_het[j] = exp(theta.sig_het*x[N_het - 1 - j] * sqrt(2.0) - 0.5*theta.sig_het*theta.sig_het);
		w_het[j] = w[j] / w_sum;
	}


	////////////////////////////
	// temporary for N_het = 1

	if (N_het == 1)
	{
		x_het[0] = 1.0;
		w_het[0] = 1.0;
	}

	for (int i = 0; i < N_age; i++)
	{
		for (int j = 0; j < N_het; j++)
		{
			x_age_het[i][j] = age_bite[i] * x_het[j];
			w_age_het[i][j] = age_demog[i] * w_het[j];
		}
	}

	////////////////////////////////
	// Boundaries of heterogeneity compartments

	x_het_bounds[0] = 0.0;

	for (int j = 1; j < N_het; j++)
	{
		x_het_bounds[j] = exp(0.5*(log(x_het[j - 1]) + log(x_het[j])));
	}

	x_het_bounds[N_het] = theta.het_max;
}


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//           //                                                                   //
//  5.2.10.  //  Setup objects for storing population information.                //
//           //                                                                   //
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

void Population::pop_setup(Params& theta)
{
	//////////////////////////////////////////////////////
	// 5.2.10.1. Set up age and heterogeneity compartments

	////////////////////////////////////////
	// 5.2.10.1.1. Bounds of age bins

	// double age_bounds_set[N_age+1] = {0.0*365.0, 20.0*365.0, 40.0*365.0, 60.0*365.0, 80.0*365.0};

	double age_bounds_set[N_age + 1] = { 0.0*365.0, 0.2*365.0, 0.4*365.0, 0.6*365.0, 0.8*365.0, 1.0*365.0,
                                         1.2*365.0, 1.4*365.0, 1.6*365.0, 1.8*365.0, 2.0*365.0,
                                         2.2*365.0, 2.4*365.0, 2.6*365.0, 2.8*365.0, 3.0*365.0,
                                         3.4*365.0, 3.8*365.0, 4.2*365.0, 4.6*365.0, 5.0*365.0,
                                         5.5*365.0, 6.0*365.0, 6.5*365.0, 7.0*365.0, 7.5*365.0, 8.0*365.0, 8.5*365.0, 9.0*365.0, 9.5*365.0, 10.0*365.0,
                                         11.0*365.0, 12.0*365.0, 13.0*365.0, 14.0*365.0, 15.0*365.0, 16.0*365.0, 17.0*365.0, 18.0*365.0, 19.0*365.0, 20.0*365.0,
                                         22.0*365.0, 24.0*365.0, 26.0*365.0, 28.0*365.0, 30.0*365.0, 32.0*365.0, 34.0*365.0, 36.0*365.0, 38.0*365.0, 40.0*365.0,
                                         45.0*365.0, 50.0*365.0, 55.0*365.0, 60.0*365.0, 65.0*365.0, 70.0*365.0, 75.0*365.0, 80.0*365.0 };

    for (int i = 0; i < (N_age + 1); i++)
	{
		age_bounds[i] = age_bounds_set[i];
	}


	////////////////////////////////////////////
	// 5.2.10.1.2. Proportion in each age bin

	for (int i = 0; i < (N_age - 1); i++)
	{
		age_demog[i] = exp(-theta.mu_H*age_bounds[i]) - exp(-theta.mu_H*age_bounds[i + 1]);
	}

	age_demog[N_age - 1] = 1.0;

	for (int i = 0; i < (N_age - 1); i++)
	{
		age_demog[N_age - 1] = age_demog[N_age - 1] - age_demog[i];
	}


	////////////////////////////////////////
	// 5.2.10.1.3. Ageing rates - formula below ensures
	//             balanced demography

	r_age[0] = theta.mu_H*(1.0 - age_demog[0]) / age_demog[0];

	for (int i = 1; i < (N_age - 1); i++)
	{
		r_age[i] = (r_age[i - 1] * age_demog[i - 1] - theta.mu_H*age_demog[i]) / age_demog[i];
	}

	r_age[N_age - 1] = 0.0;


	////////////////////////////////////////
	// 5.2.10.1.4. Age-dependent mosquito biting rates

	for (int i = 0; i < N_age; i++)
	{
		age_mids[i] = 0.5*(age_bounds[i] + age_bounds[i + 1]);
	}

	for (int i = 0; i < N_age; i++)
	{
		age_bite[i] = 1.0 - theta.rho_age*exp(-age_mids[i] / theta.age_0);
	}

	P_age_bite = exp(-t_step / theta.age_0);


	///////////////////////////////////////
	// Ensure total bites are normalised

	double omega_age = 0.0;

	for (int i = 0; i < N_age; i++)
	{
		omega_age = omega_age + age_demog[i] * age_bite[i];
	}

	omega_age = 1 / omega_age;


	for (int i = 0; i < N_age; i++)
	{
		age_bite[i] = omega_age * age_bite[i];
	}


	//////////////////////////////////////////////////////
	// 5.2.10.2. Find age category closest to 20 yr old woman

	index_age_20 = 0;

	double age_diff = (age_mids[0] - 20.0*365.0)*(age_mids[0] - 20.0*365.0);

	for (int i = 1; i < N_age; i++)
	{
		if ((age_mids[i] - 20.0*365.0)*(age_mids[i] - 20.0*365.0) < age_diff)
		{
			age_diff = (age_mids[i] - 20.0*365.0)*(age_mids[i] - 20.0*365.0);
			index_age_20 = i;
		}
	}


	//////////////////////////////////////////////////////
	// 5.2.10.3. Find indexes for age limits of pregnancy

	index_preg_age_low = 0;

	age_diff = (age_mids[0] - theta.preg_age_low)*(age_mids[0] - theta.preg_age_low);

	for (int i = 1; i < N_age; i++)
	{
		if ((age_mids[i] - theta.preg_age_low)*(age_mids[i] - theta.preg_age_low) < age_diff)
		{
			age_diff = (age_mids[i] - theta.preg_age_low)*(age_mids[i] - theta.preg_age_low);
			index_preg_age_low = i;
		}
	}


	index_preg_age_high = 0;

	age_diff = (age_mids[0] - theta.preg_age_high)*(age_mids[0] - theta.preg_age_high);

	for (int i = 1; i < N_age; i++)
	{
		if ((age_mids[i] - theta.preg_age_high)*(age_mids[i] - theta.preg_age_high) < age_diff)
		{
			age_diff = (age_mids[i] - theta.preg_age_high)*(age_mids[i] - theta.preg_age_high);
			index_preg_age_high = i;
		}
	}


	//////////////////////////////////////////////////////
	// 5.2.10.4. Find indexes for age limits of occupational exposure

	index_occup_age_low = 0;

	age_diff = (age_mids[0] - theta.occup_age_low)*(age_mids[0] - theta.occup_age_low);

	for (int i = 1; i < N_age; i++)
	{
		if ((age_mids[i] - theta.occup_age_low)*(age_mids[i] - theta.occup_age_low) < age_diff)
		{
			age_diff = (age_mids[i] - theta.occup_age_low)*(age_mids[i] - theta.occup_age_low);
			index_occup_age_low = i;
		}
	}


	index_occup_age_high = 0;

	age_diff = (age_mids[0] - theta.occup_age_high)*(age_mids[0] - theta.occup_age_high);

	for (int i = 1; i < N_age; i++)
	{
		if ((age_mids[i] - theta.occup_age_high)*(age_mids[i] - theta.occup_age_high) < age_diff)
		{
			age_diff = (age_mids[i] - theta.occup_age_high)*(age_mids[i] - theta.occup_age_high);
			index_occup_age_high = i;
		}
	}


	////////////////////////////////////////
	// 5.2.10.5.   Heterogeneity in expsoure and age demographics.
	//             Weights for heterogeneity calculate using Gaussian-Hemite quadrature.
	//             The function also fills out the N_age*N_het matrix

	gauher(theta);

	for (int i = 0; i < N_age; i++)
	{
		for (int j = 0; j < N_het; j++)
		{
			w_gen_age_het[0][0][i][j] = 0.5*(1.0 - theta.P_occup)*age_demog[i]*w_het[j];
			w_gen_age_het[0][1][i][j] = 0.5*theta.P_occup*age_demog[i]*w_het[j];
			w_gen_age_het[0][2][i][j] = 0.5*age_demog[i]*w_het[j];

			x_gen_age_het[0][0][i][j] = age_bite[i]*x_het[j];
			x_gen_age_het[0][1][i][j] = age_bite[i]*x_het[j];
			x_gen_age_het[0][2][i][j] = age_bite[i]*x_het[j];
		}
	}

	for (int i = 0; i < N_age; i++)
	{
		for (int j = 0; j < N_het; j++)
		{
			w_gen_age_het[1][0][i][j] = 0.0;
			w_gen_age_het[1][1][i][j] = 0.0;
			w_gen_age_het[1][2][i][j] = 0.0;

			x_gen_age_het[1][0][i][j] = 0.0;
			x_gen_age_het[1][1][i][j] = 0.0;
			x_gen_age_het[1][2][i][j] = 0.0;

			if ((i > index_occup_age_low) && (i <= index_occup_age_high))
			{
				w_gen_age_het[1][1][i][j] = 0.5*theta.P_occup*age_demog[i]*w_het[j];

				x_gen_age_het[1][1][i][j] = age_bite[i]*x_het[j];
			}
		}
	}



	for (int j = 0; j < N_het; j++)
	{
		w_gen_het_at_birth[0][j] = 0.5*(1.0 - theta.P_occup)*w_het[j];
		w_gen_het_at_birth[1][j] = 0.5*theta.P_occup*w_het[j];
		w_gen_het_at_birth[2][j] = 0.5*w_het[j];
	}

	double denom_gij[N_mosq];

	for (int v = 0; v < N_mosq; v++)
	{
		denom_gij[v] = 0.0;

		for (int g = 0; g < N_gen; g++)
		{
			for (int i = 0; i < N_age; i++)
			{
				for (int j = 0; j < N_het; j++)
				{
					denom_gij[v] = denom_gij[v] + w_gen_age_het[v][g][i][j] * x_gen_age_het[v][g][i][j];
				}
			}
		}
	}

	for (int v = 0; v < N_mosq; v++)
	{
		for (int g = 0; g < N_gen; g++)
		{
			for (int i = 0; i < N_age; i++)
			{
				for (int j = 0; j < N_het; j++)
				{
					x_gen_age_het[v][g][i][j] = x_gen_age_het[v][g][i][j] / denom_gij[v];
				}
			}
		}
	}


	///////////////////////////////////////////////////
	// 5.2.10.6. Proportion of total population with
	//           occupational exposure

	denom_w_occup = 0.0;

	for (int g = 0; g < N_gen; g++)
	{
		for (int i = 0; i < N_age; i++)
		{
			for (int j = 0; j < N_het; j++)
			{
				denom_w_occup = denom_w_occup + w_gen_age_het[1][g][i][j];
			}
		}
	}


	//////////////////////////////////////////////////////////////
	// 5.2.10.7. Proportion with low G6PD activity

	theta.low_G6PD_activity[0] = (1.0 - theta.G6PD_prev)*0.5*(1.0 + erf( (3.0 - theta.mu_G6PD_nor) /(theta.sig_G6PD_nor*sqrt(2.0)) ) ) +
			                             theta.G6PD_prev*0.5*(1.0 + erf((3.0 - theta.mu_G6PD_def) / (theta.sig_G6PD_def*sqrt(2.0))));
	theta.low_G6PD_activity[1] = (1.0 - theta.G6PD_prev)*0.5*(1.0 + erf((3.0 - theta.mu_G6PD_nor) / (theta.sig_G6PD_nor*sqrt(2.0)))) +
			                             theta.G6PD_prev*0.5*(1.0 + erf((3.0 - theta.mu_G6PD_def) / (theta.sig_G6PD_def*sqrt(2.0))));
	theta.low_G6PD_activity[2] = (1.0 - theta.G6PD_prev)*(1.0 - theta.G6PD_prev)*0.5*(1.0 + erf((3.0 - theta.mu_G6PD_nor) / (theta.sig_G6PD_nor*sqrt(2.0)))) +
			                         2.0*(1.0 - theta.G6PD_prev)*theta.G6PD_prev*0.5*(1.0 + erf((3.0 - theta.mu_G6PD_het) / (theta.sig_G6PD_het*sqrt(2.0)))) +
		 	                                     theta.G6PD_prev*theta.G6PD_prev*0.5*(1.0 + erf((3.0 - theta.mu_G6PD_def) / (theta.sig_G6PD_def*sqrt(2.0))));



	//////////////////////////////////////////////////////////////
	// 5.2.10.8. Proportion receiving effective treatment

	treat_PQeff.resize(N_gen);
	for (int g = 0; g < N_gen; g++)
	{
		treat_PQeff[g].resize(N_age);
	}

	//////////////////////////////////////////////////////////////
	// 5.2.10.8.1. Case-management with blood-stage drugs

	if (theta.CM_regimen == 0)
	{
		for (int g = 0; g < N_gen; g++)
		{
			for (int i = 0; i < N_age; i++)
			{
				treat_PQeff[g][i] = 0.0;
			}
		}
	}


	//////////////////////////////////////////////////////////////
	// 5.2.10.8.2. Case-management with primaquine

	if (theta.CM_regimen == 1)
	{
		for (int g = 0; g < N_gen; g++)
		{
			for (int i = 0; i < N_age; i++)
			{
				treat_PQeff[g][i] = (1.0 - theta.low_G6PD_activity[g])*theta.CM_PQ_adhere*theta.CM_PQ_eff;
			}
		}

		// Account for low CYP2D6 activity if primaquine

		for (int g = 0; g < N_gen; g++)
		{
			for (int i = 0; i < N_age; i++)
			{
				treat_PQeff[g][i] = treat_PQeff[g][i] * (1.0 - theta.CYP2D6_prev);
			}
		}

		// Infants under the age of 6 months

		for (int g = 0; g < N_gen; g++)
		{
			for (int i = 0; i < N_age; i++)
			{
				if (age_bounds[1 + i] < 0.5*365.0)
				{
					treat_PQeff[g][i] = 0.0;
				}

				if ((age_bounds[i] < 0.5*365.0) && (age_bounds[1 + i] > 0.5*365.0))
				{
					treat_PQeff[g][i] = treat_PQeff[g][i] * (0.5*365.0 - age_bounds[i]) / (age_bounds[1 + i] - age_bounds[i]);
				}
			}
		}

		// Pregnant women

		for (int i = 0; i < N_age; i++)
		{
			if ((i >= index_preg_age_low) && (i <= index_preg_age_high))
			{
				treat_PQeff[2][i] = (1.0 - theta.P_preg)*treat_PQeff[2][i];
			}
		}
	}


	//////////////////////////////////////////////////////////////
	// 5.2.10.8.3. Case-management with tafenoquine

	if (theta.CM_regimen == 2)
	{
		for (int g = 0; g < N_gen; g++)
		{
			for (int i = 0; i < N_age; i++)
			{
				treat_PQeff[g][i] = (1.0 - theta.low_G6PD_activity[g])*theta.CM_PQ_adhere*theta.CM_PQ_eff;
			}
		}


		// Infants under the age of 6 months

		for (int g = 0; g < N_gen; g++)
		{
			for (int i = 0; i < N_age; i++)
			{
				if (age_bounds[1 + i] < 0.5*365.0)
				{
					treat_PQeff[g][i] = 0.0;
				}

				if ((age_bounds[i] < 0.5*365.0) && (age_bounds[1 + i] > 0.5*365.0))
				{
					treat_PQeff[g][i] = treat_PQeff[g][i] * (0.5*365.0 - age_bounds[i]) / (age_bounds[1 + i] - age_bounds[i]);
				}
			}
		}

		// Pregnant women

		for (int i = 0; i < N_age; i++)
		{
			if ((i >= index_preg_age_low) && (i <= index_preg_age_high))
			{
				treat_PQeff[2][i] = (1.0 - theta.P_preg)*treat_PQeff[2][i];
			}
		}
	}



}



////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//           //                                                                   //
//  5.2.11.  //  This calculates the non-seasonal equilibrium solution from a     //
//           //  comparable deterministic model. This will give an approximately  //
//           //  correct solution for a seasonal setting.                         //
//           //                                                                   //
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

void Population::pop_at_equil(Params& theta)
{
	cout << "Calculating equilibrium solution.....   " << endl;

	if (equil_setup_count > 0)
	{
		cout << "Repeat number = " << equil_setup_count << endl;
		cout << "(due to PQ/TQ case management at initialisation)" << endl;
	}

	//////////////////////////////////////////////////////////////
	// 5.2.11.1. Object for storing equilibrium solution

	yH_eq.resize(N_gen);
	for (int g = 0; g < N_gen; g++)
	{
		yH_eq[g].resize(N_age);
		for (int i = 0; i < N_age; i++)
		{
			yH_eq[g][i].resize(N_het);
			for (int j = 0; j < N_het; j++)
			{
				yH_eq[g][i][j].resize(K_max + 1);
				for (int k = 0; k < (K_max + 1); k++)
				{
					yH_eq[g][i][j][k].resize(N_H_comp);
				}
			}
		}
	}

	lam_eq.resize(N_gen);
	for (int g = 0; g < N_gen; g++)
	{
		lam_eq[g].resize(N_age);
		for (int i = 0; i < N_age; i++)
		{
			lam_eq[g][i].resize(N_het);
		}
	}

	A_par_eq.resize(N_gen);
	for (int g = 0; g < N_gen; g++)
	{
		A_par_eq[g].resize(N_age);
		for (int i = 0; i < N_age; i++)
		{
			A_par_eq[g][i].resize(N_het);
			for (int j = 0; j < N_het; j++)
			{
				A_par_eq[g][i][j].resize(K_max + 1);
			}
		}
	}

	A_par_eq_mean.resize(N_gen);
	for (int g = 0; g < N_gen; g++)
	{
		A_par_eq_mean[g].resize(N_age);
		for (int i = 0; i < N_age; i++)
		{
			A_par_eq_mean[g][i].resize(N_het);
		}
	}

	A_clin_eq.resize(N_gen);
	for (int g = 0; g < N_gen; g++)
	{
		A_clin_eq[g].resize(N_age);
		for (int i = 0; i < N_age; i++)
		{
			A_clin_eq[g][i].resize(N_het);
			for (int j = 0; j < N_het; j++)
			{
				A_clin_eq[g][i][j].resize(K_max + 1);
			}
		}
	}

	A_clin_eq_mean.resize(N_gen);
	for (int g = 0; g < N_gen; g++)
	{
		A_clin_eq_mean[g].resize(N_age);
		for (int i = 0; i < N_age; i++)
		{
			A_clin_eq_mean[g][i].resize(N_het);
		}
	}

	phi_LM_eq.resize(N_gen);
	for (int g = 0; g < N_gen; g++)
	{
		phi_LM_eq[g].resize(N_age);
		for (int i = 0; i < N_age; i++)
		{
			phi_LM_eq[g][i].resize(N_het);
			for (int j = 0; j < N_het; j++)
			{
				phi_LM_eq[g][i][j].resize(K_max + 1);
			}
		}
	}

	phi_D_eq.resize(N_gen);
	for (int g = 0; g < N_gen; g++)
	{
		phi_D_eq[g].resize(N_age);
		for (int i = 0; i < N_age; i++)
		{
			phi_D_eq[g][i].resize(N_het);
			for (int j = 0; j < N_het; j++)
			{
				phi_D_eq[g][i][j].resize(K_max + 1);
			}
		}
	}

	r_PCR_eq.resize(N_gen);
	for (int g = 0; g < N_gen; g++)
	{
		r_PCR_eq[g].resize(N_age);
		for (int i = 0; i < N_age; i++)
		{
			r_PCR_eq[g][i].resize(N_het);
			for (int j = 0; j < N_het; j++)
			{
				r_PCR_eq[g][i][j].resize(K_max + 1);
			}
		}
	}


	vector<vector<double>> MM;
	MM.resize(N_H_comp * (K_max + 1));
	for (int l = 0; l < N_H_comp * (K_max + 1); l++)
	{
		MM[l].resize(N_H_comp * (K_max + 1));
	}

	vector<double> bb(N_H_comp * (K_max + 1));
	vector<double> xx(N_H_comp * (K_max + 1));


	///////////////////////////////////////////////////////////
	// 5.2.11.2. PQ transmission matrices

	theta_HPZzero.resize(N_H_comp);
	for (int p = 0; p < (N_H_comp); p++)
	{
		theta_HPZzero[p].resize(N_gen);
		for (int g = 0; g < N_gen; g++)
		{
			theta_HPZzero[p][g].resize(N_age);
			for (int i = 0; i < N_age; i++)
			{
				theta_HPZzero[p][g][i].resize(N_het);
				for (int j = 0; j < N_het; j++)
				{
					theta_HPZzero[p][g][i][j].resize(K_max + 1);
				}
			}
		}
	}

	theta_HPZzero_weight.resize(N_gen);
	for (int g = 0; g < N_gen; g++)
	{
		theta_HPZzero_weight[g].resize(N_age);
		for (int i = 0; i < N_age; i++)
		{
			theta_HPZzero_weight[g][i].resize(N_het);
			for (int j = 0; j < N_het; j++)
			{
				theta_HPZzero_weight[g][i][j].resize(K_max + 1);
			}
		}
	}


	///////////////////////////////////////
	// Set to zero on first initialisation

	if ( equil_setup_count == 0 )
	{
		for (int p = 0; p < (N_H_comp); p++)
		{
			for (int g = 0; g < N_gen; g++)
			{
				for (int i = 0; i < N_age; i++)
				{
					for (int j = 0; j < N_het; j++)
					{
						for (int k = 0; k < (K_max + 1); k++)
						{
							theta_HPZzero[p][g][i][j][k] = 0.0;
						}
					}
				}
			}
		}

		for (int g = 0; g < N_gen; g++)
		{
			for (int i = 0; i < N_age; i++)
			{
				for (int j = 0; j < N_het; j++)
				{
					for (int k = 0; k < (K_max + 1); k++)
					{
						theta_HPZzero_weight[g][i][j][k] = 0.0;
					}
				}
			}
		}

	}


	//////////////////////////////////////////////////////////////
	// 5.2.11.3. Equilibrium force of infection
	//
	//   Only the contribution from mosquito bites

	for (int g = 0; g < N_gen; g++)
	{
		for (int i = 0; i < N_age; i++)
		{
			for (int j = 0; j < N_het; j++)
			{
				lam_eq[g][i][j] = theta.EIR_dom_equil*theta.bb*x_gen_age_het[0][g][i][j] + theta.EIR_occ_equil*theta.bb*x_gen_age_het[1][g][i][j];
			}
		}
	}


	//////////////////////////////////////////////////////////////
	// 5.2.11.4. Equilibrium number of batches of relapses

	vector<vector<vector<vector<double>>>> HH_eq;
	HH_eq.resize(N_gen);
	for (int g = 0; g < N_gen; g++)
	{
		HH_eq[g].resize(N_age);
		for (int i = 0; i < N_age; i++)
		{
			HH_eq[g][i].resize(N_het);
			for (int j = 0; j < N_het; j++)
			{
				HH_eq[g][i][j].resize(K_max + 1);
			}
		}
	}

	vector<double> HH_bb(K_max + 1);
	vector<double> HH_xx(K_max + 1);

	vector<vector<double>> HH_mat;
	HH_mat.resize(K_max + 1);
	for (int k = 0; k < (K_max + 1); k++)
	{
		HH_mat[k].resize(K_max + 1);
	}

	double HH_denom;

	for (int g = 0; g < N_gen; g++)
	{
		for (int j = 0; j < N_het; j++)
		{
			/////////////////////////////
			// Youngest age category

			for (int k = 0; k < (K_max + 1); k++)
			{
				HH_bb[k] = 0.0;
			}
			HH_bb[0] = w_gen_het_at_birth[g][j] * theta.mu_H;

			// Baseline mmodel for hypnozoite batches

			for (int k1 = 0; k1 < (K_max + 1); k1++)
			{
				for (int k2 = 0; k2 < (K_max + 1); k2++)
				{
					HH_mat[k1][k2] = lam_eq[g][0][j] * theta.H_MAT[k1][k2] + theta.gamma_L*theta.L_MAT[k1][k2] - theta.mu_H*theta.D_MAT[k1][k2] - r_age[0] * theta.D_MAT[k1][k2];
				}
			}


			// Accounting for transitions to zero due to PQ/TQ treatment

			for (int k2 = 0; k2 < (K_max + 1); k2++)
			{
				HH_mat[0][k2] = HH_mat[0][k2] + theta_HPZzero_weight[g][0][j][k2] * lam_eq[g][0][j];
			}
			for (int k2 = 0; k2 < (K_max); k2++)
			{
				HH_mat[k2 + 1][k2] = HH_mat[k2 + 1][k2] - theta_HPZzero_weight[g][0][j][k2]*lam_eq[g][0][j];
			}
			HH_mat[K_max][K_max] = HH_mat[K_max][K_max] - theta_HPZzero_weight[g][0][j][K_max] * lam_eq[g][0][j];

			for (int k2 = 0; k2 < (K_max + 1); k2++)
			{
				HH_mat[0][k2] = HH_mat[0][k2] + ((double)(k2))*theta_HPZzero_weight[g][0][j][k2] * theta.ff;
			}
			for (int k2 = 0; k2 < (K_max + 1); k2++)
			{
				HH_mat[k2][k2] = HH_mat[k2][k2] - ((double)(k2))*theta_HPZzero_weight[g][0][j][k2] * theta.ff;
			}


			// Apply linear algebra step

			inv_MM_bb(HH_mat, HH_bb, HH_xx, K_max + 1);

			for (int k = 0; k < (K_max + 1); k++)
			{
				HH_eq[g][0][j][k] = -HH_xx[k];
			}


			///////////////////////////////////
			// Older age categories

			for (int i = 1; i < N_age; i++)
			{
				for (int k = 0; k < (K_max + 1); k++)
				{
					HH_bb[k] = -r_age[i - 1] * HH_xx[k];
				}


				// Baseline mmodel for hypnozoite batches

				for (int k1 = 0; k1 < (K_max + 1); k1++)
				{
					for (int k2 = 0; k2 < (K_max + 1); k2++)
					{
						HH_mat[k1][k2] = lam_eq[g][i][j] * theta.H_MAT[k1][k2] + theta.gamma_L*theta.L_MAT[k1][k2] - theta.mu_H*theta.D_MAT[k1][k2] - r_age[i] * theta.D_MAT[k1][k2];
					}
				}


				// Accounting for transitions to zero due to PQ/TQ treatment

				for (int k2 = 0; k2 < (K_max + 1); k2++)
				{
					HH_mat[0][k2] = HH_mat[0][k2] + theta_HPZzero_weight[g][i][j][k2] * lam_eq[g][i][j];
				}
				for (int k2 = 0; k2 < (K_max); k2++)
				{
					HH_mat[k2 + 1][k2] = HH_mat[k2 + 1][k2] - theta_HPZzero_weight[g][i][j][k2]*lam_eq[g][i][j];
				}
				HH_mat[K_max][K_max] = HH_mat[K_max][K_max] - theta_HPZzero_weight[g][i][j][K_max]*lam_eq[g][i][j];

				for (int k2 = 0; k2 < (K_max + 1); k2++)
				{
					HH_mat[0][k2] = HH_mat[0][k2] + ((double)(k2))*theta_HPZzero_weight[g][i][j][k2] * theta.ff;
				}
				for (int k2 = 0; k2 < (K_max + 1); k2++)
				{
					HH_mat[k2][k2] = HH_mat[k2][k2] - ((double)(k2))*theta_HPZzero_weight[g][i][j][k2] * theta.ff;
				}


				// Apply linear algebra step

				inv_MM_bb(HH_mat, HH_bb, HH_xx, K_max + 1);

				for (int k = 0; k < (K_max + 1); k++)
				{
					HH_eq[g][i][j][k] = -HH_xx[k];
				}
			}
		}
	}

	//////////////////////////////////////////////////////////////
	// 5.2.11.5. Equilibrium levels of immunity

	double w_HH[K_max + 1];
	vector<double> A_par_vec(K_max + 1);
	vector<double> A_clin_vec(K_max + 1);

	vector<double> ODE_eq_vec(K_max + 1);

	vector<vector<double>> ODE_eq_MAT;
	ODE_eq_MAT.resize(K_max + 1);
	for (int k = 0; k < (K_max + 1); k++)
	{
		ODE_eq_MAT[k].resize(K_max + 1);
	}


	double G_VEC[K_max + 1];
	double LAM_MAT[K_max + 1][K_max + 1];
	double GAM_MAT[K_max + 1][K_max + 1];


	//////////////////////////////////////////////////////////////
	// 5.2.11.5.1. ANTI-PARASITE IMMUNITY

	for (int g = 0; g < N_gen; g++)
	{
		for (int j = 0; j < N_het; j++)
		{
			/////////////////////////////
			//                         //
			//  Youngest age category  //
			//                         //
			/////////////////////////////

			/////////////////////////////
			// Vector part of ODE

			G_VEC[0] = 0.0;
			for (int k = 1; k < (K_max + 1); k++)
			{
				G_VEC[k] = lam_eq[g][0][j] * HH_eq[g][0][j][k - 1] / HH_eq[g][0][j][k] + theta.ff*((double)k);
			}
			G_VEC[K_max] = G_VEC[K_max] + lam_eq[g][0][j];

			for (int k = 0; k < (K_max + 1); k++)
			{
				G_VEC[k] = G_VEC[k] / (G_VEC[k] * theta.u_par + 1.0);
			}


			for (int k = 0; k < (K_max + 1); k++)
			{
				ODE_eq_vec[k] = G_VEC[k];
			}


			/////////////////////////////
			// Matrix part of ODE

			for (int k1 = 0; k1 < (K_max + 1); k1++)
			{
				for (int k2 = 0; k2 < (K_max + 1); k2++)
				{
					LAM_MAT[k1][k2] = 0.0;
					GAM_MAT[k1][k2] = 0.0;
				}
			}

			for (int k = 0; k < K_max; k++)
			{
				LAM_MAT[k][k] = -1.0;
				LAM_MAT[k + 1][k] = HH_eq[g][0][j][k] / HH_eq[g][0][j][k + 1];
			}

			for (int k = 1; k < (K_max + 1); k++)
			{
				GAM_MAT[k][k] = -((double)k);
				GAM_MAT[k - 1][k] = ((double)k)*HH_eq[g][0][j][k] / HH_eq[g][0][j][k - 1];
			}

			for (int k1 = 0; k1 < (K_max + 1); k1++)
			{
				for (int k2 = 0; k2 < (K_max + 1); k2++)
				{
					ODE_eq_MAT[k1][k2] = -(lam_eq[g][0][j] * LAM_MAT[k1][k2] + theta.gamma_L*GAM_MAT[k1][k2] -
						theta.r_par*theta.D_MAT[k1][k2] - theta.mu_H*theta.D_MAT[k1][k2] - r_age[0] * theta.D_MAT[k1][k2]);
				}
			}


			/////////////////////////////
			// Solve matrix equation

			inv_MM_bb(ODE_eq_MAT, ODE_eq_vec, A_par_vec, K_max + 1);

			for (int k = 0; k < (K_max + 1); k++)
			{
				A_par_eq[g][0][j][k] = A_par_vec[k];
			}


			/////////////////////////////
			//                         //
			//  Older age categories   //
			//                         //
			/////////////////////////////

			for (int i = 1; i < N_age; i++)
			{
				/////////////////////////////
				// Vector part of ODE

				G_VEC[0] = 0.0;
				for (int k = 1; k < (K_max + 1); k++)
				{
					G_VEC[k] = lam_eq[g][i][j] * HH_eq[g][i][j][k - 1] / HH_eq[g][i][j][k] + theta.ff*((double)k);
				}
				G_VEC[K_max] = G_VEC[K_max] + lam_eq[g][i][j];

				for (int k = 0; k < (K_max + 1); k++)
				{
					G_VEC[k] = G_VEC[k] / (G_VEC[k] * theta.u_par + 1.0);
				}


				for (int k = 0; k < (K_max + 1); k++)
				{
					ODE_eq_vec[k] = G_VEC[k] + r_age[i - 1] * A_par_eq[g][i - 1][j][k] * HH_eq[g][i - 1][j][k] / HH_eq[g][i][j][k];
				}


				/////////////////////////////
				// Matrix part of ODE

				for (int k1 = 0; k1 < (K_max + 1); k1++)
				{
					for (int k2 = 0; k2 < (K_max + 1); k2++)
					{
						LAM_MAT[k1][k2] = 0.0;
						GAM_MAT[k1][k2] = 0.0;
					}
				}

				for (int k = 0; k < K_max; k++)
				{
					LAM_MAT[k][k] = -1.0;
					LAM_MAT[k + 1][k] = HH_eq[g][i][j][k] / HH_eq[g][i][j][k + 1];
				}

				for (int k = 1; k < (K_max + 1); k++)
				{
					GAM_MAT[k][k] = -((double)k);
					GAM_MAT[k - 1][k] = ((double)k)*HH_eq[g][i][j][k] / HH_eq[g][i][j][k - 1];
				}

				for (int k1 = 0; k1 < (K_max + 1); k1++)
				{
					for (int k2 = 0; k2 < (K_max + 1); k2++)
					{
						ODE_eq_MAT[k1][k2] = -(lam_eq[g][i][j] * LAM_MAT[k1][k2] + theta.gamma_L*GAM_MAT[k1][k2] -
							theta.r_par*theta.D_MAT[k1][k2] - theta.mu_H*theta.D_MAT[k1][k2] - r_age[i] * theta.D_MAT[k1][k2]);
					}
				}


				/////////////////////////////
				// Solve matrix equation

				inv_MM_bb(ODE_eq_MAT, ODE_eq_vec, A_par_vec, K_max + 1);

				for (int k = 0; k < (K_max + 1); k++)
				{
					A_par_eq[g][i][j][k] = A_par_vec[k];
				}

			}
		}
	}

	//////////////////////////////////////////////////////////////
	// 5.2.11.5.2. MEAN ANTI-PARASITE IMMUNITY

	for (int g = 0; g < N_gen; g++)
	{
		for (int i = 0; i < N_age; i++)
		{
			for (int j = 0; j < N_het; j++)
			{
				/////////////////////////////
				// Set up weights and rates

				HH_denom = 0.0;

				for (int k = 0; k < (K_max + 1); k++)
				{
					w_HH[k] = HH_eq[g][i][j][k];
					HH_denom = HH_denom + HH_eq[g][i][j][k];
				}

				for (int k = 0; k < (K_max + 1); k++)
				{
					w_HH[k] = w_HH[k] / HH_denom;
				}

				/////////////////////////////
				// Average level of immunity

				A_par_eq_mean[g][i][j] = 0.0;

				for (int k = 0; k < (K_max + 1); k++)
				{
					A_par_eq_mean[g][i][j] = A_par_eq_mean[g][i][j] + A_par_eq[g][i][j][k] * w_HH[k];
				}
			}
		}
	}

	//////////////////////////////////////////////////////////////
	// 5.2.11.5.3. ANTI-CLINICAL IMMUNITY

	for (int g = 0; g < N_gen; g++)
	{
		for (int j = 0; j < N_het; j++)
		{
			/////////////////////////////
			//                         //
			//  Youngest age category  //
			//                         //
			/////////////////////////////

			/////////////////////////////
			// Vector part of ODE

			G_VEC[0] = 0.0;
			for (int k = 1; k < (K_max + 1); k++)
			{
				G_VEC[k] = lam_eq[g][0][j] * (HH_eq[g][0][j][k - 1] / HH_eq[g][0][j][k]) + theta.ff*((double)k);
			}
			G_VEC[K_max] = G_VEC[K_max] + lam_eq[g][0][j];

			for (int k = 0; k < (K_max + 1); k++)
			{
				G_VEC[k] = G_VEC[k] / (G_VEC[k] * theta.u_clin + 1.0);
			}


			for (int k = 0; k < (K_max + 1); k++)
			{
				ODE_eq_vec[k] = G_VEC[k];
			}


			/////////////////////////////
			// Matrix part of ODE

			for (int k1 = 0; k1 < (K_max + 1); k1++)
			{
				for (int k2 = 0; k2 < (K_max + 1); k2++)
				{
					LAM_MAT[k1][k2] = 0.0;
					GAM_MAT[k1][k2] = 0.0;
				}
			}

			for (int k = 0; k < K_max; k++)
			{
				LAM_MAT[k][k] = -1.0;
				LAM_MAT[k + 1][k] = HH_eq[g][0][j][k] / HH_eq[g][0][j][k + 1];
			}

			for (int k = 1; k < (K_max + 1); k++)
			{
				GAM_MAT[k][k] = -((double)k);
				GAM_MAT[k - 1][k] = ((double)k)*HH_eq[g][0][j][k] / HH_eq[g][0][j][k - 1];
			}

			for (int k1 = 0; k1 < (K_max + 1); k1++)
			{
				for (int k2 = 0; k2 < (K_max + 1); k2++)
				{
					ODE_eq_MAT[k1][k2] = -(lam_eq[g][0][j] * LAM_MAT[k1][k2] + theta.gamma_L*GAM_MAT[k1][k2] -
						theta.r_clin*theta.D_MAT[k1][k2] - theta.mu_H*theta.D_MAT[k1][k2] - r_age[0] * theta.D_MAT[k1][k2]);
				}
			}


			/////////////////////////////
			// Solve matrix equation

			inv_MM_bb(ODE_eq_MAT, ODE_eq_vec, A_clin_vec, K_max + 1);

			for (int k = 0; k < (K_max + 1); k++)
			{
				A_clin_eq[g][0][j][k] = A_clin_vec[k];
			}


			/////////////////////////////
			//                         //
			//  Older age categories   //
			//                         //
			/////////////////////////////

			for (int i = 1; i < N_age; i++)
			{
				/////////////////////////////
				// Vector part of ODE

				G_VEC[0] = 0.0;
				for (int k = 1; k < (K_max + 1); k++)
				{
					G_VEC[k] = lam_eq[g][i][j] * HH_eq[g][i][j][k - 1] / HH_eq[g][i][j][k] + theta.ff*((double)k);
				}
				G_VEC[K_max] = G_VEC[K_max] + lam_eq[g][i][j];

				for (int k = 0; k < (K_max + 1); k++)
				{
					G_VEC[k] = G_VEC[k] / (G_VEC[k] * theta.u_clin + 1.0);
				}


				for (int k = 0; k < (K_max + 1); k++)
				{
					ODE_eq_vec[k] = G_VEC[k] + r_age[i - 1] * A_clin_eq[g][i - 1][j][k] * HH_eq[g][i - 1][j][k] / HH_eq[g][i][j][k];
				}


				/////////////////////////////
				// Matrix part of ODE

				for (int k1 = 0; k1 < (K_max + 1); k1++)
				{
					for (int k2 = 0; k2 < (K_max + 1); k2++)
					{
						LAM_MAT[k1][k2] = 0.0;
						GAM_MAT[k1][k2] = 0.0;
					}
				}

				for (int k = 0; k < K_max; k++)
				{
					LAM_MAT[k][k] = -1.0;
					LAM_MAT[k + 1][k] = HH_eq[g][i][j][k] / HH_eq[g][i][j][k + 1];
				}

				for (int k = 1; k < (K_max + 1); k++)
				{
					GAM_MAT[k][k] = -((double)k);
					GAM_MAT[k - 1][k] = ((double)k)*HH_eq[g][i][j][k] / HH_eq[g][i][j][k - 1];
				}

				for (int k1 = 0; k1 < (K_max + 1); k1++)
				{
					for (int k2 = 0; k2 < (K_max + 1); k2++)
					{
						ODE_eq_MAT[k1][k2] = -(lam_eq[g][i][j] * LAM_MAT[k1][k2] + theta.gamma_L*GAM_MAT[k1][k2] -
							theta.r_clin*theta.D_MAT[k1][k2] - theta.mu_H*theta.D_MAT[k1][k2] - r_age[i] * theta.D_MAT[k1][k2]);
					}
				}


				/////////////////////////////
				// Solve matrix equation

				inv_MM_bb(ODE_eq_MAT, ODE_eq_vec, A_clin_vec, K_max + 1);

				for (int k = 0; k < (K_max + 1); k++)
				{
					A_clin_eq[g][i][j][k] = A_clin_vec[k];
				}

			}
		}
	}

	//////////////////////////////////////////////////////////////
	// 5.2.11.5.4. MEAN ANTI-CLINICAL IMMUNITY

	for (int g = 0; g < N_gen; g++)
	{
		for (int i = 0; i < N_age; i++)
		{
			for (int j = 0; j < N_het; j++)
			{
				/////////////////////////////
				// Set up weights and rates

				HH_denom = 0.0;

				for (int k = 0; k < (K_max + 1); k++)
				{
					w_HH[k] = HH_eq[g][i][j][k];
					HH_denom = HH_denom + HH_eq[g][i][j][k];
				}

				for (int k = 0; k < (K_max + 1); k++)
				{
					w_HH[k] = w_HH[k] / HH_denom;
				}

				/////////////////////////////
				// Average level of immunity

				A_clin_eq_mean[g][i][j] = 0.0;

				for (int k = 0; k < (K_max + 1); k++)
				{
					A_clin_eq_mean[g][i][j] = A_clin_eq_mean[g][i][j] + A_clin_eq[g][i][j][k] * w_HH[k];
				}
			}
		}
	}

	//////////////////////////////////////////////////////////////
	// 5.2.11.5.5. ADD IN MATERNAL IMMUNITY

	for (int g = 0; g < N_gen; g++)
	{
		for (int j = 0; j < N_het; j++)
		{
			for (int i = 0; i < N_age; i++)
			{
				for (int k = 0; k < (K_max + 1); k++)
				{
					A_par_eq[g][i][j][k] = A_par_eq[g][i][j][k] + A_par_eq_mean[2][index_age_20][j] * theta.P_mat*exp(-age_mids[i] / theta.d_mat);
					A_clin_eq[g][i][j][k] = A_clin_eq[g][i][j][k] + A_clin_eq_mean[2][index_age_20][j] * theta.P_mat*exp(-age_mids[i] / theta.d_mat);
				}
			}
		}
	}

	//////////////////////////////////////////////////////////////
	// 5.2.11.5.6. EFFECTS OF IMMUNITY

	for (int g = 0; g < N_gen; g++)
	{
		for (int j = 0; j < N_het; j++)
		{
			for (int i = 0; i < N_age; i++)
			{
				for (int k = 0; k < (K_max + 1); k++)
				{
					phi_LM_eq[g][i][j][k] = theta.phi_LM_min + (theta.phi_LM_max - theta.phi_LM_min) / (1.0 + pow(A_par_eq[g][i][j][k] / theta.A_LM_50pc, theta.K_LM));
					phi_D_eq[g][i][j][k] = theta.phi_D_min + (theta.phi_D_max - theta.phi_D_min) / (1.0 + pow(A_clin_eq[g][i][j][k] / theta.A_D_50pc, theta.K_D));
					r_PCR_eq[g][i][j][k] = 1.0 / (theta.d_PCR_min + (theta.d_PCR_max - theta.d_PCR_min) / (1.0 + pow(A_par_eq[g][i][j][k] / theta.A_PCR_50pc, theta.K_PCR)));
				}
			}
		}
	}


	///////////////////////////////////////////////
	// 5.2.11.7. Equilibrium states (applying linear algebra steps)

	for (int g = 0; g < N_gen; g++)
	{
		for (int j = 0; j < N_het; j++)
		{
			///////////////////////////////////////////////
			// 5.2.11.7.1. Youngest age category

			MM_ij(g, 0, j, theta, r_age, MM,
				lam_eq, phi_LM_eq, phi_D_eq, r_PCR_eq, treat_PQeff);

			for (int k = 0; k < N_H_comp*(K_max + 1); k++)
			{
				bb[k] = 0.0;
			}

			bb[0] = -w_gen_het_at_birth[g][j] * theta.mu_H;

			inv_MM_bb(MM, bb, xx, N_H_comp*(K_max + 1));


			/////////////////////////////////////
			// Fill out equilibrium levels

			for (int c = 0; c < N_H_comp; c++)
			{
				for (int k = 0; k < (K_max + 1); k++)
				{
					yH_eq[g][0][j][k][c] = xx[c*(K_max + 1) + k];
				}
			}


			///////////////////////////////////////////////
			// 5.2.11.7.2. Older age categories

			for (int i = 1; i < N_age; i++)
			{
				MM_ij(g, i, j, theta, r_age, MM,
					lam_eq, phi_LM_eq, phi_D_eq, r_PCR_eq, treat_PQeff);

				for (int c = 0; c < N_H_comp*(K_max + 1); c++)
				{
					bb[c] = -r_age[i - 1] * xx[c];
				}

				inv_MM_bb(MM, bb, xx, N_H_comp * (K_max + 1));


				/////////////////////////////////////
				// Fill out equilibrium levels

				for (int c = 0; c < N_H_comp; c++)
				{
					for (int k = 0; k < (K_max + 1); k++)
					{
						yH_eq[g][i][j][k][c] = xx[c*(K_max + 1) + k];
					}
				}
			}
		}
	}

	/*
	cout << "Testing yH sum......" << endl;
	double yH_sum = 0.0;
	for (int g = 0; g < N_gen; g++)
	{
	for (int i = 0; i < N_age; i++)
	{
	for (int j = 0; j < N_het; j++)
	{
	for (int k = 0; k < (K_max + 1); k++)
	{
	for (int c = 0; c < N_H_comp; c++)
	{
	yH_sum = yH_sum + yH_eq[g][i][j][k][c];
	}
	}
	}
	}
	}
	cout << "yH_sum = " << yH_sum << endl;
	system("PAUSE");
	*/


	//////////////////////////////////////////////////////////////
	// 5.2.11.6. Proportion of transitionts to zero HPZ

	for (int g = 0; g < N_gen; g++)
	{
		for (int j = 0; j < N_het; j++)
		{
			for (int i = 0; i < N_age; i++)
			{
				for (int k = 0; k < (K_max + 1); k++)
				{
					theta_HPZzero[0][g][i][j][k] = phi_LM_eq[g][i][j][k] * phi_D_eq[g][i][j][k] * theta.CM_cover*treat_PQeff[g][i];
					theta_HPZzero[1][g][i][j][k] = phi_LM_eq[g][i][j][k] * phi_D_eq[g][i][j][k] * theta.CM_cover*treat_PQeff[g][i];
					theta_HPZzero[2][g][i][j][k] = phi_D_eq[g][i][j][k] * theta.CM_cover*treat_PQeff[g][i];
					theta_HPZzero[3][g][i][j][k] = 0.0;
					theta_HPZzero[4][g][i][j][k] = 0.0;
					theta_HPZzero[5][g][i][j][k] = 0.0;
				}
			}
		}
	}


	double theta_HPZzero_weight_denom;

	for (int g = 0; g < N_gen; g++)
	{
		for (int j = 0; j < N_het; j++)
		{
			for (int i = 0; i < N_age; i++)
			{
				for (int k = 0; k < (K_max + 1); k++)
				{
					theta_HPZzero_weight[g][i][j][k] = 0.0;
					theta_HPZzero_weight_denom = 0.0;

					for (int p = 0; p < N_H_comp; p++)
					{
						theta_HPZzero_weight[g][i][j][k] = theta_HPZzero_weight[g][i][j][k] + theta_HPZzero[p][g][i][j][k] * yH_eq[g][i][j][k][p];

						theta_HPZzero_weight_denom = theta_HPZzero_weight_denom + yH_eq[g][i][j][k][p];
					}

					theta_HPZzero_weight[g][i][j][k] = theta_HPZzero_weight[g][i][j][k] / theta_HPZzero_weight_denom;
				}
			}
		}
	}


	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	// 5.2.11.8. Calculate equilibrium of model in mosquitoes

	theta.lam_M[0] = 0.0;

	for (int g = 0; g < N_gen; g++)
	{
		for (int i = 0; i < N_age; i++)
		{
			for (int j = 0; j < N_het; j++)
			{
				for (int k = 0; k < (K_max + 1); k++)
				{
					theta.lam_M[0] = theta.lam_M[0] + x_gen_age_het[0][g][i][j]*theta.aa[0]*( theta.c_PCR*yH_eq[g][i][j][k][1] + theta.c_LM*yH_eq[g][i][j][k][2] +
						                                                                      theta.c_D*yH_eq[g][i][j][k][3]   + theta.c_T*yH_eq[g][i][j][k][4] );
				}
			}
		}
	}

	theta.lam_M[1] = 0.0;

	for (int i = 0; i < N_age; i++)
	{
		for (int j = 0; j < N_het; j++)
		{
			for (int k = 0; k < (K_max + 1); k++)
			{
				theta.lam_M[1] = theta.lam_M[1] + x_gen_age_het[1][1][i][j]*theta.aa[1]*( theta.c_PCR*yH_eq[1][i][j][k][1] + theta.c_LM*yH_eq[1][i][j][k][2] +
					                                                                      theta.c_D*yH_eq[1][i][j][k][3]   + theta.c_T*yH_eq[1][i][j][k][4] );
			}
		}
	}


	double I_M_star[N_mosq];
	for (int v = 0; v < N_mosq; v++)
	{
		I_M_star[v] = theta.lam_M[v] * exp(-theta.mu_M[v] * theta.tau_M[v]) / (theta.lam_M[v] + theta.mu_M[v]);
	}

	theta.mm_0[0] = theta.EIR_dom_equil / (theta.aa[0] * I_M_star[0]);
	theta.mm_0[1] = theta.EIR_occ_equil / (theta.aa[1] * I_M_star[1]);



	for (int v = 0; v < N_mosq; v++)
	{
		yM[v][0] = 2.0*theta.omega_larvae[v] * theta.mu_M[v] * theta.d_L_larvae*(1.0 + theta.d_pupae*theta.mu_P)*theta.mm_0[v];
		yM[v][1] = 2.0*theta.mu_M[v] * theta.d_L_larvae*(1.0 + theta.d_pupae*theta.mu_P)*theta.mm_0[v];
		yM[v][2] = 2.0*theta.d_pupae*theta.mu_M[v] * theta.mm_0[v];
		yM[v][3] = theta.mm_0[v] * (theta.mu_M[v] / (theta.lam_M[v] + theta.mu_M[v]));
		yM[v][4] = theta.mm_0[v] * (theta.lam_M[v] / (theta.lam_M[v] + theta.mu_M[v]))*(1.0 - exp(-theta.mu_M[v] * theta.tau_M[v]));
		yM[v][5] = theta.mm_0[v] * (theta.lam_M[v] / (theta.lam_M[v] + theta.mu_M[v]))*exp(-theta.mu_M[v] * theta.tau_M[v]);

		theta.Karry[v] = theta.mm_0[v] * 2.0*theta.d_L_larvae*theta.mu_M[v] * (1.0 + theta.d_pupae*theta.mu_P)*theta.gamma_larvae*(theta.omega_larvae[v] + 1.0) /
			(theta.omega_larvae[v] / (theta.mu_L0*theta.d_E_larvae) - 1.0 / (theta.mu_L0*theta.d_L_larvae) - 1.0);   // Larval carry capacity


		if (theta.Karry[v] < 1.0e-10) { theta.Karry[v] = 1.0e-10; } //
	}


	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	// 5.2.11.9. Proportion in each age and heterogeneity stratified category

	//////////////////////////////////////////
	// Fill out vector of lagged lam_M*S_M

	theta.lam_S_M_track.resize(N_mosq);

	for (int v = 0; v < N_mosq; v++)
	{
		for (int k = 0; k < theta.M_track; k++)
		{
			theta.lam_S_M_track[v].push_back(theta.lam_M[v] * yM[v][3]);
		}
	}


	//////////////////////////////////////////////////////////////
	// Calculate probability for each compartment

	yH_eq_cumsum.resize(N_gen);
	for (int g = 0; g < N_gen; g++)
	{
		yH_eq_cumsum[g].resize(N_age);
		for (int i = 0; i < N_age; i++)
		{
			yH_eq_cumsum[g][i].resize(N_het);
			for (int j = 0; j < N_het; j++)
			{
				yH_eq_cumsum[g][i][j].resize(N_H_comp*(K_max + 1));
			}
		}
	}

	for (int g = 0; g < N_gen; g++)
	{
		for (int i = 0; i < N_age; i++)
		{
			for (int j = 0; j < N_het; j++)
			{
				yH_eq_cumsum[g][i][j][0] = yH_eq[g][i][j][0][0];

				for (int k = 1; k < (K_max + 1); k++)
				{
					yH_eq_cumsum[g][i][j][k] = yH_eq_cumsum[g][i][j][k - 1] + yH_eq[g][i][j][k][0];
				}

				for (int c = 1; c < N_H_comp; c++)
				{
					yH_eq_cumsum[g][i][j][c*(K_max + 1)] = yH_eq_cumsum[g][i][j][c*(K_max + 1) - 1] + yH_eq[g][i][j][0][c];

					for (int k = 1; k < (K_max + 1); k++)
					{
						yH_eq_cumsum[g][i][j][c*(K_max + 1) + k] = yH_eq_cumsum[g][i][j][c*(K_max + 1) + k - 1] + yH_eq[g][i][j][k][c];
					}
				}
			}
		}
	}

	for (int g = 0; g < N_gen; g++)
	{
		for (int i = 0; i < N_age; i++)
		{
			for (int j = 0; j < N_het; j++)
			{
				for (int k = 0; k < N_H_comp*(K_max + 1); k++)
				{
					yH_eq_cumsum[g][i][j][k] = yH_eq_cumsum[g][i][j][k] / yH_eq_cumsum[g][i][j][N_H_comp*(K_max + 1) - 1];
				}
			}
		}
	}



	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	// 5.2.11.10. Output metrics to screen

	/////////////////////////////
	// Initial API estimate

	double API_init_estimate = 0.0;

	for (int g = 0; g < N_gen; g++)
	{
		for (int i = 0; i < N_age; i++)
		{
			for (int j = 0; j < N_het; j++)
			{
				for (int k = 0; k < (K_max + 1); k++)
				{
					API_init_estimate = API_init_estimate + HH_eq[g][i][j][k] * (lam_eq[g][i][j] + theta.ff*((double)k))*phi_LM_eq[g][i][j][k] * phi_D_eq[g][i][j][k] * theta.CM_cover;
				}
			}
		}
	}

	API_init_estimate = 1000 * 365 * API_init_estimate;


	/////////////////////////////
	// Prevalence estiamtes

	double PvPR_PCR_pop = 0.0;

	for (int g = 0; g < N_gen; g++)
	{
		for (int i = 0; i < N_age; i++)
		{
			for (int j = 0; j < N_het; j++)
			{
				for (int k = 0; k < (K_max + 1); k++)
				{
					PvPR_PCR_pop = PvPR_PCR_pop + yH_eq[g][i][j][k][1] + yH_eq[g][i][j][k][2] + yH_eq[g][i][j][k][3] + yH_eq[g][i][j][k][4];
				}
			}
		}
	}

	double PvPR_PCR_occ = 0.0, occ_denom = 0.0;

	for (int i = 0; i < N_age; i++)
	{
		for (int j = 0; j < N_het; j++)
		{
			for (int k = 0; k < (K_max + 1); k++)
			{
				PvPR_PCR_occ = PvPR_PCR_occ + yH_eq[1][i][j][k][1] + yH_eq[1][i][j][k][2] + yH_eq[1][i][j][k][3] + yH_eq[1][i][j][k][4];

				occ_denom = occ_denom + HH_eq[1][i][j][k];
			}
		}
	}

	PvPR_PCR_occ = PvPR_PCR_occ / occ_denom;


	/////////////////////////////
	// Output

	cout << endl;
	cout << "Initial API (all population) = " << API_init_estimate << endl;
	cout << endl;


	cout << "Domestic transmission............" << endl;

	cout << "EIR = " << 365.0*theta.aa[0] * yM[0][5] << endl;
	cout << "I_M = " << I_M_star[0] << "\t" << yM[0][5] / theta.mm_0[0] << endl;
	cout << "mm = " << theta.mm_0[0] << endl;
	cout << "lam_M = " << theta.lam_M[0] << endl;

	cout << endl;


	cout << "Population-level PCR prevalence = " << PvPR_PCR_pop << endl;

	cout << endl;


	if (theta.risk_occup > 0.0)
	{
		cout << "Occupational transmission............" << endl;

		cout << "EIR = " << 365.0*theta.aa[1] * yM[1][5] << endl;
		cout << "I_M = " << I_M_star[1] << "\t" << yM[1][5] / theta.mm_0[1] << endl;
		cout << "mm = " << theta.mm_0[1] << endl;
		cout << "lam_M = " << theta.lam_M[1] << endl;
		cout << endl;

		cout << "Occupational PCR prevalence = " << PvPR_PCR_occ << endl;
		cout << endl;

	}
	else
	{
		cout << "Occupational risk not simulated...." << endl;
	}




}





////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//           //                                                                   //
//  5.2.12.  //  Initialise a population of individuals at equilibrium            //
//           //                                                                   //
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

void Population::ind_at_equil(Params& theta)
{
	///////////////////////////////////////////////////////////////////////////
	// 5.2.12.1. Temporary objects for setting up individuals

	double rand_comp;
	double age_start, zeta_start;

	int i_index, j_index;


	float GMN_parm[(N_int)*(N_int + 3) / 2 + 1];
	float GMN_work[N_int];
	float GMN_zero[N_int];
	float zz_GMN[N_int];

	for (int k = 0; k < N_int; k++)
	{
		GMN_zero[k] = 0.0;
	}


	///////////////////////////////////////////////////////////////////////////
	// 5.2.12.2. Loop through and create N_pop individuals

	for (int n = 0; n < N_pop; n++)
	{
		//////////////////////////////////////////////////////////////////
		// 5.2.12.2.1. Assign age and heterogeneity

		age_start = genexp(theta.age_mean);

		while (age_start > theta.age_max)
		{
			age_start = genexp(theta.age_mean);
		}

		zeta_start = exp(gennor(-0.5*theta.sig_het*theta.sig_het, theta.sig_het));

		while (zeta_start > theta.het_max)
		{
			zeta_start = exp(gennor(-0.5*theta.sig_het*theta.sig_het, theta.sig_het));
		}


		//////////////////////////////////////////////////////////////////
		// 5.2.12.2.2. Find appropriate age and het compartment

		i_index = 0;

		for (int i = 0; i < N_age; i++)
		{
			if ((age_start > age_bounds[i]) && (age_start <= age_bounds[i + 1]))
			{
				i_index = i;
			}
		}


		j_index = 0;

		for (int j = 0; j < N_het; j++)
		{
			if ((zeta_start > x_het_bounds[j]) && (zeta_start < x_het_bounds[j + 1]))
			{
				j_index = j;
			}
		}

		//////////////////////////////////////////////////////////////////
		// 5.2.12.2.3. Construct a new individual

		double q_rand;

		Individual HH(age_start, zeta_start);

		HH.S = 0;
		HH.I_PCR = 0;
		HH.I_LM = 0;
		HH.I_D = 0;
		HH.T = 0;
		HH.P = 0;


		//////////////////////////////////////////////////////////////////
		// Assign gender and G6PD status
		// 0 = male (domestic); 1 = male (occupational); 2 = female
		// For males: 0 = normal (hemizygous); 1 = deficient (hemizygous)
		// For females: 0 = normal (homozygous); 1 = deficient (homozygous); 2 = deficient (heterozygous)

		if (genunf(0.0, 1.0) < 0.5)
		{
			if (genunf(0.0, 1.0) > theta.P_occup)
			{
				// male (domestic)
				HH.gender = 0;
			}
			else {
				// male (occupational)
				HH.gender = 1;
			}
		}
		else {
			// female
			HH.gender = 2;
		}

		if ((HH.gender == 0) || (HH.gender == 1))
		{

			if (genunf(0.0, 1.0) < theta.G6PD_prev)
			{
				// males, deficient (hemizygous)

				HH.G6PD_def = 1;
				HH.G6PD_activity = gennor(theta.mu_G6PD_def, theta.sig_G6PD_def);
			}
			else {
				// males, normal (hemizygous)

				HH.G6PD_def = 0;
				HH.G6PD_activity = gennor(theta.mu_G6PD_nor, theta.sig_G6PD_nor);
			}
		}
		else {

			q_rand = genunf(0.0, 1.0);

			if (q_rand <= theta.G6PD_prev*theta.G6PD_prev)
			{
				// females, deficient (homozygous)

				HH.G6PD_def = 1;
				HH.G6PD_activity = gennor(theta.mu_G6PD_def, theta.sig_G6PD_def);
			}

			if ((q_rand > theta.G6PD_prev*theta.G6PD_prev) && (q_rand <= theta.G6PD_prev*theta.G6PD_prev + 2 * theta.G6PD_prev*(1.0 - theta.G6PD_prev)))
			{
				// females, deficient (heterozygous)

				HH.G6PD_def = 2;
				HH.G6PD_activity = gennor(theta.mu_G6PD_het, theta.sig_G6PD_het);
			}

			if (q_rand > (theta.G6PD_prev*theta.G6PD_prev + 2 * theta.G6PD_prev*(1.0 - theta.G6PD_prev)))
			{
				// females, normal (homozygous)

				HH.G6PD_def = 0;
				HH.G6PD_activity = gennor(theta.mu_G6PD_nor, theta.sig_G6PD_nor);
			}
		}


		if (genunf(0.0, 1.0) < theta.CYP2D6_prev)
		{
			HH.CYP2D6 = 1;
		}
		else {
			HH.CYP2D6 = 0;
		}


		HH.T_last_BS = 1000000.0;


		///////////////////////////////////////////
		//  5.2.12.2.4. An indicator for pregnancy appropriate age
		//                Only women between ages 18 and 40 can be pregnant.

		HH.pregnant = 0;
		HH.preg_timer = 0.0;

		if (HH.gender == 2)
		{
			if ((HH.age > theta.preg_age_low) && (HH.age <= theta.preg_age_high))
			{
				HH.preg_age = 1;

				if (genunf(0.0, 1.0) < theta.P_preg)
				{
					HH.pregnant   = 1;
					HH.preg_timer = genunf(0.0, 269.0);
				}
			}
			else {
				HH.preg_age = 0;
			}
		}


		///////////////////////////////////////////////////////////////////
		// 5.2.12.2.5. Randomly assign a state according to equilibrium probabilities

		rand_comp = genunf(0.0, 1.0);

		if (rand_comp <= yH_eq_cumsum[HH.gender][i_index][j_index][0])
		{
			HH.S = 1;
			HH.Hyp = 0;
		}

		for (int k = 1; k < (K_max + 1); k++)
		{
			if ((rand_comp > yH_eq_cumsum[HH.gender][i_index][j_index][k - 1]) && (rand_comp <= yH_eq_cumsum[HH.gender][i_index][j_index][k]))
			{
				HH.S = 1;
				HH.Hyp = k;
			}
		}

		for (int c = 1; c < N_H_comp; c++)
		{
			for (int k = 0; k < (K_max + 1); k++)
			{
				if ((rand_comp > yH_eq_cumsum[HH.gender][i_index][j_index][c*(K_max + 1) + k - 1]) && (rand_comp <= yH_eq_cumsum[HH.gender][i_index][j_index][c*(K_max + 1) + k]))
				{
					if (c == 1) { HH.I_PCR = 1; }
					if (c == 2) { HH.I_LM = 1; }
					if (c == 3) { HH.I_D = 1; }
					if (c == 4) { HH.T = 1; }
					if (c == 5) { HH.P = 1; }

					HH.Hyp = k;
				}
			}
		}

		if (rand_comp > yH_eq_cumsum[HH.gender][i_index][j_index][N_H_comp*(K_max + 1) - 1])
		{
			HH.P = 1;
			HH.Hyp = K_max;
		}


		////////////////////////////////////////////
		//  5.2.12.2.6. Initialise immunity

		HH.A_par_mat  = A_par_eq_mean[HH.gender][index_age_20][j_index] * theta.P_mat*exp(-age_mids[i_index] / theta.d_mat);
		HH.A_clin_mat = A_clin_eq_mean[HH.gender][index_age_20][j_index] * theta.P_mat*exp(-age_mids[i_index] / theta.d_mat);

		HH.A_par = A_par_eq[HH.gender][i_index][j_index][HH.Hyp]   - HH.A_par_mat;
		HH.A_clin = A_clin_eq[HH.gender][i_index][j_index][HH.Hyp] - HH.A_clin_mat;


		HH.A_par_boost = 1;
		HH.A_clin_boost = 1;

		HH.A_par_timer = -1.0;
		HH.A_clin_timer = -1.0;

		HH.AQ8_proph = 0;
		HH.AQ8_proph_timer = -1.0;


		////////////////////////////////////////////
		//  5.2.12.2.7. A vector for storing lagged force of infection

		for (int k = 0; k < theta.H_track; k++)
		{
			HH.lam_bite_track.push_back(lam_eq[HH.gender][i_index][j_index]);
		}

		for (int k = 0; k < theta.H_track; k++)
		{
			HH.lam_rel_track.push_back(HH.Hyp*theta.ff);
		}


		////////////////////////////////////////////////////////
		// 5.2.12.2.8. Give individuals their life-long intervention access score

		for (int p = 0; p < N_int; p++)
		{
			for (int q = 0; q < N_int; q++)
			{
				theta.V_int_dummy[p][q] = theta.V_int[p][q];
			}
		}

		setgmn(GMN_zero, *theta.V_int_dummy, N_int, GMN_parm);

		genmn(GMN_parm, zz_GMN, GMN_work);

		for (int k = 0; k < N_int; k++)
		{
			HH.zz_int[k] = zz_GMN[k];
		}


		////////////////////////////////////////////////////////
		// 5.2.12.2.9. Individuals begin without interventions

		HH.LLIN = 0;
		HH.IRS = 0;
		HH.IVM = 0;

		for (int v = 0; v < N_mosq; v++)
		{
			HH.w_VC[v] = 1.0;
			HH.y_VC[v] = 1.0;
			HH.z_VC[v] = 0.0;
		}


		///////////////////////////////////////////////////
		////////////////////////////////////////////////////
		// 5.2.12.2.10. Add this person to the vector of people

		people.push_back(move(HH));
	}

	///////////////////////////////////////////////////////////
	// 5.2.12.3. Proportion of bites received by each person

	pi_n.resize(N_pop);
	for (int n = 0; n < N_pop; n++)
	{
		pi_n[n].resize(N_mosq);
	}

	for (int n = 0; n < N_pop; n++)
	{
		pi_n[n][0] = people[n].zeta_het*(1.0 - theta.rho_age*exp(-people[n].age*theta.age_0_inv));

		pi_n[n][1] = 0.0;

		if (people[n].gender == 1)
		{
			if ((people[n].age > theta.occup_age_low) && ((people[n].age <= theta.occup_age_high)))
			{
				pi_n[n][1] = people[n].zeta_het*(1.0 - theta.rho_age*exp(-people[n].age*theta.age_0_inv));
			}
		}
	}


	double SIGMA_PI[N_mosq];

	for (int v = 0; v < N_mosq; v++)
	{
		SIGMA_PI[v] = 0.0;
		for (int n = 0; n < N_pop; n++)
		{
			SIGMA_PI[v] = SIGMA_PI[v] + pi_n[n][v];
		}

		for (int n = 0; n < N_pop; n++)
		{
			pi_n[n][v] = pi_n[n][v] / SIGMA_PI[v];
		}
	}

	//////////////////////////////////////////////////////////
	// 5.2.12.4. Initialise population level quantitites

	//////////////////////////////////////////////
	// 5.2.12.4.1. Vector control quantities

	for (int v = 0; v < N_mosq; v++)
	{
		SUM_pi_w[v] = 0;
		SUM_pi_z[v] = 0;


		for (int n = 0; n < N_pop; n++)
		{
			SUM_pi_w[v] = SUM_pi_w[v] + pi_n[n][v] * people[n].w_VC[v];
			SUM_pi_z[v] = SUM_pi_z[v] + pi_n[n][v] * people[n].z_VC[v];
		}
	}


	for (int v = 0; v < N_mosq; v++)
	{
		W_VC[v] = 1.0 - theta.Q_0[v] + theta.Q_0[v] * SUM_pi_w[v];
		Z_VC[v] = theta.Q_0[v] * SUM_pi_z[v];

		delta_1_VC[v] = theta.delta_1 / (1.0 - Z_VC[v]);
		delta_VC[v] = delta_1_VC[v] + theta.delta_2;

		p_1_VC[v] = exp(-theta.mu_M[v] * delta_1_VC[v]);
		mu_M_VC[v] = -log(p_1_VC[v] * theta.p_2[v]) / delta_VC[v];

		Q_VC[v] = 1.0 - (1.0 - theta.Q_0[v]) / W_VC[v];

		aa_VC[v] = Q_VC[v] / delta_VC[v];

		exp_muM_tauM_VC[v] = exp(-mu_M_VC[v] * theta.tau_M[v]);
		beta_VC[v] = theta.eps_max[v] * mu_M_VC[v] / (exp(delta_VC[v] * mu_M_VC[v]) - 1.0);
	}


	/////////////////////////////////////////////////////////////////////////
	// 5.2.12.4.2. The rate at which person n is bitten by a single mosquito

	lam_n.resize(N_pop);
	for (int n = 0; n < N_pop; n++)
	{
		lam_n[n].resize(N_mosq);
	}


	for (int n = 0; n < N_pop; n++)
	{
		for (int v = 0; v < N_mosq; v++)
		{
			lam_n[n][v] = theta.aa[v] * pi_n[n][v];
		}
	}


	/////////////////////////////////////////////////////////////////////////
	// 5.2.12.4.3. Output an overview of the initial set up

	double S_ind = 0.0, I_PCR_ind = 0.0, I_LM_ind = 0.0, I_D_ind = 0.0, T_ind = 0.0, P_ind = 0.0;
	double S_eqq = 0.0, I_PCR_eqq = 0.0, I_LM_eqq = 0.0, I_D_eqq = 0.0, T_eqq = 0.0, P_eqq = 0.0;

	for (int n = 0; n < N_pop; n++)
	{
		S_ind = S_ind + people[n].S;
		I_PCR_ind = I_PCR_ind + people[n].I_PCR;
		I_LM_ind = I_LM_ind + people[n].I_LM;
		I_D_ind = I_D_ind + people[n].I_D;
		T_ind = T_ind + people[n].T;
		P_ind = P_ind + people[n].P;
	}

	for (int g = 0; g < N_gen; g++)
	{
		for (int i = 0; i < N_age; i++)
		{
			for (int j = 0; j < N_het; j++)
			{
				for (int k = 0; k < (K_max + 1); k++)
				{
					S_eqq = S_eqq + yH_eq[g][i][j][k][0];
					I_PCR_eqq = I_PCR_eqq + yH_eq[g][i][j][k][1];
					I_LM_eqq = I_LM_eqq + yH_eq[g][i][j][k][2];
					I_D_eqq = I_D_eqq + yH_eq[g][i][j][k][3];
					T_eqq = T_eqq + yH_eq[g][i][j][k][4];
					P_eqq = P_eqq + yH_eq[g][i][j][k][5];
				}
			}
		}
	}

	cout << "Initialising population of individuals...." << endl;
	cout << endl;


	cout << "S = " << ((double)S_ind) / N_pop << "\t" << S_eqq << endl;
	cout << "I_PCR = " << ((double)I_PCR_ind) / N_pop << "\t" << I_PCR_eqq << endl;
	cout << "I_LM = " << ((double)I_LM_ind) / N_pop << "\t" << I_LM_eqq << endl;
	cout << "I_D = " << ((double)I_D_ind) / N_pop << "\t" << I_D_eqq << endl;
	cout << "T = " << ((double)T_ind) / N_pop << "\t" << T_eqq << endl;
	cout << "P = " << ((double)P_ind) / N_pop << "\t" << P_eqq << endl;

	cout << endl;
}
