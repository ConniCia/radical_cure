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

#include "Population.h"

#include <cmath>


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  6.1.1. Derivatives of mosquito ODE model                                //
//                                                                          //
//  0.01721421 = 2*pi/365                                                   //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

void mosq_derivs(const double t, double(&yM)[N_mosq][N_M_comp], double(&dyMdt)[N_mosq][N_M_comp], Params& theta, Population& POP)
{
	double Karry_seas_inv[N_mosq];

	for (int v = 0; v < N_mosq; v++)
	{
		Karry_seas_inv[v] = 1.0 / (theta.Karry[v] * (theta.dry_seas[v] + (1 - theta.dry_seas[v])*pow(0.5 + 0.5*cos(0.01721421*(t - theta.t_peak_seas[v])), theta.kappa_seas[v]) / theta.denom_seas[v]));

		//Karry_seas_inv[v] = 1.0/theta.Karry[v];

		dyMdt[v][0] = POP.beta_VC[v] * (yM[v][3] + yM[v][4] + yM[v][5]) - yM[v][0] / theta.d_E_larvae - yM[v][0] * theta.mu_E0*(1.0 + (yM[v][0] + yM[v][1])*Karry_seas_inv[v]);
		dyMdt[v][1] = yM[v][0] / theta.d_E_larvae - yM[v][1] / theta.d_L_larvae - yM[v][1] * theta.mu_L0*(1.0 + theta.gamma_larvae*(yM[v][0] + yM[v][1])*Karry_seas_inv[v]);
		dyMdt[v][2] = yM[v][1] / theta.d_L_larvae - yM[v][2] / theta.d_pupae - yM[v][2] * theta.mu_P;
		dyMdt[v][3] = 0.5*yM[v][2] / theta.d_pupae - theta.lam_M[v] * yM[v][3] - POP.mu_M_VC[v] * yM[v][3];
		dyMdt[v][4] = +theta.lam_M[v] * yM[v][3] - theta.lam_S_M_track[v][0] * POP.exp_muM_tauM_VC[v] - POP.mu_M_VC[v] * yM[v][4];
		dyMdt[v][5] = +theta.lam_S_M_track[v][0] * POP.exp_muM_tauM_VC[v] - POP.mu_M_VC[v] * yM[v][5];
	}
}

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  6.1.2. Runge-Kutta 4 step updater for mosquito model                    //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

void mosq_rk4(const double t, const double t_step_mosq, double(&yM)[N_mosq][N_M_comp], Params& theta, Population& POP)
{
	double k1_yM[N_mosq][N_M_comp], k2_yM[N_mosq][N_M_comp], k3_yM[N_mosq][N_M_comp], k4_yM[N_mosq][N_M_comp], yM_temp[N_mosq][N_M_comp];


	//////////////////////////
	// step 1

	mosq_derivs(t, yM, k1_yM, theta, POP);


	//////////////////////////
	// step 2

	for (int v = 0; v < N_mosq; v++)
	{
		for (int k = 0; k < N_M_comp; k++)
		{
			yM_temp[v][k] = yM[v][k] + 0.5*t_step_mosq*k1_yM[v][k];
		}
	}

	mosq_derivs(t + 0.5*t_step_mosq, yM_temp, k2_yM, theta, POP);


	//////////////////////////
	// step 3

	for (int v = 0; v < N_mosq; v++)
	{
		for (int k = 0; k < N_M_comp; k++)
		{
			yM_temp[v][k] = yM[v][k] + 0.5*t_step_mosq*k2_yM[v][k];
		}
	}

	mosq_derivs(t + 0.5*t_step_mosq, yM_temp, k3_yM, theta, POP);


	//////////////////////////
	// step 4

	for (int v = 0; v < N_mosq; v++)
	{
		for (int k = 0; k < N_M_comp; k++)
		{
			yM_temp[v][k] = yM[v][k] + t_step_mosq * k3_yM[v][k];
		}
	}

	mosq_derivs(t + t_step_mosq, yM_temp, k4_yM, theta, POP);


	//////////////////////////
	// output

	for (int v = 0; v < N_mosq; v++)
	{
		for (int k = 0; k < N_M_comp; k++)
		{
			yM[v][k] = yM[v][k] + t_step_mosq * k1_yM[v][k] / 6.0 + t_step_mosq * k2_yM[v][k] / 3.0 + t_step_mosq * k3_yM[v][k] / 3.0 + t_step_mosq * k4_yM[v][k] / 6.0;
		}
	}

}


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  6.1.3. Update step for mosquitoes                                        //
//                                                                           // 
//       For every human step we take mosq_steps (=10) steps for mosquitoes  //
//       The smaller step size ensures that the ODE solver works smoothly.   //
//       Especially an issue for the larval stages                           //   
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void mosquito_step(double t, Params& theta, Population& POP)
{
	//////////////////////////////////
	// Set up mosquito state vector

	double yM[N_mosq][N_M_comp];

	for (int v = 0; v < N_mosq; v++)
	{
		for (int k = 0; k < N_M_comp; k++)
		{
			yM[v][k] = POP.yM[v][k];
		}
	}


	double t_step_mosq = (double(t_step)) / (double(mosq_steps));


	//////////////////////////////////
	// Force of infection on mosquitoes

	for (int v = 0; v < N_mosq; v++)
	{
		theta.lam_M[v] = 0.0;
	}

	for (int n = 0; n < POP.N_pop; n++)
	{
		for (int v = 0; v < N_mosq; v++)
		{
			theta.lam_M[v] = theta.lam_M[v] + POP.lam_n[n][v] * (theta.c_PCR*POP.people[n].I_PCR + theta.c_LM*POP.people[n].I_LM +
				theta.c_D*POP.people[n].I_D + theta.c_T*POP.people[n].T);
		}
	}


	//////////////////////////////////////
	// Carry out the mosq_steps

	for (int j = 0; j < mosq_steps; j++)
	{
		mosq_rk4(t, t_step_mosq, yM, theta, POP);

		for (int v = 0; v < N_mosq; v++)
		{
			theta.lam_S_M_track[v].push_back(theta.lam_M[v] * yM[v][3]);
			theta.lam_S_M_track[v].erase(theta.lam_S_M_track[v].begin());
		}
	}

	for (int v = 0; v < N_mosq; v++)
	{
		for (int k = 0; k < N_M_comp; k++)
		{
			POP.yM[v][k] = yM[v][k];
		}
	}
}
