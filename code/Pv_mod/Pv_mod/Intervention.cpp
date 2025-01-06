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

#include "Intervention.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>
#include "randlib.h"


////////////////////////////////////////////////////////////
//                                                        //
//  3.2.1.  Function declarations                         //
//                                                        //
////////////////////////////////////////////////////////////

double phi_inv(double pp, double mu, double sigma);
double G6PD_SD_BioSensor(double G6PD_true);
int CH_sample(double *xx, int nn);

/////////////////////////////////////
// Read intervention data from input files

Intervention::Intervention(const char *coverage_File)
{
	/////////////////////////////////////////////////////////////////////////
	//                                                                     //
	// 3.2.2.  Read in intervention coverage                               //
	//                                                                     //
	/////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////
	// 3.2.2.1.  Read in matrix of years and coverages

	cout << "Read in intervention coverage file............" << endl;

	std::ifstream coverage_Stream(coverage_File);

	if (coverage_Stream.fail())
	{
		std::cout << "Failure reading in data." << endl;
	}


	/////////////////////////////////////////////////////////////////////////
	// 3.2.2.2. Note that the matrix we read in may have variable size.
	//          We first read in the first line of the file, and then matching
	//          the spaces between columns to get the number of interventions 
	//          rounds.
	//
	//          There's very likely a much more effective way to do this.

	string line_read;
	int N_cov_rounds = 0;

	getline(coverage_Stream, line_read);


	std::string str(line_read);
	std::string str2(" ");

	std::size_t space_find = str.find(str2);


	// Go through the string finding spaces (" ")

	while (space_find < 10000)
	{
		space_find = line_read.find(str2);

		line_read.erase(0, space_find + 1);

		N_cov_rounds = N_cov_rounds + 1;
	}

	N_cov_rounds = N_cov_rounds - 2;

	cout << "Number of intervention rounds: " << N_cov_rounds << endl;

	coverage_Stream.close();


	/////////////////////////////////////////////////////////////////////////
	// 3.2.2.3. Note that the matrix we read in may have variable size.

	std::ifstream coverage_Stream2(coverage_File);

	string discard;

	vector<vector<double>> coverage;
	coverage.resize(0);

	vector<double> cov_read(N_cov_rounds + 1);

	for (int i = 0; i < 151; i++)
	{
		coverage_Stream2 >> discard;

		for (int j = 1; j < (N_cov_rounds + 2); j++)
		{
			coverage_Stream2 >> cov_read[j - 1];
		}

		coverage.push_back(cov_read);
	}


	/////////////////////////////////////////////////////////////////////////
	// 3.2.2.4. Fill out Intervention structure

	for (int j = 0; j < N_cov_rounds; j++)
	{
		//////////////////////////////////////////////////////////////
		// Intervention 0	
		// Front-line treatment - blood-stage drugs

		if ((coverage[0][j] > -0.5) && (coverage[1][j] > -0.5))
		{
			CM0_year.push_back(     coverage[0][j]*365.0 );
			CM0_cover.push_back(    coverage[1][j] );
			CM0_CQ_eff.push_back(   coverage[2][j] );
			CM0_CQ_proph.push_back( coverage[3][j] );
		}


		//////////////////////////////////////////////////////////////
		// Intervention 1
		// Front-line treatment - primaquine

		if ((coverage[0][j] > -0.5) && (coverage[4][j] > -0.5))
		{
			CM1_year.push_back(                 coverage[0][j]*365.0 );
			CM1_cover.push_back(                coverage[4][j] );
			CM1_CQ_eff.push_back(               coverage[5][j] );
			CM1_CQ_eff_wPQ.push_back(           coverage[6][j] );
			CM1_CQ_proph.push_back(             coverage[7][j] );
			CM1_PQ_eff.push_back(               coverage[8][j] );
			CM1_PQ_proph.push_back(             coverage[9][j] );
			CM1_PQ_adhere.push_back(            coverage[10][j] );
			CM1_PQ_lowage.push_back(            coverage[11][j] );
			CM1_PQ_G6PD_risk.push_back(   (int)(coverage[12][j]) );
			CM1_PQ_CYP2D6_risk.push_back( (int)(coverage[13][j]) );
			CM1_PQ_preg_risk.push_back(   (int)(coverage[14][j]) );
			CM1_G6PD_test.push_back(      (int)(coverage[15][j]) );
		}


		//////////////////////////////////////////////////////////////
		// Intervention 2
		// Front-line treatment - tafenoquine

		if ((coverage[0][j] > -0.5) && (coverage[16][j] > -0.5))
		{
			CM2_year.push_back(                 coverage[0][j]*365.0 );
			CM2_cover.push_back(                coverage[16][j] );
			CM2_CQ_eff.push_back(               coverage[17][j] );
			CM2_CQ_eff_wPQ.push_back(           coverage[18][j] );
			CM2_CQ_proph.push_back(             coverage[19][j] );
			CM2_PQ_eff.push_back(               coverage[20][j] );
			CM2_PQ_proph.push_back(             coverage[21][j] );
			CM2_PQ_adhere.push_back(            coverage[22][j] );
			CM2_PQ_lowage.push_back(            coverage[23][j] );
			CM2_PQ_G6PD_risk.push_back(   (int)(coverage[24][j]) );
			CM2_PQ_CYP2D6_risk.push_back( (int)(coverage[25][j]) );
			CM2_PQ_preg_risk.push_back(   (int)(coverage[26][j]) );
			CM2_TQ_eff.push_back(               coverage[27][j] );
			CM2_TQ_proph.push_back(             coverage[28][j] );
			CM2_TQ_adhere.push_back(            coverage[29][j] );
			CM2_TQ_lowage.push_back(            coverage[30][j] );
			CM2_TQ_G6PD_risk.push_back(   (int)(coverage[31][j]) );
			CM2_TQ_CYP2D6_risk.push_back( (int)(coverage[32][j]) );
			CM2_TQ_preg_risk.push_back(   (int)(coverage[33][j]) );
			CM2_G6PD_test.push_back(      (int)(coverage[34][j]) );
		}


		//////////////////////////////////////////////////////////////
		// Intervention 3
		// LLINs

		if ((coverage[0][j] > -0.5) && (coverage[35][j] > -0.5))
		{
			LLIN_year.push_back(  coverage[0][j]*365.0 );
			LLIN_cover.push_back( coverage[35][j] );
		}


		//////////////////////////////////////////////////////////////
		// Intervention 4
		// IRS

		if ((coverage[0][j] > -0.5) && (coverage[36][j] > -0.5))
		{
			IRS_year.push_back(  coverage[0][j]*365.0 );
			IRS_cover.push_back( coverage[36][j] );
		}


		//////////////////////////////////////////////////////////////
		// Intervention 5
		// MDA - chloroquine

		if ((coverage[0][j] > -0.5) && (coverage[37][j] > -0.5))
		{
			MDA0_year.push_back(     coverage[0][j]*365.0 );
			MDA0_cover.push_back(    coverage[37][j] );
			MDA0_CQ_eff.push_back(   coverage[38][j] );
			MDA0_CQ_proph.push_back( coverage[39][j] );
		}


		//////////////////////////////////////////////////////////////
		// Intervention 6
		// MDA - chloroquine plus primaquine

		if ((coverage[0][j] > -0.5) && (coverage[40][j] > -0.5))
		{
			MDA1_year.push_back(                 coverage[0][j]*365.0 );
			MDA1_cover.push_back(                coverage[40][j] );
			MDA1_CQ_eff.push_back(               coverage[41][j] );
			MDA1_CQ_eff_wPQ.push_back(           coverage[42][j] );
			MDA1_CQ_proph.push_back(             coverage[43][j] );
			MDA1_PQ_eff.push_back(               coverage[44][j] );
			MDA1_PQ_proph.push_back(             coverage[45][j] );
			MDA1_PQ_adhere.push_back(            coverage[46][j] );
			MDA1_PQ_lowage.push_back(            coverage[47][j] );
			MDA1_PQ_G6PD_risk.push_back(   (int)(coverage[48][j]) );
			MDA1_PQ_CYP2D6_risk.push_back( (int)(coverage[49][j]) );
			MDA1_PQ_preg_risk.push_back(   (int)(coverage[50][j]) );
			MDA1_G6PD_test.push_back(      (int)(coverage[51][j]) );
		}


		//////////////////////////////////////////////////////////////
		// Intervention 7
		// MDA - blood-stage drugs plus tafenoquine

		if ((coverage[0][j] > -0.5) && (coverage[52][j] > -0.5))
		{
			MDA2_year.push_back(                 coverage[0][j]*365.0 );
			MDA2_cover.push_back(                coverage[52][j] );
			MDA2_CQ_eff.push_back(               coverage[53][j] );
			MDA2_CQ_eff_wPQ.push_back(           coverage[54][j] );
			MDA2_CQ_proph.push_back(             coverage[55][j] );
			MDA2_PQ_eff.push_back(               coverage[56][j] );
			MDA2_PQ_proph.push_back(             coverage[57][j] );
			MDA2_PQ_adhere.push_back(            coverage[58][j] );
			MDA2_PQ_lowage.push_back(            coverage[59][j] );
			MDA2_PQ_G6PD_risk.push_back(   (int)(coverage[60][j]) );
			MDA2_PQ_CYP2D6_risk.push_back( (int)(coverage[61][j]) );
			MDA2_PQ_preg_risk.push_back(   (int)(coverage[62][j]) );
			MDA2_TQ_eff.push_back(               coverage[63][j] );
			MDA2_TQ_proph.push_back(             coverage[64][j] );
			MDA2_TQ_adhere.push_back(            coverage[65][j] );
			MDA2_TQ_lowage.push_back(            coverage[66][j] );
			MDA2_TQ_G6PD_risk.push_back(   (int)(coverage[67][j]) );
			MDA2_TQ_CYP2D6_risk.push_back( (int)(coverage[68][j]) );
			MDA2_TQ_preg_risk.push_back(   (int)(coverage[69][j]) );
			MDA2_G6PD_test.push_back(      (int)(coverage[70][j]) );
		}


		//////////////////////////////////////////////////////////////
		// Intervention 8
		// MSAT - blood-stage drugs only

		if ((coverage[0][j] > -0.5) && (coverage[71][j] > -0.5))
		{
			MSAT0_year.push_back(     coverage[0][j]*365.0 );
			MSAT0_cover.push_back(    coverage[71][j] );
			MSAT0_RDT_PCR.push_back(  coverage[72][j] );
			MSAT0_sens.push_back(     coverage[73][j] );
			MSAT0_CQ_eff.push_back(   coverage[74][j] ); 
			MSAT0_CQ_proph.push_back( coverage[75][j] );
		}


		//////////////////////////////////////////////////////////////
		// Intervention 9
		// MSAT - chloroquine plus primaquine

		if ((coverage[0][j] > -0.5) && (coverage[76][j] > -0.5))
		{
			MSAT1_year.push_back(                 coverage[0][j]*365.0 );
			MSAT1_cover.push_back(                coverage[76][j] );
			MSAT1_RDT_PCR.push_back(              coverage[77][j] );
			MSAT1_sens.push_back(                 coverage[78][j] );
			MSAT1_CQ_eff.push_back(               coverage[79][j] );
			MSAT1_CQ_eff_wPQ.push_back(           coverage[80][j] );
			MSAT1_CQ_proph.push_back(             coverage[81][j] );
			MSAT1_PQ_eff.push_back(               coverage[82][j] );
			MSAT1_PQ_proph.push_back(             coverage[83][j] );
			MSAT1_PQ_adhere.push_back(            coverage[84][j] );
			MSAT1_PQ_lowage.push_back(            coverage[85][j] );
			MSAT1_PQ_G6PD_risk.push_back(   (int)(coverage[86][j] ));
			MSAT1_PQ_CYP2D6_risk.push_back( (int)(coverage[87][j] ));
			MSAT1_PQ_preg_risk.push_back(   (int)(coverage[88][j] ));
			MSAT1_G6PD_test.push_back(      (int)(coverage[89][j] ));
		}


		//////////////////////////////////////////////////////////////
		// Intervention 10
		// MSAT - chloroquine plus tafenoquine

		if ((coverage[0][j] > -0.5) && (coverage[90][j] > -0.5))
		{
			MSAT2_year.push_back(                 coverage[0][j]*365.0 );
			MSAT2_cover.push_back(                coverage[90][j] );
			MSAT2_RDT_PCR.push_back(              coverage[91][j] );
			MSAT2_sens.push_back(                 coverage[92][j] );
			MSAT2_CQ_eff.push_back(               coverage[93][j] );
			MSAT2_CQ_eff_wPQ.push_back(           coverage[94][j] );
			MSAT2_CQ_proph.push_back(             coverage[95][j] );
			MSAT2_PQ_eff.push_back(               coverage[96][j] );
			MSAT2_PQ_proph.push_back(             coverage[97][j] );
			MSAT2_PQ_adhere.push_back(            coverage[98][j] );
			MSAT2_PQ_lowage.push_back(            coverage[99][j] );
			MSAT2_PQ_G6PD_risk.push_back(   (int)(coverage[100][j]) );
			MSAT2_PQ_CYP2D6_risk.push_back( (int)(coverage[101][j]) );
			MSAT2_PQ_preg_risk.push_back(   (int)(coverage[102][j]) );
			MSAT2_TQ_eff.push_back(               coverage[103][j] );
			MSAT2_TQ_proph.push_back(             coverage[104][j] );
			MSAT2_TQ_adhere.push_back(            coverage[105][j] );
			MSAT2_TQ_lowage.push_back(            coverage[106][j] );
			MSAT2_TQ_G6PD_risk.push_back(   (int)(coverage[107][j]) );
			MSAT2_TQ_CYP2D6_risk.push_back( (int)(coverage[108][j]) );
			MSAT2_TQ_preg_risk.push_back(   (int)(coverage[109][j]) );
			MSAT2_G6PD_test.push_back(      (int)(coverage[110][j]) );
		}


		//////////////////////////////////////////////////////////////
		// Intervention 11
		// STAT - chloroquine plus primaquine

		if ((coverage[0][j] > -0.5) && (coverage[111][j] > -0.5))
		{
			STAT1_year.push_back(                 coverage[0][j]*365.0 );
			STAT1_cover.push_back(                coverage[111][j] );
			STAT1_sens.push_back(                 coverage[112][j] );
			STAT1_spec.push_back(                 coverage[113][j] );
			STAT1_RDT_PCR.push_back(              coverage[114][j] );
			STAT1_CQ_eff.push_back(               coverage[115][j] );
			STAT1_CQ_eff_wPQ.push_back(           coverage[116][j] );
			STAT1_CQ_proph.push_back(             coverage[117][j] );
			STAT1_PQ_eff.push_back(               coverage[118][j] );
			STAT1_PQ_proph.push_back(             coverage[119][j] );
			STAT1_PQ_adhere.push_back(            coverage[120][j] );
			STAT1_PQ_lowage.push_back(            coverage[121][j] );
			STAT1_PQ_G6PD_risk.push_back(   (int)(coverage[122][j]) );
			STAT1_PQ_CYP2D6_risk.push_back( (int)(coverage[123][j]) );
			STAT1_PQ_preg_risk.push_back(   (int)(coverage[124][j]) );
			STAT1_G6PD_test.push_back(      (int)(coverage[125][j]) );
		}


		//////////////////////////////////////////////////////////////
		// Intervention 12
		// STAT - chloroquine plus tafenoquine

		if ((coverage[0][j] > -0.5) && (coverage[126][j] > -0.5))
		{
			STAT2_year.push_back(                 coverage[0][j]*365.0 );
			STAT2_cover.push_back(                coverage[126][j] );
			STAT2_sens.push_back(                 coverage[127][j] );
			STAT2_spec.push_back(                 coverage[128][j] );
			STAT2_RDT_PCR.push_back(              coverage[129][j] );
			STAT2_CQ_eff.push_back(               coverage[130][j] );
			STAT2_CQ_eff_wPQ.push_back(           coverage[131][j] );
			STAT2_CQ_proph.push_back(             coverage[132][j] );
			STAT2_PQ_eff.push_back(               coverage[133][j] );
			STAT2_PQ_proph.push_back(             coverage[134][j] );
			STAT2_PQ_adhere.push_back(            coverage[135][j] );
			STAT2_PQ_lowage.push_back(            coverage[136][j] );
			STAT2_PQ_G6PD_risk.push_back(   (int)(coverage[137][j]) );
			STAT2_PQ_CYP2D6_risk.push_back( (int)(coverage[138][j]) );
			STAT2_PQ_preg_risk.push_back(   (int)(coverage[139][j]) );
			STAT2_TQ_eff.push_back(               coverage[140][j] );
			STAT2_TQ_proph.push_back(             coverage[141][j] );
			STAT2_TQ_adhere.push_back(            coverage[142][j] );
			STAT2_TQ_lowage.push_back(            coverage[143][j] );
			STAT2_TQ_G6PD_risk.push_back(   (int)(coverage[144][j]) );
			STAT2_TQ_CYP2D6_risk.push_back( (int)(coverage[145][j]) );
			STAT2_TQ_preg_risk.push_back(   (int)(coverage[146][j]) );
			STAT2_G6PD_test.push_back(      (int)(coverage[147][j]) );
		}


		//////////////////////////////////////////////////////////////
		// Intervention 12
		// IVM

		if ((coverage[0][j] > -0.5) && (coverage[148][j] > -0.5))
		{
			IVM_year.push_back(      coverage[0][j]*365.0 );
			IVM_cover.push_back(     coverage[148][j] );
			d_IVM_0.push_back(       coverage[149][j] );
			IVM_half_life.push_back( coverage[150][j] );
		}
	}

}


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  3.2.3.  Vector control distribution                                     //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

void Intervention::distribute(double t, Params& theta, Population& POP)
{
	double QQ;

	bool MSAT_pos;
	bool STAT_pos;


	//////////////////////////////////////////////////////////
	// Intervention 0: case management with chloroquine

	for (int m = 0; m < CM0_year.size(); m++)
	{
		if ((t > CM0_year[m] - 0.5*t_step) &&
			(t < CM0_year[m] + 0.51*t_step))
		{
			cout << "New CQ case management regimen" << endl;

			theta.CM_regimen = 0;

			theta.CM_cover    = CM0_cover[m];
			theta.CM_CQ_eff   = CM0_CQ_eff[m];
			theta.CM_CQ_proph = CM0_CQ_proph[m];

			theta.CM_CQ_coveff = CM0_cover[m]* CM0_CQ_eff[m];

			theta.r_P = 1.0/ theta.CM_CQ_proph;
		}
	}


	//////////////////////////////////////////////////////////
	// Intervention 1: case management with chloroquine plus primaquine

	for (int m = 0; m < CM1_year.size(); m++)
	{
		if ((t > CM1_year[m] - 0.5*t_step) &&
			(t < CM1_year[m] + 0.51*t_step))
		{
			cout << "New front-line PQ treatment" << endl;

			theta.CM_regimen = 1;

			theta.CM_cover          = CM1_cover[m];
			theta.CM_CQ_eff         = CM1_CQ_eff[m];
			theta.CM_CQ_eff_wPQ     = CM1_CQ_eff_wPQ[m];
			theta.CM_CQ_proph       = CM1_CQ_proph[m];
			theta.CM_PQ_eff         = CM1_PQ_eff[m];
			theta.CM_PQ_proph       = CM1_PQ_proph[m];
			theta.CM_PQ_adhere      = CM1_PQ_adhere[m];
			theta.CM_PQ_lowage      = CM1_PQ_lowage[m];
			theta.CM_PQ_G6PD_risk   = CM1_PQ_G6PD_risk[m];
			theta.CM_PQ_CYP2D6_risk = CM1_PQ_CYP2D6_risk[m];
			theta.CM_PQ_preg_risk   = CM1_PQ_preg_risk[m];
			theta.CM_G6PD_test      = CM1_G6PD_test[m];

			theta.CM_CQ_coveff = CM1_cover[m]*CM1_CQ_eff[m];
			theta.CM_PQ_coveff = CM1_cover[m]*CM1_PQ_eff[m]*CM1_PQ_adhere[m];

			theta.r_P = 1.0/theta.CM_CQ_proph;
		}
	}


	//////////////////////////////////////////////////////////
	// Intervention 2: case management with chloroquine plus tafenoquine

	for (int m = 0; m < CM2_year.size(); m++)
	{
		if ((t > CM2_year[m] - 0.5*t_step) &&
			(t < CM2_year[m] + 0.51*t_step))
		{
			cout << "New front-line TQ treatment" << endl;

			theta.CM_regimen = 2;

			theta.CM_cover          = CM2_cover[m];
			theta.CM_CQ_eff         = CM2_CQ_eff[m];
			theta.CM_CQ_eff_wPQ     = CM2_CQ_eff_wPQ[m];
			theta.CM_CQ_proph       = CM2_CQ_proph[m];
			theta.CM_PQ_eff         = CM2_PQ_eff[m];
			theta.CM_PQ_proph       = CM2_PQ_proph[m];
			theta.CM_PQ_adhere      = CM2_PQ_adhere[m];
			theta.CM_PQ_lowage      = CM2_PQ_lowage[m];
			theta.CM_PQ_G6PD_risk   = CM2_PQ_G6PD_risk[m];
			theta.CM_PQ_CYP2D6_risk = CM2_PQ_CYP2D6_risk[m];
			theta.CM_PQ_preg_risk   = CM2_PQ_preg_risk[m];
			theta.CM_TQ_eff         = CM2_TQ_eff[m];
			theta.CM_TQ_proph       = CM2_TQ_proph[m];
			theta.CM_TQ_adhere      = CM2_TQ_adhere[m];
			theta.CM_TQ_lowage      = CM2_TQ_lowage[m];
			theta.CM_TQ_G6PD_risk   = CM2_TQ_G6PD_risk[m];
			theta.CM_TQ_CYP2D6_risk = CM2_TQ_CYP2D6_risk[m];
			theta.CM_TQ_preg_risk   = CM2_TQ_preg_risk[m];
			theta.CM_G6PD_test      = CM2_G6PD_test[m];


			theta.CM_CQ_coveff = CM2_cover[m]*CM2_CQ_eff[m];
			theta.CM_PQ_coveff = CM2_cover[m]*CM2_PQ_eff[m]*CM2_PQ_adhere[m];
			theta.CM_TQ_coveff = CM2_cover[m]*CM2_TQ_eff[m]*CM2_TQ_adhere[m];

			theta.r_P = 1.0/fmax( theta.CM_CQ_proph, theta.CM_TQ_proph);
		}
	}


	//////////////////////////////////////////////////////////
	// Intervention 3: LLINS

	for (int m = 0; m < LLIN_year.size(); m++)
	{
		if ((t > LLIN_year[m] - 0.5*t_step) &&
			(t < LLIN_year[m] + 0.51*t_step))
		{
			cout << "LLIN distribution" << endl;

			try
			{
				QQ = phi_inv(LLIN_cover[m], 0.0, sqrt(1.0 + theta.sig_round_LLIN*theta.sig_round_LLIN));
			}
			catch (const char* e)
			{
				std::cerr << e << std::endl;
				exit(1);
			}

			for (int n = 0; n < POP.N_pop; n++)
			{
				if (gennor(POP.people[n].zz_int[0], theta.sig_round_LLIN) < QQ)
				{
					POP.people[n].LLIN = 1;
					POP.people[n].LLIN_age = 0.0;

					for (int v = 0; v < N_mosq; v++)
					{
						POP.people[n].d_LLIN[v] = theta.d_LLIN_0[v];
						POP.people[n].r_LLIN[v] = theta.r_LLIN_0[v];
						POP.people[n].s_LLIN[v] = 1.0 - POP.people[n].d_LLIN[v] - POP.people[n].r_LLIN[v];
					}
				}
			}
		}
	}


	//////////////////////////////////////////////////////////
	// Intervention 4: IRS

	for (int m = 0; m < IRS_year.size(); m++)
	{
		if ((t > IRS_year[m] - 0.5*t_step) &&
			(t < IRS_year[m] + 0.51*t_step))
		{
			cout << "IRS distribution" << endl;

			try
			{
				QQ = phi_inv(IRS_cover[m], 0.0, sqrt(1.0 + theta.sig_round_IRS*theta.sig_round_IRS));
			}
			catch (const char* e)
			{
				std::cerr << e << std::endl;
				exit(1);
			}

			for (int n = 0; n < POP.N_pop; n++)
			{
				if (gennor(POP.people[n].zz_int[1], theta.sig_round_IRS) < QQ)
				{
					POP.people[n].IRS = 1;
					POP.people[n].IRS_age = 0.0;

					for (int v = 0; v < N_mosq; v++)
					{
						POP.people[n].d_IRS[v] = theta.d_IRS_0[v];
						POP.people[n].r_IRS[v] = theta.r_IRS_0[v];
						POP.people[n].s_IRS[v] = 1.0 - POP.people[n].d_IRS[v] - POP.people[n].r_IRS[v];
					}
				}
			}
		}
	}


	//////////////////////////////////////////////////////////
	// Intervention 5: MDA (CQ)

	for (int m = 0; m < MDA0_year.size(); m++)
	{
		if ((t > MDA0_year[m] - 0.5*t_step) &&
			(t < MDA0_year[m] + 0.51*t_step))
		{
			cout << "MDA (CQ) distribution" << endl;

			////////////////////////////////////////////
			// Update parameters

			theta.MDA0_cover    = MDA0_cover[m];
			theta.MDA0_CQ_eff   = MDA0_CQ_eff[m];
			theta.MDA0_CQ_proph = MDA0_CQ_proph[m];


			////////////////////////////////////////////
			// Implement intervention

			try
			{
				QQ = phi_inv(theta.MDA0_cover, 0.0, sqrt(1.0 + theta.sig_round_MDA*theta.sig_round_MDA));
			}
			catch (const char* e)
			{
				std::cerr << e << std::endl;
				exit(1);
			}

			for (int n = 0; n < POP.N_pop; n++)
			{
				POP.people[n].CQ_treat  = 0;
				POP.people[n].PQ_treat  = 0;
				POP.people[n].TQ_treat  = 0;
				POP.people[n].G6PD_test = 0;


				if (gennor(POP.people[n].zz_int[2], theta.sig_round_MDA) < QQ)
				{
					POP.people[n].CQ_treat = 1;
					POP.people[n].CQ_effective = 0;

					if (genunf(0.0, 1.0) < theta.MDA0_CQ_eff)
					{
						POP.people[n].CQ_effective = 1;
					}

					/////////////////////////////////////////////////////////////////////
					// ACTION: administer blood-stage drug

					if (POP.people[n].CQ_effective == 1 )
					{
						if (POP.people[n].S     == 1) { POP.people[n].S = 0;     POP.people[n].P = 1; }
						if (POP.people[n].I_PCR == 1) { POP.people[n].I_PCR = 0; POP.people[n].P = 1; }
						if (POP.people[n].I_LM  == 1) { POP.people[n].I_LM = 0;  POP.people[n].P = 1; }
						if (POP.people[n].I_D   == 1) { POP.people[n].I_D = 0;   POP.people[n].T = 1; }
					}
					else
					{
						if (POP.people[n].S     == 1) { POP.people[n].S = 0;     POP.people[n].P = 1; }
						if (POP.people[n].I_PCR == 1) { POP.people[n].I_PCR = 0; POP.people[n].P = 1; }
						if (POP.people[n].I_LM  == 1) { POP.people[n].I_LM = 0;  POP.people[n].P = 1; }
						if (POP.people[n].I_D   == 1) { POP.people[n].I_D = 0;   POP.people[n].I_LM = 1; }
					}
				}
			}
		}
	}


	//////////////////////////////////////////////////////////
	// Intervention 6: MDA (chloroquine and PQ)

	for (int m = 0; m < MDA1_year.size(); m++)
	{
		if ((t > MDA1_year[m] - 0.5*t_step) &&
			(t < MDA1_year[m] + 0.51*t_step))
		{
			cout << "MDA (CQ+PQ) distribution" << endl;

			////////////////////////////////////////////
			// Update parameters

			theta.MDA1_cover          = MDA1_cover[m];
			theta.MDA1_CQ_eff         = MDA1_CQ_eff[m];
			theta.MDA1_CQ_eff_wPQ     = MDA1_CQ_eff_wPQ[m];
			theta.MDA1_CQ_proph       = MDA1_CQ_proph[m];
			theta.MDA1_PQ_eff         = MDA1_PQ_eff[m];
			theta.MDA1_PQ_proph       = MDA1_PQ_proph[m];
			theta.MDA1_PQ_adhere	  = MDA1_PQ_adhere[m];
			theta.MDA1_PQ_lowage      = MDA1_PQ_lowage[m];
			theta.MDA1_PQ_G6PD_risk   = MDA1_PQ_G6PD_risk[m];
			theta.MDA1_PQ_CYP2D6_risk = MDA1_PQ_CYP2D6_risk[m];
			theta.MDA1_PQ_preg_risk   = MDA1_PQ_preg_risk[m];
			theta.MDA1_G6PD_test      = MDA1_G6PD_test[m];
			
			////////////////////////////////////////////
			// Implement intervention

			try
			{
				QQ = phi_inv(theta.MDA1_cover, 0.0, sqrt(1.0 + theta.sig_round_MDA*theta.sig_round_MDA));
			}
			catch (const char* e)
			{
				std::cerr << e << std::endl;
				exit(1);
			}

			for (int n = 0; n < POP.N_pop; n++)
			{
				POP.people[n].CQ_treat  = 0;
				POP.people[n].PQ_treat  = 0;
				POP.people[n].TQ_treat  = 0;
				POP.people[n].G6PD_test = 0;


				if (gennor(POP.people[n].zz_int[3], theta.sig_round_MDA) < QQ)
				{				
					/////////////////////////////////////////////////////////////////////
					// Is PQ administered?

					POP.people[n].PQ_treat = 1;
					POP.people[n].G6PD_test = 0;

					/////////////////////////////////////////////////////////////////////
					// Exclude PQ because of young age

					if (POP.people[n].age < theta.MDA1_PQ_lowage)
					{
						POP.people[n].PQ_treat = 0;
					}

					/////////////////////////////////////////////////////////////////////
					// Exclude PQ because of pregancy

					if ((theta.MDA1_PQ_preg_risk == 1) && (POP.people[n].pregnant == 1))
					{
						POP.people[n].PQ_treat = 0;
					}

					//////////////////////////////////////////////////////////////////////////////////////////
					// Is G6PD testing administered ? If so, count test

					if ((POP.people[n].PQ_treat == 1) && (theta.MDA1_G6PD_test == 1))
					{
						POP.people[n].G6PD_test = 1;
					}

					/////////////////////////////////////////////////////////////////////
					// Exclude PQ because of G6PD deficiency if there is testing

					if ((POP.people[n].G6PD_test == 1) && (theta.MDA1_PQ_G6PD_risk == 1))
					{
						POP.people[n].G6PD_read = G6PD_SD_BioSensor( POP.people[n].G6PD_activity );

						if (POP.people[n].G6PD_read < 3.0)
						{
							POP.people[n].PQ_treat = 0;
						}
					}

					/////////////////////////////////////////////////////////////////////
					// Is PQ effective?

					POP.people[n].PQ_effective = 1;

					if (genunf(0.0, 1.0) > theta.MDA1_PQ_eff)
					{
						POP.people[n].PQ_effective = 0;
					}

					/////////////////////////////////////////////////////////////////////
					// Is PQ adhered to?

					if (genunf(0.0, 1.0) < (1 - theta.MDA1_PQ_adhere))
					{
						POP.people[n].PQ_effective = 0;
					}

					/////////////////////////////////////////////////////////////////////
					// Is PQ metabolised?

					if ((theta.MDA1_PQ_CYP2D6_risk == 1) && (POP.people[n].CYP2D6 == 1))
					{
						POP.people[n].PQ_effective = 0;
					}


					/////////////////////////////////////////////////////////////////////
					// Was there PQ overtreatment?

					if ((POP.people[n].PQ_treat == 1) && (POP.people[n].Hyp == 0))
					{
						POP.people[n].PQ_overtreat = 1;
					}

					if ((POP.people[n].PQ_treat == 1) && (POP.people[n].T_last_BS > 270.0))
					{
						POP.people[n].PQ_overtreat_9m = 1;
					}


					/////////////////////////////////////////////////////
					// Chloroquine is always administered.
					// Is chloroquine effective

					POP.people[n].CQ_treat = 1;
					POP.people[n].CQ_effective = 0;

					if (POP.people[n].PQ_treat == 0)
					{
						if (genunf(0.0, 1.0) < theta.MDA1_CQ_eff)
						{
							POP.people[n].CQ_effective = 1;
						}
					}
					if (POP.people[n].PQ_treat == 1)
					{
						if (genunf(0.0, 1.0) < theta.MDA1_CQ_eff_wPQ)
						{
							POP.people[n].CQ_effective = 1;
						}
					}


					/////////////////////////////////////////////////////////////////////
					// ACTION: administer chloroquine
									   
					if (POP.people[n].CQ_effective == 1)
					{
						if (POP.people[n].S     == 1) { POP.people[n].S = 0;     POP.people[n].P = 1; }
						if (POP.people[n].I_PCR == 1) { POP.people[n].I_PCR = 0; POP.people[n].P = 1; }
						if (POP.people[n].I_LM  == 1) { POP.people[n].I_LM = 0;  POP.people[n].P = 1; }
						if (POP.people[n].I_D   == 1) { POP.people[n].I_D = 0;   POP.people[n].T = 1; }
					}
					else
					{
						if (POP.people[n].S     == 1) { POP.people[n].S = 0;     POP.people[n].P = 1; }
						if (POP.people[n].I_PCR == 1) { POP.people[n].I_PCR = 0; POP.people[n].P = 1; }
						if (POP.people[n].I_LM  == 1) { POP.people[n].I_LM = 0;  POP.people[n].P = 1; }
						if (POP.people[n].I_D   == 1) { POP.people[n].I_D = 0;   POP.people[n].I_LM = 1; }
					}

					/////////////////////////////////////////////////////////////////////
					// ACTION: administer primaquine

					if ((POP.people[n].PQ_treat == 1) && (POP.people[n].PQ_effective == 1))
					{
						POP.people[n].Hyp = 0;                                  // Hypnozoites cleared

						POP.people[n].AQ8_proph = 1;                            // Put under prophylaxis 
						POP.people[n].AQ8_proph_timer = theta.MDA1_PQ_proph;    // Timer for prophylaxis set

						// Developing liver hepatic stages killed

						for (int z = 0; z < POP.people[n].lam_bite_track.size(); z++)
						{
							POP.people[n].lam_bite_track[z] = 0.0;
						}
						for (int z = 0; z < POP.people[n].lam_rel_track.size(); z++)
						{
							POP.people[n].lam_rel_track[z] = 0.0;
						}
					}
				}
			}
		}
	}


	//////////////////////////////////////////////////////////
	// Intervention 7: MDA (CQ and TQ)

	for (int m = 0; m < MDA2_year.size(); m++)
	{
		if ((t > MDA2_year[m] - 0.5*t_step) &&
			(t < MDA2_year[m] + 0.51*t_step))
		{
			cout << "MDA (CQ+TQ) distribution" << endl;

			////////////////////////////////////////////
			// Update parameters

			theta.MDA2_cover           = MDA2_cover[m];
			theta.MDA2_CQ_eff          = MDA2_CQ_eff[m];
			theta.MDA2_CQ_eff_wPQ      = MDA2_CQ_eff_wPQ[m];
			theta.MDA2_CQ_proph        = MDA2_CQ_proph[m];
			theta.MDA2_PQ_eff          = MDA2_PQ_eff[m];
			theta.MDA2_PQ_proph        = MDA2_PQ_proph[m];
			theta.MDA2_PQ_adhere       = MDA2_PQ_adhere[m];
			theta.MDA2_PQ_G6PD_risk    = MDA2_PQ_G6PD_risk[m];
			theta.MDA2_PQ_CYP2D6_risk  = MDA2_PQ_CYP2D6_risk[m];
			theta.MDA2_PQ_preg_risk    = MDA2_PQ_preg_risk[m];
			theta.MDA2_PQ_lowage       = MDA2_PQ_lowage[m];
			theta.MDA2_TQ_eff          = MDA2_TQ_eff[m];
			theta.MDA2_TQ_proph        = MDA2_TQ_proph[m];
			theta.MDA2_TQ_adhere	   = MDA2_TQ_adhere[m];
			theta.MDA2_TQ_G6PD_risk    = MDA2_TQ_G6PD_risk[m];
			theta.MDA2_TQ_CYP2D6_risk  = MDA2_TQ_CYP2D6_risk[m];
			theta.MDA2_TQ_preg_risk    = MDA2_TQ_preg_risk[m];
			theta.MDA2_TQ_lowage       = MDA2_TQ_lowage[m];
			theta.MDA2_G6PD_test	   = MDA2_G6PD_test[m];


			////////////////////////////////////////////
			// Implement intervention

			try
			{
				QQ = phi_inv(theta.MDA2_cover, 0.0, sqrt(1.0 + theta.sig_round_MDA*theta.sig_round_MDA));
			}
			catch (const char* e)
			{
				std::cerr << e << std::endl;
				exit(1);
			}

			for (int n = 0; n < POP.N_pop; n++)
			{
				POP.people[n].CQ_treat  = 0;
				POP.people[n].PQ_treat  = 0;
				POP.people[n].TQ_treat  = 0;
				POP.people[n].G6PD_test = 0;

				if (gennor(POP.people[n].zz_int[4], theta.sig_round_MDA) < QQ)
				{
					/////////////////////////////////////////////////////////////////////
					// Primaquine administered to people 6 months to 16 years

					if ((POP.people[n].age > theta.MDA2_PQ_lowage) && (POP.people[n].age < theta.MDA2_TQ_lowage))
					{
						POP.people[n].PQ_treat  = 1;
						POP.people[n].TQ_treat  = 0;
						POP.people[n].G6PD_test = 0;

						/////////////////////////////////////////////////////////////////////
						// Exclude PQ because of pregancy

						if ((theta.MDA2_PQ_preg_risk == 1) && (POP.people[n].pregnant == 1))
						{
							POP.people[n].PQ_treat = 0;
						}

						//////////////////////////////////////////////////////////////////////////////////////////
						// Is G6PD testing administered ? If so, count test

						if ((POP.people[n].PQ_treat == 1) && (theta.MDA2_G6PD_test == 1))
						{
							POP.people[n].G6PD_test = 1;
						}

						/////////////////////////////////////////////////////////////////////
						// Exclude PQ because of G6PD deficiency

						if ((POP.people[n].G6PD_test == 1) && (theta.MDA2_PQ_G6PD_risk == 1))
						{
							POP.people[n].G6PD_read = G6PD_SD_BioSensor(POP.people[n].G6PD_activity);

							if (POP.people[n].G6PD_read < 3.0)
							{
								POP.people[n].PQ_treat = 0;
							}
						}

						/////////////////////////////////////////////////////////////////////
						// Is PQ effective?

						POP.people[n].PQ_effective = 1;

						if (genunf(0.0, 1.0) > theta.MDA2_PQ_eff)
						{
							POP.people[n].PQ_effective = 0;
						}

						/////////////////////////////////////////////////////////////////////
						// Is PQ adhered to?

						if (genunf(0.0, 1.0) < 1 - theta.MDA2_PQ_adhere)
						{
							POP.people[n].PQ_effective = 0;
						}

						/////////////////////////////////////////////////////////////////////
						// Is PQ metabolised?

						if ((theta.MDA2_PQ_CYP2D6_risk == 1) && (POP.people[n].CYP2D6 == 1))
						{
							POP.people[n].PQ_effective = 0;
						}


						/////////////////////////////////////////////////////////////////////
						// Was there PQ overtreatment?

						if ((POP.people[n].PQ_treat == 1) && (POP.people[n].Hyp == 0))
						{
							POP.people[n].PQ_overtreat = 1;
						}

						if ((POP.people[n].PQ_treat == 1) && (POP.people[n].T_last_BS > 270.0))
						{
							POP.people[n].PQ_overtreat_9m = 1;
						}
					}


					/////////////////////////////////////////////////////////////////////
					// Tafenoquine administered to people older than 16 years

					if (POP.people[n].age >= theta.MDA2_TQ_lowage)
					{
						POP.people[n].TQ_treat  = 1;
						POP.people[n].PQ_treat  = 0;
						POP.people[n].G6PD_test = 0;

						/////////////////////////////////////////////////////////////////////
						// Exclude TQ because of pregancy

						if ((theta.MDA2_TQ_preg_risk == 1) && (POP.people[n].pregnant == 1))
						{
							POP.people[n].TQ_treat = 0;
						}

						//////////////////////////////////////////////////////////////////////////////////////////
						// Is G6PD testing administered to those >16 years and not pregnant? If so, count test

						if ((POP.people[n].TQ_treat == 1) && (theta.MDA2_G6PD_test == 1))
						{
							POP.people[n].G6PD_test = 1;
						}

						/////////////////////////////////////////////////////////////////////
						// Exclude TQ because of G6PD deficiency

						if ((POP.people[n].G6PD_test == 1) && (theta.MDA2_TQ_G6PD_risk == 1))
						{
							POP.people[n].G6PD_read = G6PD_SD_BioSensor(POP.people[n].G6PD_activity);

							if (POP.people[n].G6PD_read < 7.0)
							{
								POP.people[n].TQ_treat = 0;

								if (POP.people[n].G6PD_read > 3.0)
								{
									POP.people[n].PQ_treat = 1;
								}
							}
						}

						/////////////////////////////////////////////////////////////////////
						// Is TQ effective?

						POP.people[n].TQ_effective = 0;

						if (POP.people[n].TQ_treat == 1)
						{
							if (genunf(0.0, 1.0) < theta.MDA2_TQ_eff)
							{
								POP.people[n].TQ_effective = 1;
							}

							/////////////////////////////////////////////////////////////////////
							// Is TQ adhered to?

							if (genunf(0.0, 1.0) < 1 - theta.MDA2_TQ_adhere)
							{
								POP.people[n].TQ_effective = 0;
							}

							/////////////////////////////////////////////////////////////////////
							// Is TQ metabolised?

							if ((theta.MDA2_TQ_CYP2D6_risk == 1) && (POP.people[n].CYP2D6 == 1))
							{
								POP.people[n].TQ_effective = 0;
							}
						}

						/////////////////////////////////////////////////////////////////////
						// Is PQ effective?

						POP.people[n].PQ_effective = 0;

						if (POP.people[n].PQ_treat == 1)
						{
							if (genunf(0.0, 1.0) < theta.MDA2_PQ_eff)
							{
								POP.people[n].PQ_effective = 1;
							}

							/////////////////////////////////////////////////////////////////////
							// Is PQ adhered to?

							if (genunf(0.0, 1.0) < (1 - theta.MDA2_PQ_adhere))
							{
								POP.people[n].PQ_effective = 0;
							}

							/////////////////////////////////////////////////////////////////////
							// Is PQ metabolised?

							if ((theta.MDA2_PQ_CYP2D6_risk == 1) && (POP.people[n].CYP2D6 == 1))
							{
								POP.people[n].PQ_effective = 0;
							}
						}


						/////////////////////////////////////////////////////////////////////
						// Was there PQ overtreatment?

						if ((POP.people[n].PQ_treat == 1) && (POP.people[n].Hyp == 0))
						{
							POP.people[n].PQ_overtreat = 1;
						}

						if ((POP.people[n].PQ_treat == 1) && (POP.people[n].T_last_BS > 270.0))
						{
							POP.people[n].PQ_overtreat_9m = 1;
						}


						/////////////////////////////////////////////////////////////////////
						// Was there TQ overtreatment?

						if ((POP.people[n].TQ_treat == 1) && (POP.people[n].Hyp == 0))
						{
							POP.people[n].TQ_overtreat = 1;
						}

						if ((POP.people[n].TQ_treat == 1) && (POP.people[n].T_last_BS > 270.0))
						{
							POP.people[n].TQ_overtreat_9m = 1;
						}

					}

					/////////////////////////////////////////////////////
					// Chloroquine is always administered
					// Is blood-stage treatment effective

					POP.people[n].CQ_treat = 1;
					POP.people[n].CQ_effective = 0;
					
					if( (POP.people[n].PQ_treat == 0) && (POP.people[n].TQ_treat == 0) )
					{
						if (genunf(0.0, 1.0) < theta.MDA2_CQ_eff)
						{
							POP.people[n].CQ_effective = 1;
						}
					}

					if( POP.people[n].PQ_treat == 1 )
					{
						if (genunf(0.0, 1.0) < theta.MDA2_CQ_eff_wPQ)
						{
							POP.people[n].CQ_effective = 1;
						}
					}
					
					if (POP.people[n].TQ_treat == 1)
					{
						POP.people[n].CQ_effective = 1;
					}

					/////////////////////////////////////////////////////////////////////
					// ACTION: administer chloroquine

					POP.people[n].CQ_treat = 1;

					if (POP.people[n].CQ_effective == 1)
					{
						if (POP.people[n].S     == 1) { POP.people[n].S = 0;     POP.people[n].P = 1; }
						if (POP.people[n].I_PCR == 1) { POP.people[n].I_PCR = 0; POP.people[n].P = 1; }
						if (POP.people[n].I_LM  == 1) { POP.people[n].I_LM = 0;  POP.people[n].P = 1; }
						if (POP.people[n].I_D   == 1) { POP.people[n].I_D = 0;   POP.people[n].T = 1; }
					}
					else
					{
						if (POP.people[n].S     == 1) { POP.people[n].S = 0;     POP.people[n].P = 1; }
						if (POP.people[n].I_PCR == 1) { POP.people[n].I_PCR = 0; POP.people[n].P = 1; }
						if (POP.people[n].I_LM  == 1) { POP.people[n].I_LM = 0;  POP.people[n].P = 1; }
						if (POP.people[n].I_D   == 1) { POP.people[n].I_D = 0;   POP.people[n].I_LM = 1; }
					}

					/////////////////////////////////////////////////////////////////////
					// ACTION: administer primaquine

					if ((POP.people[n].PQ_treat == 1) && (POP.people[n].PQ_effective == 1))
					{
						POP.people[n].Hyp = 0;                                  // Hypnozoites cleared

						POP.people[n].AQ8_proph = 1;                            // Put under prophylaxis 
						POP.people[n].AQ8_proph_timer = theta.MDA2_PQ_proph;    // Timer for prophylaxis set

						// Developing liver hepatic stages killed

						for (int z = 0; z < POP.people[n].lam_bite_track.size(); z++)
						{
							POP.people[n].lam_bite_track[z] = 0.0;
						}
						for (int z = 0; z < POP.people[n].lam_rel_track.size(); z++)
						{
							POP.people[n].lam_rel_track[z] = 0.0;
						}
					}


					/////////////////////////////////////////////////////////////////////
					// ACTION: administer tafenoquine

					if ((POP.people[n].TQ_treat == 1) && (POP.people[n].TQ_effective == 1))
					{
						POP.people[n].Hyp = 0;                                  // Hypnozoites cleared

						POP.people[n].AQ8_proph = 1;                            // Put under prophylaxis 
						POP.people[n].AQ8_proph_timer = theta.MDA2_TQ_proph;    // Timer for prophylaxis set

						// Developing liver hepatic stages killed

						for (int z = 0; z < POP.people[n].lam_bite_track.size(); z++)
						{
							POP.people[n].lam_bite_track[z] = 0.0;
						}
						for (int z = 0; z < POP.people[n].lam_rel_track.size(); z++)
						{
							POP.people[n].lam_rel_track[z] = 0.0;
						}
					}
				}
			}
		}
	}


	//////////////////////////////////////////////////////////
	// Intervention 8: MSAT (chloroquine)

	for (int m = 0; m < MSAT0_year.size(); m++)
	{
		if ((t > MSAT0_year[m] - 0.5*t_step) &&
			(t < MSAT0_year[m] + 0.51*t_step))
		{
			cout << "MSAT (CQ) distribution" << endl;

			////////////////////////////////////////////
			// Update parameters

			theta.MSAT0_cover    = MSAT0_cover[m];
			theta.MSAT0_RDT_PCR  = MSAT0_RDT_PCR[m];
			theta.MSAT0_sens     = MSAT0_sens[m];
			theta.MSAT0_CQ_eff   = MSAT0_CQ_eff[m];
			theta.MSAT0_CQ_proph = MSAT0_CQ_proph[m];


			////////////////////////////////////////////
			// Implement intervention

			try
			{
				QQ = phi_inv(theta.MSAT0_cover, 0.0, sqrt(1.0 + theta.sig_round_MDA*theta.sig_round_MDA));
			}
			catch (const char* e)
			{
				std::cerr << e << std::endl;
				exit(1);
			}

			for (int n = 0; n < POP.N_pop; n++)
			{
				POP.people[n].CQ_treat  = 0;
				POP.people[n].PQ_treat  = 0;
				POP.people[n].TQ_treat  = 0;
				POP.people[n].G6PD_test = 0;

				if (gennor(POP.people[n].zz_int[5], theta.sig_round_MDA) < QQ)
				{
					/////////////////////////////////////////////////////
					// Blood-stage treatment is only administered to MSAT positive individuals

					MSAT_pos = 0;

					////////////////////////////////////////////////
					// Diagnosis by RDT, assumed the same as LM 	

					if (theta.MSAT0_RDT_PCR == 1)
					{
						if ((POP.people[n].I_LM == 1) || (POP.people[n].I_D == 1) || (POP.people[n].T == 1))
						{
							if (genunf(0.0, 1.0) < theta.MSAT0_sens)
							{
								MSAT_pos = 1;
							}
						}
					}


					////////////////////////////////////////////////
					// Diagnosis by PCR 	

					if (theta.MSAT0_RDT_PCR == 2)
					{
						if ((POP.people[n].I_PCR == 1) || (POP.people[n].I_LM == 1) || (POP.people[n].I_D == 1) || (POP.people[n].T == 1))
						{
							if (genunf(0.0, 1.0) < theta.MSAT0_sens)
							{
								MSAT_pos = 1;
							}
						}
					}


					/////////////////////////////////////////////////////////////////////
					// ACTION: administer chloroquine

					if (MSAT_pos == 1)
					{
						/////////////////////////////////////////////////////
						// Is blood-stage treatment effective

						POP.people[n].CQ_effective = 0;

						if (genunf(0.0, 1.0) < theta.MSAT0_CQ_eff)
						{
							POP.people[n].CQ_effective = 1;
						}

						POP.people[n].CQ_treat = 1;

						if (POP.people[n].CQ_effective == 1)
						{
							if (POP.people[n].S     == 1) { POP.people[n].S = 0;     POP.people[n].P = 1; }
							if (POP.people[n].I_PCR == 1) { POP.people[n].I_PCR = 0; POP.people[n].P = 1; }
							if (POP.people[n].I_LM  == 1) { POP.people[n].I_LM = 0;  POP.people[n].P = 1; }
							if (POP.people[n].I_D   == 1) { POP.people[n].I_D = 0;   POP.people[n].T = 1; }
						}
						else
						{
							if (POP.people[n].S     == 1) { POP.people[n].S = 0;     POP.people[n].P = 1; }
							if (POP.people[n].I_PCR == 1) { POP.people[n].I_PCR = 0; POP.people[n].P = 1; }
							if (POP.people[n].I_LM  == 1) { POP.people[n].I_LM = 0;  POP.people[n].P = 1; }
							if (POP.people[n].I_D   == 1) { POP.people[n].I_D = 0;   POP.people[n].I_LM = 1; }
						}
					}

				}
			}
		}
	}


	//////////////////////////////////////////////////////////
	// Intervention 9: MSAT (chloroquine and PQ)

	for (int m = 0; m < MSAT1_year.size(); m++)
	{
		if ((t > MSAT1_year[m] - 0.5*t_step) &&
			(t < MSAT1_year[m] + 0.51*t_step))
		{
			cout << "MSAT (CQ+PQ) distribution" << endl;

			////////////////////////////////////////////
			// Update parameters

			theta.MSAT1_cover          = MSAT1_cover[m];
			theta.MSAT1_RDT_PCR        = MSAT1_RDT_PCR[m];
			theta.MSAT1_sens           = MSAT1_sens[m];
			theta.MSAT1_CQ_eff         = MSAT1_CQ_eff[m];
			theta.MSAT1_CQ_eff_wPQ     = MSAT1_CQ_eff_wPQ[m];
			theta.MSAT1_CQ_proph       = MSAT1_CQ_proph[m];
			theta.MSAT1_PQ_eff         = MSAT1_PQ_eff[m];
			theta.MSAT1_PQ_proph       = MSAT1_PQ_proph[m];
			theta.MSAT1_PQ_adhere	   = MSAT1_PQ_adhere[m];
			theta.MSAT1_PQ_lowage      = MSAT1_PQ_lowage[m];
			theta.MSAT1_PQ_G6PD_risk   = MSAT1_PQ_G6PD_risk[m];
			theta.MSAT1_PQ_CYP2D6_risk = MSAT1_PQ_CYP2D6_risk[m];
			theta.MSAT1_PQ_preg_risk   = MSAT1_PQ_preg_risk[m];
			theta.MSAT1_G6PD_test      = MSAT1_G6PD_test[m];


			////////////////////////////////////////////
			// Implement intervention

			try
			{
				QQ = phi_inv(theta.MSAT1_cover, 0.0, sqrt(1.0 + theta.sig_round_MDA*theta.sig_round_MDA));
			}
			catch (const char* e)
			{
				std::cerr << e << std::endl;
				exit(1);
			}

			for (int n = 0; n < POP.N_pop; n++)
			{
				POP.people[n].CQ_treat  = 0;
				POP.people[n].PQ_treat  = 0;
				POP.people[n].TQ_treat  = 0;
				POP.people[n].G6PD_test = 0;

				if (gennor(POP.people[n].zz_int[6], theta.sig_round_MDA) < QQ)
				{
					/////////////////////////////////////////////////////
					// Treatment administered only to positive individuals

					MSAT_pos = 0;

					////////////////////////////////////////////////
					// Diagnosis by RDT, assumed the same as LM 	

					if (theta.MSAT1_RDT_PCR == 1)
					{
						if ((POP.people[n].I_LM == 1) || (POP.people[n].I_D == 1) || (POP.people[n].T == 1))
						{
							if (genunf(0.0, 1.0) < theta.MSAT1_sens)
							{
								MSAT_pos = 1;
							}
						}
					}


					////////////////////////////////////////////////
					// Diagnosis by PCR 	

					if (theta.MSAT1_RDT_PCR == 2)
					{
						if ((POP.people[n].I_PCR == 1) || (POP.people[n].I_LM == 1) || (POP.people[n].I_D == 1) || (POP.people[n].T == 1))
						{
							if (genunf(0.0, 1.0) < theta.MSAT1_sens)
							{
								MSAT_pos = 1;
							}
						}
					}


					if (MSAT_pos == 1)
					{
						/////////////////////////////////////////////////////////////////////
						// Is PQ administered?

						POP.people[n].PQ_treat = 1;
						POP.people[n].G6PD_test = 0;
						
						/////////////////////////////////////////////////////////////////////
						// Exclude PQ because of young age

						if (POP.people[n].age < theta.MSAT1_PQ_lowage)
						{
							POP.people[n].PQ_treat = 0;
						}

						/////////////////////////////////////////////////////////////////////
						// Exclude PQ because of pregancy

						if ((theta.MSAT1_PQ_preg_risk == 1) && (POP.people[n].pregnant == 1))
						{
							POP.people[n].PQ_treat = 0;
						}

						//////////////////////////////////////////////////////////////////////////////////////////
						// Is G6PD testing administered ? If so, count test

						if ((POP.people[n].PQ_treat == 1) && (theta.MSAT1_G6PD_test == 1))
						{
							POP.people[n].G6PD_test = 1;
						}

						/////////////////////////////////////////////////////////////////////
						// Exclude PQ because of G6PD deficiency if there is testing

						if ((POP.people[n].G6PD_test == 1) && (theta.MSAT1_PQ_G6PD_risk == 1))
						{
							POP.people[n].G6PD_read = G6PD_SD_BioSensor(POP.people[n].G6PD_activity);

							if (POP.people[n].G6PD_read < 3.0)
							{
								POP.people[n].PQ_treat = 0;
							}
						}

						/////////////////////////////////////////////////////////////////////
						// Is PQ effective?

						POP.people[n].PQ_effective = 1;

						if (genunf(0.0, 1.0) > theta.MSAT1_PQ_eff)
						{
							POP.people[n].PQ_effective = 0;
						}

						/////////////////////////////////////////////////////////////////////
						// Is PQ adhered to?

						if (genunf(0.0, 1.0) < 1 - theta.MSAT1_PQ_adhere)
						{
							POP.people[n].PQ_effective = 0;
						}

						/////////////////////////////////////////////////////////////////////
						// Is PQ metabolised?

						if ((theta.MSAT1_PQ_CYP2D6_risk == 1) && (POP.people[n].CYP2D6 == 1))
						{
							POP.people[n].PQ_effective = 0;
						}


						/////////////////////////////////////////////////////////////////////
						// Was there PQ overtreatment?

						if ((POP.people[n].PQ_treat == 1) && (POP.people[n].Hyp == 0))
						{
							POP.people[n].PQ_overtreat = 1;
						}

						if ((POP.people[n].PQ_treat == 1) && (POP.people[n].T_last_BS > 270.0))
						{
							POP.people[n].PQ_overtreat_9m = 1;
						}


						/////////////////////////////////////////////////////
						// Is chloroquine effective

						POP.people[n].CQ_effective = 0;

						if( POP.people[n].PQ_treat == 0 )
						{
							if (genunf(0.0, 1.0) < theta.MSAT1_CQ_eff)
							{
								POP.people[n].CQ_effective = 1;
							}
						}
						if (POP.people[n].PQ_treat == 1)
						{
							if (genunf(0.0, 1.0) < theta.MSAT1_CQ_eff_wPQ)
							{
								POP.people[n].CQ_effective = 1;
							}
						}


						/////////////////////////////////////////////////////////////////////
						// ACTION: administer blood-stage drug

						POP.people[n].CQ_treat = 1;

						if (POP.people[n].CQ_effective == 1)
						{
							if (POP.people[n].S     == 1) { POP.people[n].S = 0;     POP.people[n].P = 1; }
							if (POP.people[n].I_PCR == 1) { POP.people[n].I_PCR = 0; POP.people[n].P = 1; }
							if (POP.people[n].I_LM  == 1) { POP.people[n].I_LM = 0;  POP.people[n].P = 1; }
							if (POP.people[n].I_D   == 1) { POP.people[n].I_D = 0;   POP.people[n].T = 1; }
						}
						else
						{
							if (POP.people[n].S     == 1) { POP.people[n].S = 0;     POP.people[n].P = 1; }
							if (POP.people[n].I_PCR == 1) { POP.people[n].I_PCR = 0; POP.people[n].P = 1; }
							if (POP.people[n].I_LM  == 1) { POP.people[n].I_LM = 0;  POP.people[n].P = 1; }
							if (POP.people[n].I_D   == 1) { POP.people[n].I_D = 0;   POP.people[n].I_LM = 1; }
						}


						/////////////////////////////////////////////////////////////////////
						// ACTION: administer primaquine

						if ((POP.people[n].PQ_treat == 1) && (POP.people[n].PQ_effective == 1))
						{
							POP.people[n].Hyp = 0;                                  // Hypnozoites cleared

							POP.people[n].AQ8_proph = 1;                            // Put under prophylaxis 
							POP.people[n].AQ8_proph_timer = theta.MSAT1_PQ_proph;    // Timer for prophylaxis set

							// Developing liver hepatic stages killed

							for (int z = 0; z < POP.people[n].lam_bite_track.size(); z++)
							{
								POP.people[n].lam_bite_track[z] = 0.0;
							}
							for (int z = 0; z < POP.people[n].lam_rel_track.size(); z++)
							{
								POP.people[n].lam_rel_track[z] = 0.0;
							}

						}

					}
				}
			}
		}
	}


	//////////////////////////////////////////////////////////
	// Intervention 10: MSAT (chloroquine and TQ)

	for (int m = 0; m < MSAT2_year.size(); m++)
	{
		if ((t > MSAT2_year[m] - 0.5*t_step) &&
			(t < MSAT2_year[m] + 0.51*t_step))
		{
			cout << "MSAT (CQ+TQ) distribution" << endl;

			////////////////////////////////////////////
			// Update parameters

			theta.MSAT2_cover          = MSAT2_cover[m];
			theta.MSAT2_RDT_PCR        = MSAT2_RDT_PCR[m];
			theta.MSAT2_sens           = MSAT2_sens[m];
			theta.MSAT2_CQ_eff         = MSAT2_CQ_eff[m];
			theta.MSAT2_CQ_eff_wPQ     = MSAT2_CQ_eff_wPQ[m];
			theta.MSAT2_CQ_proph       = MSAT2_CQ_proph[m];
			theta.MSAT2_PQ_eff         = MSAT2_PQ_eff[m];
			theta.MSAT2_PQ_proph       = MSAT2_PQ_proph[m];
			theta.MSAT2_PQ_adhere      = MSAT2_PQ_adhere[m];
			theta.MSAT2_PQ_lowage      = MSAT2_PQ_lowage[m];
			theta.MSAT2_PQ_G6PD_risk   = MSAT2_PQ_G6PD_risk[m];
			theta.MSAT2_PQ_CYP2D6_risk = MSAT2_PQ_CYP2D6_risk[m];
			theta.MSAT2_PQ_preg_risk   = MSAT2_PQ_preg_risk[m];
			theta.MSAT2_TQ_eff         = MSAT2_TQ_eff[m];
			theta.MSAT2_TQ_proph       = MSAT2_TQ_proph[m];
			theta.MSAT2_TQ_adhere      = MSAT2_TQ_adhere[m];
			theta.MSAT2_TQ_lowage      = MSAT2_TQ_lowage[m];
			theta.MSAT2_TQ_G6PD_risk   = MSAT2_TQ_G6PD_risk[m];
			theta.MSAT2_TQ_CYP2D6_risk = MSAT2_TQ_CYP2D6_risk[m];
			theta.MSAT2_TQ_preg_risk   = MSAT2_TQ_preg_risk[m];
			theta.MSAT2_G6PD_test      = MSAT2_G6PD_test[m];


			////////////////////////////////////////////
			// Implement intervention

			try
			{
				QQ = phi_inv(theta.MSAT2_cover, 0.0, sqrt(1.0 + theta.sig_round_MDA*theta.sig_round_MDA));
			}
			catch (const char* e)
			{
				std::cerr << e << std::endl;
				exit(1);
			}

			for (int n = 0; n < POP.N_pop; n++)
			{
				POP.people[n].CQ_treat = 0;
				POP.people[n].PQ_treat = 0;
				POP.people[n].TQ_treat = 0;
				POP.people[n].G6PD_test = 0;

				if (gennor(POP.people[n].zz_int[7], theta.sig_round_MDA) < QQ)
				{
					/////////////////////////////////////////////////////
					// Blood-stage treatment is always administered

					MSAT_pos = 0;


					////////////////////////////////////////////////
					// Diagnosis by RDT, assumed the same as LM 	

					if (theta.MSAT2_RDT_PCR == 1)
					{
						if ((POP.people[n].I_LM == 1) || (POP.people[n].I_D == 1) || (POP.people[n].T == 1))
						{
							if (genunf(0.0, 1.0) < theta.MSAT2_sens)
							{
								MSAT_pos = 1;
							}
						}
					}


					////////////////////////////////////////////////
					// Diagnosis by PCR 	

					if (theta.MSAT2_RDT_PCR == 2)
					{
						if ((POP.people[n].I_PCR == 1) || (POP.people[n].I_LM == 1) || (POP.people[n].I_D == 1) || (POP.people[n].T == 1))
						{
							if (genunf(0.0, 1.0) < theta.MSAT2_sens)
							{
								MSAT_pos = 1;
							}
						}
					}


					if (MSAT_pos == 1)
					{
						/////////////////////////////////////////////////////////////////////
						// Primaquine administered to people 6 months to 16 years

						if ((POP.people[n].age > theta.MSAT2_PQ_lowage) && (POP.people[n].age < theta.MSAT2_TQ_lowage))
						{
							POP.people[n].PQ_treat = 1;

							/////////////////////////////////////////////////////////////////////
							// Exclude PQ because of pregancy

							if ((theta.MSAT2_PQ_preg_risk == 1) && (POP.people[n].pregnant == 1))
							{
								POP.people[n].PQ_treat = 0;
							}

							//////////////////////////////////////////////////////////////////////////////////////////
							// Is G6PD testing administered ? If so, count test

							if ((POP.people[n].PQ_treat == 1) && (theta.MSAT2_G6PD_test == 1))
							{
								POP.people[n].G6PD_test = 1;
							}

							/////////////////////////////////////////////////////////////////////
							// Exclude PQ because of G6PD deficiency if tested

							if ((POP.people[n].G6PD_test == 1) && (theta.MSAT2_PQ_G6PD_risk == 1))
							{
								POP.people[n].G6PD_read = G6PD_SD_BioSensor(POP.people[n].G6PD_activity);

								if (POP.people[n].G6PD_read < 3.0)
								{
									POP.people[n].PQ_treat = 0;
								}
							}

							/////////////////////////////////////////////////////////////////////
							// Is PQ effective?

							POP.people[n].PQ_effective = 0;

							if (genunf(0.0, 1.0) < theta.MSAT2_PQ_eff)
							{
								POP.people[n].PQ_effective = 1;
							}

							/////////////////////////////////////////////////////////////////////
							// Is PQ adhered to?

							if (genunf(0.0, 1.0) < (1 - theta.MSAT2_PQ_adhere))
							{
								POP.people[n].PQ_effective = 0;
							}

							/////////////////////////////////////////////////////////////////////
							// Is PQ metabolised?

							if ((theta.MSAT2_PQ_CYP2D6_risk == 1) && (POP.people[n].CYP2D6 == 1))
							{
								POP.people[n].PQ_effective = 0;
							}

							/////////////////////////////////////////////////////////////////////////////////
							// If PQ is administered and effective, set the number of HPZ batches to zero

							if ((POP.people[n].PQ_treat == 1) && (POP.people[n].PQ_effective == 1))
							{
								POP.people[n].Hyp = 0;

								POP.people[n].AQ8_proph = 1;
								POP.people[n].AQ8_proph_timer = theta.MSAT2_PQ_proph;
							}

							/////////////////////////////////////////////////////////////////////
							// Was there PQ overtreatment?

							if ((POP.people[n].PQ_treat == 1) && (POP.people[n].Hyp == 0))
							{
								POP.people[n].PQ_overtreat = 1;
							}

							if ((POP.people[n].PQ_treat == 1) && (POP.people[n].T_last_BS > 270.0))
							{
								POP.people[n].PQ_overtreat_9m = 1;
							}

						}


						/////////////////////////////////////////////////////////////////////
						// Tafenoquine administered to people older than 16 years

						if (POP.people[n].age >= theta.MSAT2_TQ_lowage)
						{
							POP.people[n].TQ_treat  = 1;
							POP.people[n].PQ_treat  = 0;
							POP.people[n].G6PD_test = 0;

							/////////////////////////////////////////////////////////////////////
							// Exclude TQ because of pregancy

							if ((theta.MSAT2_TQ_preg_risk == 1) && (POP.people[n].pregnant == 1))
							{
								POP.people[n].TQ_treat = 0;
							}

							//////////////////////////////////////////////////////////////////////////////////////////
							// Is G6PD testing administered to those >16 years and not pregnant? If so, count test

							if ((POP.people[n].TQ_treat == 1) && (theta.MSAT2_G6PD_test == 1))
							{
								POP.people[n].G6PD_test = 1;
							}

							/////////////////////////////////////////////////////////////////////
							// Exclude TQ because of G6PD deficiency

							if ((POP.people[n].G6PD_test == 1) && (theta.MSAT2_TQ_G6PD_risk == 1))
							{
								POP.people[n].G6PD_read = G6PD_SD_BioSensor(POP.people[n].G6PD_activity);

								if (POP.people[n].G6PD_read < 7.0)
								{
									POP.people[n].TQ_treat = 0;

									if (POP.people[n].G6PD_read > 3.0)
									{
										POP.people[n].PQ_treat = 1;
									}
								}
							}

							/////////////////////////////////////////////////////////////////////
							// Is TQ (or PQ) effective?

							POP.people[n].TQ_effective = 0;

							if (POP.people[n].TQ_treat == 1)
							{
								if (genunf(0.0, 1.0) < theta.MSAT2_TQ_eff)
								{
									POP.people[n].TQ_effective = 1;
								}

								/////////////////////////////////////////////////////////////////////
								// Is TQ adhered to?

								if (genunf(0.0, 1.0) < (1 - theta.MSAT2_TQ_adhere))
								{
									POP.people[n].TQ_effective = 0;
								}

								/////////////////////////////////////////////////////////////////////
								// Is TQ metabolised?

								if ((theta.MSAT2_TQ_CYP2D6_risk == 1) && (POP.people[n].CYP2D6 == 1))
								{
									POP.people[n].TQ_effective = 0;
								}
							}

							if (POP.people[n].PQ_treat == 1)
							{
								if (genunf(0.0, 1.0) < theta.MSAT2_PQ_eff)
								{
									POP.people[n].PQ_effective = 1;
								}

								/////////////////////////////////////////////////////////////////////
								// Is PQ adhered to?

								if (genunf(0.0, 1.0) < (1 - theta.MSAT2_PQ_adhere))
								{
									POP.people[n].PQ_effective = 0;
								}

								/////////////////////////////////////////////////////////////////////
								// Is PQ metabolised?

								if ((theta.MSAT2_PQ_CYP2D6_risk == 1) && (POP.people[n].CYP2D6 == 1))
								{
									POP.people[n].PQ_effective = 0;
								}
							}

							/////////////////////////////////////////////////////////////////////
							// Was there PQ overtreatment?

							if ((POP.people[n].PQ_treat == 1) && (POP.people[n].Hyp == 0))
							{
								POP.people[n].PQ_overtreat = 1;
							}

							if ((POP.people[n].PQ_treat == 1) && (POP.people[n].T_last_BS > 270.0))
							{
								POP.people[n].PQ_overtreat_9m = 1;
							}

							
							/////////////////////////////////////////////////////////////////////
							// Was there TQ overtreatment?

							if ((POP.people[n].TQ_treat == 1) && (POP.people[n].Hyp == 0))
							{
								POP.people[n].TQ_overtreat = 1;
							}

							if ((POP.people[n].TQ_treat == 1) && (POP.people[n].T_last_BS > 270.0))
							{
								POP.people[n].TQ_overtreat_9m = 1;
							}
						}

						/////////////////////////////////////////////////////
						// Is chloroquine effective

						POP.people[n].CQ_effective = 0;

						if ((POP.people[n].PQ_treat == 0) & (POP.people[n].TQ_treat == 0))
						{
							if (genunf(0.0, 1.0) < theta.MSAT2_CQ_eff)
							{
								POP.people[n].CQ_effective = 1;
							}
						}

						if( POP.people[n].PQ_treat == 1 )
						{
							if (genunf(0.0, 1.0) < theta.MSAT2_CQ_eff_wPQ)
							{
								POP.people[n].CQ_effective = 1;
							}
						}

						if( POP.people[n].TQ_treat == 1 )
						{
							POP.people[n].CQ_effective = 1;
						}

						////////////////////////////////////////////////////////////
						// ACTION: administer chloroquine

						POP.people[n].CQ_treat = 1;

						if (POP.people[n].CQ_effective == 1)
						{
							if (POP.people[n].S     == 1) { POP.people[n].S = 0;     POP.people[n].P = 1; }
							if (POP.people[n].I_PCR == 1) { POP.people[n].I_PCR = 0; POP.people[n].P = 1; }
							if (POP.people[n].I_LM  == 1) { POP.people[n].I_LM = 0;  POP.people[n].P = 1; }
							if (POP.people[n].I_D   == 1) { POP.people[n].I_D = 0;   POP.people[n].T = 1; }
						}
						else
						{
							if (POP.people[n].S     == 1) { POP.people[n].S = 0;     POP.people[n].P = 1; }
							if (POP.people[n].I_PCR == 1) { POP.people[n].I_PCR = 0; POP.people[n].P = 1; }
							if (POP.people[n].I_LM  == 1) { POP.people[n].I_LM = 0;  POP.people[n].P = 1; }
							if (POP.people[n].I_D   == 1) { POP.people[n].I_D = 0;   POP.people[n].I_LM = 1; }
						}


						/////////////////////////////////////////////////////////////////////
						// ACTION: administer primaquine

						if ((POP.people[n].PQ_treat == 1) && (POP.people[n].PQ_effective == 1))
						{
							POP.people[n].Hyp = 0;                                  // Hypnozoites cleared

							POP.people[n].AQ8_proph = 1;                            // Put under prophylaxis 
							POP.people[n].AQ8_proph_timer = theta.MSAT2_PQ_proph;    // Timer for prophylaxis set

							// Developing liver hepatic stages killed

							for (int z = 0; z < POP.people[n].lam_bite_track.size(); z++)
							{
								POP.people[n].lam_bite_track[z] = 0.0;
							}
							for (int z = 0; z < POP.people[n].lam_rel_track.size(); z++)
							{
								POP.people[n].lam_rel_track[z] = 0.0;
							}
						}


						/////////////////////////////////////////////////////////////////////
						// ACTION: administer tafenoquine

						if ((POP.people[n].TQ_treat == 1) && (POP.people[n].TQ_effective == 1))
						{
							POP.people[n].Hyp = 0;                                  // Hypnozoites cleared

							POP.people[n].AQ8_proph = 1;                            // Put under prophylaxis 
							POP.people[n].AQ8_proph_timer = theta.MSAT2_TQ_proph;    // Timer for prophylaxis set

							// Developing liver hepatic stages killed

							for (int z = 0; z < POP.people[n].lam_bite_track.size(); z++)
							{
								POP.people[n].lam_bite_track[z] = 0.0;
							}
							for (int z = 0; z < POP.people[n].lam_rel_track.size(); z++)
							{
								POP.people[n].lam_rel_track[z] = 0.0;
							}

						}
						
					}

				}
			}
		}
	}


	//////////////////////////////////////////////////////////
	// Intervention 11: STAT (CQ and PQ)

	for (int m = 0; m < STAT1_year.size(); m++)
	{
		if ((t > STAT1_year[m] - 0.5*t_step) &&
			(t < STAT1_year[m] + 0.51*t_step))
		{
			cout << "STAT (CQ+PQ) distribution" << endl;

			////////////////////////////////////////////
			// Update parameters

			theta.STAT1_cover          = STAT1_cover[m];
			theta.STAT1_sens           = STAT1_sens[m];
			theta.STAT1_spec           = STAT1_spec[m];
			theta.STAT1_RDT_PCR        = STAT1_RDT_PCR[m];
			theta.STAT1_CQ_eff         = STAT1_CQ_eff[m];
			theta.STAT1_CQ_eff_wPQ     = STAT1_CQ_eff_wPQ[m];
			theta.STAT1_CQ_proph       = STAT1_CQ_proph[m];
			theta.STAT1_PQ_eff         = STAT1_PQ_eff[m];
			theta.STAT1_PQ_proph       = STAT1_PQ_proph[m];
			theta.STAT1_PQ_adhere      = STAT1_PQ_adhere[m];
			theta.STAT1_PQ_lowage      = STAT1_PQ_lowage[m];
			theta.STAT1_PQ_G6PD_risk   = STAT1_PQ_G6PD_risk[m];
			theta.STAT1_PQ_CYP2D6_risk = STAT1_PQ_CYP2D6_risk[m];
			theta.STAT1_PQ_preg_risk   = STAT1_PQ_preg_risk[m];
			theta.STAT1_G6PD_test      = STAT1_G6PD_test[m];


			////////////////////////////////////////////
			// Implement intervention

			try
			{
				QQ = phi_inv(theta.STAT1_cover, 0.0, sqrt(1.0 + theta.sig_round_MDA*theta.sig_round_MDA));
			}
			catch (const char* e)
			{
				std::cerr << e << std::endl;
				exit(1);
			}

			for (int n = 0; n < POP.N_pop; n++)
			{
				POP.people[n].CQ_treat  = 0;
				POP.people[n].PQ_treat  = 0;
				POP.people[n].TQ_treat  = 0;
				POP.people[n].G6PD_test = 0;

				if (gennor(POP.people[n].zz_int[8], theta.sig_round_MDA) < QQ)
				{
					/////////////////////////////////////////////////////////////////////
					// SSAT screening for blood-stage infection in the last 9 months
					//
					// There are two options here. Option 1 define over-treatment on the 
					// basis of blood-stage infection with the last 9 month. Option 2
					// defines over-treatment on the basis of presence of hypnozoites.

					STAT_pos = 0;

					// OPTION 1
					
					if( (POP.people[n].T_last_BS <= 270.0) && (genunf(0.0, 1.0) < theta.STAT1_sens) )
					{
						STAT_pos = 1;
					}

					if( (POP.people[n].T_last_BS > 270.0) && (genunf(0.0, 1.0) > theta.STAT1_spec) )
					{
						STAT_pos = 1;
					}
					

					// OPTION 2
					/*
					if ((POP.people[n].Hyp > 0) && (genunf(0.0, 1.0) < theta.STAT1_sens))
					{
						STAT_pos = 1;
					}

					if ((POP.people[n].Hyp == 0) && (genunf(0.0, 1.0) > theta.STAT1_spec))
					{
						STAT_pos = 1;
					}
					*/


					/////////////////////////////////
					// Additional detection of BS infections

					// RDT or LM
					if (theta.STAT1_RDT_PCR == 1)
					{
						if (POP.people[n].I_LM == 1)
						{
							STAT_pos = 1;
						}
					}

					// molecular diagnosis (PCR)
					if (theta.STAT1_RDT_PCR == 2)
					{
						if( (POP.people[n].I_PCR == 1) || (POP.people[n].I_LM == 1) )
						{
							STAT_pos = 1;
						}
					}

					// clinically ill or recently treated people always test positive

					if ((POP.people[n].I_D == 1) || (POP.people[n].T == 1) || (POP.people[n].P == 1))
					{
						STAT_pos = 1;
					}


					/////////////////////////////////
					// Proceed if diagnosed positive

					if (STAT_pos == 1)
					{
						POP.people[n].PQ_treat = 1;


						/////////////////////////////////////////////////////////////////////
						// Exclude PQ because of young age

						if (POP.people[n].age < theta.STAT1_PQ_lowage)
						{
							POP.people[n].PQ_treat = 0;
						}

						/////////////////////////////////////////////////////////////////////
						// Exclude PQ because of pregancy

						if ((theta.STAT1_PQ_preg_risk == 1) && (POP.people[n].pregnant == 1))
						{
							POP.people[n].PQ_treat = 0;
						}

						//////////////////////////////////////////////////////////////////////////////////////////
						// Is G6PD testing administered to those >6 months and not pregnant? If so, count test

						if ((POP.people[n].PQ_treat == 1) && (theta.STAT1_G6PD_test == 1))
						{
							POP.people[n].G6PD_test = 1;
						}

						/////////////////////////////////////////////////////////////////////
						// Exclude PQ because of G6PD deficiency

						if ((POP.people[n].G6PD_test == 1) && (theta.STAT1_PQ_G6PD_risk == 1))
						{
							POP.people[n].G6PD_read = G6PD_SD_BioSensor(POP.people[n].G6PD_activity);

							if (POP.people[n].G6PD_read < 3.0)
							{
								POP.people[n].PQ_treat = 0;
							}
						}

						/////////////////////////////////////////////////////////////////////
						// Is PQ effective?

						POP.people[n].PQ_effective = 0;

						if (genunf(0.0, 1.0) < theta.STAT1_PQ_eff)
						{
							POP.people[n].PQ_effective = 1;
						}

						/////////////////////////////////////////////////////////////////////
						// Is PQ adhered to?

						if (genunf(0.0, 1.0) < (1 - theta.STAT1_PQ_adhere))
						{
							POP.people[n].PQ_effective = 0;
						}

						/////////////////////////////////////////////////////////////////////
						// Is PQ metabolised?

						if ((theta.STAT1_PQ_CYP2D6_risk == 1) && (POP.people[n].CYP2D6 == 1))
						{
							POP.people[n].PQ_effective = 0;
						}


						/////////////////////////////////////////////////////////////////////
						// Was there PQ overtreatment?

						if ((POP.people[n].PQ_treat == 1) && (POP.people[n].Hyp == 0))
						{
							POP.people[n].PQ_overtreat = 1;
						}

						if ((POP.people[n].PQ_treat == 1) && (POP.people[n].T_last_BS > 270.0))
						{
							POP.people[n].PQ_overtreat_9m = 1;
						}


						/////////////////////////////////////////////////////
						// Is chloroquine effective

						POP.people[n].CQ_effective = 0;

						if (POP.people[n].PQ_treat == 0)
						{
							if (genunf(0.0, 1.0) < theta.STAT1_CQ_eff)
							{
								POP.people[n].CQ_effective = 1;
							}
						}

						if (POP.people[n].PQ_treat == 1)
						{
							if (genunf(0.0, 1.0) < theta.STAT1_CQ_eff_wPQ)
							{
								POP.people[n].CQ_effective = 1;
							}
						}


						/////////////////////////////////////////////////////////////////////
						// ACTION: administer chloroquine

						POP.people[n].CQ_treat = 1;

						if (POP.people[n].CQ_effective == 1)
						{
							if (POP.people[n].S     == 1) { POP.people[n].S     = 0; POP.people[n].P = 1; }
							if (POP.people[n].I_PCR == 1) { POP.people[n].I_PCR = 0; POP.people[n].P = 1; }
							if (POP.people[n].I_LM  == 1) { POP.people[n].I_LM  = 0; POP.people[n].P = 1; }
							if (POP.people[n].I_D   == 1) { POP.people[n].I_D   = 0; POP.people[n].T = 1; }
						}
						else
						{
							if (POP.people[n].S     == 1) { POP.people[n].S     = 0; POP.people[n].P    = 1; }
							if (POP.people[n].I_PCR == 1) { POP.people[n].I_PCR = 0; POP.people[n].P    = 1; }
							if (POP.people[n].I_LM  == 1) { POP.people[n].I_LM  = 0; POP.people[n].P    = 1; }
							if (POP.people[n].I_D   == 1) { POP.people[n].I_D   = 0; POP.people[n].I_LM = 1; }
						}


						/////////////////////////////////////////////////////////////////////
						// ACTION: administer primaquine

						if ((POP.people[n].PQ_treat == 1) && (POP.people[n].PQ_effective == 1))
						{
							POP.people[n].Hyp = 0;                                  // Hypnozoites cleared

							POP.people[n].AQ8_proph = 1;                            // Put under prophylaxis 
							POP.people[n].AQ8_proph_timer = theta.STAT1_PQ_proph;    // Timer for prophylaxis set

							// Developing liver hepatic stages killed

							for (int z = 0; z < POP.people[n].lam_bite_track.size(); z++)
							{
								POP.people[n].lam_bite_track[z] = 0.0;
							}
							for (int z = 0; z < POP.people[n].lam_rel_track.size(); z++)
							{
								POP.people[n].lam_rel_track[z] = 0.0;
							}

						}
					}
				}
			}
		}
	}


	//////////////////////////////////////////////////////////
	// Intervention 12: STAT (chloroquine and TQ)

	for (int m = 0; m < STAT2_year.size(); m++)
	{
		if ((t > STAT2_year[m] - 0.5*t_step) &&
			(t < STAT2_year[m] + 0.51*t_step))
		{
			cout << "STAT (CQ+TQ) distribution" << endl;

			////////////////////////////////////////////
			// Update parameters

			theta.STAT2_cover          = STAT2_cover[m];
			theta.STAT2_sens           = STAT2_sens[m];
			theta.STAT2_spec           = STAT2_spec[m];
			theta.STAT2_RDT_PCR        = STAT2_RDT_PCR[m];
			theta.STAT2_CQ_eff         = STAT2_CQ_eff[m];
			theta.STAT2_CQ_eff_wPQ     = STAT2_CQ_eff_wPQ[m];
			theta.STAT2_CQ_proph       = STAT2_CQ_proph[m];
			theta.STAT2_PQ_eff         = STAT2_PQ_eff[m];
			theta.STAT2_PQ_proph       = STAT2_PQ_proph[m];
			theta.STAT2_PQ_adhere      = STAT2_PQ_adhere[m];
			theta.STAT2_PQ_lowage      = STAT2_PQ_lowage[m];
			theta.STAT2_PQ_G6PD_risk   = STAT2_PQ_G6PD_risk[m];
			theta.STAT2_PQ_CYP2D6_risk = STAT2_PQ_CYP2D6_risk[m];
			theta.STAT2_PQ_preg_risk   = STAT2_PQ_preg_risk[m];
			theta.STAT2_TQ_eff         = STAT2_TQ_eff[m];
			theta.STAT2_TQ_proph       = STAT2_TQ_proph[m];
			theta.STAT2_TQ_adhere      = STAT2_TQ_adhere[m];
			theta.STAT2_TQ_lowage      = STAT2_TQ_lowage[m];
			theta.STAT2_TQ_G6PD_risk   = STAT2_TQ_G6PD_risk[m];
			theta.STAT2_TQ_CYP2D6_risk = STAT2_TQ_CYP2D6_risk[m];
			theta.STAT2_TQ_preg_risk   = STAT2_TQ_preg_risk[m];
			theta.STAT2_G6PD_test      = STAT2_G6PD_test[m];


			////////////////////////////////////////////
			// Implement intervention

			try
			{
				QQ = phi_inv(theta.STAT2_cover, 0.0, sqrt(1.0 + theta.sig_round_MDA*theta.sig_round_MDA));
			}
			catch (const char* e)
			{
				std::cerr << e << std::endl;
				exit(1);
			}

			for (int n = 0; n < POP.N_pop; n++)
			{
				POP.people[n].CQ_treat  = 0;
				POP.people[n].PQ_treat  = 0;
				POP.people[n].TQ_treat  = 0;
				POP.people[n].G6PD_test = 0;

				if (gennor(POP.people[n].zz_int[9], theta.sig_round_MDA) < QQ)
				{
					/////////////////////////////////////////////////////////////////////
					// SSAT screening for blood-stage infection in the last 9 months
					//
					// There are two options here. Option 1 define over-treatment on the 
					// basis of blood-stage infection with the last 9 month. Option 2
					// defines over-treatment on the basis of presence of hypnozoites.

					STAT_pos = 0;

					// OPTION 1
					
					if( (POP.people[n].T_last_BS <= 270.0) && (genunf(0.0, 1.0) < theta.STAT2_sens) )
					{
						STAT_pos = 1;
					}

					if( (POP.people[n].T_last_BS > 270.0) && (genunf(0.0, 1.0) > theta.STAT2_spec) )
					{
						STAT_pos = 1;
					}
					

					// OPTION 2
					/*
					if ((POP.people[n].Hyp > 0) && (genunf(0.0, 1.0) < theta.STAT2_sens))
					{
						STAT_pos = 1;
					}

					if ((POP.people[n].Hyp == 0) && (genunf(0.0, 1.0) > theta.STAT2_spec))
					{
						STAT_pos = 1;
					}
					*/


					/////////////////////////////////
					// Additional detection of BS infections

					// RDT or LM
					if (theta.STAT2_RDT_PCR == 1)
					{
						if (POP.people[n].I_LM == 1)
						{
							STAT_pos = 1;
						}
					}

					// molecular diagnosis (PCR)
					if (theta.STAT2_RDT_PCR == 2)
					{
						if ((POP.people[n].I_PCR == 1) || (POP.people[n].I_LM == 1))
						{
							STAT_pos = 1;
						}
					}

					// clinically ill or recently treated people always test positive

					if ((POP.people[n].I_D == 1) || (POP.people[n].T == 1) || (POP.people[n].P == 1))
					{
						STAT_pos = 1;
					}


					if (STAT_pos == 1)
					{
						POP.people[n].PQ_treat  = 0;
						POP.people[n].TQ_treat  = 0;
						POP.people[n].G6PD_test = 0;


						/////////////////////////////////////////////////////////////////////
						// Primaquine administered to people 6 months to 16 years

						if ((POP.people[n].age > theta.STAT2_PQ_lowage) && (POP.people[n].age < theta.STAT2_TQ_lowage))
						{
							POP.people[n].PQ_treat = 1;

							/////////////////////////////////////////////////////////////////////
							// Exclude PQ because of pregancy

							if ((theta.STAT2_PQ_preg_risk == 1) && (POP.people[n].pregnant == 1))
							{
								POP.people[n].PQ_treat = 0;
							}

							//////////////////////////////////////////////////////////////////////////////////////////
							// Is G6PD testing administered ? If so, count test

							if ((POP.people[n].PQ_treat == 1) && (theta.STAT2_G6PD_test == 1))
							{
								POP.people[n].G6PD_test = 1;
							}

							/////////////////////////////////////////////////////////////////////
							// Exclude PQ because of G6PD deficiency

							if ((POP.people[n].G6PD_test == 1) && (theta.STAT2_PQ_G6PD_risk == 1))
							{
								POP.people[n].G6PD_read = G6PD_SD_BioSensor(POP.people[n].G6PD_activity);

								if (POP.people[n].G6PD_read < 3.0)
								{
									POP.people[n].PQ_treat = 0;
								}
							}

							/////////////////////////////////////////////////////////////////////
							// Is PQ effective?

							POP.people[n].PQ_effective = 0;

							if (genunf(0.0, 1.0) < theta.STAT2_PQ_eff)
							{
								POP.people[n].PQ_effective = 1;
							}

							/////////////////////////////////////////////////////////////////////
							// Is PQ adhered to?

							if (genunf(0.0, 1.0) < (1 - theta.STAT2_PQ_adhere))
							{
								POP.people[n].PQ_effective = 0;
							}

							/////////////////////////////////////////////////////////////////////
							// Is PQ metabolised?

							if ((theta.STAT2_PQ_CYP2D6_risk == 1) && (POP.people[n].CYP2D6 == 1))
							{
								POP.people[n].PQ_effective = 0;
							}

							/////////////////////////////////////////////////////////////////////
							// Was there PQ overtreatment?

							if ((POP.people[n].PQ_treat == 1) && (POP.people[n].Hyp == 0))
							{
								POP.people[n].PQ_overtreat = 1;
							}

							if ((POP.people[n].PQ_treat == 1) && (POP.people[n].T_last_BS > 270.0))
							{
								POP.people[n].PQ_overtreat_9m = 1;
							}
						}


						/////////////////////////////////////////////////////////////////////
						// Tafenoquine administered to people older than 16 years

						if (POP.people[n].age >= theta.STAT2_TQ_lowage)
						{
							POP.people[n].TQ_treat  = 1;
							POP.people[n].PQ_treat  = 0;
							POP.people[n].G6PD_test = 0;

							/////////////////////////////////////////////////////////////////////
							// Exclude TQ because of pregancy

							if ((theta.STAT2_TQ_preg_risk == 1) && (POP.people[n].pregnant == 1))
							{
								POP.people[n].TQ_treat = 0;
							}

							//////////////////////////////////////////////////////////////////////////////////////////
							// Is G6PD testing administered ? If so, count test

							if ((POP.people[n].TQ_treat == 1) && (theta.STAT2_G6PD_test == 1))
							{
								POP.people[n].G6PD_test = 1;
							}

							/////////////////////////////////////////////////////////////////////
							// Exclude TQ because of G6PD deficiency

							if ((POP.people[n].G6PD_test == 1) && (theta.STAT2_TQ_G6PD_risk == 1))
							{
								POP.people[n].G6PD_read = G6PD_SD_BioSensor(POP.people[n].G6PD_activity);

								if (POP.people[n].G6PD_read < 7.0)
								{
									POP.people[n].TQ_treat = 0;

									if (POP.people[n].G6PD_read > 3.0)
									{
										POP.people[n].PQ_treat = 1;
									}
								}
							}


							/////////////////////////////////////////////////////////////////////
							// Is TQ and PQ effective?

							POP.people[n].TQ_effective = 0;

							if (POP.people[n].TQ_treat == 1)
							{
								if (genunf(0.0, 1.0) < theta.STAT2_TQ_eff)
								{
									POP.people[n].TQ_effective = 1;
								}

								/////////////////////////////////////////////////////////////////////
								// Is TQ adhered to?

								if (genunf(0.0, 1.0) < (1 - theta.STAT2_TQ_adhere))
								{
									POP.people[n].PQ_effective = 0;
								}

								/////////////////////////////////////////////////////////////////////
								// Is TQ metabolised?

								if ((theta.STAT2_TQ_CYP2D6_risk == 1) && (POP.people[n].CYP2D6 == 1))
								{
									POP.people[n].TQ_effective = 0;
								}
							}

							if (POP.people[n].PQ_treat == 1)
							{
								if (genunf(0.0, 1.0) < theta.STAT2_PQ_eff)
								{
									POP.people[n].PQ_effective = 1;
								}

								/////////////////////////////////////////////////////////////////////
								// Is PQ adhered to?

								if (genunf(0.0, 1.0) < (1 - theta.STAT2_PQ_adhere))
								{
									POP.people[n].PQ_effective = 0;
								}

								/////////////////////////////////////////////////////////////////////
								// Is PQ metabolised?

								if ((theta.STAT2_PQ_CYP2D6_risk == 1) && (POP.people[n].CYP2D6 == 1))
								{
									POP.people[n].PQ_effective = 0;
								}
							}

							/////////////////////////////////////////////////////////////////////
							// Was there PQ overtreatment?

							if ((POP.people[n].PQ_treat == 1) && (POP.people[n].Hyp == 0))
							{
								POP.people[n].PQ_overtreat = 1;
							}

							if ((POP.people[n].PQ_treat == 1) && (POP.people[n].T_last_BS > 270.0))
							{
								POP.people[n].PQ_overtreat_9m = 1;
							}

							/////////////////////////////////////////////////////////////////////
							// Was there TQ overtreatment?

							if ((POP.people[n].TQ_treat == 1) && (POP.people[n].Hyp == 0))
							{
								POP.people[n].TQ_overtreat = 1;
							}

							if ((POP.people[n].TQ_treat == 1) && (POP.people[n].T_last_BS > 270.0))
							{
								POP.people[n].TQ_overtreat_9m = 1;
							}

						}


						/////////////////////////////////////////////////////
						// Is chloroquine effective

						POP.people[n].CQ_effective = 0;

						if ((POP.people[n].PQ_treat == 0) && (POP.people[n].TQ_treat == 0))
						{
							if (genunf(0.0, 1.0) < theta.STAT2_CQ_eff)
							{
								POP.people[n].CQ_effective = 1;
							}
						}

						if( POP.people[n].PQ_treat == 1 )
						{
							if (genunf(0.0, 1.0) < theta.STAT2_CQ_eff_wPQ)
							{
								POP.people[n].CQ_effective = 1;
							}
						}
						
						if (POP.people[n].TQ_treat == 1)
						{
							POP.people[n].CQ_effective = 1;
						}

						/////////////////////////////////////////////////////////////////////
						// ACTION: administer chloroquine

						POP.people[n].CQ_treat = 1;

						if (POP.people[n].CQ_effective == 1)
						{
							if (POP.people[n].S     == 1) { POP.people[n].S     = 0; POP.people[n].P = 1; }
							if (POP.people[n].I_PCR == 1) { POP.people[n].I_PCR = 0; POP.people[n].P = 1; }
							if (POP.people[n].I_LM  == 1) { POP.people[n].I_LM  = 0; POP.people[n].P = 1; }
							if (POP.people[n].I_D   == 1) { POP.people[n].I_D   = 0; POP.people[n].T = 1; }
						}
						else
						{
							if (POP.people[n].S     == 1) { POP.people[n].S     = 0; POP.people[n].P    = 1; }
							if (POP.people[n].I_PCR == 1) { POP.people[n].I_PCR = 0; POP.people[n].P    = 1; }
							if (POP.people[n].I_LM  == 1) { POP.people[n].I_LM  = 0; POP.people[n].P    = 1; }
							if (POP.people[n].I_D   == 1) { POP.people[n].I_D   = 0; POP.people[n].I_LM = 1; }
						}

						/////////////////////////////////////////////////////////////////////
						// ACTION: administer primaquine

						if ((POP.people[n].PQ_treat == 1) && (POP.people[n].PQ_effective == 1))
						{
							POP.people[n].Hyp = 0;                                  // Hypnozoites cleared

							POP.people[n].AQ8_proph = 1;                            // Put under prophylaxis 
							POP.people[n].AQ8_proph_timer = theta.STAT2_PQ_proph;    // Timer for prophylaxis set

							// Developing liver hepatic stages killed

							for (int z = 0; z < POP.people[n].lam_bite_track.size(); z++)
							{
								POP.people[n].lam_bite_track[z] = 0.0;
							}
							for (int z = 0; z < POP.people[n].lam_rel_track.size(); z++)
							{
								POP.people[n].lam_rel_track[z] = 0.0;
							}

						}

						/////////////////////////////////////////////////////////////////////
						// ACTION: administer tafenoquine

						if ((POP.people[n].TQ_treat == 1) && (POP.people[n].TQ_effective == 1))
						{
							POP.people[n].Hyp = 0;                                  // Hypnozoites cleared

							POP.people[n].AQ8_proph = 1;                            // Put under prophylaxis 
							POP.people[n].AQ8_proph_timer = theta.STAT2_TQ_proph;    // Timer for prophylaxis set

							// Developing liver hepatic stages killed

							for (int z = 0; z < POP.people[n].lam_bite_track.size(); z++)
							{
								POP.people[n].lam_bite_track[z] = 0.0;
							}
							for (int z = 0; z < POP.people[n].lam_rel_track.size(); z++)
							{
								POP.people[n].lam_rel_track[z] = 0.0;
							}
						}

					}

				}
			}
		}
	}


	//////////////////////////////////////////////////////////
	// Intervention 13: IVM

	for (int m = 0; m < IVM_year.size(); m++)
	{
		if ((t > (IVM_year[m] - 0.5*t_step)) &&
			(t < (IVM_year[m] + 0.51*t_step)))
		{

			try
			{
				QQ = phi_inv(IVM_cover[m], 0.0, sqrt(1.0 + theta.sig_round_MDA*theta.sig_round_MDA));
			}
			catch (const char* e)
			{
				std::cerr << e << std::endl;
				exit(1);
			}
			theta.P_IVM_decay = exp( -t_step*CONST_LOG_2/IVM_half_life[m] );

			cout << IVM_cover[m] << "\t" << theta.P_IVM_decay << endl;

			for (int n = 0; n < POP.N_pop; n++)
			{
				if (gennor(POP.people[n].zz_int[10], theta.sig_round_MDA) < QQ)
				{
					POP.people[n].IVM = 1;
					POP.people[n].IVM_age = 0.0;

					POP.people[n].d_IVM = d_IVM_0[m];
					POP.people[n].s_IVM = 1.0 - POP.people[n].d_IVM;
				}
			}
		}
	}
}


///////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////// 
//         //                                                                            // 
//  3.2.4. //  Inverse of the cumulative normal distribution function required           //
//         //  for implementing correlated intervention coverage.                        //
/////////////                                                                            //
/////////////  The code is based on the following website                                //
/////////////  http://www.johndcook.com/blog/cpp_phi_inverse/                            //
/////////////  which is based the algorithm in Abramowitz and Stegun formula 26.2.23.    // 
/////////////                                                                            //
/////////////  The absolute value of the error should be less than 4.5 e-4 and it tests  //
/////////////  out nicely in R.                                                          //
/////////////                                                                            //
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

double phi_inv(double pp, double mu, double sigma)
{
	if (pp < 0.0 || pp > 1.0)
	{
		throw("bad vlaue of pp (coverage) in phi_inv");
	}
	else if (pp == 0.0)
	{
		return -std::numeric_limits<double>::infinity();
	}
	else if (pp == 1.0)
	{
		return std::numeric_limits<double>::infinity();
	}

	double cc[] = { 2.515517, 0.802853, 0.010328 };
	double dd[] = { 1.432788, 0.189269, 0.001308 };

	double tt, temp;

	if (pp < 0.5)
	{
		tt = sqrt(-2.0*log(pp));

		temp = -(tt - ((cc[2] * tt + cc[1])*tt + cc[0]) / (((dd[2] * tt + dd[1])*tt + dd[0])*tt + 1.0));

	}
	else {

		tt = sqrt(-2.0*log(1.0 - pp));

		temp = (tt - ((cc[2] * tt + cc[1])*tt + cc[0]) / (((dd[2] * tt + dd[1])*tt + dd[0])*tt + 1.0));
	}

	return mu + sigma * temp;
}
