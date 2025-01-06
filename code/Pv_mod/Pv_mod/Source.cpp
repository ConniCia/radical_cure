/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
///                                                                       ///
///  Individual-based Plasmodium vivax transmission model.                ///
///                                                                       ///
///  Dr Michael White                                                     ///
///  Institut Pasteur                                                     ///
///  michael.white@pasteur.fr                                             ///
///                                                                       ///
///  With contributions from Dr Thomas Obadia, Dr Narimane Nekkab         ///
///  and Diggory Hardy                                                    ///
///                                                                       ///
///  Please feel free to use and modify if you wish. However,             ///
///  please provide appropriate acknowledgement and get in touch          ///
///  if you have any questions. This is not necessarily the               ///
///  final, canonical version of the code - contact me to get that.       ///
///                                                                       ///
///  There is a prize of a pint for reporting any serious bugs or         ///
///  finding something new that results in >20% speed up.                 ///
///  First prize to this goes to Diggory who achieved >80% speeed up,     ///
///  which will earn him 4 pints or his chosen equivalent.                ///
///                                                                       ///
///  Model code is split up into multiple files as follows:               ///
///                                                                       ///
///  1.1.  Source.cpp                                                     ///
///        This file                                                      ///
///                                                                       ///
///  2.1.  Params.h                                                       ///
///  2.2.  Params.cpp                                                     ///
///        The Params structure stores input parameters and has           ///
///        associated code for reading parameters from input files.       ///
///                                                                       ///
///  3.1.  Intervention.h                                                 ///
///  3.2.  Intervention.cpp                                               ///
///        The Intervention struct stores intervention parameters.        ///
///                                                                       ///
///  4.1.  Individual.h                                                   ///
///  4.2   Individual.cpp                                                 ///
///        A class is created which stores all the information of a       ///
///        single individual.                                             ///
///        Details of the stochastic individual-based model for each      ///
///        person. Transitions occur with a fixed time step according to  ///
///        compting hazards                                               ///
///                                                                       ///
///  5.1.  Population.h                                                   ///
///  5.2.  Population.cpp                                                 ///
///        A structure called Population stores all individuals.          ///
///        This set of functions calculates the equilibrium set up of the ///
///        population. It is only called once while the population is     ///
///        being initialised.                                             ///
///                                                                       ///
///  6.1.  Mosquito.cpp                                                   ///
///        Mosquitoes are simulated using a deterministic ODE solver.     ///
///                                                                       ///
///  7.1.  Simulation.h                                                   ///
///  7.2.  Simulation.cpp                                                 ///
///        Creation of class for storing simulation output, running       ///
///        and writing output                                             ///
///                                                                       ///
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

#include "Simulation.h"

#include <iostream>
#include <cmath>
#include <time.h>
#include "randlib.h"


////////////////////////////////////////////
//                                        //
//  1.1.1. Initiate main object           //
//                                        //
////////////////////////////////////////////

int main(int argc, char** argv)
{
	setall(std::stol(argv[1]), 7);

	clock_t clock_time;
	clock_time = clock();


	////////////////////////////////////////////
	//                                        //
	//  1.1.2. Read in file names             //
	//                                        //
	////////////////////////////////////////////

	// do we have the correct command line?
	if (argc != 5 + N_mosq)
	{
		std::cout << "Incorrect command line.\n";
		return 0;
	}

	const char* parameter_File = argv[2];

	const char* mosquito_File[N_mosq];
	for (int v = 0; v < N_mosq; v++)
	{
		mosquito_File[v] = argv[3 + v];
	}

	const char* coverage_File = argv[3 + N_mosq];
	const char* output_File   = argv[4 + N_mosq];


	////////////////////////////////////////////
	//                                        //
	// 1.1.3. Initialise objects              //
	//                                        //
	////////////////////////////////////////////

	Population POP;
	Params PV_MOD_PAR;


	///////////////////////////////////////////////
	//                                           //
	// 1.1.4. Read in model parameters           //
	//        and create an intervention object  //
	//                                           //
	///////////////////////////////////////////////

	SimTimes times = PV_MOD_PAR.read(parameter_File, mosquito_File);
	POP.N_pop = PV_MOD_PAR.N_pop;

	Intervention INTVEN(coverage_File);


	///////////////////////////////////////////////////////////////////////////
	//                                                                       //
	// 1.1.5. Initialise Population of individuals                           //
	//        Note that they begin with exponential age distribution         //
	//        and susceptible without immunity                               //
	//                                                                       //
	///////////////////////////////////////////////////////////////////////////

	cout << "Initialise population of individuals for simulation at equilbirium EIR of " << 365.0*PV_MOD_PAR.EIR_dom_equil << endl;
	cout << endl;


	////////////////////////////////
	// Create population objects

	POP.pop_setup(PV_MOD_PAR);

	/////////////////////////////////////
	// Calculate population equilibrium

	POP.equil_setup_count = 0;

	POP.pop_at_equil(PV_MOD_PAR);


	/////////////////////////////////////
	// Additional initialisations for PQ or TQ
	// case management at beginning

	if( (PV_MOD_PAR.CM_regimen == 1) || (PV_MOD_PAR.CM_regimen == 2) )
	{
		POP.equil_setup_count = 1;

		for (int eq = 0; eq < N_eq_setup; eq++)
		{
			POP.pop_at_equil(PV_MOD_PAR);

			POP.equil_setup_count = POP.equil_setup_count + 1;
		}
	}


	////////////////////////////////////////////
	// Initialise population of individuals at equilibrium

	POP.ind_at_equil(PV_MOD_PAR);


	cout << "Population of size " << POP.N_pop << " initialised!" << endl;
	cout << endl;


	/////////////////////////////////////////////////////////////////////////
	//                                                                     //
	// 1.1.6. Create Simulation object                                     //
	//                                                                     //
	/////////////////////////////////////////////////////////////////////////

	Simulation SIM(times);


	//////////////////////////////////////////////////////
	//                                                  //
	// 1.1.7. Begin stochastic simulations              //
	//                                                  //
	//////////////////////////////////////////////////////

	cout << "Starting model simulations......." << endl;

	SIM.run(PV_MOD_PAR, POP, INTVEN);

	cout << "Model simulations completed....." << endl;
	cout << endl;


	//////////////////////////////////////////////////////
	//                                                  //
	// 1.1.8. Output to file                            //
	//                                                  //
	//////////////////////////////////////////////////////

	SIM.write_output(output_File);


	cout << "Time taken: " << ((double)clock() - clock_time) / ((double)CLOCKS_PER_SEC) << " seconds" << endl;


	return 0;
}
