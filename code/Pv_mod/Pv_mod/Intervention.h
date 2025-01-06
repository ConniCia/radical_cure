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
///  A class for storing all intervention parameters and coverages        ///
///                                                                       ///
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

// Include guard
#ifndef PVIVAX_MODEL_INTERVENTION
#define PVIVAX_MODEL_INTERVENTION

#include "Population.h"


/////////////////////////////////////////////////////////////////////////////////////////
//                                                                                     //
// 3.1.1.  Define a structure for details of interventions                             //
//                                                                                     //
/////////////////////////////////////////////////////////////////////////////////////////

class Intervention
{
public:
	//////////////////////////////////////////////////////////////////////////
	//  Class constructors and destructors
	//////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////
	// Constructor: read intervention data from input files
	Intervention(const char *coverage_File);


	////////////////////////////////////////////////////
	// Copy and move constructors

	// Delete unwanted copy constructors
	Intervention(Intervention&) = delete;
	void operator= (Intervention&) = delete;
	
	// Allow default move constructors
	Intervention(Intervention&&) = default;
	Intervention& operator= (Intervention&&) = default;


	//////////////////////////////////////////////////////////////////////////
	//  Class member functions
	//////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////
	// Distribute interventions
	
	void distribute(double t, Params& theta, Population& POP);


private:
	//////////////////////////////////////////////////////////////////////////
	//  Data
	//////////////////////////////////////////////////////////////////////////

	////////////////////////////////
	// Intervention 0 
	// Treatment regimen = 0; chloroquine

	vector<double> CM0_year;        // start of case management regimen 0
	vector<double> CM0_cover;       // proportion of symptomatic cases treated
	vector<double> CM0_CQ_eff;      // blood-stage treatment efficacy
	vector<double> CM0_CQ_proph;    // blood-stage treatment prophylaxis


	////////////////////////////////
	// Intervention 1 
	// Treatment regimen = 1; primaquine

	vector<double> CM1_year;                 // start of case management regimen 1
	vector<double> CM1_cover;                // proportion of symptomatic cases treated
	vector<double> CM1_CQ_eff;               // chloroquine efficacy (on its own)
	vector<double> CM1_CQ_eff_wPQ;           // chloroquine efficacy (co-administered with PQ)
	vector<double> CM1_CQ_proph;             // duration of prophylaxis of chloroquine
	vector<double> CM1_PQ_eff;               // primaquine efficacy
	vector<double> CM1_PQ_proph;             // duration of primaquine prophylaxis (number of days + 1)
	vector<double> CM1_PQ_adhere;            // adherence to full primaquine regimen
	vector<double> CM1_PQ_lowage;            // youngest age for primaquine treatment regimen
	vector<int>    CM1_PQ_G6PD_risk;         // Risk in G6PD-deficient individuals (0 = no risk; 1 = risk)
	vector<int>	   CM1_PQ_CYP2D6_risk;       // Risk of not working in low CYP2D6 metabolizers  (0 = no risk; 1 = risk, i.e. doesn't work)
	vector<int>    CM1_PQ_preg_risk;         // Risk in pregnant women  (0 = no risk; 1 = risk)
	vector<int>    CM1_G6PD_test;		  	 // Is G6PD testing implemented in case management (0 = no; 1 = yes)


	////////////////////////////////
	// Intervention 2 
	// Treatment regimen = 2; tafenoquine

	vector<double> CM2_year;                 // start of blood-stage treatment regimen
	vector<double> CM2_cover;                // proportion of symptomatic cases treated
	vector<double> CM2_CQ_eff;               // chloroquine efficacy (on its own)
	vector<double> CM2_CQ_eff_wPQ;           // chloroquine efficacy (co-administered with PQ)
	vector<double> CM2_CQ_proph;             // duration of prophylaxis of chloroquine
	vector<double> CM2_PQ_eff;               // primaquine efficacy
	vector<double> CM2_PQ_proph;             // duration of primaquine prophylaxis (number of days + 1)
	vector<double> CM2_PQ_adhere;            // adherence to full primaquine regimen
	vector<double> CM2_PQ_lowage;            // youngest age for primaquine treatment regimen
	vector<int>	   CM2_PQ_G6PD_risk;         // Risk in G6PD-deficient individuals (0 = no risk; 1 = risk)
	vector<int>	   CM2_PQ_CYP2D6_risk;       // Risk of not working in low CYP2D6 metabolizers  (0 = no risk; 1 = risk, i.e. doesn't work)
	vector<int>	   CM2_PQ_preg_risk;         // Risk in pregnant women  (0 = no risk; 1 = risk)
	vector<double> CM2_TQ_eff;               // tafenoquine efficacy
	vector<double> CM2_TQ_proph;             // duration of tafenoquine prophylaxis (number of days + 1)
	vector<double> CM2_TQ_adhere;            // adherence to full tafenoquine regimen
	vector<double> CM2_TQ_lowage;            // youngest age for tafenoquine treatment regimen
	vector<int>	   CM2_TQ_G6PD_risk;         // Risk in G6PD-deficient individuals (0 = no risk; 1 = risk)
	vector<int>	   CM2_TQ_CYP2D6_risk;       // Risk of not working in low CYP2D6 metabolizers  (0 = no risk; 1 = risk, i.e. doesn't work)
	vector<int>	   CM2_TQ_preg_risk;         // Risk in pregnant women  (0 = no risk; 1 = risk)
	vector<int>    CM2_G6PD_test;		     // Is G6PD testing implemented in case management (0 = no; 1 = yes)


	////////////////////////////////
	// Intervention 3 
	// LLINs

	vector<double> LLIN_year;      
	vector<double> LLIN_cover;


	////////////////////////////////
	// Intervention 4 
	// IRS

	vector<double> IRS_year;
	vector<double> IRS_cover;


	////////////////////////////////
	// Intervention 5 parameters
	// MDA CQ parameters 

	vector<double> MDA0_year;             // Year of MDA
	vector<double> MDA0_cover;            // Coverage 
	vector<double> MDA0_CQ_eff;           // Efficacy of chloroquine
	vector<double> MDA0_CQ_proph;         // Duration of blood - stage prophylaxis


	////////////////////////////////
	// Intervention 6
	// MDA PQ parameters 

	vector<double> MDA1_year;             // Year of MDA   
	vector<double> MDA1_cover;            // Coverage
	vector<double> MDA1_CQ_eff;           // Efficacy of chloroquine (on its own)
	vector<double> MDA1_CQ_eff_wPQ;       // Efficacy of chloroquine (co-administerd with PQ)
	vector<double> MDA1_CQ_proph;         // Duration of chloroquine prophylaxis
	vector<double> MDA1_PQ_eff;           // Efficacy of primaquine
	vector<double> MDA1_PQ_proph;         // Duration of primaquine prophylaxis
	vector<double> MDA1_PQ_adhere;        // Primaquine adherence
	vector<double> MDA1_PQ_lowage;        // Lower age limit for primaquine treatment (in days)
	vector<int>    MDA1_PQ_G6PD_risk;     // Risk in G6PD - deficient individuals
	vector<int>    MDA1_PQ_CYP2D6_risk;   // Risk of not working in low CYP2D6 metabolizers
	vector<int>    MDA1_PQ_preg_risk;     // Risk in pregnant women
	vector<int>    MDA1_G6PD_test;	      // G6PD testing for PQ eligible cases  (0 = no test; 1 = test)


	////////////////////////////////
	// Intervention 7
	// MDA TQ parameters 

	vector<double> MDA2_year;             // Year of MDA   
	vector<double> MDA2_cover;            // Coverage
	vector<double> MDA2_CQ_eff;           // Efficacy of chloroquine (on its own)
	vector<double> MDA2_CQ_eff_wPQ;       // Efficacy of chloroquine (co-adminsitered with PQ)
	vector<double> MDA2_CQ_proph;         // Duration of chloroquine prophylaxis
	vector<double> MDA2_PQ_eff;           // Efficacy of primaquine
	vector<double> MDA2_PQ_proph;         // Duration of primaquine prophylaxis
	vector<double> MDA2_PQ_adhere;        // Primaquine adherence
	vector<double> MDA2_PQ_lowage;        // Lower age limit for primaquine treatment (in days)
	vector<int>    MDA2_PQ_G6PD_risk;     // Risk in G6PD - deficient individuals
	vector<int>    MDA2_PQ_CYP2D6_risk;   // Risk of not working in low CYP2D6 metabolizers
	vector<int>    MDA2_PQ_preg_risk;     // Risk in pregnant women
	vector<double> MDA2_TQ_eff;           // Efficacy of tafenoquine
	vector<double> MDA2_TQ_proph;         // Duration of tafenoquine prophylaxis
	vector<double> MDA2_TQ_adhere;        // Tafenoquine adherence
	vector<double> MDA2_TQ_lowage;        // Lower age limit for tafenoquine treatment (in days)
	vector<int>    MDA2_TQ_G6PD_risk;     // Risk in G6PD - deficient individuals
	vector<int>    MDA2_TQ_CYP2D6_risk;   // Risk of not working in low CYP2D6 metabolizers
	vector<int>    MDA2_TQ_preg_risk;     // Risk in pregnant women
	vector<int>    MDA2_G6PD_test;		  // G6PD testing for PQ eligible cases  (0 = no test; 1 = test)


	///////////////////////////////
	// Intervention 8
	// MSAT with CQ parameters 

	vector<double> MSAT0_year;            // Year of MSAT
	vector<double> MSAT0_cover;           // Coverage 
	vector<double> MSAT0_RDT_PCR;         // What diagnostic tool : 1 for RDT(= LM); 2 for PCR
	vector<double> MSAT0_sens;            // Sensitivity of diagnostic tool
	vector<double> MSAT0_CQ_eff;          // Efficacy of chloroquine
	vector<double> MSAT0_CQ_proph;        // Duration of chloroquine prophylaxis


	///////////////////////////////
	// Intervention 9
	// MSAT with PQ parameters 

	vector<double> MSAT1_year;             // Year of MSAT
	vector<double> MSAT1_cover;            // Coverage
	vector<double> MSAT1_RDT_PCR;          // What diagnostic tool : 1 for RDT(= LM); 2 for PCR
	vector<double> MSAT1_sens;             // Sensitivity of diagnostic tool
	vector<double> MSAT1_CQ_eff;           // Efficacy of chloroquine (on its own)
	vector<double> MSAT1_CQ_eff_wPQ;       // Efficacy of chloroquine (co-administered with PQ)
	vector<double> MSAT1_CQ_proph;         // Duration of chloroquine prophylaxis
	vector<double> MSAT1_PQ_eff;           // Efficacy of primaquine
	vector<double> MSAT1_PQ_proph;         // Duration of primaquine prophylaxis
	vector<double> MSAT1_PQ_adhere;        // Primaquine adherence
	vector<double> MSAT1_PQ_lowage;        // Lower age limit for primaquine treatment (in days)
	vector<int>    MSAT1_PQ_G6PD_risk;     // Risk in G6PD - deficient individuals
	vector<int>    MSAT1_PQ_CYP2D6_risk;   // Risk of not working in low CYP2D6 metabolizers
	vector<int>    MSAT1_PQ_preg_risk;     // Risk in pregnant women
	vector<int>    MSAT1_G6PD_test;	       // G6PD testing for PQ eligible cases  (0 = no test; 1 = test)


	///////////////////////////////
	// Intervention 10 parameters
	// MSAT with TQ 

	vector<double> MSAT2_year;             // Year of MSAT
	vector<double> MSAT2_cover;            // Coverage
	vector<double> MSAT2_RDT_PCR;          // What diagnostic tool : 1 for RDT(= LM); 2 for PCR
	vector<double> MSAT2_sens;             // Sensitivity of diagnostic tool
	vector<double> MSAT2_CQ_eff;           // Efficacy of chloroquine (on its own)
	vector<double> MSAT2_CQ_eff_wPQ;       // Efficacy of chloroquine (co-administerd with PQ)
	vector<double> MSAT2_CQ_proph;         // Duration of chloroquine prophylaxis
	vector<double> MSAT2_PQ_eff;           // Efficacy of primaquine
	vector<double> MSAT2_PQ_proph;         // Duration of primaquine prophylaxis
	vector<double> MSAT2_PQ_adhere;        // Primaquine adherence
	vector<double> MSAT2_PQ_lowage;        // Lower age limit for primaquine treatment (in days)
	vector<int>    MSAT2_PQ_G6PD_risk;     // Risk in G6PD - deficient individuals
	vector<int>    MSAT2_PQ_CYP2D6_risk;   // Risk of not working in low CYP2D6 metabolizers
	vector<int>    MSAT2_PQ_preg_risk;     // Risk in pregnant women
	vector<double> MSAT2_TQ_eff;           // Efficacy of tafenoquine
	vector<double> MSAT2_TQ_proph;         // Duration of tafenoquine prophylaxis
	vector<double> MSAT2_TQ_adhere;        // Tafenoquine adherence
	vector<double> MSAT2_TQ_lowage;        // Lower age limit for tafenoquine treatment (in days)
	vector<int>    MSAT2_TQ_G6PD_risk;     // Risk in G6PD - deficient individuals
	vector<int>    MSAT2_TQ_CYP2D6_risk;   // Risk of not working in low CYP2D6 metabolizers
	vector<int>    MSAT2_TQ_preg_risk;     // Risk in pregnant women
	vector<int>    MSAT2_G6PD_test;	       // G6PD testing for PQ eligible cases  (0 = no test; 1 = test)


	///////////////////////////////
	// Intervention 11 parameters
	// STAT with PQ 

	vector<double> STAT1_year;             // Year of STAT
	vector<double> STAT1_cover;            // Coverage
	vector<double> STAT1_sens;             // Sensitivity of diagnostic tool
	vector<double> STAT1_spec;             // Sensitivity of diagnostic tool
	vector<double> STAT1_RDT_PCR;          // What diagnostic tool (extra parasitological testing)): 0 = nothing; 1 for RDT(= LM); 2 for PCR
	vector<double> STAT1_CQ_eff;           // Efficacy of chloroquine (on its own)
	vector<double> STAT1_CQ_eff_wPQ;       // Efficacy of chloroquine (co-administered with PQ)
	vector<double> STAT1_CQ_proph;         // Duration of chloroquine prophylaxis
	vector<double> STAT1_PQ_eff;           // Efficacy of primaquine
	vector<double> STAT1_PQ_proph;         // Duration of primaquine prophylaxis
	vector<double> STAT1_PQ_adhere;        // Primaquine adherence
	vector<double> STAT1_PQ_lowage;        // Lower age limit for primaquine treatment (in days)
	vector<int>    STAT1_PQ_G6PD_risk;     // Risk in G6PD - deficient individuals
	vector<int>    STAT1_PQ_CYP2D6_risk;   // Risk of not working in low CYP2D6 metabolizers
	vector<int>    STAT1_PQ_preg_risk;     // Risk in pregnant women
	vector<int>    STAT1_G6PD_test;	       // G6PD testing for PQ eligible cases  (0 = no test; 1 = test)


	///////////////////////////////
	// Intervention 12 parameters
	// STAT with TQ 

	vector<double> STAT2_year;             // Year of STAT
	vector<double> STAT2_cover;            // Coverage
	vector<double> STAT2_sens;             // Sensitivity of diagnostic tool
	vector<double> STAT2_spec;             // Sensitivity of diagnostic tool
	vector<double> STAT2_RDT_PCR;          // What diagnostic tool (extra parasitological testing)): 0 = nothing; 1 for RDT(= LM); 2 for PCR
	vector<double> STAT2_CQ_eff;           // Efficacy of chloroquine (on its own)
	vector<double> STAT2_CQ_eff_wPQ;       // Efficacy of chloroquine (co-administered with PQ)
	vector<double> STAT2_CQ_proph;         // Duration of chloroquine prophylaxis
	vector<double> STAT2_PQ_eff;           // Efficacy of primaquine
	vector<double> STAT2_PQ_proph;         // Duration of primaquine prophylaxis
	vector<double> STAT2_PQ_adhere;        // Primaquine adherence
	vector<double> STAT2_PQ_lowage;        // Lower age limit for primaquine treatment (in days)
	vector<int>    STAT2_PQ_G6PD_risk;     // Risk in G6PD - deficient individuals
	vector<int>    STAT2_PQ_CYP2D6_risk;   // Risk of not working in low CYP2D6 metabolizers
	vector<int>    STAT2_PQ_preg_risk;     // Risk in pregnant women
	vector<double> STAT2_TQ_eff;           // Efficacy of tafenoquine
	vector<double> STAT2_TQ_proph;         // Duration of tafenoquine prophylaxis
	vector<double> STAT2_TQ_adhere;        // Tafenoquine adherence
	vector<double> STAT2_TQ_lowage;        // Lower age limit for tafenoquine treatment (in days)
	vector<int>    STAT2_TQ_G6PD_risk;     // Risk in G6PD - deficient individuals
	vector<int>    STAT2_TQ_CYP2D6_risk;   // Risk of not working in low CYP2D6 metabolizers
	vector<int>    STAT2_TQ_preg_risk;     // Risk in pregnant women
	vector<int>    STAT2_G6PD_test;	       // G6PD testing for PQ eligible cases  (0 = no test; 1 = test)


	/////////////////////////////////////////////////////////
	// Intervention 13 parameters
	// Ivermectin 

	vector<double> IVM_year;
	vector<double> IVM_cover;
	vector<double> d_IVM_0;          // Probability mosquito dies during feeding attempt
	vector<double> IVM_half_life;    // Half-life of loss of LLINs

};

#endif
