/*==============================================================================
Project: LiFe
Theme: Apex Method
Module: Problem-Parameters.h (Problem Parameters)
Prefix: PP
Author: Leonid B. Sokolinsky
This source code has been produced with using BSF-skeleton
==============================================================================*/
//-------------------------- Compilation Modes -----------------------
//#define PP_DEBUG
//#define PP_PATH "Problems/"
#define PP_PATH "C:/EXTRACTED/"
#define PP_WITHOUT_APEX
//=========================== Problem Parameters =========================
/**#define PP_PROBLEM_NAME		"afiro"		 //========
#define PP_M 27		// Number of equations (number of rows in *.mtx)
#define PP_N 51		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define EPS							1E-6
#define PP_EPS_DIR_LENGTH			EPS		// Minimal Length of Direction Vector
#define PP_EPS_OBJ					1E-3	// |F(u) - PP_EXACT_OBJ_VALUE| <= PP_EPS_OBJ
#define PP_EPS_SHIFT				EPS		// Minimal shift to stop motion
#define PP_EPS_ZERO_COMPARE			EPS		// Zero comparison precision
#define PP_EPS_ZERO_DIR				1E-5	// Accuracy of zero value for direction vector coordinates for direction vector coordinates
#define PP_EXACT_OBJ_VALUE		464.7531	// Exact maximum value of objective function
#define PP_INFINITY					1E+308	// High bound in *_hi.mtx
#define PP_SIGMA_TO_APEX			5000
#define PP_LOW_COST_PERCENTILE		0		// Level below which value of low cost variable falls (must be in [0,1])
#define PP_GAP						1E-4	// Maximum gap from polytope surface (not too small!!!)
#define PP_OBJECTIVE_VECTOR_LENGTH	100		// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		20		// Start length of shift vector
#define PP_MODE_BLOCK_HCV_VARIABLE	false	// block high cos variables
#define PP_MODE_USE_LCV_VARIABLE	false	// utilize low cost and zero variables
//------------------------------------------------------------------/**/

/**#define PP_PROBLEM_NAME		"adlittle"	//========
#define PP_M 56		// Number of equations (number of rows in *.mtx)
#define PP_N 138	// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define EPS							1E-2
#define PP_EPS_DIR_LENGTH			1E-3	// Minimal Length of Direction Vector
#define PP_EPS_OBJ					250		// |F(u) - PP_EXACT_OBJ_VALUE| <= PP_EPS_OBJ
#define PP_EPS_SHIFT				1E-4	// Minimal shift to stop motion
#define PP_EPS_ZERO_COMPARE			EPS		// Zero comparison precision
#define PP_EPS_ZERO_DIR				1E-5	// Accuracy of zero value for direction vector coordinates
#define PP_EXACT_OBJ_VALUE		-225494.963	// Exact maximum value of objective function
#define PP_INFINITY					1E+308	// High bound in *_hi.mtx
#define PP_SIGMA_TO_APEX			1E+4
#define PP_GAP						0.2		// Minimum distance to polytope (not too small!!!)
#define PP_LOW_COST_PERCENTILE		2
#define PP_OBJECTIVE_VECTOR_LENGTH	0.1		// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		1		// Start length of shift vector
#define PP_MODE_BLOCK_HCV_VARIABLE	false	// block high cos variables
#define PP_MODE_USE_LCV_VARIABLE	false	// utilize low cost and zero variables
//------------------------------------------------------------------/**/

/**#define PP_PROBLEM_NAME		"blend"	//==========================================
#define PP_M 74		// Number of equations (number of rows in *.mtx)
#define PP_N 114	// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define EPS							1E-3
#define PP_EPS_DIR_LENGTH			EPS		// Minimal Length of Direction Vector
#define PP_EPS_OBJ					1E-3	// |F(u) - PP_EXACT_OBJ_VALUE| <= PP_EPS_OBJ
#define PP_EPS_SHIFT				1E-5	// Minimal shift to stop motion
#define PP_EPS_ZERO_COMPARE			EPS		// Zero comparison precision
#define PP_EPS_ZERO_DIR				1E-5	// Accuracy of zero value for direction vector coordinates
#define PP_EXACT_OBJ_VALUE			30.812	// Exact maximum value of objective function
#define PP_INFINITY					1E+308	// Highest bound in *_hi.mtx
#define PP_SIGMA_TO_APEX			125
#define PP_GAP						0.2		// Minimum distance to polytope (not too small!!!)
#define PP_LOW_COST_PERCENTILE		0.02
#define PP_OBJECTIVE_VECTOR_LENGTH	1		// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		1		// Start length of shift vector
#define PP_MODE_BLOCK_HCV_VARIABLE	false	// block high cos variables
#define PP_MODE_USE_LCV_VARIABLE	false	// utilize low cost and zero variables
//------------------------------------------------------------------/**/

/**#define PP_PROBLEM_NAME		"fit1d"	//==========================================
#define PP_M 24		// Number of equations (number of rows in *.mtx)
#define PP_N 1049	// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define EPS							1E-4
#define PP_EPS_DIR_LENGTH			EPS		// Minimal Length of Direction Vector
#define PP_EPS_OBJ					1E-2	// |F(u) - PP_EXACT_OBJ_VALUE| <= PP_EPS_OBJ
#define PP_EPS_SHIFT				EPS		// Minimal shift to stop motion
#define PP_EPS_ZERO_COMPARE			EPS		// Zero comparison precision
#define PP_EPS_ZERO_DIR				EPS		// Accuracy of zero value for direction vector coordinates
#define PP_EXACT_OBJ_VALUE			9146.378// Exact maximum value of objective function
#define PP_INFINITY					1E+308	// Highest bound in *_hi.mtx
#define PP_SIGMA_TO_APEX			250
#define PP_GAP						1E-2	// Minimum distance to polytope (not too small!!!)
#define PP_LOW_COST_PERCENTILE		0
#define PP_OBJECTIVE_VECTOR_LENGTH	0.1		// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		0.1		// Start length of shift vector
#define PP_MODE_BLOCK_HCV_VARIABLE	false	// block high cos variables
#define PP_MODE_USE_LCV_VARIABLE	false	// utilize low cost and zero variables
//------------------------------------------------------------------/**/

/**#define PP_PROBLEM_NAME		"kb2"	//==========================================
#define PP_M 43	// Number of equations (number of rows in *.mtx)
#define PP_N 68	// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define PP_EPS_ZERO_COMPARE			1E-3	// Zero comparison precision
#define PP_EPS_DIR_LENGTH			1E-3	// Minimal Length of Direction Vector
#define PP_EPS_OBJ					1E-2	// |F(u) - PP_EXACT_OBJ_VALUE| <= PP_EPS_OBJ
#define PP_EPS_SHIFT				1E-4	// Minimal shift to stop motion
#define PP_EPS_ZERO_DIR				1E-5	// Accuracy of zero value for direction vector coordinates
#define PP_EXACT_OBJ_VALUE			1749.9	// <-- Exact maximum value of objective function
#define PP_INFINITY					1E+308	// Highest bound in *_hi.mtx
#define PP_SIGMA_TO_APEX			1E+7
#define PP_GAP						0.2		// Minimum distance to polytope (not too small!!!)
#define PP_LOW_COST_PERCENTILE		0.1
#define PP_OBJECTIVE_VECTOR_LENGTH	0.01	// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		1		// Start length of shift vector
#define PP_MODE_BLOCK_HCV_VARIABLE	false	// block high cos variables
#define PP_MODE_USE_LCV_VARIABLE	false	// utilize low cost and zero variables
//------------------------------------------------------------------/**/

/**#define PP_PROBLEM_NAME		"recipe"	//==========================================
#define PP_M 91		// Number of equations (number of rows in *.mtx)
#define PP_N 204	// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define PP_EPS_ZERO_COMPARE			1E-5	// Zero comparison precision
#define PP_EPS_DIR_LENGTH			1E-5	// Minimal Length of Direction Vector
#define PP_EPS_OBJ					1E-2	// |F(u) - PP_EXACT_OBJ_VALUE| <= PP_EPS_OBJ
#define PP_EPS_SHIFT				1E-5	// Minimal shift to stop motion
#define PP_EPS_ZERO_DIR				1E-5	// Accuracy of zero value for direction vector coordinates
#define PP_EXACT_OBJ_VALUE			266.61	// <-- Exact maximum value of objective function
#define PP_INFINITY					1E+308	// Highest bound in *_hi.mtx
#define PP_SIGMA_TO_APEX			1E+7
#define PP_GAP						0.2		// Minimum distance to polytope (not too small!!!)
#define PP_LOW_COST_PERCENTILE		0.1
#define PP_OBJECTIVE_VECTOR_LENGTH	0.1		// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		1		// Start length of shift vector
#define PP_MODE_BLOCK_HCV_VARIABLE	false	// block high cos variables
#define PP_MODE_USE_LCV_VARIABLE	true	// utilize low cost and zero variables
//------------------------------------------------------------------/**/

/**#define PP_PROBLEM_NAME		"sc50a"	//==========================================
#define PP_M 50	// Number of equations (number of rows in *.mtx)
#define PP_N 78	// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define PP_EPS_ZERO_COMPARE			1E-8	// Zero comparison precision
#define PP_EPS_DIR_LENGTH			1E-7	// Minimal Length of Direction Vector
#define PP_EPS_OBJ					1E-2	// |F(u) - PP_EXACT_OBJ_VALUE| <= PP_EPS_OBJ
#define PP_EPS_SHIFT				1E-8	// Minimal shift to stop motion
#define PP_EPS_ZERO_DIR				1E-8	// Accuracy of zero value for direction vector coordinates
#define PP_EXACT_OBJ_VALUE			64.575	// <-- Exact maximum value of objective function
#define PP_INFINITY					1E+308	// Highest bound in *_hi.mtx
#define PP_SIGMA_TO_APEX			5E+4
#define PP_GAP						0.2		// Minimum distance to polytope (not too small!!!)
#define PP_LOW_COST_PERCENTILE		0
#define PP_OBJECTIVE_VECTOR_LENGTH	0.01	// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		1		// Start length of shift vector
#define PP_MODE_BLOCK_HCV_VARIABLE	false	// block high cos variables
#define PP_MODE_USE_LCV_VARIABLE	false	// utilize low cost and zero variables
//------------------------------------------------------------------/**/

/**#define PP_PROBLEM_NAME		"sc50b"	//==========================================
#define PP_M 50		// Number of equations (number of rows in *.mtx)
#define PP_N 78	// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define PP_EPS_ZERO_COMPARE			1E-8	// Zero comparison precision
#define PP_EPS_DIR_LENGTH			1E-7	// Minimal Length of Direction Vector
#define PP_EPS_OBJ					1E-2	// |F(u) - PP_EXACT_OBJ_VALUE| <= PP_EPS_OBJ
#define PP_EPS_SHIFT				1E-8	// Minimal shift to stop motion
#define PP_EPS_ZERO_DIR				1E-8	// Accuracy of zero value for direction vector coordinates
#define PP_EXACT_OBJ_VALUE			70		// <-- Exact maximum value of objective function
#define PP_INFINITY					1E+308	// Highest bound in *_hi.mtx
#define PP_SIGMA_TO_APEX			1E+6
#define PP_GAP						0.2		// Minimum distance to polytope (not too small!!!)
#define PP_LOW_COST_PERCENTILE		0
#define PP_OBJECTIVE_VECTOR_LENGTH	0.01	// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		1		// Start length of shift vector
#define PP_MODE_BLOCK_HCV_VARIABLE	false	// block high cos variables
#define PP_MODE_USE_LCV_VARIABLE	false	// utilize low cost and zero variables
//------------------------------------------------------------------/**/

/**#define PP_PROBLEM_NAME		"sc105"	//==========================================
#define PP_M 105	// Number of equations (number of rows in *.mtx)
#define PP_N 163	// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define PP_EPS_ZERO_COMPARE			1E-8	// Zero comparison precision
#define PP_EPS_DIR_LENGTH			1E-7	// Minimal Length of Direction Vector
#define PP_EPS_OBJ					0.364	// |F(u) - PP_EXACT_OBJ_VALUE| <= PP_EPS_OBJ
#define PP_EPS_SHIFT				1E-8	// Minimal shift to stop motion
#define PP_EPS_ZERO_DIR				1E-8	// Accuracy of zero value for direction vector coordinates
#define PP_EXACT_OBJ_VALUE			52.202	// <------------ Exact maximum value of objective function
#define PP_INFINITY					1E+308	// Highest bound in *_hi.mtx
#define PP_SIGMA_TO_APEX			1E+5
#define PP_GAP						0.0001	// Minimum distance to polytope (not too small!!!)
#define PP_LOW_COST_PERCENTILE		0
#define PP_OBJECTIVE_VECTOR_LENGTH	0.0001	// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		1		// Start length of shift vector
#define PP_MODE_BLOCK_HCV_VARIABLE	false	// block high cos variables
#define PP_MODE_USE_LCV_VARIABLE	false	// utilize low cost and zero variables
//------------------------------------------------------------------/**/

/**#define PP_PROBLEM_NAME		"share2b" //==========================================
#define PP_M 96		// Number of equations (number of rows in *.mtx)
#define PP_N 162	// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define EPS							1E-3
#define PP_EPS_ZERO_COMPARE			EPS		// Zero comparison precision
#define PP_EPS_DIR_LENGTH			1E-5	// Minimal Length of Direction Vector
#define PP_EPS_OBJ					1E-2	// |F(u) - PP_EXACT_OBJ_VALUE| <= PP_EPS_OBJ
#define PP_EPS_SHIFT				1E-5	// Minimal shift to stop motion
#define PP_EPS_ZERO_DIR				1E-5	// Accuracy of zero value for direction vector coordinates
#define PP_EXACT_OBJ_VALUE			415.73	// <-- Exact maximum value of objective function
#define PP_INFINITY					1E+308	// Highest bound in *_hi.mtx
#define PP_SIGMA_TO_APEX			100
#define PP_GAP						0.2		// Minimum distance to polytope (not too small!!!)
#define PP_LOW_COST_PERCENTILE		0
#define PP_OBJECTIVE_VECTOR_LENGTH	0.01	// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		1		// Start length of shift vector
#define PP_MODE_BLOCK_HCV_VARIABLE	false	// block high cos variables
#define PP_MODE_USE_LCV_VARIABLE	false	// utilize low cost and zero variables
//------------------------------------------------------------------/**/

/**#define PP_PROBLEM_NAME		"simple1" //==========================================
#define PP_M 4		// Number of equations (number of rows in *.mtx)
#define PP_N 7		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define EPS							1E-8
#define PP_EPS_DIR_LENGTH			EPS		// Minimal Length of Direction Vector
#define PP_EPS_OBJ					EPS		// |F(u) - PP_EXACT_OBJ_VALUE| <= PP_EPS_OBJ
#define PP_EPS_SHIFT				EPS		// Minimal shift to stop motion
#define PP_EPS_ZERO_COMPARE			EPS		// Zero comparison precision
#define PP_EPS_ZERO_DIR				EPS		// Accuracy of zero value for direction vector coordinates for direction vector coordinates
#define PP_EXACT_OBJ_VALUE			55000	// Exact maximum value of objective function
#define PP_INFINITY					1E+308	// Highest bound in *_hi.mtx
#define PP_SIGMA_TO_APEX			5000
#define PP_LOW_COST_PERCENTILE		0		// Level below which value of low cost variable falls (must be in [0,1])
#define PP_GAP						1E-2	// Maximum gap from polytope surface (not too small!!!)
#define PP_OBJECTIVE_VECTOR_LENGTH	0.1		// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		5		// Start length of shift vector
#define PP_MODE_BLOCK_HCV_VARIABLE	true	// block high cos variables
#define PP_MODE_USE_LCV_VARIABLE	false	// utilize low cost and zero variables
//------------------------------------------------------------------/**/

/**#define PP_PROBLEM_NAME		"simple1.1" //==========================================
#define PP_M 4		// Number of equations (number of rows in *.mtx)
#define PP_N 7		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define EPS							1E-8
#define PP_EPS_DIR_LENGTH			EPS		// Minimal Length of Direction Vector
#define PP_EPS_OBJ					EPS		// |F(u) - PP_EXACT_OBJ_VALUE| <= PP_EPS_OBJ
#define PP_EPS_SHIFT				EPS		// Minimal shift to stop motion
#define PP_EPS_ZERO_COMPARE			EPS		// Zero comparison precision
#define PP_EPS_ZERO_DIR				EPS		// Accuracy of zero value for direction vector coordinates for direction vector coordinates
#define PP_EXACT_OBJ_VALUE			40000	// Exact maximum value of objective function
#define PP_INFINITY					1E+308	// Highest bound in *_hi.mtx
#define PP_SIGMA_TO_APEX			5000
#define PP_LOW_COST_PERCENTILE		0		// Level below which value of low cost variable falls (must be in [0,1])
#define PP_GAP						1E-2	// Maximum gap from polytope surface (not too small!!!)
#define PP_OBJECTIVE_VECTOR_LENGTH	0.1		// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		5		// Start length of shift vector
#define PP_MODE_BLOCK_HCV_VARIABLE	true	// block high cos variables
#define PP_MODE_USE_LCV_VARIABLE	false	// utilize low cost and zero variables
//------------------------------------------------------------------/**/

/**#define PP_PROBLEM_NAME		"simple2" //==========================================
#define PP_M 5		// Number of equations (number of rows in *.mtx)
#define PP_N 8		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define EPS							1E-8
#define PP_EPS_DIR_LENGTH			EPS		// Minimal Length of Direction Vector
#define PP_EPS_OBJ					EPS		// |F(u) - PP_EXACT_OBJ_VALUE| <= PP_EPS_OBJ
#define PP_EPS_SHIFT				EPS		// Minimal shift to stop motion
#define PP_EPS_ZERO_COMPARE			EPS		// Zero comparison precision
#define PP_EPS_ZERO_DIR				EPS		// Accuracy of zero value for direction vector coordinates for direction vector coordinates
#define PP_EXACT_OBJ_VALUE			59000	// Exact maximum value of objective function
#define PP_INFINITY					1E+308	// Highest bound in *_hi.mtx
#define PP_SIGMA_TO_APEX			5000
#define PP_LOW_COST_PERCENTILE		0		// Level below which value of low cost variable falls (must be in [0,1])
#define PP_GAP						1E-2	// Maximum gap from polytope surface (not too small!!!)
#define PP_OBJECTIVE_VECTOR_LENGTH	0.1		// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		5		// Start length of shift vector
#define PP_MODE_BLOCK_HCV_VARIABLE	true	// block high cos variables
#define PP_MODE_USE_LCV_VARIABLE	false	// utilize low cost and zero variables
//------------------------------------------------------------------/**/

/**#define PP_PROBLEM_NAME		"simple3" //==========================================
#define PP_M 6		// Number of equations (number of rows in *.mtx)
#define PP_N 8		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define EPS							1E-8
#define PP_EPS_DIR_LENGTH			EPS		// Minimal Length of Direction Vector
#define PP_EPS_OBJ					EPS		// |F(u) - PP_EXACT_OBJ_VALUE| <= PP_EPS_OBJ
#define PP_EPS_SHIFT				EPS		// Minimal shift to stop motion
#define PP_EPS_ZERO_COMPARE			EPS		// Zero comparison precision
#define PP_EPS_ZERO_DIR				EPS		// Accuracy of zero value for direction vector coordinates for direction vector coordinates
#define PP_EXACT_OBJ_VALUE			60000	// Exact maximum value of objective function
#define PP_INFINITY					1E+308	// Highest bound in *_hi.mtx
#define PP_SIGMA_TO_APEX			5000
#define PP_LOW_COST_PERCENTILE		0		// Level below which value of low cost variable falls (must be in [0,1])
#define PP_GAP						1E-2	// Maximum gap from polytope surface (not too small!!!)
#define PP_OBJECTIVE_VECTOR_LENGTH	0.1		// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		5		// Start length of shift vector
#define PP_MODE_BLOCK_HCV_VARIABLE	true	// block high cos variables
#define PP_MODE_USE_LCV_VARIABLE	false	// utilize low cost and zero variables
//------------------------------------------------------------------/**/

/**#define PP_PROBLEM_NAME		"simple1min" //==========================================
#define PP_M 5		// Number of equations (number of rows in *.mtx)
#define PP_N 8		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define EPS							1E-8
#define PP_EPS_DIR_LENGTH			EPS		// Minimal Length of Direction Vector
#define PP_EPS_OBJ					1E-1	// |F(u) - PP_EXACT_OBJ_VALUE| <= PP_EPS_OBJ
#define PP_EPS_SHIFT				EPS		// Minimal shift to stop motion
#define PP_EPS_ZERO_COMPARE			EPS		// Zero comparison precision
#define PP_EPS_ZERO_DIR				1E-5	// Accuracy of zero value for direction vector coordinates for direction vector coordinates
#define PP_EXACT_OBJ_VALUE			-5000	// Exact maximum value of objective function
#define PP_INFINITY					1E+308	// Highest bound in *_hi.mtx
#define PP_SIGMA_TO_APEX			10
#define PP_LOW_COST_PERCENTILE		0		// Level below which value of low cost variable falls (must be in [0,1])
#define PP_GAP						1E-2	// Maximum gap from polytope surface (not too small!!!)
#define PP_OBJECTIVE_VECTOR_LENGTH	0.1		// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		5		// Start length of shift vector
#define PP_MODE_BLOCK_HCV_VARIABLE	true	// block high cos variables
#define PP_MODE_USE_LCV_VARIABLE	false	// utilize low cost and zero variables
//------------------------------------------------------------------/**/

/**#define PP_PROBLEM_NAME		"simple_zcv" //==========================================
#define PP_M 4		// Number of equations (number of rows in *.mtx)
#define PP_N 7		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define EPS							1E-8
#define PP_EPS_DIR_LENGTH			EPS		// Minimal Length of Direction Vector
#define PP_EPS_OBJ					1E-1	// |F(u) - PP_EXACT_OBJ_VALUE| <= PP_EPS_OBJ
#define PP_EPS_SHIFT				EPS		// Minimal shift to stop motion
#define PP_EPS_ZERO_COMPARE			EPS		// Zero comparison precision
#define PP_EPS_ZERO_DIR				1E-5	// Accuracy of zero value for determining direction vector coordinates
#define PP_EXACT_OBJ_VALUE			50000	// Exact maximum value of objective function
#define PP_INFINITY					1E+308	// Highest bound in *_hi.mtx
#define PP_SIGMA_TO_APEX			5000
#define PP_LOW_COST_PERCENTILE		0		// Level below which value of low cost variable falls (must be in [0,1])
#define PP_GAP						1E-2	// Maximum gap from polytope surface (not too small!!!)
#define PP_OBJECTIVE_VECTOR_LENGTH	1		// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		5		// Start length of shift vector
#define PP_MODE_BLOCK_HCV_VARIABLE	false	// block high cos variables
#define PP_MODE_USE_LCV_VARIABLE	true	// utilize low cost and zero variables
//------------------------------------------------------------------/**/

/**#define PP_PROBLEM_NAME		"simple_lcv" //==========================================
#define PP_M 4		// Number of equations (number of rows in *.mtx)
#define PP_N 7		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define EPS							1E-8
#define PP_EPS_DIR_LENGTH			EPS		// Minimal Length of Direction Vector
#define PP_EPS_OBJ					1E-1	// |F(u) - PP_EXACT_OBJ_VALUE| <= PP_EPS_OBJ
#define PP_EPS_SHIFT				EPS		// Minimal shift to stop motion
#define PP_EPS_ZERO_COMPARE			EPS		// Zero comparison precision
#define PP_EPS_ZERO_DIR				1E-5	// Accuracy of zero value for determining direction vector coordinates
#define PP_EXACT_OBJ_VALUE			50000.2	// Exact maximum value of objective function
#define PP_INFINITY					1E+308	// Highest bound in *_hi.mtx
#define PP_SIGMA_TO_APEX			5000
#define PP_LOW_COST_PERCENTILE		0.01	// Level below which value of low cost variable falls (must be in [0,1])
#define PP_GAP						1E-2	// Maximum gap from polytope surface (not too small!!!)
#define PP_OBJECTIVE_VECTOR_LENGTH	1		// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		5		// Start length of shift vector
#define PP_MODE_BLOCK_HCV_VARIABLE	false	// block high cos variables
#define PP_MODE_USE_LCV_VARIABLE	true	// utilize low cost and zero variables
//------------------------------------------------------------------/**/

/**#define PP_PROBLEM_NAME		"simple_lcv_neg" //==========================================
#define PP_M 4		// Number of equations (number of rows in *.mtx)
#define PP_N 7		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define EPS							1E-8
#define PP_EPS_DIR_LENGTH			EPS		// Minimal Length of Direction Vector
#define PP_EPS_OBJ					1E-1	// |F(u) - PP_EXACT_OBJ_VALUE| <= PP_EPS_OBJ
#define PP_EPS_SHIFT				EPS		// Minimal shift to stop motion
#define PP_EPS_ZERO_COMPARE			EPS		// Zero comparison precision
#define PP_EPS_ZERO_DIR				1E-5	// Accuracy of zero value for determining direction vector coordinates
#define PP_EXACT_OBJ_VALUE			49999.8	// Exact maximum value of objective function
#define PP_INFINITY					1E+308	// Highest bound in *_hi.mtx
#define PP_SIGMA_TO_APEX			5000
#define PP_LOW_COST_PERCENTILE		0.01	// Level below which value of low cost variable falls (must be in [0,1])
#define PP_GAP						1E-2	// Maximum gap from polytope surface (not too small!!!)
#define PP_OBJECTIVE_VECTOR_LENGTH	1		// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		5		// Start length of shift vector
#define PP_MODE_BLOCK_HCV_VARIABLE	false	// block high cos variables
#define PP_MODE_USE_LCV_VARIABLE	true	// utilize low cost and zero variables
//------------------------------------------------------------------/**/
/**/#define PP_PROBLEM_NAME		"dataset" //==========================================
#define PP_M 11		// Number of equations (number of rows in *.mtx)
#define PP_N 14		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define EPS							1E-8
#define PP_EPS_DIR_LENGTH			EPS		// Minimal Length of Direction Vector
#define PP_EPS_OBJ					1E-1	// |F(u) - PP_EXACT_OBJ_VALUE| <= PP_EPS_OBJ
#define PP_EPS_SHIFT				EPS		// Minimal shift to stop motion
#define PP_EPS_ZERO_COMPARE			EPS		// Zero comparison precision
#define PP_EPS_ZERO_DIR				1E-5	// Accuracy of zero value for determining direction vector coordinates
#define PP_EPS_DISTANCE				1.0		// Minimal distance between points in a trace
#define PP_EXACT_OBJ_VALUE			50000	// Exact maximum value of objective function
#define PP_INFINITY					1E+308	// Highest bound in *_hi.mtx
#define PP_SIGMA_TO_APEX			0	// Distance from apex base to apex point
#define PP_LOW_COST_PERCENTILE		0.01	// Percentile for low cost variable (must be in [0,1])
#define PP_GAP						1E-2	// Maximum gap from polytope surface (not too small!!!)
#define PP_OBJECTIVE_VECTOR_LENGTH	0.1		// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		5		// Start length of shift vector
#define PP_MODE_BLOCK_HCV_VARIABLE	true	// block high cost variables
#define PP_MODE_USE_LCV_VARIABLE	false	// utilize low cost and zero variables
//------------------------------------------------------------------/**/
//================================ Common Paramrters ===========================
#define PP_MM (PP_MODE_BLOCK_HCV_VARIABLE ? 2*PP_M+3*PP_N-1 : 2*PP_M+2*PP_N)
#define PP_MAX_NUM_SHIFTS_SAME_LENGTH	5 // Maximal number of shifts with the same length
#define PP_MAX_ITER_COUNT				40000 // Maximal count of iterations
#define PP_ADD_FLAG						PP_N
//#define PP_RANDOM_X0
//-------------------------- Input/Outpoot Parameters ---------------------------
#define PP_OUTPUT_LIMIT	8	// Number of Elements to output
#define PP_TRACE_LIMIT 50	// Maximum points in a problem trace
#define PP_MATRIX_OUTPUT	// To output Matrix
#define PP_SETW 18
//
//------------------------- Matrix format ----------------
#define PP_MTX_PREFIX			"lp_"
#define PP_MTX_POSTFIX_A		"_A.mtx.pack"
#define PP_MTX_POSTFIX_B		"_b.mtx.pack"
#define PP_MTX_POSTFIX_LO		"_lo.mtx.pack"
#define PP_MTX_POSTFIX_HI		"_hi.mtx.pack"
#define PP_MTX_POSTFIX_C		"_c.mtx.pack"
#define PP_MTX_POSTFIX_X0		"_x0.mtx.pack" // Starting point on polytope
#define PP_MTX_POSTFIX_TR		"_tr.mtx.pack"
#define PP_MTX_POSTFIX_SO		"_so.mtx" // Solution
//
//-------------------------- Jobs  -----------------------
#define PP_JOB_PSEUDOPOJECTION	0 
#define PP_JOB_CHECK			1		
//-------------------------- Process states --------------------------
#define PP_STATE_START						0
#define PP_STATE_FIND_INTERIOR_POINT		1
#define PP_STATE_FIND_INITIAL_APPROXIMATION	2
#define PP_STATE_DETERMINE_DIRECTION		3
#define PP_STATE_MOVING_ALONG_SURFACE		4
#define PP_STATE_LANDING					5

//------------- Vector_ProjectOnHalfspace exit codes -------------
#define PP_EXITCODE_DEGENERATE_INEQUALITY		0
#define PP_EXITCODE_POINT_BELONGS_TO_HALFSPACE	1
#define PP_EXITCODE_SUCCESSFUL_PROJECTING		2