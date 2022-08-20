/*==============================================================================
Project: LiFe
Theme: Apex Method (No MPI)
Module: Problem-Parameters.h (Problem Parameters)
Prefix: PP
Author: Leonid B. Sokolinsky
This source code has been produced with using BSF-skeleton
==============================================================================*/
//-------------------------- Compilation Modes -----------------------
#define PP_DEBUG
#define PP_PATH "D:/OneDrive/_private/Программы/BSF/BSF Target/Problems/"
//#define PP_PATH ""
//=========================== Problem Parameters =========================
/**#define PP_PROBLEM_NAME		"afiro"		 //========
#define PP_M 27		// Number of equations (number of rows in *.mtx)
#define PP_N 51		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define EPS							1E-6
#define PP_EPS_DIR_LENGTH			EPS		// Minimal Length of Direction Vector
#define PP_EPS_OBJECTIVE			1E-3		// < |F(u)-F(w)|
#define PP_EPS_SHIFT				EPS		// Minimal shift to stop motion
#define PP_EPS_ZERO_COMPARE			EPS		// Zero comparison precision
#define PP_EPS_ZERO_DIR				1E-5		// Accuracy of zero value for direction vector coordinates for direction vector coordinates
#define PP_EXACT_OBJ_VALUE		464.7531	// Exact maximum value of objective function
#define PP_INFINITY					1E+308	// High bound in *_hi.mtx
#define PP_DISTANCE_TO_APEX			5000
#define PP_LOW_COST_PERCENTILE		0		// Level below which value of low cost variable falls (must be in [0,1])
#define PP_GAP						1E-4	// Maximum gap from polytope surface (not too small!!!)
#define PP_OBJECTIVE_VECTOR_LENGTH	100		// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		20		// Start length of shift vector
#define PP_BLOCK_OBJ_VARIABLE		false	// straightens the trace, but can lead to an incorrect solution
//------------------------------------------------------------------/**/

/**#define PP_PROBLEM_NAME		"adlittle"	//========
#define PP_M 56		// Number of equations (number of rows in *.mtx)
#define PP_N 138	// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define EPS							1E-8
#define PP_EPS_DIR_LENGTH			EPS		// Minimal Length of Direction Vector
#define PP_EPS_OBJECTIVE			1E-3	// < |F(u)-F(w)|
#define PP_EPS_SHIFT				EPS		// Minimal shift to stop motion
#define PP_EPS_ZERO_COMPARE			EPS		// Zero comparison precision
#define PP_EPS_ZERO_DIR				1E-7	// Accuracy of zero value for direction vector coordinates for direction vector coordinates
#define PP_EXACT_OBJ_VALUE		-225494.963	// Exact maximum value of objective function
#define PP_INFINITY					1E+308	// High bound in *_hi.mtx
//#define PP_DISTANCE_TO_APEX			200000
#define PP_DISTANCE_TO_APEX			0
#define PP_GAP						1E-3	// Maximum gap from polytope surface (not too small!!!)
#define PP_LOW_COST_PERCENTILE		0.001	// Level below which value of low cost variable falls (must be in [0,1])
#define PP_OBJECTIVE_VECTOR_LENGTH	1E-4	// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		1		// Start length of shift vector
#define PP_BLOCK_OBJ_VARIABLE		false	// straightens the trace, but can lead to an incorrect solution
//------------------------------------------------------------------/**/

/**/#define PP_PROBLEM_NAME		"blend"	//==========================================
#define PP_M 74		// Number of equations (number of rows in *.mtx)
#define PP_N 114	// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define EPS							1E-8
#define PP_EPS_DIR_LENGTH			EPS		// Minimal Length of Direction Vector
#define PP_EPS_OBJECTIVE			1E-7		// < |F(u)-F(w)|
#define PP_EPS_SHIFT				EPS		// Minimal shift to stop motion
#define PP_EPS_ZERO_COMPARE			EPS		// Zero comparison precision
#define PP_EPS_ZERO_DIR				1E-7	// Accuracy of zero value for direction vector coordinates
#define PP_EXACT_OBJ_VALUE			30.812	// Exact maximum value of objective function
#define PP_INFINITY					1E+308	// Highest bound in *_hi.mtx
//#define PP_DISTANCE_TO_APEX		300
#define PP_DISTANCE_TO_APEX			0
#define PP_GAP						5E-3	// Minimum distance to polytope (not too small!!!)
#define PP_LOW_COST_PERCENTILE		0.01
#define PP_OBJECTIVE_VECTOR_LENGTH	1E-4	// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		1		// Start length of shift vector
#define PP_BLOCK_OBJ_VARIABLE		false	// straightens the trace, but can lead to an incorrect solution
//------------------------------------------------------------------/**/

/**#define PP_PROBLEM_NAME		"fit1d"	// -9146.378  //==========================================
#define PP_M 24		// Number of equations (number of rows in *.mtx)
#define PP_N 1049	// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define PP_EPS_DIR_LENGTH			1E-4	// Minimal Length of Direction Vector
#define PP_EPS_OBJECTIVE			1E-5	// < |F(u)-F(w)|
#define PP_EPS_SHIFT				1E-6	// Minimal shift to stop motion
#define PP_EPS_ZERO_DIR				1E-4	// Accuracy of zero value for direction vector coordinates
#define PP_EPS_ZERO_COMPARE			1E-3	// Comparison precision
#define PP_EXACT_OBJ_VALUE			-9146.378	// Exact maximum value of objective function
#define PP_GAP_MIN					1E-2	// Minimum distance to polytope (not too small!!!)
#define PP_GAP_MAX					1E-2		// Maximum distance from polytope
#define PP_OBJECTIVE_VECTOR_LENGTH	1E-1	// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		1		// Start length of shift vector
#define PP_BLOCK_OBJ_VARIABLE		true	// straightens the trace, but can lead to an incorrect solution
//------------------------------------------------------------------/**/

/**#define PP_PROBLEM_NAME		"kb2"	// -1749.9  //==========================================
#define PP_M 43	// Number of equations (number of rows in *.mtx)
#define PP_N 68	// Number of variables (number of cols in *.mtx)/**/

/**#define PP_PROBLEM_NAME		"recipe"	// -266.61  //==========================================
#define PP_M 91		// Number of equations (number of rows in *.mtx)
#define PP_N 204	// Number of variables (number of cols in *.mtx)/**/

/**#define PP_PROBLEM_NAME		"sc50a"	// -64.575
#define PP_M 50		// Number of equations (number of rows in *.mtx)
#define PP_N 78	// Number of variables (number of cols in *.mtx)/**/

/**#define PP_PROBLEM_NAME		"sc50b"	//==========================================
#define PP_M 50		// Number of equations (number of rows in *.mtx)
#define PP_N 78	// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define PP_GAP_MIN					1E-2	// Minimal distance to polytope (not too small!!!)
#define PP_EPS_SHIFT				1E-4	// Minimal shift to stop motion
#define PP_EPS_ZERO_DIR					1E-6	// Accuracy of zero value for direction vector coordinates
#define PP_EPS_ZERO_COMPARE				1E-6	// Comparison precision
#define PP_EPS_DIR_LENGTH					1E-6	// Minimal Length of Direction Vector
#define PP_EPS_OBJECTIVE			1E-4	// < |F(u)-F(w)|
#define PP_EXACT_OBJ_VALUE			-70		// Exact maximum value of objective function
#define PP_OBJECTIVE_VECTOR_LENGTH	45		// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		5		// Start length of shift vector
//------------------------------------------------------------------/**/

/**#define PP_PROBLEM_NAME		"share2b" // -415.73  //==========================================
#define PP_M 96		// Number of equations (number of rows in *.mtx)
#define PP_N 162	// Number of variables (number of cols in *.mtx)/**/

/**#define PP_PROBLEM_NAME		"simple1" //==========================================
#define PP_M 4		// Number of equations (number of rows in *.mtx)
#define PP_N 7		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define EPS							1E-8
#define PP_EPS_DIR_LENGTH			EPS		// Minimal Length of Direction Vector
#define PP_EPS_OBJECTIVE			EPS		// < |F(u)-F(w)|
#define PP_EPS_SHIFT				EPS		// Minimal shift to stop motion
#define PP_EPS_ZERO_COMPARE			EPS		// Zero comparison precision
#define PP_EPS_ZERO_DIR				EPS		// Accuracy of zero value for direction vector coordinates for direction vector coordinates
#define PP_EXACT_OBJ_VALUE			55000	// Exact maximum value of objective function
#define PP_INFINITY					1E+308	// Highest bound in *_hi.mtx
#define PP_DISTANCE_TO_APEX			5000
#define PP_LOW_COST_PERCENTILE		0		// Level below which value of low cost variable falls (must be in [0,1])
#define PP_GAP						1E-2	// Maximum gap from polytope surface (not too small!!!)
#define PP_OBJECTIVE_VECTOR_LENGTH	0.1		// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		5		// Start length of shift vector
#define PP_BLOCK_OBJ_VARIABLE		true	// straightens the trace, but can lead to an incorrect solution
//------------------------------------------------------------------/**/

/**#define PP_PROBLEM_NAME		"simple1.1" //==========================================
#define PP_M 4		// Number of equations (number of rows in *.mtx)
#define PP_N 7		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define EPS							1E-8
#define PP_EPS_DIR_LENGTH			EPS		// Minimal Length of Direction Vector
#define PP_EPS_OBJECTIVE			EPS		// < |F(u)-F(w)|
#define PP_EPS_SHIFT				EPS		// Minimal shift to stop motion
#define PP_EPS_ZERO_COMPARE			EPS		// Zero comparison precision
#define PP_EPS_ZERO_DIR				EPS		// Accuracy of zero value for direction vector coordinates for direction vector coordinates
#define PP_EXACT_OBJ_VALUE			40000	// Exact maximum value of objective function
#define PP_INFINITY					1E+308	// Highest bound in *_hi.mtx
#define PP_DISTANCE_TO_APEX			5000
#define PP_LOW_COST_PERCENTILE		0		// Level below which value of low cost variable falls (must be in [0,1])
#define PP_GAP						1E-2	// Maximum gap from polytope surface (not too small!!!)
#define PP_OBJECTIVE_VECTOR_LENGTH	0.1		// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		5		// Start length of shift vector
#define PP_BLOCK_OBJ_VARIABLE		true	// straightens the trace, but can lead to an incorrect solution
//------------------------------------------------------------------/**/

/**#define PP_PROBLEM_NAME		"simple2" //==========================================
#define PP_M 5		// Number of equations (number of rows in *.mtx)
#define PP_N 8		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define EPS							1E-8
#define PP_EPS_DIR_LENGTH			EPS		// Minimal Length of Direction Vector
#define PP_EPS_OBJECTIVE			EPS		// < |F(u)-F(w)|
#define PP_EPS_SHIFT				EPS		// Minimal shift to stop motion
#define PP_EPS_ZERO_COMPARE			EPS		// Zero comparison precision
#define PP_EPS_ZERO_DIR				EPS		// Accuracy of zero value for direction vector coordinates for direction vector coordinates
#define PP_EXACT_OBJ_VALUE			59000	// Exact maximum value of objective function
#define PP_INFINITY					1E+308	// Highest bound in *_hi.mtx
#define PP_DISTANCE_TO_APEX			5000
#define PP_LOW_COST_PERCENTILE		0		// Level below which value of low cost variable falls (must be in [0,1])
#define PP_GAP						1E-2	// Maximum gap from polytope surface (not too small!!!)
#define PP_OBJECTIVE_VECTOR_LENGTH	0.1		// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		5		// Start length of shift vector
#define PP_BLOCK_OBJ_VARIABLE		true	// straightens the trace, but can lead to an incorrect solution
//------------------------------------------------------------------/**/

/**#define PP_PROBLEM_NAME		"simple3" //==========================================
#define PP_M 6		// Number of equations (number of rows in *.mtx)
#define PP_N 8		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define EPS							1E-8
#define PP_EPS_DIR_LENGTH			EPS		// Minimal Length of Direction Vector
#define PP_EPS_OBJECTIVE			EPS		// < |F(u)-F(w)|
#define PP_EPS_SHIFT				EPS		// Minimal shift to stop motion
#define PP_EPS_ZERO_COMPARE			EPS		// Zero comparison precision
#define PP_EPS_ZERO_DIR				EPS		// Accuracy of zero value for direction vector coordinates for direction vector coordinates
#define PP_EXACT_OBJ_VALUE			60000	// Exact maximum value of objective function
#define PP_INFINITY					1E+308	// Highest bound in *_hi.mtx
#define PP_DISTANCE_TO_APEX			5000
#define PP_LOW_COST_PERCENTILE		0		// Level below which value of low cost variable falls (must be in [0,1])
#define PP_GAP						1E-2	// Maximum gap from polytope surface (not too small!!!)
#define PP_OBJECTIVE_VECTOR_LENGTH	0.1		// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		5		// Start length of shift vector
#define PP_BLOCK_OBJ_VARIABLE		true	// straightens the trace, but can lead to an incorrect solution
//------------------------------------------------------------------/**/

/**#define PP_PROBLEM_NAME		"simple1min" //==========================================
#define PP_M 5		// Number of equations (number of rows in *.mtx)
#define PP_N 8		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define EPS							1E-8
#define PP_EPS_DIR_LENGTH			EPS		// Minimal Length of Direction Vector
#define PP_EPS_OBJECTIVE			1E-1		// < |F(u)-F(w)|
#define PP_EPS_SHIFT				EPS		// Minimal shift to stop motion
#define PP_EPS_ZERO_COMPARE			EPS		// Zero comparison precision
#define PP_EPS_ZERO_DIR				1E-5		// Accuracy of zero value for direction vector coordinates for direction vector coordinates
#define PP_EXACT_OBJ_VALUE			-5000	// Exact maximum value of objective function
#define PP_INFINITY					1E+308	// Highest bound in *_hi.mtx
#define PP_DISTANCE_TO_APEX			10
#define PP_LOW_COST_PERCENTILE		0		// Level below which value of low cost variable falls (must be in [0,1])
#define PP_GAP						1E-2	// Maximum gap from polytope surface (not too small!!!)
#define PP_OBJECTIVE_VECTOR_LENGTH	0.1		// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		5		// Start length of shift vector
#define PP_BLOCK_OBJ_VARIABLE		true	// straightens the trace, but can lead to an incorrect solution
//------------------------------------------------------------------/**/

/**#define PP_PROBLEM_NAME		"simple1free" //==========================================
#define PP_M 4		// Number of equations (number of rows in *.mtx)
#define PP_N 7		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define EPS							1E-8
#define PP_EPS_DIR_LENGTH			EPS		// Minimal Length of Direction Vector
#define PP_EPS_OBJECTIVE			5E-1	// < |F(u)-F(w)|
#define PP_EPS_SHIFT				EPS		// Minimal shift to stop motion
#define PP_EPS_ZERO_COMPARE			EPS		// Zero comparison precision
#define PP_EPS_ZERO_DIR				1E-5	// Accuracy of zero value for determining direction vector coordinates
#define PP_EXACT_OBJ_VALUE			50000	// Exact maximum value of objective function
#define PP_INFINITY					1E+308	// Highest bound in *_hi.mtx
#define PP_DISTANCE_TO_APEX			5000
#define PP_LOW_COST_PERCENTILE		0		// Level below which value of low cost variable falls (must be in [0,1])
#define PP_GAP						1E-2	// Maximum gap from polytope surface (not too small!!!)
#define PP_OBJECTIVE_VECTOR_LENGTH	1		// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		5		// Start length of shift vector
#define PP_BLOCK_OBJ_VARIABLE		true	// straightens the trace, but can lead to an incorrect solution
//------------------------------------------------------------------/**/
//================================ Common Paramrters ===========================
#define PP_MM (PP_BLOCK_OBJ_VARIABLE ? 2*PP_M+3*PP_N-1 : 2*PP_M+2*PP_N)
#define PP_MAX_NUM_SHIFTS_SAME_LENGTH	5 // Maximal number of shifts with the same length
#define PP_MAX_ITER_COUNT				10000000000 // Maximal count of iterations
#define PP_ADD_FLAG						PP_N
//-------------------------- Input/Outpoot Parameters ---------------------------
#define PP_OUTPUT_LIMIT	8	// Number of Elements to output
//#define PP_MATRIX_OUTPUT	// To output Matrix
#define PP_SETW 16
#define PP_TRACE_FILE "trace.txt"
//------------------------- Matrix format ----------------
#define PP_MTX_PREFIX			"lp_"
#define PP_MTX_POSTFIX_A		".mtx"
#define PP_MTX_POSTFIX_B		"_b.mtx"
#define PP_MTX_POSTFIX_LO		"_lo.mtx"
#define PP_MTX_POSTFIX_HI		"_hi.mtx"
#define PP_MTX_POSTFIX_C		"_c.mtx"
#define PP_MTX_POSTFIX_X0		"_x0.mtx" // Starting point on polytope
#define PP_MTX_POSTFIX_SO		"_so.mtx" // Solution
#define PP_MTX_POSTFIX_TR		"_tr.mtx" // Trace
//-------------------------- Jobs  -----------------------
#define PP_JOB_PSEUDOPOJECTION	0 
#define PP_JOB_CHECK			1		
#define PP_JOB_CHECK_S			2		
//-------------------------- Process states --------------------------
#define PP_STATE_START					0
#define PP_STATE_DETERMINE_DIRECTION	1
#define MOVING_ALONG_SURFACE			2
#define PP_STATE_LANDING				3
#define PP_STATE_FIND_FEASIBLE_POINT	4

//------------- Vector_ProjectOnHalfspace exit codes -------------
#define PP_EXITCODE_DEGENERATE_INEQUALITY		0
#define PP_EXITCODE_POINT_BELONGS_TO_HALFSPACE	1
#define PP_EXITCODE_SUCCESSFUL_PROJECTING		2