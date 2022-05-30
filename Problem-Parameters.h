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
#define MTX_FORMAT
#define PP_PATH "C:/TEMP/"
//#define PP_PATH "Problems/"
//#define PP_PAUSE
//#define PP_MAJOR_COORDINATES_CAN_NOT_DECREASE	// straightens the trace, but can lead to an incorrect solution
//=========================== Problem Parameters =========================
#ifdef MTX_FORMAT

/**#define PP_MTX_PROBLEM_NAME		"afiro"		//  464.7531
#define PP_M 27		// Number of equations (number of rows in *.mtx)
#define PP_N 51		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define PP_EPS_RELAX				1E-10	// Precision for relaxation
#define PP_EPS_IN					1E-2	// Minimal distance to polytope (not too small!!!)
#define PP_EPS_SHIFT				1E-4	// Minimal shift
#define PP_EPS_ZERO					1E-6	// Accuracy of zero value
#define PP_EPS_COMPARE				1E-6	// Comparison precision
#define PP_EPS_DIR					1E-6	// Minimal Length of Direction Vector
#define PP_EPS_OBJECTIVE			1E-2	// < |F(u)-F(u)|
#define PP_LAMBDA					1.99	// 0<PP_LAMBDA<2 relaxation speed
#define PP_OBJECTIVE_VECTOR_LENGTH	50		// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		5	// Start length of shift vector
//------------------------------------------------------------------/**/
/**/

/**/#define PP_MTX_PROBLEM_NAME		"adlittle"	// -225494.963
#define PP_M 56		// Number of equations (number of rows in *.mtx)
#define PP_N 138	// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define PP_EPS_RELAX				1E-10	// Precision for relaxation
//#define PP_EPS_IN					5E-1	// Minimal distance to polytope (not too small!!!)
#define PP_EPS_IN					1E-1	// Minimal distance to polytope (not too small!!!)
#define PP_EPS_SHIFT				1E-4	// Minimal shift
#define PP_EPS_ZERO					1E-6	// Accuracy of zero value
#define PP_EPS_COMPARE				1E-6	// Comparison precision
#define PP_EPS_DIR					1E-6	// Minimal Length of Direction Vector
#define PP_EPS_OBJECTIVE			1E-2	// < |F(u)-F(u)|
#define PP_LAMBDA					1.99	// 0<PP_LAMBDA<2 relaxation speed
#define PP_OBJECTIVE_VECTOR_LENGTH	0.1		// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		1		// Start length of shift vector
//------------------------------------------------------------------/**/
/**/

/**#define PP_MTX_PROBLEM_NAME		"blend"	// -30.812
#define PP_M 74		// Number of equations (number of rows in *.mtx)
#define PP_N 114	// Number of variables (number of cols in *.mtx)/**/

/**#define PP_MTX_PROBLEM_NAME		"fit1d"	// -9146.378
#define PP_M 24		// Number of equations (number of rows in *.mtx)
#define PP_N 1049	// Number of variables (number of cols in *.mtx)/**/

/**#define PP_MTX_PROBLEM_NAME		"kb2"	// -1749.9
#define PP_M 43	// Number of equations (number of rows in *.mtx)
#define PP_N 68	// Number of variables (number of cols in *.mtx)/**/

/**#define PP_MTX_PROBLEM_NAME		"recipe"	// -266.61
#define PP_M 91		// Number of equations (number of rows in *.mtx)
#define PP_N 204	// Number of variables (number of cols in *.mtx)/**/

/**#define PP_MTX_PROBLEM_NAME		"sc50a"	// -64.575
#define PP_M 50		// Number of equations (number of rows in *.mtx)
#define PP_N 78	// Number of variables (number of cols in *.mtx)/**/

/**#define PP_MTX_PROBLEM_NAME		"sc50b"	// -70
#define PP_M 50		// Number of equations (number of rows in *.mtx)
#define PP_N 78	// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define PP_EPS_RELAX				1E-11	// Precision for relaxation
#define PP_EPS_IN					1E-2	// Minimal distance to polytope (not too small!!!)
#define PP_EPS_SHIFT				1E-4	// Minimal shift
#define PP_EPS_ZERO					1E-6	// Accuracy of zero value
#define PP_EPS_COMPARE				1E-6	// Comparison precision
#define PP_EPS_DIR					1E-6	// Minimal Length of Direction Vector
#define PP_EPS_OBJECTIVE			1E-4	// < |F(u)-F(u)|
#define PP_LAMBDA					1.99	// 0<PP_LAMBDA<2 relaxation speed
#define PP_OBJECTIVE_VECTOR_LENGTH	45		// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		5	// Start length of shift vector
//------------------------------------------------------------------/**/

/**#define PP_MTX_PROBLEM_NAME		"share2b" // -415.73
#define PP_M 96		// Number of equations (number of rows in *.mtx)
#define PP_N 162	// Number of variables (number of cols in *.mtx)/**/

/**#define PP_MTX_PROBLEM_NAME		"simple1"
#define PP_M 4		// Number of equations (number of rows in *.mtx)
#define PP_N 7		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------
#define PP_EPS_RELAX				1E-10	// Precision for relaxation
#define PP_EPS_IN					1E-2	// Minimal distance to polytope (not too small!!!)
#define PP_EPS_SHIFT				1E-4	// Minimal shift
#define PP_EPS_ZERO					1E-8	// Accuracy of zero value
#define PP_EPS_COMPARE				1E-4	// Comparison precision
#define PP_EPS_DIR					1E-6	// Minimal Length of Direction Vector
#define PP_EPS_OBJECTIVE			1E-3	// < |F(u)-F(u)|
#define PP_LAMBDA					1		// 0<PP_LAMBDA<2 relaxation speed
#define PP_OBJECTIVE_VECTOR_LENGTH	1		// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		5		// Start length of shift vector
//------------------------------------------------------------------/**/

#define PP_MM (2*PP_M+2*PP_N)	// Maximal number of inequalities
#else
#define PP_MTX_PROBLEM_NAME		"lpp"
#define PP_M 7		// Number of inequalities
#define PP_N 3		// Space dimension
#ifdef  PP_MAJOR_COORDINATES_CAN_NOT_DECREASE
#define PP_MM (PP_M + PP_N - 1)	// Maximal number of inequalities including additional
#else
#define PP_MM PP_M			// Maximal number of inequalities including additional
#endif
#endif

//--------------------------------------
#define PP_SF 200					// Scale factor
#define PP_EPS_OBJECTIVE_VECTOR_LENGTH	0.01// Length of Objective Vector
#define PP_MAX_NUM_SHIFTS_SAME_LENGTH	5 // Maximal number of shifts with the same length
#define PP_MAX_ITER_COUNT				10000000000 // Maximal count of iterations

//-------------------------- Input/Outpoot Parameters ---------------------------
#define PP_OUTPUT_LIMIT	8	// Number of Elements to output
//#define PP_MATRIX_OUTPUT	// To output Matrix
#define PP_SETW 16
//------------------------- LPP format ----------------
#define PP_LPP_FILE "lpp.txt"
//------------------------- Matrix format ----------------
#define PP_MTX_PREFIX			"lp_"
#define PP_MTX_POSTFIX_A		".mtx"
#define PP_MTX_POSTFIX_B		"_b.mtx"
#define PP_MTX_POSTFIX_LO		"_lo.mtx"
#define PP_MTX_POSTFIX_HI		"_hi.mtx"
#define PP_MTX_POSTFIX_C		"_c.mtx"
#define PP_MTX_POSTFIX_X0		"_x0.mtx" // Starting point on polytope
#define PP_MTX_POSTFIX_SO		"_so.mtx" // Solution
//-------------------------- Jobs  -----------------------
#define PP_JOB_PSEUDOPOJECTION	0 
#define PP_JOB_CHECK			1		
#define PP_JOB_CHECK_S			2		
//-------------------------- Process states --------------------------
#define PP_STATE_START					0
#define PP_STATE_DETERMINE_DIRECTION	1
#define PP_STATE_MOVE_AND_CHECK			2
#define PP_STATE_LANDING				3
#define PP_MOVE_INSIDE_POLYTOPE			4
#define PP_STATE_FIND_START_POINT		5