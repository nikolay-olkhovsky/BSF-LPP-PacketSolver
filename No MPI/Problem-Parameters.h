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
#define MTX_FORMAT
//#define PP_PAUSE
//#define PP_MAJOR_COORDINATES_CAN_NOT_DECREASE // straightens the trace, but can lead to an incorrect solution
//=========================== Problem Parameters =========================
#define PP_M 8		// simple							// Maximal number of inequalities
#define PP_N 7		// simple							// Maximal Space Dimension
//#define PP_M 27		// afiro							// Maximal number of inequalities
//#define PP_N 51		// afiro							// Maximal Space Dimension
//#define PP_M 56		// adlittle							// Maximal number of inequalities
//#define PP_N 138		// adlittle							// Maximal Space Dimension
//--------------------------------------

#ifdef MTX_FORMAT
#define PP_MM (2*PP_M+PP_N)	// Maximal Total number of inequalities
#else
#define PP_MM PP_M			// Maximal Total number of inequalities
#endif

#ifdef PP_MAJOR_COORDINATES_CAN_NOT_DECREASE
#define PP_MAX_MM (PP_M + PP_N - 1)	// Maximal number of inequalities including additional
#else
#define PP_MAX_MM PP_MM			// Maximal number of inequalities including additional
#endif

//--------------------------------------
#define PP_SF 200													// Scale factor

#define PP_EPS_RELAX		1E-11									// Precision for relaxation
#define PP_EPS_IN			1E-2									// Not too small!
#define PP_EPS_SHIFT		1E-3
#define PP_EPS_ZERO			1E-8
#define PP_EPS_DIR			1E-6
#define PP_EPS_OBJECTIVE	1E-2

#define PP_MAX_NUM_SHIFTS_SAME_LENGTH 5						// Maximal number of shifts with the same length

#define PP_MAX_ITER_COUNT	1000000000						// Maximal count of iterations
//#define PP_OBJECTIVE_VECTOR_LENGTH ((PT_float_T)PP_SF/200)	// Length of Objective Vector
#define PP_OBJECTIVE_VECTOR_LENGTH 1
//#define PP_START_SHIFT_LENGTH (PP_SF/20)					// Start length of shift vector
#define PP_START_SHIFT_LENGTH 10					// Start length of shift vector

//-------------------------- Input/Outpoot Parameters ---------------------------
#define PP_OUTPUT_LIMIT	8	// Number of Elements to output
#define PP_MATRIX_OUTPUT	// To output Matrix
#define PP_SETW 12
#define PP_PATH "C:/TEMP/"
//#define PP_PATH ""
#define PP_SOLUTION_FILE "solution.txt"
#define PP_TRACE_FILE "trace.txt"
//------------------------- LPP format ----------------
#define PP_LPP_FILE "lpp.txt"
//------------------------- Matrix format ----------------
#define PP_MTX_PROBLEM_NAME		"simple"
//#define PP_MTX_PROBLEM_NAME		"afiro"
//#define PP_MTX_PROBLEM_NAME		"adlittle"
#define PP_MTX_PREFIX			"lp_"
#define PP_MTX_POSTFIX_A		".mtx"
#define PP_MTX_POSTFIX_B		"_b.mtx"
#define PP_MTX_POSTFIX_LO		"_lo.mtx"
#define PP_MTX_POSTFIX_HI		"_hi.mtx"
#define PP_MTX_POSTFIX_C		"_c.mtx"
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