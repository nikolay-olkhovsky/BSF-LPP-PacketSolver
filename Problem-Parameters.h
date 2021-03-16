/*==============================================================================
Project: LiFe
Theme: Apex Method (Predictor + Corrector)
Module: Problem-Parameters.h (Problem Parameters)
Prefix: PP
Author(s): Leonid B. Sokolinsky, Irina M. Sokolinskaya
This source code has been produced with using BSF-skeleton
==============================================================================*/
//-------------------------- Compilation Modes -----------------------
#define PP_DEBUG
//=========================== Problem Parameters =========================
#define PP_N 4				// n - Dimension of Space
#define PP_NUM_OF_NATURAL_INEQUALITIES (PP_N)				// Number of natural inequalities
#define PP_M (2*PP_N + PP_NUM_OF_NATURAL_INEQUALITIES + 1)	// Total number of inequalities of given system
#define PP_MM (PP_M + 2 * (PP_N - 2))	// Maximal number of inequalities including additional
#define PP_SF 200			// Scale factor
#define PP_DIST_TO_APEX (PP_N*PP_SF*100000)// Distance to Apex Point 

#define PP_EPS_RELAX		1E-11			// Precision for relaxation
#define PP_EPS_SHIFT		1E-10
#define PP_EPS_IN			1E-1
#define PP_EPS_ZERO			1E-5
#define PP_EPS_DIR			1E-3
#define PP_EPS_OBJECTIVE	1E-1
#define PP_MAX_NUM_SHIFTS_SAME_LENGTH 5	// Maximal number of shifts with the same length
#define PP_MAX_NUM_SEQ_SHIFTS 500		// Maximal number of sequential shifts

#define PP_MAX_ITER_COUNT	10000000					// Maximal count of iterations
#define PP_OBJECTIVE_VECTOR_LENGTH ((double)PP_SF/200)	// Length of Objective Vector
#define PP_START_SHIFT_LENGTH (PP_SF/20)				// Start length of shift vector

//-------------------------- Outpoot Parameters ---------------------------
#define PP_OUTPUT_LIMIT	8	// Number of Elements to output
#define PP_MATRIX_OUTPUT	// To output Matrix
//#define PP_SETW 14
#define PP_SETW 20
//-------------------------- Jobs  -----------------------
#define PP_JOB_PSEUDOPOJECTION	0 
#define PP_JOB_CHECK			1		
#define PP_JOB_HYPERPLANE_TOUCH	2
//-------------------------- Process states --------------------------
#define PP_STATE_START						0
#define PP_STATE_FIND_BEGINNING_OF_PATH		1
#define PP_STATE_DETERMINE_DIRECTION		2
#define PP_STATE_MOVE_AND_CHECK				3
//-------------------------- Input ----------------------
#define PP_LPP_FILE "lpp.txt"
// Input data file in the following format:
/*
------------ begin of file -------------
m n
A_11 A_12 ... A_1n b_1
A_21 A_22 ... A_2n b_2
...
A_m1 A_m2 ... A_mn b_m
c_1 c_2 ... c_n
------------ end of file----------------
*/

//-------------------------- Output ---------------------------------------------
#define PP_SOLUTION_FILE "solution.txt" // Solution of LPP	|
// Output data file in the following format:
/*
------------ begin of file -------------
n
x_1 x_2 ... x_n
------------ end of file----------------
*/