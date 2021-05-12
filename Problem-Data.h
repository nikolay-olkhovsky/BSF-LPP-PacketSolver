/*==============================================================================
Project: LiFe
Theme: Apex Method
Module: Problem-Data.h (Problem Data)
Prefix: PD
Author: Leonid B. Sokolinsky
This source code has been produced with using BSF-skeleton
==============================================================================*/
#include "Problem-Types.h"		// Problem Parameters 
using namespace std;
//========================== Problem variables ====================================
static int		PD_m;						// Upper bound of valuable inequalities
static int		PD_numShiftsSameLength;		// Number of shifts with the same length
static int		PD_numSeqShifts;			// Number of sequential shifts
static int		PD_numDetDir;				// Number of sequential states "Determine Direction"
static bool		PD_newInequations;
static bool		PD_pointIn;
static double	PD_shiftLength;				// Shift length
static int		PD_state;
//========================== Problem structures ====================================
static PT_matrix_T PD_A;
static PT_column_T PD_b;					// Column of the constant terms of the system Ax <= PD_b
static PT_vector_T PD_c;					// Objective Function Coefficients
static PT_vector_T PD_apex;					// Apex Point
static PT_vector_T PD_basePoint;			// Base point on Polytope
static PT_vector_T PD_direction;			// Unit vector to set shift direction
static PT_vector_T PD_objectiveUnitVector;	// c/=||c||
static PT_vector_T PD_relaxationVector;
//========================== INput/Output ====================================
static string PD_lppFile; /* LPP file in the following format:
------------ begin of file -------------
m n
A_11 A_12 ... A_1n b_1
A_21 A_22 ... A_2n b_2
...
A_m1 A_m2 ... A_mn b_m
c_1 c_2 ... c_n
------------ end of file----------------*/

static string PD_solutionFile; /* LPP solution file in the following format:
------------ begin of file -------------
n
x_1 x_2 ... x_n
------------ end of file----------------*/

//=========================== Debug Data =====================================
#ifdef PP_DEBUG
static bool PD_pointTouch[PP_M];
#endif // PP_DEBUG