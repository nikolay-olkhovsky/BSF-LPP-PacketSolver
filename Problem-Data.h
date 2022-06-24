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
static int		PD_n;						// Space dimension
static int		PD_numShiftsSameLength;		// Number of shifts with the same length
static int		PD_numDetDir;				// Number of sequential states "Determine Direction"
static bool		PD_newInequalities;
static double	PD_ObjectiveVectorLength;	// Length of objective vector
static bool		PD_pointIn;
static double	PD_shiftLength;				// Shift length
static int		PD_state;
//========================== Problem structures ====================================
static PT_matrix_T PD_A;					// Matrix of coefficients of inequalities 
static PT_column_T PD_b;					// Column of the constant terms of the system Ax <= PD_b
static PT_vector_T PD_c;					// Objective Function Coefficients
static PT_vector_T PD_basePoint;			// Base point on Polytope
static PT_vector_T PD_direction;			// Unit vector to set shift direction
static PT_vector_T PD_hi;					// Higher bound
static PT_vector_T PD_lo;					// Lower bound
static PT_vector_T PD_objectiveUnitVector;	// = c/||c||
static PT_vector_T PD_objectiveVector;		// = PD_objectiveUnitVector * PP_OBJECTIVE_VECTOR_LENGTH
static int PD_objI[PP_N];					// Index of objective variables in absolute descending order
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

//
//

//========================== Matrix format ====================================
// nor - number of rows; noc - number of columns; non - number of non-zero values
static string PD_MTX_File_A; /* File of matrix A (only non-zero elements):
------------ begin of file -------------
nor noc non // nor=m (number of inequalities); noc=n (dimension); non - number of non-zero values
i_1 j_1 A_{{i_1}{j_1}} // i_1 - index of row; j_1 - index of column
...
i_k j_k A_{{i_k}{j_k}}
------------ end of file----------------*/
static string PD_MTX_File_b; /* File of column b:
------------ begin of file -------------
nor noc // nor=m (number of inequalities); noc=1
b_1
...
b_{nor}
------------ end of file----------------*/
static string PD_MTX_File_lo; /* File of variable lower bounds:
------------ begin of file -------------
nor noc // nor=n (dimension); noc=1
lo_1
...
lo_{nor}
------------ end of file----------------*/
static string PD_MTX_File_hi; /* File of variable higher bounds:
------------ begin of file -------------
nor noc // nor=n (dimension); noc=1
lo_1
...
lo_{nor}
------------ end of file----------------*/
static string PD_MTX_File_c; /* File of variable higher bounds:
------------ begin of file -------------
nor noc // nor=n (dimension); noc=1
c_1
...
c_{nor}
------------ end of file----------------*/
static string PD_MTX_File_x0;/* File of starting point in the following format:
------------ begin of file -------------
nor noc // nor=n (dimension); noc=1
x_1 x_2 ... x_n
------------ end of file----------------*/
static string PD_MTX_File_so;/* Solution file in the following format:
------------ begin of file -------------
nor noc // nor=n (dimension); noc=1
x_1 x_2 ... x_n
------------ end of file----------------*/