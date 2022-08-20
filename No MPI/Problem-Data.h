/*==============================================================================
Project: LiFe
Theme: Apex Method (No MPI)
Module: Problem-Data.h (Problem Data)
Prefix: PP
Author: Leonid B. Sokolinsky
This source code has been produced with using BSF-skeleton
==============================================================================*/
#include "Problem-Types.h"		// Problem Parameters 
using namespace std;
//========================== Problem variables ====================================
static int		PD_m;						// Current number of inequalities
static int		PD_n;						// Curtrent space dimension
static int		PD_numShiftsSameLength;		// Number of shifts with the same length
static int		PD_numDetDir;				// Number of sequential states "Determine Direction"
static bool		PD_pointIn;					// Point is inside polytope
static double	PD_shiftLength;				// Shift length
static int		PD_state;					// State of Job Dispatcher (see PC_bsf_JobDispatcher)
static int		PD_indexToBlock;			// Index of variable to block
static int		PD_firstZcvI;				// Index of first zero-costs variable
static int		PD_firstLcvI;				// Index of first low-costs variable
static double	PD_minCostPercentile;		// Minimum cost percentile
static bool		PD_utilizeLowCostVariables;
static double	PD_baseObjValue;
//========================== Problem structures ====================================
static PT_matrix_T PD_A;					// Matrix of coefficients of inequalities
static PT_column_T PD_b;					// Column of the constant terms of the system Ax <= PD_b
static PT_vector_T PD_c;					// Objective Function Coefficients
static PT_vector_T PD_apexPoint;			// Apex point
static PT_vector_T PD_basePoint;			// Base point on Polytope
static PT_vector_T PD_direction;			// Unit vector to set shift direction
static PT_vector_T PD_hi;					// Higher bound
static PT_vector_T PD_lo;					// Lower bound
static PT_vector_T PD_unitObjVector;		// = c/||c||
static PT_vector_T PD_objVector;			// = PD_unitObjectiveVector * PP_OBJECTIVE_VECTOR_LENGTH
static int PD_objI[PP_N];					// Index of objective variables in absolute descending order
//========================== Input/Output ====================================
static string PD_problemName;
static string PD_solutionFile;
static string PD_MTX_File_tr;
static FILE*  PD_traceStream;
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
x_1
...
x_{nor}
------------ end of file----------------*/
static string PD_MTX_File_so;/* Solution file in the following format:
------------ begin of file -------------
nor noc // nor=n (dimension); noc=1
x_1
...
x_{nor}
------------ end of file----------------*/