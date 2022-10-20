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
static double	PD_objF_u;
static double	PD_objF_w;
static int		PD_problemCounter;			// Index of the problem in the package
static int		PD_problemsNumber;			// Number of problems in the package
static PT_float_T PD_problemTrace[PP_TRACE_LIMIT][PP_N];	// Trace points for a problem in the package
static int		PD_traceIndex;				// Index of current trace point
//========================== Variables from functions ====================================
static PT_vector_T relaxationVector;		// used in PC_bsf_ProcessResults
static PT_float_T fvA[PP_MM];				// used in Conversion
static bool Flag[PP_N];						// used in Conversion
static int fvEqI;							// used in Conversion
static bool single;							// used in Conversion
static bool PD_delete[PP_MM];				// used in Conversion
//========================== Problem structures ====================================
static PT_matrix_T PD_A;					// Matrix of coefficients of inequalities 
static PT_column_T PD_b;					// Column of the constant terms of the system Ax <= PD_b
static PT_vector_T PD_c;					// Objective Function Coefficients
static PT_vector_T PD_apexPoint;			// Apex point
static PT_vector_T PD_x0;					// Starting point (defauil is 0)
static PT_vector_T PD_u;					// Base point on Polytope
static PT_vector_T PD_direction;			// Unit vector to set shift direction
static PT_vector_T PD_hi;					// Higher bound
static PT_vector_T PD_lo;					// Lower bound
static PT_vector_T PD_unitObjVector;		// = c/||c||
static PT_vector_T PD_objVector;			// = PD_unitObjectiveVector * PP_OBJECTIVE_VECTOR_LENGTH
static int PD_objI[PP_N];					// Index of objective variables in absolute descending order
//========================== Input/Output ====================================
static string PD_problemName;
//
//
// 
//========================== Matrix format ====================================
// nor - number of rows; noc - number of columns; non - number of non-zero values
static FILE* PD_stream_A;
static string PD_MTX_File_A; /* File of matrix A (only non-zero elements):
------------ begin of file -------------
nor noc non // nor=m (number of inequalities); noc=n (dimension); non - number of non-zero values
i_1 j_1 A_{{i_1}{j_1}} // i_1 - index of row; j_1 - index of column
...
i_k j_k A_{{i_k}{j_k}}
------------ end of file----------------*/
static FILE* PD_stream_b;
static string PD_MTX_File_b; /* File of column b:
------------ begin of file -------------
nor noc // nor=m (number of inequalities); noc=1
b_1
...
b_{nor}
------------ end of file----------------*/
static FILE* PD_stream_lo;
static string PD_MTX_File_lo; /* File of variable lower bounds:
------------ begin of file -------------
nor noc // nor=n (dimension); noc=1
lo_1
...
lo_{nor}
------------ end of file----------------*/
static FILE* PD_stream_hi;
static string PD_MTX_File_hi; /* File of variable higher bounds:
------------ begin of file -------------
nor noc // nor=n (dimension); noc=1
lo_1
...
lo_{nor}
------------ end of file----------------*/
static FILE* PD_stream_c;
static string PD_MTX_File_c; /* File of variable higher bounds:
------------ begin of file -------------
nor noc // nor=n (dimension); noc=1
c_1
...
c_{nor}
------------ end of file----------------*/
static FILE* PD_stream_x0;
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
static FILE* PD_stream_tr;
static string PD_MTX_File_tr;