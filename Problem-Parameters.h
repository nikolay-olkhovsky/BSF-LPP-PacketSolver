/*==============================================================================
Project: LiFe
Theme: Apex Method (Predictor + Corrector)
Module: Problem-Parameters.h (Problem Parameters)
Prefix: PP
Author(s): Leonid B. Sokolinsky, Irina M. Sokolinskaya
This source code has been produced with using BSF-skeleton
==============================================================================*/
#pragma once
//=========================== Problem Parameters =========================
#define PP_N 8				// Dimension of Space
#define PP_M (2*PP_N+2)		// Number of inequations
#define PP_SF	200			// Scale factor

#define PP_DIST_TO_APEX (PP_N*PP_SF*100)// Distance to Apex Point 
#define PP_EPS_RELAX	1E-9			// Precision for relaxation
#define PP_EPS_SHIFT	1E-6
#define PP_EPS_ZERO		1E-5
#define PP_MAX_ITER_COUNT	1000000		// Maximal count of iterations
#define PP_OBJECTIVE_VECTOR_LENGTH 1	// Length of Objective Vector
#define PP_START_SHIFT_LENGTH 10		// Start length of shift vector

//-------------------------- Outpoot Parameters ---------------------------
#define PP_OUTPUT_LIMIT	8	// Number of Elements to output
#define PP_MATRIX_OUTPUT	// To output Matrix
//-------------------------- Job Cases -----------------------
#define PP_CASE_APEX_POJECTION			0 // 0. Apex Pseudo-pojection
#define PP_CASE_MOVEMENT_ON_POLYTOPE	1 // 1. Movement on Polytope
#define PP_CASE_POINT_PSEUDOPROJECTION	2 // 2. Point Pseudo-projection
//-------------------------- Compilation Modes -----------------------
//#define PP_DEBUG
