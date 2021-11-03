/*==============================================================================
Project: LiFe
Theme: Apex Method
Module: Problem-bsfTypes.h (Predefined BSF Problem Types)
Prefix: PT_bsf
Author: Leonid B. Sokolinsky
This source code is a part of BSF Skeleton
==============================================================================*/
#pragma once
#include "Problem-Types.h"		// Problem Types 
//=========================== BSF Types =========================
struct PT_bsf_parameter_T {		// Type of Parameter for workers (current approximation)
	PT_vector_T x;	// Current approximation
	int i;			// Index of new inequality
	PT_vector_T a;	// Left part of new inequality
	PT_float_T b;	// Right part of new inequality
};

struct PT_bsf_mapElem_T {		// Type of map-list elements
	PT_float_T* a;				// Pointer to row: a[0],...,a[n-1]
	PT_float_T* b;				// Pointer to constant term b: a_0x_0+...+a_{n-1}x_{n-1} \leq b
};

struct PT_bsf_reduceElem_T {	// Type of reduce-list elements for Job 0 (default)	
	PT_vector_T projection;	// Point projection onto hyperplane
};

struct PT_bsf_reduceElem_T_1 {	// Type of reduce-list elements for Job 1
	bool pointIn;
};

struct PT_bsf_reduceElem_T_2 {	// Type of reduce-list elements for Job 2
	bool pointIn;
};

struct PT_bsf_reduceElem_T_3 {	// Type of reduce-list elements for Job 3
// Not used
};