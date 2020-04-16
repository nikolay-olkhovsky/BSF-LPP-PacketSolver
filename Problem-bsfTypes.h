/*==============================================================================
Project: Bulk Synchronous Farm (BSF)
Theme: BSF Skeleton
Module: Problem-bsfTypes.h (Predefined BSF Problem Types)
Prefix: PT_bsf
Author: Leonid B. Sokolinsky
This source code is a part of BSF Skeleton
==============================================================================*/
#pragma once
#include "Problem-Types.h"		// Problem Types 
//=========================== BSF Types =========================
struct PT_bsf_parameter_T {				// Parameter for workers (current approximation)
	PT_vector_T x;		// Current approximation
};

struct PT_bsf_mapElem_T {	// Element of the map list
	int inequalityNo;
};

struct PT_bsf_reduceElem_T {	// Element of reduce list	
	PT_vector_T projection;	// Point projection onto hyperplane
};

struct PT_bsf_reduceElem_T_1 {
	bool pointIn;
};

struct PT_bsf_reduceElem_T_2 {
	PT_vector_T projection;	// Point projection onto hyperplane
};

struct PT_bsf_reduceElem_T_3 {
// Not used
};