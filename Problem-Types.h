/*==============================================================================
Project: LiFe
Theme: Apex Method
Module: Problem-Types.h (BSF Types)
Prefix: PT
Author: Leonid B. Sokolinsky
This source code has been produced with using BSF-skeleton
==============================================================================*/			
#pragma once
#include "Problem-Include.h"		// Problem "Include" Files
#include "Problem-Parameters.h"		// Problem Parameters 
//=========================== Problem Types =========================
typedef double PT_float_T;
typedef PT_float_T PT_vector_T[PP_N];			// Vector of dimension n
typedef PT_float_T	PT_matrix_T[PP_MM][PP_N + 1];
typedef PT_float_T	PT_column_T[PP_MM];