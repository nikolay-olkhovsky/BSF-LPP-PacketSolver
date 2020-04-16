/*==============================================================================
Project: LiFe
Theme: Apex Method (Predictor + Corrector)
Module: Problem-Types.h (BSF Types)
Prefix: PT
Author(s): Leonid B. Sokolinsky, Irina M. Sokolinskaya
This source code has been produced with using BSF-skeleton
==============================================================================*/			
#pragma once
#include "Problem-Include.h"		// Problem "Include" Files
#include "Problem-Parameters.h"		// Problem Parameters 
//=========================== Problem Types =========================
typedef double PT_vector_T[PP_N];			// Vector of dimensional n
typedef double PT_matrix_T[PP_M][PP_N];		// Matrix m x n
typedef double PT_column_T[PP_M];			// Column of m elements