/*==============================================================================
Project: LiFe
Theme: Apex Method
Module: Problem-Forwards.h (Problem Function Forwards)
Author: Leonid B. Sokolinsky
This source code has been produced with using BSF-skeleton
==============================================================================*/
#include "Problem-bsfTypes.h"
#include "Problem-Types.h"
//====================== Problem Functions ===========================
PT_float_T	Distance(PT_vector_T x, PT_vector_T y);
bool		GetDirection(PT_vector_T startPoint, PT_vector_T endPoint, PT_vector_T unitVector);
PT_float_T	ObjectiveF(PT_vector_T x);
void		ObjectiveUnitVector(PT_vector_T objectiveUnitVector);
bool		PointInHalfspace(PT_vector_T point, PT_vector_T a, PT_float_T b);
bool		PointInHalfspace_s(PT_vector_T point, PT_vector_T a, PT_float_T b);
bool		PointInPolytope_s(PT_vector_T x);
void		ProblemOutput(double elapsedTime);
bool		SaveSolution(PT_vector_T x, const char* filename);
void		Shift(PT_vector_T basePoint, PT_vector_T direction, double siftLength, PT_vector_T endPoint);


void		Vector_Addition(PT_vector_T x, PT_vector_T y, PT_vector_T z);
void		Vector_Copy(PT_vector_T fromPoint, PT_vector_T toPoint);
void		Vector_DivideByNumber(PT_vector_T x, double r, PT_vector_T y);
PT_float_T	Vector_DotProductSquare(PT_vector_T x, PT_vector_T y);
PT_float_T	Vector_NormSquare(PT_vector_T x);
void		Vector_MinusEquals(PT_vector_T equalPoint, PT_vector_T minusVector);
void		Vector_MultiplyByNumber(PT_vector_T x, double r, PT_vector_T y);
void		Vector_PlusEquals(PT_vector_T equalVector, PT_vector_T plusVector);
bool		Vector_ProjectOnHalfspace(PT_vector_T point, PT_vector_T a, PT_float_T b, PT_vector_T projection);
void		Vector_Round(PT_vector_T x);
void		Vector_Relaxation(PT_vector_T sumOfProjections, int numberOfProjections, PT_vector_T PD_relaxationVector);
void		Vector_Subtraction(PT_vector_T x, PT_vector_T y, PT_vector_T z);
void		Vector_Unit(PT_vector_T vector, PT_vector_T unitVector);

//====================== Macros ================================
#define PF_MIN(x,y) (x<y?x:y)
#define PF_MAX(x,y) (x>y?x:y)