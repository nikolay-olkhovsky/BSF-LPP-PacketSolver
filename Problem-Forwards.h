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
void		ApexPoint(PT_vector_T innerPont, PT_vector_T apexPoint);
bool		Conversion();
void		DetermineDirection(PT_bsf_parameter_T* parameter, bool* exit, bool* repeat);
PT_float_T	Distance(PT_vector_T x, PT_vector_T y);
bool		GetDirection(PT_vector_T startPoint, PT_vector_T endPoint, PT_vector_T unitVector);
bool		OpenDataFiles();
bool		CloseDataFiles();
void		WriteTrace(PT_vector_T x);
bool		LoadMatrixFormat();
void		MakeObjVector(PT_vector_T c, PT_vector_T objVector);
void		MovingOnSurface(PT_vector_T ptr_unitVectorToSurface, PT_vector_T basePoint, PT_vector_T x, bool* goOn);
PT_float_T	ObjF(PT_vector_T x);
void		UnitObjVector(PT_vector_T objUnitVector);
//
bool		PointInHalfspace(PT_vector_T point, PT_vector_T a, PT_float_T b);
bool		PointInHalfspace_s(PT_vector_T point, PT_vector_T a, PT_float_T b);
bool		PointInPolytope_s(PT_vector_T x);
void		ProblemOutput(double elapsedTime);
double		ProblemScale();
bool		SavePoint(PT_vector_T x, const char* filename, double elapsedTime);
void		Shift(PT_vector_T basePoint, PT_vector_T direction, PT_float_T siftLength, PT_vector_T endPoint);
void		ShrinkUnitVector(PT_vector_T objUnitVector, int shrinkBound);
void		SkipComments(FILE* stream);
void		SortObjVarI();
void		Vector_Addition(PT_vector_T x, PT_vector_T y, PT_vector_T z);
void		Vector_Copy(PT_vector_T fromPoint, PT_vector_T toPoint);
void		Vector_DivideByNumber(PT_vector_T x, double r, PT_vector_T y);
void		Vector_DivideEquals(PT_vector_T x, double r);
PT_float_T	Vector_DotProduct(PT_vector_T x, PT_vector_T y);
void		Vector_EpsZero(PT_vector_T x);
PT_float_T	Vector_Norm(PT_vector_T x);
PT_float_T	Vector_NormSquare(PT_vector_T x);
void		Vector_MinusEquals(PT_vector_T equalPoint, PT_vector_T minusVector);
void		Vector_MultiplyByNumber(PT_vector_T x, double r, PT_vector_T y);
void		Vector_MultiplyEquals(PT_vector_T x, double r);
void		Vector_PlusEquals(PT_vector_T equalVector, PT_vector_T plusVector);
void		Vector_ProjectOnHalfspace(PT_vector_T point, PT_vector_T a, PT_float_T b, PT_vector_T projection, int* exitCode);
void		Vector_ResetToZero(PT_vector_T x);
void		Vector_Round(PT_vector_T x);
void		Vector_Subtraction(PT_vector_T x, PT_vector_T y, PT_vector_T z);
void		Vector_Unit(PT_vector_T vector);
//
//====================== Macros ================================
#define PF_MIN(x,y) (x<y?x:y)
#define PF_MAX(x,y) (x>y?x:y)
#define PF_SIGN(x) (x>0?1:-1)