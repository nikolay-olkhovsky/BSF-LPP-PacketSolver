/*==============================================================================
Project: LiFe
Theme: Apex Method (Predictor + Corrector)
Module: Problem-Forwards.h (Problem Function Forwards)
Author(s): Leonid B. Sokolinsky, Irina M. Sokolinskaya
This source code has been produced with using BSF-skeleton
==============================================================================*/
#include "Problem-bsfTypes.h"
#include "Problem-Types.h"
//====================== Problem Functions ===========================
static void		GetDirection(PT_vector_T startPoint, PT_vector_T endPoint, PT_vector_T unitVector);
static double	ObjectiveF(PT_vector_T x);
static bool		PointIn(PT_vector_T point, int halfSpaceNo);
static void		Projection(PT_vector_T point, int hyperplaneNo, double normSquare, PT_vector_T projection);
static void		Shift(PT_vector_T basePoint, PT_vector_T direction, double siftLength, PT_vector_T endPoint);
static void		Vector_Addition(PT_vector_T x, PT_vector_T y, PT_vector_T z);
static void		Vector_Copy(PT_vector_T fromPoint, PT_vector_T toPoint);
static void		Vector_DivideByNumber(PT_vector_T x, double r, PT_vector_T y);
static double	Vector_DotProductSquare(PT_vector_T x, PT_vector_T y);
static void		Vector_From_x_to_y(PT_vector_T x, PT_vector_T y, double vectorLength, PT_vector_T vectorFrom_x_to_y);
static double	Vector_NormSquare(PT_vector_T x);
static void		Vector_MinusEquals(PT_vector_T equalPoint, PT_vector_T minusVector);
static void		Vector_MultiplyByNumber(PT_vector_T x, double r, PT_vector_T y);
static void		Vector_PlusEquals(PT_vector_T equalVector, PT_vector_T plusVector);
static void		Vector_Round(PT_vector_T x);
static void		Vector_Relaxation(PT_vector_T sumOfProjections, int numberOfProjections, PT_vector_T relaxationVector);
static void		Vector_Subtraction(PT_vector_T x, PT_vector_T y, PT_vector_T z);
static void		Vector_Unit(PT_vector_T vector, PT_vector_T unitVector);
//====================== Macros ================================
#define PF_MIN(x,y) (x<y?x:y)
#define PF_MAX(x,y) (x>y?x:y)