/*==============================================================================
Project: LiFe
Theme: Apex Method (Predictor + Corrector)
Module: Problem-bsfCode.cpp (Implementation of the Problem)
Prefix: PI
Author(s): Leonid B. Sokolinsky, Irina M. Sokolinskaya
This source code has been produced with using BSF-skeleton
==============================================================================*/
#include "Problem-Data.h"			// Problem Types 
#include "Problem-Forwards.h"		// Problem Function Forwards
#include "Problem-bsfParameters.h"	// BSF-skeleton parameters
#include "BSF-SkeletonVariables.h"	// Skeleton Variables
using namespace std;

void PC_bsf_SetInitParameter(PT_bsf_parameter_T* parameter) {
	for (int j = 0; j < PP_N; j++) // Generating initial approximation
		parameter->x[j] = PD_apex[j];
};

void PC_bsf_Init(bool* success) {

	// Generating Objective Function Coefficients
	for (int j = 0; j < PP_N; j++)
		PD_c[j] = ((double)PP_N - j)*10;

	// Generating a point inside the polytope
	for (int j = 0; j < PP_N; j++)
		PD_basePoint[j] = (double)PP_SF / 2;

	// Generating Objective Unit Vector
	double c_norm = sqrt(Vector_NormSquare(PD_c));
	Vector_DivideByNumber(PD_c,c_norm ,PD_objectiveUnitVector);

	// Generating Coordinates of Apex Point
	Vector_MultiplyByNumber(PD_objectiveUnitVector, PP_DIST_TO_APEX, PD_apex);
	Vector_PlusEquals(PD_apex, PD_basePoint);

	if (PD_apex[PP_N - 1] < PP_SF * 1.1) {
		if (BSF_sv_mpiRank == 0)
			cout << "PC_bsf_Init.Error: PD_apex[PP_N - 1] = " << PD_apex[PP_N - 1] << " < PP_SF * 1.1 = " << PP_SF * 1.1 << endl;
		*success = false;
		return;
	}

	// Generating A & b
	for (int i = 0; i < PP_N; i++) {
		for (int j = 0; j < PP_N; j++) {
			PD_A[i][j] = 0;
			/* debug */// { if (i >= PP_M) exit(13); if (j >= PP_N) exit(14); }
		}
		PD_A[i][i] = 1;
		PD_b[i] = PP_SF;
	};
	for (int j = 0; j < PP_N; j++)
		PD_A[PP_N][j] = 1;
	PD_b[PP_N] = (double)PP_SF * ((double)PP_N - 1) + (double)PP_SF / 2;
	for (int j = 0; j < PP_N; j++)
		PD_A[PP_N + 1][j] = -1;
	PD_b[PP_N + 1] = -PP_SF / 2;
	for (int i = PP_N + 2; i < PP_M; i++) {
		for (int j = 0; j < PP_N; j++)
			PD_A[i][j] = 0;
		PD_A[i][i - PP_N - 2] = -1;
		/* debug */// { if (i >= PP_M) exit(13); if (i - PP_N - 2 >= PP_N) exit(14); }
		PD_b[i] = 0;
	};

	// Calculating a_i norm squares
	for (int i = 0; i < PP_M; i++) {
		PD_normSquare_a[i] = 0;
		for (int j = 0; j < PP_N; j++)
			PD_normSquare_a[i] += PD_A[i][j] * PD_A[i][j];
	};
	
	/* debug *//* if (BSF_sv_mpiRank == 0) {
		cout << "----------PC_bsf_Init-------------" << endl;
		for (int i = 0; i < PP_M; i++)
			cout << "InequalityNo = " << i << "\tNorm Square = " << PD_normSquare_a[i] << endl;
		//system("pause");
	};/* end debug */
}

void PC_bsf_SetListSize(int* listSize) {
	*listSize = PP_M;
}

void PC_bsf_SetMapListElem(PT_bsf_mapElem_T* elem, int i) {
	elem->inequalityNo = i;
}

// 0. Apex Pseudo-pojection
void PC_bsf_MapF(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T* reduceElem, int* success // 1 - reduceElem was produced successfully; 0 - otherwise
){
	double factor;

	factor = (PD_b[mapElem->inequalityNo] - Vector_DotProductSquare(BSF_sv_parameter.x, PD_A[mapElem->inequalityNo])) / PD_normSquare_a[mapElem->inequalityNo];

	if (factor > 0)
		*success = false;
	else
		for (int j = 0; j < PP_N; j++) {
			reduceElem->projection[j] = factor * PD_A[mapElem->inequalityNo][j];
			/* debug */// { if (mapElem->inequalityNo >= PP_M || mapElem->inequalityNo < 0) exit(13); }
		}

	/* debug *//* if (BSF_sv_mpiRank == 0) {
		cout << "Hyperplane No = " << mapElem->inequalityNo << "\tProjection: ";
		if (factor >= 0)
			cout << "\tsuccess = false" << endl;
		else {
			cout << "\tProjection Vector:";
			for (int j = 0; j < PP_N; j++)
				cout << setw(20) << reduceElem->projection[j];
		};
		cout << endl;
		system("pause");
	};/* end debug */
};

void PC_bsf_MapF_3(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_3* reduceElem, int* success) {
	// not used
};

// 1. Movement on Polytope
void PC_bsf_MapF_1(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_1* reduceElem,int* success) {
	if (PointIn(BSF_sv_parameter.x, mapElem->inequalityNo))
		reduceElem->pointIn = true;
	else
		reduceElem->pointIn = false;
}

// 2. Point Pseudo-projection
void PC_bsf_MapF_2(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_2* reduceElem,int* success) {
	double factor;

	/* debug *//*  if (BSF_sv_mpiRank == 0) {
		cout << "PC_bsf_MapF_2: parameter->x:";
		for (int j = 0; j < PP_N; j++)
			cout << setw(20) << BSF_sv_parameter.x[j];
		cout << endl;
		//system("pause");
	} /* end debug */

	factor = (PD_b[mapElem->inequalityNo] - Vector_DotProductSquare(BSF_sv_parameter.x, PD_A[mapElem->inequalityNo])) / PD_normSquare_a[mapElem->inequalityNo];

	if (factor > 0)
		*success = false;
	else
		for (int j = 0; j < PP_N; j++) {
			reduceElem->projection[j] = factor * PD_A[mapElem->inequalityNo][j];
		}

	/* debug *//* if (BSF_sv_mpiRank == 0) {
		cout << BSF_sv_mpiRank << "------------ > PC_bsf_ProcessResults_2------------ >" << endl;
		cout << "Hyperplane No = " << mapElem->inequalityNo << "\tProjection: ";
		if (factor >= 0)
			cout << "\tsuccess = false" << endl;
		else {
			cout << "\tProjection Vector: ";
			for (int j = 0; j < PP_N; j++)
				cout << setw(12) << reduceElem->projection[j];
		};
		cout << endl;
		//system("pause");
	};/* end debug */
}

// 0. Apex Pseudo-pojection
void PC_bsf_ReduceF(PT_bsf_reduceElem_T* x, PT_bsf_reduceElem_T* y, PT_bsf_reduceElem_T* z) { // z = x + y
	Vector_Addition(x->projection, y->projection, z->projection); 
}

void PC_bsf_ReduceF_3(PT_bsf_reduceElem_T_3* x, PT_bsf_reduceElem_T_3* y, PT_bsf_reduceElem_T_3* z) {
	// not used
}

//0. Apex pseudo-projection onto Politope
void PC_bsf_ProcessResults(
	PT_bsf_reduceElem_T* reduceResult,
	int reduceCounter, // Number of successfully produced Elrments of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* nextJob,
	bool* exit // "true" if Stopping Criterion is satisfied, and "false" otherwise
) {
#ifdef PP_DEBUG
	if (PD_caseApexPojection) {
		cout << "------------> 0: APEX POJECTION  -------------->" << endl;
		PD_caseApexPojection = false;
	}
	/* debug */// if (BSF_sv_iterCounter > 20000) *exit = true; return;
#endif // PP_DEBUG

	/* debug *//* for (int j = 0; j < PP_N; j++)
		parameter->x[j] = 200;
	parameter->x[PP_N - 1] = 0;

	cout << "parameter->x:";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++)
		cout << setw(20) << parameter->x[j];
	if (PP_OUTPUT_LIMIT < PP_N) cout << "	..." << setw(20) << parameter->x[PP_N - 1];
	cout << endl;

	PT_vector_T objectiveVector;
	Vector_MultiplyByNumber(PD_objectiveUnitVector, PP_OBJECTIVE_VECTOR_LENGTH, objectiveVector);
	Vector_Copy(parameter->x, PD_basePoint);
	Vector_PlusEquals(parameter->x, objectiveVector);

	*jobCase = PP_CASE_POINT_PSEUDOPROJECTION;
	PD_casePointProjection = true;

	cout << "x for ps-pr.:";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++)
		cout << setw(20) << parameter->x[j];
	if (PP_OUTPUT_LIMIT < PP_N) cout << "	..." << setw(20) << parameter->x[PP_N - 1];
	cout << endl;
	//system("pause");
	return;
	/* end debug */

#ifdef PP_MAX_ITER_COUNT
	if (BSF_sv_iterCounter > PP_MAX_ITER_COUNT) {
		cout << "-------------> PC_bsf_ProcessResults: Acceptable maximum number of iterations is exceeded: PP_MAX_ITER_COUNT = "
			<< PP_MAX_ITER_COUNT << endl;
		*exit = true;
		return;
	};
#endif // PP_MAX_ITER_COUNT

	PT_vector_T relaxationVector;
	Vector_Relaxation(reduceResult->projection, reduceCounter, relaxationVector);
	Vector_PlusEquals(parameter->x, relaxationVector);

	if (Vector_NormSquare(relaxationVector) < PP_EPS_RELAX * PP_EPS_RELAX)
	{
		Vector_Round(parameter->x);

#ifdef PP_DEBUG
		PD_casePointProjection = true;
		cout << "parameter->x:";
		for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++)
			cout << setw(20) << parameter->x[j];
		if (PP_OUTPUT_LIMIT < PP_N) cout << "	..." << setw(20) << parameter->x[PP_N - 1];
		cout << endl;
		//system("pause");
#endif // PP_DEBUG

		/* debug */// *exit = true; return;

		PT_vector_T objectiveVector;
		Vector_MultiplyByNumber(PD_objectiveUnitVector, PP_OBJECTIVE_VECTOR_LENGTH, objectiveVector);
		Vector_Copy(parameter->x, PD_basePoint);
		Vector_PlusEquals(parameter->x, objectiveVector);

		*nextJob = PP_CASE_POINT_PSEUDOPROJECTION;

#ifdef PP_DEBUG
		PD_casePointProjection = true;
		cout << "x for ps-pr.:";
		for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++)
			cout << setw(20) << parameter->x[j];
		if (PP_OUTPUT_LIMIT < PP_N) cout << "	..." << setw(20) << parameter->x[PP_N - 1];
		cout << endl;
		//system("pause");
#endif // PP_DEBUG
	};
}

// 1. Movement on Polytope
void PC_bsf_ReduceF_1(PT_bsf_reduceElem_T_1* x, PT_bsf_reduceElem_T_1* y, PT_bsf_reduceElem_T_1* z) {
	z->pointIn = x->pointIn && y->pointIn;
}

// 2. Point Pseudo-projection
void PC_bsf_ReduceF_2(PT_bsf_reduceElem_T_2* x, PT_bsf_reduceElem_T_2* y, PT_bsf_reduceElem_T_2* z) {
	Vector_Addition(x->projection, y->projection, z->projection);
}

void PC_bsf_ProcessResults_3(
	PT_bsf_reduceElem_T_3* reduceResult,
	int reduceCounter, // Number of successfully produced Elrments of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* nextJob,
	bool* exit // "true" if Stopping Criterion is satisfied, and "false" otherwise
) {
	// not used
}

// 1. Movement on Polytope  ========================================================
void PC_bsf_ProcessResults_1(
	PT_bsf_reduceElem_T_1* reduceResult,
	int reduceCounter, // Number of successfully produced Elrments of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* nextJob,
	bool* exit // "true" if Stopping Criterion is satisfied, and "false" otherwise
) {
#ifdef PP_DEBUG
	if (PD_caseMovementOnPlytope) {
		cout << "=================> 1: MOVEMENT ON PLYTOPE  =================>" << endl;
		PD_caseMovementOnPlytope = false;
		/* debug */// if (BSF_sv_iterCounter > 20000) *exit = true; return;
	}
#endif // PP_DEBUG
#ifdef PP_MAX_ITER_COUNT
	if (BSF_sv_iterCounter > PP_MAX_ITER_COUNT) {
		cout << "-------------> Acceptable maximum number of iterations is exceeded: PP_MAX_ITER_COUNT = " << PP_MAX_ITER_COUNT << endl;
		*exit = true;
		return;
	}
#endif // PP_MAX_ITER_COUNT

	if (reduceResult->pointIn) {
		Vector_Copy(parameter->x, PD_basePoint);
		Shift(PD_basePoint, PD_direction, PD_siftLength, parameter->x);
		return;
	}

	if (PD_siftLength >= PP_EPS_SHIFT) {
		PD_siftLength /= 2;
		Shift(PD_basePoint, PD_direction, PD_siftLength, parameter->x);
		return;
	}

	PT_vector_T objectiveVector;
	Vector_MultiplyByNumber(PD_objectiveUnitVector, PP_OBJECTIVE_VECTOR_LENGTH, objectiveVector);
	Vector_Addition(PD_basePoint, objectiveVector, parameter->x);
#ifdef PP_DEBUG
	PD_casePointProjection = true;
	cout << "Point for pp:";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++)
		cout << setw(20) << parameter->x[j];
	if (PP_OUTPUT_LIMIT < PP_N) cout << "	..." << setw(20) << parameter->x[PP_N - 1];
	cout << endl;
#endif // PP_DEBUG

	* nextJob = PP_CASE_POINT_PSEUDOPROJECTION;
	Vector_Round(PD_basePoint);

#ifdef PP_DEBUG
	PD_casePointProjection = true;
	cout << "PD_basePoint:";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++)
		cout << setw(20) << PD_basePoint[j];
	if(PP_OUTPUT_LIMIT < PP_N) cout << "	..." << setw(20) << PD_basePoint[PP_N - 1];
	cout << endl;

	/* debug */// *exit = true;

#endif // PP_DEBUG
}

// 2. Point Pseudo-projection  ========================================================
void PC_bsf_ProcessResults_2(
	PT_bsf_reduceElem_T_2* reduceResult,
	int reduceCounter, // Number of successfully produced Elrments of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* nextJob,
	bool* exit // "true" if Stopping Criterion is satisfied, and "false" otherwise
) {
#ifdef PP_DEBUG
	if (PD_casePointProjection) {
		cout << "=================> 2: POINT PSEUDO-PROJECTION  =================>" << endl;
		PD_casePointProjection = false;
	}
	/* debug */// if (BSF_sv_iterCounter > 20000) *exit = true; return;
#endif // PP_DEBUG

#ifdef PP_MAX_ITER_COUNT
	if (BSF_sv_iterCounter > PP_MAX_ITER_COUNT) {
		cout << "-------------> Acceptable maximum number of iterations is exceeded: PP_MAX_ITER_COUNT = " << PP_MAX_ITER_COUNT << endl;
		*exit = true;
		return;
	}
#endif // PP_MAX_ITER_COUNT

	if (reduceCounter == 0) {
		cout << "PC_bsf_ProcessResults_2: Error: reduceCounter == 0 " << endl;
		*exit = true;
		return; 
	}


	/* debug *//* for (int j = 0; j < PP_N -1 ; j++) {
		parameter->x[j] = 200;
		PD_direction[j] = 0;
	}
	parameter->x[PP_N - 1] = 0;
	PD_direction[PP_N - 1] = 1;

	cout << "parameter->x:";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++)
		cout << setw(20) << parameter->x[j];
	if (PP_OUTPUT_LIMIT < PP_N) cout << "	..." << setw(20) << parameter->x[PP_N - 1];
	cout << endl;
	cout << "PD_direction:";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++)
		cout << setw(20) << PD_direction[j];
	if (PP_OUTPUT_LIMIT < PP_N) cout << "	..." << setw(20) << PD_direction[PP_N - 1];
	cout << endl;

	PD_siftLength = PP_START_SHIFT_LENGTH;
	Shift(PD_basePoint, PD_direction, PD_siftLength, parameter->x);

	*jobCase = PP_CASE_MOVEMENT_ON_POLYTOPE;
	PD_caseMovementOnPlytope = true;

	//system("pause");
	return;
	/* end debug */

	PT_vector_T relaxationVector;
	Vector_Relaxation(reduceResult->projection, reduceCounter, relaxationVector);
	Vector_PlusEquals(parameter->x, relaxationVector);

	if (Vector_NormSquare(relaxationVector) >= PP_EPS_RELAX * PP_EPS_RELAX) 
		return; 

	Vector_Round(parameter->x);
	GetDirection(PD_basePoint, parameter->x, PD_direction);

	/* debug *//* cout << "parameter->x:";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++)
		cout << setw(20) << parameter->x[j];
	if (PP_OUTPUT_LIMIT < PP_N) cout << "	..." << setw(20) << parameter->x[PP_N - 1];
	cout << endl;
	/* end debug */

	if (ObjectiveF(parameter->x) > ObjectiveF(PD_basePoint) + PP_EPS_ZERO) {
#ifdef PP_DEBUG
		cout << "Objective value in base point = " << ObjectiveF(PD_basePoint) << " < Objective value in new point = " << ObjectiveF(parameter->x) << endl;
#endif // PP_DEBUG

		PD_siftLength = PP_START_SHIFT_LENGTH;
		Shift(PD_basePoint, PD_direction, PD_siftLength, parameter->x);

		*nextJob = PP_CASE_MOVEMENT_ON_POLYTOPE;
		PD_caseMovementOnPlytope = true;
		return;
	}

#ifdef PP_DEBUG
	cout << "Objective value in base point = " << ObjectiveF(PD_basePoint) << " >= Objective value in new point = " << ObjectiveF(parameter->x) << endl;
#endif // PP_DEBUG

	*exit = true; return;
}

void PC_bsf_ParametersOutput(PT_bsf_parameter_T parameter) {
	cout << "=================================================== Quest & Target ====================================================" << endl;
	cout << "Number of Workers: " << BSF_sv_numOfWorkers << endl;
#ifdef PP_BSF_OMP
#ifdef PP_BSF_NUM_THREADS
	cout << "Number of Threads: " << PP_BSF_NUM_THREADS << endl;
#else
	cout << "Number of Threads: " << omp_get_num_procs() << endl;
#endif // PP_BSF_NUM_THREADS
#else
	cout << "OpenMP is turned off!" << endl;
#endif // PP_BSF_OMP
#ifdef PP_BSF_FRAGMENTED_MAP_LIST
	cout << "Map List is Fragmented." << endl;
#else
	cout << "Map List is not Fragmented." << endl;
#endif // PP_BSF_NUM_THREADS
	cout << "Dimension: N = " << PP_N << endl;
	cout << "Number of Constraints: M = " << PP_M << endl;
	cout << "Scale Factor: SF = " << PP_SF << endl;
	cout << "Distance to Apex = " << PP_DIST_TO_APEX << endl;
	cout << "Eps_Relax = " << PP_EPS_RELAX << endl;
	//cout << "Eps_In = " << PP_EPS_IN << endl;
	cout << "Eps_Shift = " << PP_EPS_SHIFT << endl;
	cout << "Length of Objective Vector: " << PP_OBJECTIVE_VECTOR_LENGTH << endl;
	cout << "Start length of shift vector: " << PP_START_SHIFT_LENGTH << endl;
#ifdef PP_MATRIX_OUTPUT
	cout << "------- Matrix A & Column b -------" << endl;
	for (int i = 0; i < PP_M; i++) {
		cout << i << ")";
		for (int j = 0; j < PP_N; j++)
			cout << setw(5) << PD_A[i][j];
		cout << "\t<=" << setw(5) << PD_b[i] << endl;
	};
#endif // PP_MATRIX_OUTPUT
	cout << "Objective Function:"; for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(6) << PD_c[j]; cout << (PP_OUTPUT_LIMIT < PP_N ? "	..." : "") << endl;
	cout << "Interior point:"; 
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(6) << PD_basePoint[j]; 
	if (PP_OUTPUT_LIMIT < PP_N) cout << "	..." << setw(12) << PD_basePoint[PP_N - 1];
	cout << endl;
	cout << "Apex Point:"; 
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(12) << round(PD_apex[j]);
	if (PP_OUTPUT_LIMIT < PP_N) cout << "	..." << setw(12) << round(PD_apex[PP_N - 1]);
	cout << endl;
	cout << "-------------------------------------------" << endl;
};

void PC_bsf_CopyParameter(PT_bsf_parameter_T parameterIn, PT_bsf_parameter_T* parameterOutP) {
	for (int i = 0; i < PP_N; i++)
		parameterOutP->x[i] = parameterIn.x[i];
};

// 0. Apex Pseudo-pojection  ==========================================================
void PC_bsf_IterOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob) {
	cout << "------------------ 0. Apex Pseudo-pojection. Iter # " << BSF_sv_iterCounter << " -------------------" << endl;
	/* debug */// cout << "Elapsed time: " << round(elapsedTime) << endl;
	cout << "Approximat. :"; 
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(20) << parameter.x[j]; 
	if (PP_OUTPUT_LIMIT < PP_N) cout << "	..." << setw(20) << parameter.x[PP_N - 1];
	cout << endl;
};

void PC_bsf_IterOutput_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	// not used
}

// 1. Movement on Polytope
void PC_bsf_IterOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	cout << "------------------ 1. Movement on Polytope. Iter # " << BSF_sv_iterCounter << " ------------------" << endl;
	/* debug */// cout << "Elapsed time: " << round(elapsedTime) << endl;
	cout << "PD_basePoint:";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++)
		cout << setw(20) << PD_basePoint[j];
	if(PP_OUTPUT_LIMIT < PP_N) cout << "	..." << setw(20) << PD_basePoint[PP_N - 1];
	cout << endl;
	cout << "PD_direction:";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++)
		cout << setw(20) << PD_direction[j];
	if(PP_OUTPUT_LIMIT < PP_N) cout << "	..." << setw(20) << PD_direction[PP_N - 1];
	cout << endl;
	cout << "Sift Length = " << PD_siftLength << endl;
};

// 2. Point Pseudo-projection
void PC_bsf_IterOutput_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	cout << "------------------ 2. Point Pseudo-projection. Iter # " << BSF_sv_iterCounter << " ------------------" << endl;
	/* debug */// cout << "Elapsed time: " << round(elapsedTime) << endl;
	cout << "parameter.x :";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++)
		cout << setw(20) << parameter.x[j];
	if(PP_OUTPUT_LIMIT < PP_N) cout << "	..." << setw(20) << parameter.x[PP_N - 1];
	cout << endl;
}

// 0. Apex Pseudo-pojection
void PC_bsf_ProblemOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	cout << "=============================================" << endl;
	cout << "Time: " << t << endl;
	cout << "Iterations: " << BSF_sv_iterCounter << endl;
	cout << "Solution:"; for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(20) << parameter.x[j]; 
	if (PP_OUTPUT_LIMIT < PP_N) cout << "	..." << setw(20) << parameter.x[PP_N - 1];
	cout << endl;
}

void PC_bsf_ProblemOutput_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,double t) {
	// not used
}

// 1. Movement on Polytope
void PC_bsf_ProblemOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,double t) {
	cout << "Time: " << t << endl;
	cout << "Iterations: " << BSF_sv_iterCounter << endl;
	cout << "Solution: ";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(5) << round(PD_basePoint[j]);
	if (PP_OUTPUT_LIMIT < PP_N) cout << "	..." << setw(20) << round(PD_basePoint[PP_N - 1]);
}

// 2. Point Pseudo-projection
void PC_bsf_ProblemOutput_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,double t) {
	cout << "=============================================" << endl;
	cout << "Time: " << t << endl;
	cout << "Iterations: " << BSF_sv_iterCounter << endl;
	cout << "Solution: "; 
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(5) << round(PD_basePoint[j]); 
	if (PP_OUTPUT_LIMIT < PP_N) cout << "	..." << setw(20) << round(PD_basePoint[PP_N - 1]);
}


//----------------------- Assigning Values to BSF-skeleton Variables (Do not modify!) -----------------------
void PC_bsfAssignAddressOffset(int value) { BSF_sv_addressOffset = value; };
void PC_bsfAssignIterCounter(int value) { BSF_sv_iterCounter = value; };
void PC_bsfAssignJobCase(int value) { BSF_sv_jobCase = value; };
void PC_bsfAssignMpiRank(int value) { BSF_sv_mpiRank = value; };
void PC_bsfAssignNumberInSublist(int value) { BSF_sv_numberInSublist = value; };
void PC_bsfAssignNumOfWorkers(int value) { BSF_sv_numOfWorkers = value; };
void PC_bsfAssignParameter(PT_bsf_parameter_T parameter) { PC_bsf_CopyParameter(parameter, &BSF_sv_parameter); }
void PC_bsfAssignSublistLength(int value) { BSF_sv_sublistLength = value; };

//---------------------------------- Problem functions -------------------------
static double Vector_DotProductSquare(PT_vector_T x, PT_vector_T y) {
	double sum = 0;
	for (int j = 0; j < PP_N; j++)
		sum += x[j] * y[j];
	return sum;
}

static void Vector_Relaxation(PT_vector_T sumOfProjections, int numberOfProjections, PT_vector_T relaxationVector) {
	for (int j = 0; j < PP_N; j++) 
		relaxationVector[j] = sumOfProjections[j] / (double)numberOfProjections;
}

static double Vector_NormSquare(PT_vector_T x) {
	double sum = 0;

	for (int j = 0; j < PP_N; j++)
		sum += x[j] * x[j];
	return sum;
}

static void Projection(PT_vector_T point, int hyperplaneNo, double normSquare, PT_vector_T projection) {
	double factor;

	factor = (PD_b[hyperplaneNo] - Vector_DotProductSquare(point, PD_A[hyperplaneNo])) / normSquare;
	for (int j = 0; j < PP_N; j++)
		projection[j] = point[j] + factor * PD_A[hyperplaneNo][j];

	/* debug *//* if (BSF_sv_mpiRank == 0) {
		cout << "---------- Projection -------------" << endl;
		cout << "point: ";
		for (int j = 0; j < PP_N; j++)
			cout << setw(20) << point[j];
		cout << "\thyperplaneNo: " << hyperplaneNo << "\tPD_A[hyperplaneNo]: ";
		for (int j = 0; j < PP_N; j++)
			cout << setw(2) << PD_A[hyperplaneNo][j];
		cout << " <= " << PD_b[hyperplaneNo] << endl;
		cout << "Vector_DotProductSquare(point, PD_A[hyperplaneNo]) = " << Vector_DotProductSquare(point, PD_A[hyperplaneNo]) << endl;
		cout << "normSquare = " << normSquare << endl;
		cout << "factor = " << factor << endl;
	}/* end debug */
}

static bool PointIn(PT_vector_T point, int halfSpaceNo) { // If the point belonges to the Halfspace 
	double sum = 0;

	for (int j = 0; j < PP_N; j++)
		sum += PD_A[halfSpaceNo][j] * point[j];
	//if (sum > PD_b[halfSpaceNo] + PP_EPS_IN)
	if (sum > PD_b[halfSpaceNo])
			return false;
	else
		return true;
}

static void Vector_From_x_to_y(PT_vector_T x, PT_vector_T y, double vectorLength, PT_vector_T vectorFrom_x_to_y) {
	double s;

	s = 0;
	for (int j = 0; j < PP_N; j++) {
		vectorFrom_x_to_y[j] = y[j] - x[j];
		s += vectorFrom_x_to_y[j] * vectorFrom_x_to_y[j];
	};
	s = vectorLength/sqrt(s);
	for (int j = 0; j < PP_N; j++)
		vectorFrom_x_to_y[j] *= s;
}

static void Shift(PT_vector_T basePoint, PT_vector_T direction, double siftLength, PT_vector_T endPoint) {
	for (int j = 0; j < PP_N; j++)
		endPoint[j] = basePoint[j] + direction[j] * siftLength;
}

// Calculating unit vector of direction from startPoint to endPoint
static void GetDirection(PT_vector_T startPoint, PT_vector_T endPoint, PT_vector_T unitVector) { 
	for (int j = 0; j < PP_N; j++)
		unitVector[j] = endPoint[j] - startPoint[j];
	double normOfUnitVector = sqrt(Vector_NormSquare(unitVector));
	for (int j = 0; j < PP_N; j++)
		unitVector[j] /= normOfUnitVector;
	for (int j = 0; j < PP_N; j++)
		if (fabs(unitVector[j]) < PP_EPS_ZERO)
			unitVector[j] = 0;
}

static void Vector_Copy(PT_vector_T fromPoint, PT_vector_T toPoint) { // toPoint = fromPoint
	for (int j = 0; j < PP_N; j++)
		toPoint[j] = fromPoint[j];
}

static void Vector_PlusEquals(PT_vector_T equalVector, PT_vector_T plusVector) { // equalVector += plusVector
	for (int j = 0; j < PP_N; j++)
		equalVector[j] += plusVector[j];
}

static void Vector_MinusEquals(PT_vector_T equalPoint, PT_vector_T minusVector) { // equalPoint += minusVector
	for (int j = 0; j < PP_N; j++)
		equalPoint[j] -= minusVector[j];
};

static void Vector_Addition(PT_vector_T x, PT_vector_T y, PT_vector_T z) {  // z = x + y
	for (int j = 0; j < PP_N; j++)
		z[j] = x[j] + y[j];
}

static void Vector_Subtraction(PT_vector_T x, PT_vector_T y, PT_vector_T z) {  // z = x - y
	for (int j = 0; j < PP_N; j++)
		z[j] = x[j] - y[j];
}

static void Vector_MultiplyByNumber(PT_vector_T x, double r, PT_vector_T y) {  // y = r*x
	for (int j = 0; j < PP_N; j++)
		y[j] = x[j] * r;
}

static void Vector_DivideByNumber(PT_vector_T x, double r, PT_vector_T y) {  // x = x/r
	for (int j = 0; j < PP_N; j++)
		y[j] = x[j] / r;
}

static void Vector_Round(PT_vector_T x) {
	double floorValue;
	double fractionalPart;
	double sign;
	double absValue;

	for (int j = 0; j < PP_N; j++) {
		if (fabs(x[j]) < PP_EPS_ZERO) {
			x[j] = 0;
			continue;
		}
		absValue = fabs(x[j]);
		sign = x[j] > 0 ? 1 : -1;
		floorValue = floor(absValue);
		fractionalPart = absValue - floorValue;
		if (1 - fractionalPart < PP_EPS_ZERO) {
			x[j] = sign * (floorValue + 1);
			continue;
		}
		if (fractionalPart < PP_EPS_ZERO)
			x[j] = sign * floorValue;
	}
}

// Calculating unit vector
static void Vector_Unit(PT_vector_T vector, PT_vector_T unitVector) {
	double normOfVector = sqrt(Vector_NormSquare(vector));
	for (int j = 0; j < PP_N; j++)
		unitVector[j] = unitVector[j] / normOfVector;
};

static double ObjectiveF(PT_vector_T x) {
	double s = 0;
	for (int j = 0; j < PP_N; j++)
		s += PD_c[j] * x[j];
	return s;
}