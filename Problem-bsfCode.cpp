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

#ifdef PP_DEBUG
#define PP_EPS_ON		5E-1
#endif // PP_DEBUG

void PC_bsf_SetInitParameter(PT_bsf_parameter_T* parameter) {
	for (int j = 0; j < PP_N; j++) // Generating initial approximation
		parameter->x[j] = PD_apex[j];
};

void PC_bsf_Init(bool* success) {
	PD_state = PP_STATE_START;
	PD_m = PP_M;
	
	if (PP_EPS_IN < 0.1) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "PP_EPS_IN must be equal or greater than 0.1!\n";
		*success = false; return;
	}

	// ------------- Load LPP data -------------------
	FILE* stream;
	float buf;
	int m, n;
	stream = fopen(PP_LPP_FILE, "r");
	if (stream == NULL) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Failure of opening file '" << PP_LPP_FILE << "'.\n";
		*success = false; return;
	}

	fscanf(stream, "%d%d", &m, &n);
	if (n != PP_N || m != PP_M) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Error in input data '" << PP_LPP_FILE << "': PP_N != n and/or PP_M != m (PP_N = "
			<< PP_N << ", n = " << n << "; PP_M = " << PP_M << ", m = " << m << ").\n";
		*success = false; return;
	}

	for (int i = 0; i < PP_M; i++) {
		for (int j = 0; j < PP_N; j++) {
			fscanf(stream, "%f", &buf);
			PD_A[i][j] = buf;
		}
		fscanf(stream, "%f", &buf);
		PD_b[i] = buf;
	}

	for (int j = 0; j < PP_N; j++) {
		fscanf(stream, "%f", &buf);
		PD_c[j] = buf;
	}
	fclose(stream);

	// Generating a point inside the polytope
	for (int j = 0; j < PP_N; j++)
		PD_basePoint[j] = (double)PP_SF / 2;

	// Generating Coordinates of Apex Point
	ObjectiveUnitVector(PD_objectiveUnitVector);
	Vector_MultiplyByNumber(PD_objectiveUnitVector, PP_DIST_TO_APEX, PD_apex);
	Vector_PlusEquals(PD_apex, PD_basePoint);

	if (PD_apex[PP_N - 1] < PP_SF * 1.1) {
		if (BSF_sv_mpiRank == 0)
			cout << "PC_bsf_Init.Error: PD_apex[PP_N - 1] = " << PD_apex[PP_N - 1] << " < PP_SF * 1.1 = " << PP_SF * 1.1 << endl;
		*success = false;
		return;
	}
}

void PC_bsf_SetListSize(int* listSize) {
	*listSize = PP_MM;
}

void PC_bsf_SetMapListElem(PT_bsf_mapElem_T* elem, int i) {
	elem->a = PD_A[i];
	elem->b = &(PD_b[i]);
}

// 0. Pseudo-pojection
void PC_bsf_MapF(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T* reduceElem, int* success // 1 - reduceElem was produced successfully; 0 - otherwise
){
	*success = Vector_ProjectOnHalfspace(BSF_sv_parameter.x, mapElem->a, *mapElem->b, reduceElem->projection);
}

// 1. Movement on Polytope
void PC_bsf_MapF_1(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_1* reduceElem, int* success) {
	if (PointInHalfspace(BSF_sv_parameter.x, mapElem->a, *mapElem->b))
		reduceElem->pointIn = true;
	else
		reduceElem->pointIn = false;
}

// 2. ...
void PC_bsf_MapF_2(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_2* reduceElem, int* success) {
	// not used
}

void PC_bsf_MapF_3(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_3* reduceElem, int* success) {
	// not used
}

// 0. Pseudo-pojection
void PC_bsf_ReduceF(PT_bsf_reduceElem_T* x, PT_bsf_reduceElem_T* y, PT_bsf_reduceElem_T* z) { // z = x + y
	Vector_Addition(x->projection, y->projection, z->projection); 
}

// 1. Check
void PC_bsf_ReduceF_1(PT_bsf_reduceElem_T_1* x, PT_bsf_reduceElem_T_1* y, PT_bsf_reduceElem_T_1* z) {
	z->pointIn = x->pointIn && y->pointIn;
}

// 2. On
void PC_bsf_ReduceF_2(PT_bsf_reduceElem_T_2* x, PT_bsf_reduceElem_T_2* y, PT_bsf_reduceElem_T_2* z) {
	z->on = x->on + y->on;
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
#ifdef PP_MAX_ITER_COUNT
	if (BSF_sv_iterCounter > PP_MAX_ITER_COUNT) {
		cout << "-------------> PC_bsf_ProcessResults: Acceptable maximum number of iterations is exceeded: PP_MAX_ITER_COUNT = "
			<< PP_MAX_ITER_COUNT << endl;
		*exit = true;
		return;
	};
#endif // PP_MAX_ITER_COUNT

	Vector_Relaxation(reduceResult->projection, reduceCounter, PD_relaxationVector);
	Vector_PlusEquals(parameter->x, PD_relaxationVector);
	//Vector_Round(parameter->x);
}

// 1. Movement on Polytope  ========================================================
void PC_bsf_ProcessResults_1(
	PT_bsf_reduceElem_T_1* reduceResult,
	int reduceCounter, // Number of successfully produced Elrments of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* nextJob,
	bool* exit // "true" if Stopping Criterion is satisfied, and "false" otherwise
) {
#ifdef PP_MAX_ITER_COUNT
	if (BSF_sv_iterCounter > PP_MAX_ITER_COUNT) {
		cout << "-------------> Acceptable maximum number of iterations is exceeded: PP_MAX_ITER_COUNT = " << PP_MAX_ITER_COUNT << endl;
		*exit = true;
		return;
	}
#endif // PP_MAX_ITER_COUNT

	PD_pointIn = reduceResult->pointIn;
}

// 2. 
void PC_bsf_ProcessResults_2(
	PT_bsf_reduceElem_T_2* reduceResult,
	int reduceCounter, // Number of successfully produced Elrments of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* nextJob,
	bool* exit // "true" if Stopping Criterion is satisfied, and "false" otherwise
) {
	// not used
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

void PC_bsf_JobDispatcher(
	PT_bsf_parameter_T* parameter, // Current Approximation
	int *job,
	bool *exit
) {
	PT_vector_T objectiveVector;

	parameter->i = 0;

	switch (PD_state) {
	case PP_STATE_START://-------------------------- Start -----------------------------
		// Preparations for finding beginnig of path
		PD_state = PP_STATE_FIND_BEGINNING_OF_PATH;
		*job = PP_JOB_PSEUDOPOJECTION;
		break;
	case PP_STATE_FIND_BEGINNING_OF_PATH://------------------ u: Find Beginning of Path -----------------------------
		if (Vector_NormSquare(PD_relaxationVector) >= PP_EPS_RELAX * PP_EPS_RELAX)
			return;

#ifdef PP_DEBUG
		cout << "\t\t\tu = ";
		for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++)
			cout << setw(PP_SETW) << parameter->x[j];
		if (PP_OUTPUT_LIMIT < PP_N) cout << "	..." << setw(PP_SETW) << parameter->x[PP_N - 1];
		cout << "\tF(u) = " << ObjectiveF(parameter->x);
		cout << endl;

		/*PointTouch(parameter->x);
		cout << "--> Hyperplanes toushing u:\t";
		for (int i = 0; i < PP_M; i++)
			if (PD_pointTouch[i])
				cout << i << "\t";
		cout << endl;/**/

		//if (SaveSolution(parameter->x, PP_SOLUTION_FILE))
		//	cout << "\nPoint u is saved into the file '" << PP_SOLUTION_FILE << "'." << endl;
		//system("pause");

#endif // PP_DEBUG

		// Preparations for determining direction
		ObjectiveUnitVector(PD_objectiveUnitVector);
		Vector_MultiplyByNumber(PD_objectiveUnitVector, PP_OBJECTIVE_VECTOR_LENGTH, objectiveVector);
		Vector_Copy(parameter->x, PD_basePoint);
		Vector_PlusEquals(parameter->x, objectiveVector);
/*#ifdef PP_DEBUG
		cout << "\t\t\tv = ";
		for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++)
			cout << setw(PP_SETW) << parameter->x[j];
		if (PP_OUTPUT_LIMIT < PP_N) cout << "	..." << setw(PP_SETW) << parameter->x[PP_N - 1];
		cout << "\tF(v) = " << ObjectiveF(parameter->x);
		cout << endl;
		//system("pause");
#endif // PP_DEBUG /**/
		* job = PP_JOB_PSEUDOPOJECTION;
		PD_state = PP_STATE_DETERMINE_DIRECTION;
		PD_numDetDir = 0;
		break;

	case PP_STATE_DETERMINE_DIRECTION://------------------------- Determine Direction -----------------------------
		if (Vector_NormSquare(PD_relaxationVector) >= PP_EPS_RELAX * PP_EPS_RELAX)
			return;

#ifdef PP_DEBUG
		cout << "\t\t\tw = ";
		for (int j = 0; j < PP_N; j++)
			cout << setw(PP_SETW) << parameter->x[j];
		//cout << "\tF(w) = " << ObjectiveF(parameter->x);
		cout << endl;/**/

		/*PointTouch(parameter->x);
		cout << "--> Hyperplanes toushing w:\t";
		for (int i = 0; i < PP_M; i++)
			if (PD_pointTouch[i])
				cout << i << "\t";
		cout << endl;
		//system("pause");/**/
#endif // PP_DEBUG

		if (fabs(ObjectiveF(parameter->x) - ObjectiveF(PD_basePoint)) < PP_EPS_ZERO) {
			*exit = true;
#ifdef PP_DEBUG
			cout << setw(PP_SETW) << "F(u) = " << ObjectiveF(PD_basePoint) << " == F(w) = " << ObjectiveF(parameter->x) << "\n";
			//system("pause");
#endif // PP_DEBUG
			return;
		}

		if (ObjectiveF(parameter->x) <= ObjectiveF(PD_basePoint) - PP_EPS_OBJECTIVE) {
			*exit = true;
#ifdef PP_DEBUG
			cout << setw(PP_SETW) << "F(u) = " << ObjectiveF(PD_basePoint) << " >= F(w) = " << ObjectiveF(parameter->x) << "\n";
			//system("pause");
#endif // PP_DEBUG
			return;
		}

		PD_numDetDir++;
		if (PD_numDetDir > 2) {
			*exit = true;
#ifdef PP_DEBUG
			cout << "Number of Squential States 'Determine Direction' is greater than 2!\n";
			//system("pause");
#endif // PP_DEBUG
			return;
		}

		// Preparations for motion
		PD_shiftLength = PP_START_SHIFT_LENGTH;
		if (!GetDirection(PD_basePoint, parameter->x, PD_direction)) {
			*exit = true;
#ifdef PP_DEBUG
			cout << "Direction is too small!\n";
			//system("pause");
#endif // PP_DEBUG
			return;
		}
#ifdef PP_DEBUG //--------------------------------------//
		cout << "\t\t\tD = ";							//
		for (int j = 0; j < PP_N; j++)					//
			cout << setw(PP_SETW) << PD_direction[j];	//
		cout << endl;									//
		//system("pause");								//
		cout << "-----------------------------------\n";//
#endif // PP_DEBUG -------------------------------------//

		PD_newInequations = false;
		for (int j = 0; j < PP_N - 2; j++) {
			if (PD_direction[j] > 0)
				break;
			if (PD_direction[j] == 0) 
				continue;
			PD_direction[j] = 0;
			parameter->i = PD_m;
			Vector_Copy(PD_A[PD_m], parameter->a);
			parameter->a[j] = 1;
			parameter->b = PD_basePoint[j];
			PD_m += 2;
			PD_newInequations = true;
			break;
		}
		if (PD_newInequations) {
			PD_state = PP_STATE_FIND_BEGINNING_OF_PATH;
			*job = PP_JOB_PSEUDOPOJECTION;
#ifdef PP_DEBUG //------------------------------------------//
			cout << "\t\t\tD = ";							//
			for (int j = 0; j < PP_N; j++)					//
				cout << setw(PP_SETW) << PD_direction[j];	//
			cout << endl;									//
			//system("pause");								//
			cout << "-----------------------------------\n";//
#endif // PP_DEBUG -----------------------------------------//
			return;
		}

		PD_numSeqShifts = 0;
		PD_numShiftsSameLength = 0;
		Shift(PD_basePoint, PD_direction, PD_shiftLength, parameter->x);
		*job = PP_JOB_CHECK;
		PD_state = PP_STATE_MOVE_AND_CHECK;
		return;

	case PP_STATE_MOVE_AND_CHECK://-------------------------- t: Move and check -----------------------------
		PD_numSeqShifts++;

		if (PD_pointIn) {
			PD_numDetDir = 0;
			PD_numShiftsSameLength++;
			if (PD_numShiftsSameLength > PP_MAX_NUM_SHIFTS_SAME_LENGTH) {
				PD_shiftLength *= 2;
				PD_numShiftsSameLength = 0;
			}
#ifdef PP_DEBUG
			cout << "Sift = " << setw(PP_SETW) << PD_shiftLength << "\tt = ";
			for (int j = 0; j < PP_N; j++)
				cout << setw(PP_SETW) << parameter->x[j];
			cout << "\tF(t) = " << setw(PP_SETW) << ObjectiveF(parameter->x);
			cout << endl;
#endif // PP_DEBUG /**/
			Vector_Copy(parameter->x, PD_basePoint);
			Shift(PD_basePoint, PD_direction, PD_shiftLength, parameter->x);
			return;
		}

		if (PD_shiftLength >= PP_EPS_SHIFT && PD_numSeqShifts < PP_MAX_NUM_SEQ_SHIFTS) {
			PD_shiftLength /= 2;
			PD_numShiftsSameLength = 0;
			Shift(PD_basePoint, PD_direction, PD_shiftLength, parameter->x);
			return;
		}

#ifdef PP_DEBUG
		cout << "\n\t\t\tu = ";
		for (int j = 0; j < PP_N; j++)
			cout << setw(PP_SETW) << PD_basePoint[j];
		cout << "\tF(u) = " << ObjectiveF(PD_basePoint);
		cout << endl;
		if (SaveSolution(PD_basePoint, PP_SOLUTION_FILE))
			cout << "Point u is saved into the file '" << PP_SOLUTION_FILE << "'." << endl;/**/
		//system("pause");

		/*PointTouch(PD_basePoint);
		cout << "--> Hyperplanes toushing t:\t";
		for (int i = 0; i < PP_M; i++)
			if (PD_pointTouch[i])
				cout << i << "\t";
		cout << endl;/**/
#endif // PP_DEBUG

		// Preparations for determining direction
		Vector_MultiplyByNumber(PD_objectiveUnitVector, PP_OBJECTIVE_VECTOR_LENGTH, objectiveVector);
		Vector_Copy(PD_basePoint, parameter->x);
		Vector_PlusEquals(parameter->x, objectiveVector);
#ifdef PP_DEBUG
		/*cout << "\t\t\tv = ";
		for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++)
			cout << setw(PP_SETW) << parameter->x[j];
		if (PP_OUTPUT_LIMIT < PP_N) cout << "	..." << setw(PP_SETW) << parameter->x[PP_N - 1];
		//cout << "\tF(v) = " << ObjectiveF(parameter->x);
		cout << endl;/**/
		//system("pause");
#endif // PP_DEBUG
		* job = PP_JOB_PSEUDOPOJECTION;
		PD_state = PP_STATE_DETERMINE_DIRECTION;
		break;
	default://------------------------------------- default -----------------------------------
		cout << "PC_bsf_JobDispatcher: Undefined state!" << endl;
		break;
	}
}

void PC_bsf_ParametersOutput(PT_bsf_parameter_T parameter) {
	cout << "=================================================== Target ====================================================" << endl;
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
	cout << "Eps Relax:\t" << PP_EPS_RELAX << endl;
	cout << "Eps Shift:\t" << PP_EPS_SHIFT << endl;
	cout << "Eps Zero:\t" << PP_EPS_ZERO << endl;
	cout << "Eps Minimal Length of Direction Vector:\t" << PP_EPS_DIR << endl;
	cout << "Length of Objective Vector: " << PP_OBJECTIVE_VECTOR_LENGTH << endl;
	cout << "Start length of shift vector: " << PP_START_SHIFT_LENGTH << endl;

#ifdef PP_MATRIX_OUTPUT
	cout << "------- Matrix PD_A & Column PD_b -------" << endl;
	for (int i = 0; i < PP_M; i++) {
		cout << i << ")";
		for (int j = 0; j < PP_N; j++)
			cout << setw(PP_SETW) << PD_A[i][j];
		cout << "\t<=" << setw(PP_SETW) << PD_b[i] << endl;
	}
#endif // PP_MATRIX_OUTPUT

	cout << "Objective Function:\t"; 
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(PP_SETW) << PD_c[j];
	cout << (PP_OUTPUT_LIMIT < PP_N ? "	..." : "") << endl;
	cout << "Interior point:\t\t"; 
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(PP_SETW) << PD_basePoint[j];
	if (PP_OUTPUT_LIMIT < PP_N) cout << "	..." << setw(PP_SETW) << PD_basePoint[PP_N - 1];
	cout << endl;
	cout << "Apex Point:\t"; 
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(PP_SETW) << round(PD_apex[j]);
	if (PP_OUTPUT_LIMIT < PP_N) cout << "	..." << setw(PP_SETW) << PD_apex[PP_N - 1];
	cout << endl;
	cout << "-------------------------------------------" << endl;
}

void PC_bsf_CopyParameter(PT_bsf_parameter_T parameterIn, PT_bsf_parameter_T* parameterOutP) {
	for (int j = 0; j < PP_N; j++) { 
		parameterOutP->x[j] = parameterIn.x[j];
	}
	if (parameterIn.i > 0) {
		for (int j = 0; j < PP_N; j++) {
			PD_A[parameterIn.i][j] = parameterIn.a[j];
			PD_A[parameterIn.i + 1][j] = -parameterIn.a[j];
		}
		PD_b[parameterIn.i] = parameterIn.b;
		PD_b[parameterIn.i + 1] = -parameterIn.b;
	}
}

// 0. Apex Pseudo-pojection  ==========================================================
void PC_bsf_IterOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob) {
	cout << "------------------ 0. Apex Pseudo-pojection. Iter # " << BSF_sv_iterCounter << " -------------------" << endl;
	/* debug */// cout << "Elapsed time: " << round(elapsedTime) << endl;
	cout << "Approximat. :"; 
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(PP_SETW) << parameter.x[j];
	if (PP_OUTPUT_LIMIT < PP_N) cout << "	..." << setw(PP_SETW) << parameter.x[PP_N - 1];
	cout << endl;
}

// 1. Movement on Polytope
void PC_bsf_IterOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	cout << "------------------ 1. Movement on Polytope. Iter # " << BSF_sv_iterCounter << " ------------------" << endl;
	/* debug */// cout << "Elapsed time: " << round(elapsedTime) << endl;
	cout << "PD_basePoint:";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++)
		cout << setw(PP_SETW) << PD_basePoint[j];
	if(PP_OUTPUT_LIMIT < PP_N) cout << "	..." << setw(PP_SETW) << PD_basePoint[PP_N - 1];
	cout << endl;
	cout << "PD_direction:";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++)
		cout << setw(PP_SETW) << PD_direction[j];
	if(PP_OUTPUT_LIMIT < PP_N) cout << "	..." << setw(PP_SETW) << PD_direction[PP_N - 1];
	cout << endl;
	cout << "Sift Length = " << PD_shiftLength << endl;
};

// 2. ...
void PC_bsf_IterOutput_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	// not used
}


void PC_bsf_IterOutput_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	// not used
}

// 0. Apex Pseudo-pojection
void PC_bsf_ProblemOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	ProblemOutput(t);
}

// 1. Movement on Polytope
void PC_bsf_ProblemOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,double t) {
	ProblemOutput(t);
}

// 2. Point Pseudo-projection
void PC_bsf_ProblemOutput_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,double t) {
	ProblemOutput(t);
}

void PC_bsf_ProblemOutput_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,double t) {
	// not used
}

//----------------------- Assigning Values to BSF-skeleton Variables (Do not modify!) -----------------------
void PC_bsfAssignAddressOffset(int value) { BSF_sv_addressOffset = value; };
void PC_bsfAssignIterCounter(int value) { BSF_sv_iterCounter = value; };
void PC_bsfAssignJobCase(int value) { BSF_sv_jobCase = value; };
void PC_bsfAssignMpiMaster(int value) { BSF_sv_mpiMaster = value; };
void PC_bsfAssignMpiRank(int value) { BSF_sv_mpiRank = value; };
void PC_bsfAssignNumberInSublist(int value) { BSF_sv_numberInSublist = value; };
void PC_bsfAssignNumOfWorkers(int value) { BSF_sv_numOfWorkers = value; };
void PC_bsfAssignParameter(PT_bsf_parameter_T parameter) { PC_bsf_CopyParameter(parameter, &BSF_sv_parameter); }
void PC_bsfAssignSublistLength(int value) { BSF_sv_sublistLength = value; };

//---------------------------------- Problem functions -------------------------
inline double Vector_DotProductSquare(PT_vector_T x, PT_vector_T y) {
	double sum = 0;
	for (int j = 0; j < PP_N; j++) {
		sum += x[j] * y[j];
	}
	return sum;
}

inline void Vector_Relaxation(PT_vector_T sumOfProjections, int numberOfProjections, PT_vector_T relaxationVector) {
	for (int j = 0; j < PP_N; j++) {
		relaxationVector[j] = sumOfProjections[j] / (double)numberOfProjections;
	}
}

inline double Vector_NormSquare(PT_vector_T x) { 
	double sum = 0;

	for (int j = 0; j < PP_N; j++) {
		sum += x[j] * x[j];
	}
	return sum;
}

inline bool PointInHalfspace // If the point belongs to the Halfspace with prescigion of PP_EPS_ZERO
(PT_vector_T point, PT_vector_T a, PT_float_T b) {
	return Vector_DotProductSquare(a, point) <= b + PP_EPS_IN;
}

inline bool PointInHalfspace_s(PT_vector_T point, PT_vector_T a, PT_float_T b) { // If the point exactly belongs to the Half-space 
	return Vector_DotProductSquare(a, point) <= b;
}

inline void Shift(PT_vector_T basePoint, PT_vector_T direction, double siftLength, PT_vector_T endPoint) {
	for (int j = 0; j < PP_N; j++) { 
		endPoint[j] = basePoint[j] + direction[j] * siftLength;
	}
}

// Calculating unit vector of direction from startPoint to endPoint
inline bool GetDirection(PT_vector_T startPoint, PT_vector_T endPoint, PT_vector_T unitVector) {
	for (int j = 0; j < PP_N; j++) { 
		unitVector[j] = endPoint[j] - startPoint[j];
	}
	double normOfUnitVector = sqrt(Vector_NormSquare(unitVector));
/*#ifdef PP_DEBUG
	cout << setw(PP_SETW) << "||w-u|| = " << normOfUnitVector << endl;
#endif // PP_DEBUG /**/
	
	if (normOfUnitVector < PP_EPS_DIR)
		return false;

	for (int j = 0; j < PP_N; j++) { 
		unitVector[j] /= normOfUnitVector;
	}

	for (int j = 0; j < PP_N; j++) { 
		if (fabs(unitVector[j]) < PP_EPS_ZERO)
			unitVector[j] = 0;
	}

	return true;
}

inline void Vector_Copy(PT_vector_T fromPoint, PT_vector_T toPoint) { // toPoint = fromPoint
	for (int j = 0; j < PP_N; j++) { 
		toPoint[j] = fromPoint[j];
	}
}

inline void Vector_PlusEquals(PT_vector_T equalVector, PT_vector_T plusVector) { // equalVector += plusVector
	for (int j = 0; j < PP_N; j++) {
		equalVector[j] += plusVector[j];
	}
}

inline void Vector_MinusEquals(PT_vector_T equalPoint, PT_vector_T minusVector) { // equalPoint += minusVector
	for (int j = 0; j < PP_N; j++) { 
		equalPoint[j] -= minusVector[j];
	}
};

inline void Vector_Addition(PT_vector_T x, PT_vector_T y, PT_vector_T z) {  // z = x + y
	for (int j = 0; j < PP_N; j++) {
		z[j] = x[j] + y[j];
	}
}

inline void Vector_Subtraction(PT_vector_T x, PT_vector_T y, PT_vector_T z) {  // z = x - y
	for (int j = 0; j < PP_N; j++) { 
		z[j] = x[j] - y[j];
	}
}

inline void Vector_MultiplyByNumber(PT_vector_T x, double r, PT_vector_T y) {  // y = r*x
	for (int j = 0; j < PP_N; j++) { 
		y[j] = x[j] * r;
	}
}

inline void Vector_DivideByNumber(PT_vector_T x, double r, PT_vector_T y) {  // x = x/r
	for (int j = 0; j < PP_N; j++) { 
		y[j] = x[j] / r;
	}
}

inline void Vector_Round(PT_vector_T x) {
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
inline void Vector_Unit(PT_vector_T vector, PT_vector_T unitVector) {
	double normOfVector = sqrt(Vector_NormSquare(vector));
	for (int j = 0; j < PP_N; j++) { 
		unitVector[j] = unitVector[j] / normOfVector;
	}
};

inline double ObjectiveF(PT_vector_T x) {
	double s = 0;
	for (int j = 0; j < PP_N; j++)
		s += PD_c[j] * x[j];
	return s;
}

static bool SaveSolution(PT_vector_T x, const char* filename) {
	FILE* stream;

	stream = fopen(filename, "w");
	if (stream == NULL) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Failure of opening file '" << filename << "'.\n";
		return false;
	}

	fprintf(stream, "%d\n", PP_N);

	for (int j = 0; j < PP_N; j++)
		fprintf(stream, "%f\t", x[j]);

	fclose(stream);
	return true;
}

// Point projection onto Half-space <a,x> <= b
inline bool // true if the point does not belong to the half-space and false otherwise 
Vector_ProjectOnHalfspace(PT_vector_T point, PT_vector_T a, PT_float_T b, PT_vector_T projection) {
	double factor;
	double aNormSquare = Vector_NormSquare(a);

	if(aNormSquare < PP_EPS_ZERO)
		return false;

	factor = (b - Vector_DotProductSquare(point, a)) / aNormSquare;

	if (factor > PP_EPS_ZERO)
		return false;

	for (int j = 0; j < PP_N; j++) { 
		projection[j] = factor * a[j];
	}

	return true;
}

inline void ProblemOutput(double elapsedTime) {
	cout << "=============================================" << endl;
	cout << "Elapsed time: " << elapsedTime << endl;
	cout << "Iterations: " << BSF_sv_iterCounter << endl;
	if (SaveSolution(PD_basePoint, PP_SOLUTION_FILE))
		cout << "Solution is saved into the file '" << PP_SOLUTION_FILE << "'." << endl;
	cout << "Solution: ";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(PP_SETW) << PD_basePoint[j];
	if (PP_OUTPUT_LIMIT < PP_N) cout << "	..." << setw(PP_SETW) << PD_basePoint[PP_N - 1];
	cout << endl;
}

inline PT_float_T Distance(PT_vector_T x, PT_vector_T y) {
	PT_vector_T z;
	Vector_Subtraction(x, y, z);
	return sqrt(Vector_NormSquare(z));
}

inline void ObjectiveUnitVector(PT_vector_T objectiveUnitVector) { // Calculating Objective Unit Vector
	double c_norm = sqrt(Vector_NormSquare(PD_c));
	Vector_DivideByNumber(PD_c, c_norm, objectiveUnitVector);
}

#ifdef PP_DEBUG //=========================== Debug Functions =====================================
static void PointTouch(PT_vector_T point) { // Mark the hyperplanes that touch the point
	for (int i = 0; i < PP_M; i++)
		if (PointOnHyperplane(point, PD_A[i], PD_b[i]))
			PD_pointTouch[i] = true;
		else
			PD_pointTouch[i] = false;
}

inline bool PointOnHyperplane(PT_vector_T point, PT_vector_T a, PT_float_T b) { // If the point belongs to the hyperplane with prescigion of PP_EPS_ZERO
	return fabs(Vector_DotProductSquare(a, point) - b) < PP_EPS_ON;
}
#endif // PP_DEBUG ===================================================================================