/*==============================================================================
Project: LiFe
Theme: Apex Method (No MPI)
Module: Problem-bsfCode.cpp (Implementation of the Problem)
Prefix: PI
Author: Leonid B. Sokolinsky
This source code has been produced with using BSF-skeleton
==============================================================================*/
#include "Problem-Data.h"			// Problem Types 
#include "Problem-Forwards.h"		// Problem Function Forwards
#include "Problem-bsfParameters.h"	// BSF-skeleton parameters
#include "BSF-SkeletonVariables.h"	// Skeleton Variables
using namespace std;

void PC_bsf_SetInitParameter(PT_bsf_parameter_T* parameter) {
	for (int j = 0; j < PD_n; j++) // Generating initial approximation
		parameter->x[j] = PD_basePoint[j];
};

void PC_bsf_Init(bool* success) {
	PD_state = PP_STATE_START;

#ifdef MTX_FORMAT
	* success = LoadMatrixFormat();
#else
	* success = LoadLppFormat();
#endif // MTX_FORMAT

	if (*success == false)
		return;

	// Open trace file
	const char* traceFile = PD_traceFile.c_str();
	PD_traceStream = fopen(traceFile, "w");
	if (PD_traceStream == NULL) {
		cout << "Failure of opening file '" << PD_traceFile << "'.\n";
		*success = false;
		system("pause");
		return;
	}

	// Generating Coordinates of starting point
	for (int j = 0; j < PD_n; j++)
		PD_basePoint[j] = 0;

	PD_m = PP_MM;

	ObjectiveUnitVector(PD_objectiveUnitVector);
	Vector_MultiplyByNumber(PD_objectiveUnitVector, PP_OBJECTIVE_VECTOR_LENGTH, PD_objectiveVector);

}

void PC_bsf_SetListSize(int* listSize) {
	*listSize = PD_m;
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

// 1. Approximate check
void PC_bsf_MapF_1(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_1* reduceElem, int* success) {
	reduceElem->pointIn = PointInHalfspace(BSF_sv_parameter.x, mapElem->a, *mapElem->b);
}

// 2. Exact check
void PC_bsf_MapF_2(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_2* reduceElem, int* success) {
	reduceElem->pointIn = PointInHalfspace_s(BSF_sv_parameter.x, mapElem->a, *mapElem->b);
}

void PC_bsf_MapF_3(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_3* reduceElem, int* success) {
	// not used
}

// 0. Pseudo-pojection
void PC_bsf_ReduceF(PT_bsf_reduceElem_T* x, PT_bsf_reduceElem_T* y, PT_bsf_reduceElem_T* z) { // z = x + y
	Vector_Addition(x->projection, y->projection, z->projection); 
}

// 1. Approximate check
void PC_bsf_ReduceF_1(PT_bsf_reduceElem_T_1* x, PT_bsf_reduceElem_T_1* y, PT_bsf_reduceElem_T_1* z) {
	z->pointIn = x->pointIn && y->pointIn;
}

// 2. Exact check
void PC_bsf_ReduceF_2(PT_bsf_reduceElem_T_2* x, PT_bsf_reduceElem_T_2* y, PT_bsf_reduceElem_T_2* z) {
	z->pointIn = x->pointIn && y->pointIn;
}

void PC_bsf_ReduceF_3(PT_bsf_reduceElem_T_3* x, PT_bsf_reduceElem_T_3* y, PT_bsf_reduceElem_T_3* z) {
	// not used
}

//0. Start
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
#ifdef PP_MAX_ITER_COUNT
	if (BSF_sv_iterCounter > PP_MAX_ITER_COUNT) {
		cout << "-------------> Acceptable maximum number of iterations is exceeded: PP_MAX_ITER_COUNT = " << PP_MAX_ITER_COUNT << endl;
		*exit = true;
		return;
	}
#endif // PP_MAX_ITER_COUNT

	PD_pointIn = reduceResult->pointIn;
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
	const char* traceFile = PD_traceFile.c_str();

	parameter->i = 0;

	switch (PD_state) {
	case PP_STATE_START://-------------------------- Start -----------------------------
		if (PointInPolytope_s(PD_basePoint)) { 
			Vector_Copy(PD_objectiveUnitVector, PD_direction);
			PD_shiftLength = PP_START_SHIFT_LENGTH;
			// Preparations for moving inside the polytope
			PD_numShiftsSameLength = 0;
			Shift(PD_basePoint, PD_direction, PD_shiftLength, parameter->x);
			*job = PP_JOB_CHECK_S;
			PD_state = PP_MOVE_INSIDE_POLYTOPE;
		}
		else {
			// Preparations for finding a start point
			PD_state = PP_STATE_FIND_START_POINT;
			*job = PP_JOB_PSEUDOPOJECTION;
#ifdef PP_DEBUG
			cout << "--------- Pseudoprojecting ------------\n";
#ifdef PP_PAUSE
			system("pause");
#endif // PP_PAUSE
#endif // PP_DEBUG
		}
		break;
	case PP_MOVE_INSIDE_POLYTOPE: //--------------------- Moving inside the polytope ----------------------------
		if (PD_pointIn) {
			PD_numDetDir = 0;
			PD_numShiftsSameLength++;
			if (PD_numShiftsSameLength > PP_MAX_NUM_SHIFTS_SAME_LENGTH) {
				PD_shiftLength *= 2;
				PD_numShiftsSameLength = 0;
			}
#ifdef PP_DEBUG
			cout << "Sift = " << setw(PP_SETW) << PD_shiftLength << "\tt = ";
			for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
				cout << setw(PP_SETW) << parameter->x[j];
			if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
			cout << "\tF(t) = " << setw(PP_SETW) << ObjectiveF(parameter->x);
			cout << endl;
#endif // PP_DEBUG
			Vector_Copy(parameter->x, PD_basePoint);
			Shift(PD_basePoint, PD_direction, PD_shiftLength, parameter->x);
			return;
		}

		if (PD_shiftLength >= PP_EPS_RELAX) {
			PD_shiftLength /= 2;
			PD_numShiftsSameLength = 0;
			Shift(PD_basePoint, PD_direction, PD_shiftLength, parameter->x);
			return;
		}
		WriteTrace(parameter->x);

		// Preparations for determining direction
		Vector_Copy(parameter->x, PD_basePoint);
		Vector_PlusEquals(parameter->x, PD_objectiveVector);
		*job = PP_JOB_PSEUDOPOJECTION;
		PD_state = PP_STATE_DETERMINE_DIRECTION;
#ifdef PP_DEBUG
		cout << "--------- Determine Direction ------------\n";
#ifdef PP_PAUSE
		system("pause");
#endif // PP_PAUSE
#endif // PP_DEBUG
		PD_numDetDir = 0;
		break;
	case PP_STATE_DETERMINE_DIRECTION://------------------------- Determine Direction -----------------------------
		if (Vector_NormSquare(PD_relaxationVector) >= PP_EPS_RELAX * PP_EPS_RELAX)
			return;

		for (int j = 0; j < PD_n; j++) {
			if (fabs(PD_basePoint[j]) < PP_EPS_ZERO)
				PD_basePoint[j] = 0;
			if (fabs(parameter->x[j]) < PP_EPS_ZERO)
				parameter->x[j] = 0;
		}

#ifdef PP_DEBUG
		cout << "\t\t\tu = ";
		for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
			cout << setw(PP_SETW) << PD_basePoint[j];
		if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
		cout << "\tF(u) = " << ObjectiveF(PD_basePoint);
		cout << endl;/**/
		cout << "\t\t\tw = ";
		for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
			cout << setw(PP_SETW) << parameter->x[j];
		if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
		cout << "\tF(w) = " << ObjectiveF(parameter->x);
		cout << endl;/**/
#endif // PP_DEBUG

		if (fabs(ObjectiveF(parameter->x) - ObjectiveF(PD_basePoint)) < PP_EPS_OBJECTIVE) {
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
		for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)					//
			cout << setw(PP_SETW) << PD_direction[j];	//
		if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
		cout << endl;									//
#endif // PP_DEBUG -------------------------------------//

#ifdef PP_MAJOR_COORDINATES_CAN_NOT_DECREASE
		PD_newInequations = false;
		for (int j = 0; j < PD_n - 1; j++) {
			if (PD_direction[j] > 0)
				break;
			if (PD_direction[j] == 0)
				continue;


			PD_A[PD_m][j] = -1;
			PD_b[PD_m] = -PD_basePoint[j];
			PD_m++;
			PD_newInequations = true;
			break;
		}
		if (PD_newInequations) {
			Vector_Copy(PD_basePoint, parameter->x);
			Vector_PlusEquals(parameter->x, PD_objectiveVector);
			return;
		}
#endif
		PD_numShiftsSameLength = 0;
		Shift(PD_basePoint, PD_direction, PD_shiftLength, parameter->x);
		*job = PP_JOB_CHECK;
		PD_state = PP_STATE_MOVE_AND_CHECK;
#ifdef PP_DEBUG
		cout << "--------- Movement on surface ------------\n";
#ifdef PP_PAUSE
		system("pause");
#endif // PP_PAUSE
#endif // PP_DEBUG
		break;

	case PP_STATE_MOVE_AND_CHECK://-------------------------- t: Move and check -----------------------------
		if (PD_pointIn) {
			PD_numDetDir = 0;
			PD_numShiftsSameLength++;
			if (PD_numShiftsSameLength > PP_MAX_NUM_SHIFTS_SAME_LENGTH) {
				PD_shiftLength *= 2;
				PD_numShiftsSameLength = 0;
			}
#ifdef PP_DEBUG
			cout << "Sift = " << setw(PP_SETW) << PD_shiftLength << "\tt = ";
			for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
				cout << setw(PP_SETW) << parameter->x[j];
			if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
			cout << "\tF(t) = " << setw(PP_SETW) << ObjectiveF(parameter->x);
			cout << endl;
#endif // PP_DEBUG /**/
			Vector_Copy(parameter->x, PD_basePoint);
			Shift(PD_basePoint, PD_direction, PD_shiftLength, parameter->x);
			return;
		}

		if (PD_shiftLength >= PP_EPS_SHIFT) {
			PD_shiftLength /= 2;
			PD_numShiftsSameLength = 0;
			Shift(PD_basePoint, PD_direction, PD_shiftLength, parameter->x);
			return;
		}

		// Preparations for landing
		Vector_Copy(PD_basePoint, parameter->x);
		*job = PP_JOB_PSEUDOPOJECTION;
		PD_state = PP_STATE_LANDING;
#ifdef PP_DEBUG
		cout << "--------- Landing ------------\n";
#endif // PP_DEBUG
		break;
	case PP_STATE_LANDING://-------------------------- Landing -----------------------------
		if (Vector_NormSquare(PD_relaxationVector) >= PP_EPS_RELAX * PP_EPS_RELAX)
			return;

		Vector_Copy(parameter->x, PD_basePoint);
		WriteTrace(PD_basePoint);

		// Preparations for determining direction
		Vector_PlusEquals(parameter->x, PD_objectiveVector);

		* job = PP_JOB_PSEUDOPOJECTION;
		PD_state = PP_STATE_DETERMINE_DIRECTION;
#ifdef PP_DEBUG
		cout << "--------- Determining direction ------------\n";
#ifdef PP_PAUSE
		system("pause");
#endif // PP_PAUSE
#endif // PP_DEBUG
		break;
	case PP_STATE_FIND_START_POINT://-------------------------- Finding a start point -----------------------------
		if (Vector_NormSquare(PD_relaxationVector) >= PP_EPS_RELAX * PP_EPS_RELAX)
			return;
		// Preparations for moving inside the polytope
		Vector_Copy(parameter->x, PD_basePoint);
		Vector_Copy(PD_objectiveUnitVector, PD_direction);
		PD_shiftLength = PP_START_SHIFT_LENGTH;
		PD_numShiftsSameLength = 0;
		Shift(PD_basePoint, PD_direction, PD_shiftLength, parameter->x);
		if (PointInPolytope_s(parameter->x)) {
			*job = PP_JOB_CHECK_S;
			PD_state = PP_MOVE_INSIDE_POLYTOPE;
#ifdef PP_DEBUG
			cout << "--------- Moving inside polytope ------------\n";
#endif // PP_DEBUG
		}
		else {
			WriteTrace(PD_basePoint);
			// Preparations for determining direction
			Vector_Copy(PD_basePoint, parameter->x);
			Vector_PlusEquals(parameter->x, PD_objectiveVector);
			*job = PP_JOB_PSEUDOPOJECTION;
			PD_state = PP_STATE_DETERMINE_DIRECTION;
			PD_numDetDir = 0;
#ifdef PP_DEBUG
			cout << "--------- Determining direction ------------\n";
#ifdef PP_PAUSE
			system("pause");
#endif // PP_PAUSE
#endif // PP_DEBUG
		}
		break;
	default://------------------------------------- default -----------------------------------
		cout << "PC_bsf_JobDispatcher: Undefined state!" << endl;
		break;
	}
}

void PC_bsf_ParametersOutput(PT_bsf_parameter_T parameter) {
	cout << "=================================================== Target ====================================================" << endl;
	//

#ifdef PP_BSF_OMP
#ifdef PP_BSF_NUM_THREADS
	cout << "Number of Threads: " << PP_BSF_NUM_THREADS << endl;
#else
	cout << "Number of Threads: " << omp_get_num_procs() << endl;
#endif 
#else
	cout << "OpenMP is turned off!" << endl;
#endif 

#ifdef PP_BSF_FRAGMENTED_MAP_LIST
	cout << "Map List is Fragmented." << endl;
#else
	cout << "Map List is not Fragmented." << endl;
#endif

	cout << "Dimension: N = " << PD_n << endl;
	cout << "Number of Constraints: M = " << PD_m_init << endl;
	cout << "Scale Factor: SF = " << PP_SF << endl;
	cout << "Eps Relax:\t" << PP_EPS_RELAX << endl;
	cout << "Eps Shift:\t" << PP_EPS_SHIFT << endl;
	cout << "Eps Zero:\t" << PP_EPS_ZERO << endl;
	cout << "Eps Minimal Length of Direction Vector:\t" << PP_EPS_DIR << endl;
	cout << "Length of Objective Vector: " << PP_OBJECTIVE_VECTOR_LENGTH << endl;
	cout << "Start length of shift vector: " << PP_START_SHIFT_LENGTH << endl;

#ifdef PP_MATRIX_OUTPUT
	cout << "------- Matrix PD_A & Column PD_b -------" << endl;
	for (int i = 0; i < PD_m_init; i++) {
		cout << i << ")";
		for (int j = 0; j < PD_n; j++)
			cout << setw(PP_SETW) << PD_A[i][j];
		cout << "\t<=" << setw(PP_SETW) << PD_b[i] << endl;
	}
#endif // PP_MATRIX_OUTPUT

	cout << "Objective Function:\t";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++) cout << setw(PP_SETW) << PD_c[j];
	cout << (PP_OUTPUT_LIMIT < PD_n ? "	..." : "") << endl;
	cout << "Starting point:\t\t";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++) cout << setw(PP_SETW) << PD_basePoint[j];
	if (PP_OUTPUT_LIMIT < PD_n)
		cout << "	...";
	cout << endl;
	cout << "-------------------------------------------" << endl;
}

void PC_bsf_CopyParameter(PT_bsf_parameter_T parameterIn, PT_bsf_parameter_T* parameterOutP) {
	for (int j = 0; j < PD_n; j++) { 
		parameterOutP->x[j] = parameterIn.x[j];
	}
	if (parameterIn.i > 0) {
		for (int j = 0; j < PD_n; j++) 
			PD_A[parameterIn.i][j] = parameterIn.a[j];
		PD_b[parameterIn.i] = parameterIn.b;
	}
}

// 0. Start
void PC_bsf_IterOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob) {
	cout << "------------------ 0. Start. Iter # " << BSF_sv_iterCounter << " -------------------" << endl;
	/* debug */// cout << "Elapsed time: " << round(elapsedTime) << endl;
	cout << "Approximat. :"; 
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++) cout << setw(PP_SETW) << parameter.x[j];
	if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
	cout << endl;
}

// 1. Movement on Polytope
void PC_bsf_IterOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	cout << "------------------ 1. Movement on Polytope. Iter # " << BSF_sv_iterCounter << " ------------------" << endl;
	/* debug */// cout << "Elapsed time: " << round(elapsedTime) << endl;
	cout << "PD_basePoint:";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
		cout << setw(PP_SETW) << PD_basePoint[j];
	if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
	cout << endl;
	cout << "PD_direction:";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
		cout << setw(PP_SETW) << PD_direction[j];
	if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
	cout << endl;
	cout << "Sift Length = " << PD_shiftLength << endl;
};

// 2. Movement iside Polytope
void PC_bsf_IterOutput_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	cout << "------------------ 2. Movement inside Polytope. Iter # " << BSF_sv_iterCounter << " ------------------" << endl;
	/* debug */// cout << "Elapsed time: " << round(elapsedTime) << endl;
	cout << "PD_basePoint:";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
		cout << setw(PP_SETW) << PD_basePoint[j];
	if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
	cout << endl;
	cout << "PD_direction:";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
		cout << setw(PP_SETW) << PD_direction[j];
	if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
	cout << endl;
	cout << "Sift Length = " << PD_shiftLength << endl;
}


void PC_bsf_IterOutput_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	// not used
}

// 0. Start
void PC_bsf_ProblemOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	ProblemOutput(t);
}

// 1. Movement on Polytope
void PC_bsf_ProblemOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,double t) {
	ProblemOutput(t);
}

// 2. Movement inside Polytope
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
	for (int j = 0; j < PD_n; j++) 
		sum += x[j] * y[j];
	return sum;
}

inline void Vector_Relaxation(PT_vector_T sumOfProjections, int numberOfProjections, PT_vector_T relaxationVector) {
	for (int j = 0; j < PD_n; j++) {
		relaxationVector[j] = sumOfProjections[j] / (double)numberOfProjections;
	}
}

inline double Vector_NormSquare(PT_vector_T x) { 
	double sum = 0;

	for (int j = 0; j < PD_n; j++) {
		sum += x[j] * x[j];
	}
	return sum;
}

inline bool PointInHalfspace // If the point belongs to the Halfspace with prescigion of PP_EPS_IN
(PT_vector_T point, PT_vector_T a, PT_float_T b) {
	return Vector_DotProductSquare(a, point) <= b + PP_EPS_IN;
}

inline bool PointInHalfspace_s(PT_vector_T point, PT_vector_T a, PT_float_T b) { // If the point belongs to the Halfspace
	return Vector_DotProductSquare(a, point) <= b;
}

inline bool PointInPolytope_s(PT_vector_T x) { // If the point belongs to the polytope
	for (int i = 0; i < PD_m; i++)
		if (!PointInHalfspace_s(x, PD_A[i], PD_b[i]))
			return false;
	return true;
}

inline void Shift(PT_vector_T basePoint, PT_vector_T direction, double siftLength, PT_vector_T endPoint) {
	for (int j = 0; j < PD_n; j++) { 
		endPoint[j] = basePoint[j] + direction[j] * siftLength;
	}
}

// Calculating unit vector of direction from startPoint to endPoint
inline bool GetDirection(PT_vector_T startPoint, PT_vector_T endPoint, PT_vector_T unitVector) {
	for (int j = 0; j < PD_n; j++) { 
		unitVector[j] = endPoint[j] - startPoint[j];
	}
	double normOfUnitVector = sqrt(Vector_NormSquare(unitVector));
/*#ifdef PP_DEBUG
	cout << setw(PP_SETW) << "||w-u|| = " << normOfUnitVector << endl;
#endif // PP_DEBUG /**/
	
	if (normOfUnitVector < PP_EPS_DIR)
		return false;

	for (int j = 0; j < PD_n; j++) { 
		unitVector[j] /= normOfUnitVector;
	}

	for (int j = 0; j < PD_n; j++) { 
		if (fabs(unitVector[j]) < PP_EPS_ZERO)
			unitVector[j] = 0;
	}

	return true;
}

inline void Vector_Copy(PT_vector_T fromPoint, PT_vector_T toPoint) { // toPoint = fromPoint
	for (int j = 0; j < PD_n; j++) {
		toPoint[j] = fromPoint[j];
	}
}

inline void Vector_PlusEquals(PT_vector_T equalVector, PT_vector_T plusVector) { // equalVector += plusVector
	for (int j = 0; j < PD_n; j++) {
		equalVector[j] += plusVector[j];
	}
}

inline void Vector_MinusEquals(PT_vector_T equalPoint, PT_vector_T minusVector) { // equalPoint += minusVector
	for (int j = 0; j < PD_n; j++) { 
		equalPoint[j] -= minusVector[j];
	}
};

inline void Vector_Addition(PT_vector_T x, PT_vector_T y, PT_vector_T z) {  // z = x + y
	for (int j = 0; j < PD_n; j++) {
		z[j] = x[j] + y[j];
	}
}

inline void Vector_Subtraction(PT_vector_T x, PT_vector_T y, PT_vector_T z) {  // z = x - y
	for (int j = 0; j < PD_n; j++) { 
		z[j] = x[j] - y[j];
	}
}

inline void Vector_MultiplyByNumber(PT_vector_T x, double r, PT_vector_T y) {  // y = r*x
	for (int j = 0; j < PD_n; j++) { 
		y[j] = x[j] * r;
	}
}

inline void Vector_DivideByNumber(PT_vector_T x, double r, PT_vector_T y) {  // x = x/r
	for (int j = 0; j < PD_n; j++) { 
		y[j] = x[j] / r;
	}
}

inline void Vector_Round(PT_vector_T x) {
	double floorValue;
	double fractionalPart;
	double sign;
	double absValue;

	for (int j = 0; j < PD_n; j++) { 
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
	for (int j = 0; j < PD_n; j++) { 
		unitVector[j] = unitVector[j] / normOfVector;
	}
};

inline double ObjectiveF(PT_vector_T x) {
	double s = 0;
	for (int j = 0; j < PD_n; j++)
		s += PD_c[j] * x[j];
	return s;
}

static bool LoadLppFormat() {
	//
	//
	const char* lppFile = PD_lppFile.c_str();

	FILE* stream;
	float buf;
	stream = fopen(lppFile, "r");

	if (stream == NULL) {
		cout << "Failure of opening file '" << lppFile << "'.\n";
		system("pause");
		return false;
	}

	if (fscanf(stream, "%d%d", &PD_m_init, &PD_n) == 0) {
		cout << "Unexpected end of file" << endl;
		system("pause");
		return false;
	}

	if (PD_n > PP_N) {
		cout << "Invalid input data: Space dimension n = " << PD_n << " must be < " << PP_N + 1 << "\n";
		system("pause");
		return false;
	}


	for (int i = 0; i < PD_m_init; i++) {
		for (int j = 0; j < PD_n; j++) {
			if (fscanf(stream, "%f", &buf) == 0) { cout << "Unexpected end of file" << endl; return false; }
			PD_A[i][j] = buf;
		}
		if (fscanf(stream, "%f", &buf) == 0) { cout << "Unexpected end of file" << endl; return false; }
		PD_b[i] = buf;
	}

	for (int j = 0; j < PD_n; j++) {
		if (fscanf(stream, "%f", &buf) == 0) { cout << "Unexpected end of file" << endl; return false; }
		PD_c[j] = buf;
	}

	fclose(stream);
	return true;
}

static bool LoadMatrixFormat() {
	int nor, // Number of matrix rows
		noc, // Number of matrix columns
		non; // Number of non-zero elements
	int m_cur;	// Current number of inequalities
	int noe;	// Number of equations
	const char* mtxFile;
	FILE* stream;// Input stream
	float buf;

	//--------------- Reading A ------------------
	PD_MTX_File_A = PP_PATH;
	PD_MTX_File_A += PP_MTX_PREFIX;
	PD_MTX_File_A += PP_MTX_PROBLEM_NAME;
	PD_MTX_File_A += PP_MTX_POSTFIX_A;
	mtxFile = PD_MTX_File_A.c_str();
	stream = fopen(mtxFile, "r+b");

	if (stream == NULL) {
		cout << "Failure of opening file '" << mtxFile << "'.\n";
		//
		return false;
	}

	SkipComments(stream);
	if (fscanf(stream, "%d%d%d", &nor, &noc, &non) < 3) {
		cout << "Unexpected end of file" << endl;
		//
		return false;
	}
	m_cur = noe = nor;
	PD_n = noc;
	PD_m_init = 2 * nor + noc;

	if (PD_n > PP_N) {
		cout << "Invalid input data: Space dimension n = " << PD_n << " must be < " << PP_N + 1 << "\n";
		//
		return false;
	}

	if (PD_m_init > PP_MM) {
		cout << "Invalid input data: number of inequalities m = " << PD_m_init << " must be < " << PP_MM + 1 << "\n";
		//
		return false;
	}

	for (int k = 0; k < non; k++) {
		int i, j; 

		if (fscanf(stream, "%d%d%f", &i, &j, &buf) < 3) {
			cout << "Unexpected end of file'" << mtxFile << "'." << endl;	return false;
		}

		i -= 1;
		j -= 1;
		if (i < 0) { cout << "Negative row index in'" << mtxFile << "'.\n" << endl;	return false; }
		if (j < 0) { cout << "Negative column index in'" << mtxFile << "'.\n" << endl;	return false; }
		PD_A[i][j] = buf;
		PD_A[i + m_cur][j] = -buf;
	}
	
	m_cur += nor;

	fclose(stream);

	//--------------- Reading b ------------------
	PD_MTX_File_b = PP_PATH;
	PD_MTX_File_b += PP_MTX_PREFIX;
	PD_MTX_File_b += PP_MTX_PROBLEM_NAME;
	PD_MTX_File_b += PP_MTX_POSTFIX_B;
	mtxFile = PD_MTX_File_b.c_str();
	stream = fopen(mtxFile, "r+b");

	if (stream == NULL) {
		cout << "Failure of opening file '" << mtxFile << "'.\n";
		return false;
	}

	SkipComments(stream);
	if (fscanf(stream, "%d%d", &nor, &noc) <2) {
		cout << "Unexpected end of file'" << mtxFile << "'." << endl;
		return false;
	}
	if (noe != nor) { cout << "Incorrect number of rows in'" << mtxFile << "'.\n"; return false; }
	if (noc != 1) { cout << "Incorrect number of columnws in'" << mtxFile << "'.\n"; return false; }

	for (int i = 0; i < noe; i++) {
		if (fscanf(stream, "%f", &buf) < 1) { 
			cout << "Unexpected end of file '" << mtxFile << "'." << endl; return false; 
		}
		PD_b[i] = buf;
		PD_b[i + noe] = -buf;
	}
	fclose(stream);

	//--------------- Reading lo ------------------
	PD_MTX_File_lo = PP_PATH;
	PD_MTX_File_lo += PP_MTX_PREFIX;
	PD_MTX_File_lo += PP_MTX_PROBLEM_NAME;
	PD_MTX_File_lo += PP_MTX_POSTFIX_LO; 
	mtxFile = PD_MTX_File_lo.c_str();
	stream = fopen(mtxFile, "r+b");

	if (stream == NULL) {
		cout << "Failure of opening file '" << mtxFile << "'.\n";
		return false;
	}

	SkipComments(stream);
	if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
		cout << "Unexpected end of file'" << mtxFile << "'." << endl;
		return false;
	}
	if (nor != PD_n) { cout << "Incorrect number of rows in'" << mtxFile << "'.\n"; return false; }
	if (noc != 1) { cout << "Incorrect number of columnws in'" << mtxFile << "'.\n"; return false; }

	for (int j = 0; j < PD_n; j++) {
		if (fscanf(stream, "%f", &buf) < 1) { cout << "Unexpected end of file '" << mtxFile << "'." << endl; return false; }
		if (buf != 0) { cout << "Non-zero lower bound in'" << mtxFile << "'.\n"; return false; }
		PD_A[m_cur + j][j] = -1;
	}
	fclose(stream);

	//--------------- Reading hi ------------------
	char s[6];
	PD_MTX_File_hi = PP_PATH;
	PD_MTX_File_hi += PP_MTX_PREFIX;
	PD_MTX_File_hi += PP_MTX_PROBLEM_NAME;
	PD_MTX_File_hi += PP_MTX_POSTFIX_HI;
	mtxFile = PD_MTX_File_hi.c_str();
	stream = fopen(mtxFile, "r+b");

	if (stream == NULL) {
		cout << "Failure of opening file '" << mtxFile << "'.\n";
		return false;
	}

	SkipComments(stream);
	if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
		cout << "Unexpected end of file'" << mtxFile << "'." << endl;
		return false;
	}
	if (nor != PD_n) { cout << "Incorrect number of rows in'" << mtxFile << "'.\n"; return false; }
	if (noc != 1) { cout << "Incorrect number of columnws in'" << mtxFile << "'.\n"; return false; }

	for (int j = 0; j < PD_n; j++) {
		if (fscanf(stream, "%s", s) < 1) { cout << "Unexpected end of file '" << mtxFile << "'." << endl; return false; }
		if (s[0] != '1' && s[1] != 'e' && s[2] != '3' && s[3] != '0' && s[4] != '8') { 
			cout << "Unexpected higher bound in'" << mtxFile << "'.\n"; return false; 
		}
	}
	fclose(stream);

	//--------------- Reading c ------------------
	PD_MTX_File_c = PP_PATH;
	PD_MTX_File_c += PP_MTX_PREFIX;
	PD_MTX_File_c += PP_MTX_PROBLEM_NAME;
	PD_MTX_File_c += PP_MTX_POSTFIX_C;
	mtxFile = PD_MTX_File_c.c_str();
	stream = fopen(mtxFile, "r+b");

	if (stream == NULL) {
		cout << "Failure of opening file '" << mtxFile << "'.\n";
		return false;
	}

	SkipComments(stream);
	if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
		cout << "Unexpected end of file'" << mtxFile << "'." << endl;
		return false;
	}
	if (nor != PD_n) { cout << "Incorrect number of rows in'" << mtxFile << "'.\n"; return false; }
	if (noc != 1) { cout << "Incorrect number of columnws in'" << mtxFile << "'.\n"; return false; }


	for (int j = 0; j < PD_n; j++) {
		if (fscanf(stream, "%f", &buf) < 0) { cout << "Unexpected end of file" << endl; return false; }
		PD_c[j] = -buf;
	}
	fclose(stream);
	return true;
}

static bool SaveSolution(PT_vector_T x, const char* filename) {
	FILE* stream;

	stream = fopen(filename, "w");
	if (stream == NULL) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Failure of opening file '" << filename << "'.\n";
		return false;
	}

	fprintf(stream, "%d\n", PD_n);

	for (int j = 0; j < PD_n; j++)
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

	for (int j = 0; j < PD_n; j++) { 
		projection[j] = factor * a[j];
	}

	return true;
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

inline void ProblemOutput(double elapsedTime) {
	cout << "=============================================" << endl;
	cout << "Elapsed time: " << elapsedTime << endl;
	cout << "Iterations: " << BSF_sv_iterCounter << endl;
	//cout << "Optimal objective value: " << setprecision(PP_SETW) << setw(PP_SETW) << ObjectiveF(PD_basePoint) << endl;
	cout << "Optimal objective value: " << setw(PP_SETW) << ObjectiveF(PD_basePoint) << endl;
	cout << "=============================================" << endl;
	fclose(PD_traceStream);
	const char* solutionFile = PD_solutionFile.c_str();
	if (SaveSolution(PD_basePoint, solutionFile))
		cout << "Solution is saved into the file '" << solutionFile << "'." << endl;
	cout << "Solution: ";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++) cout << setw(PP_SETW) << PD_basePoint[j];
	if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
	cout << endl;
}

inline void SkipComments(FILE* stream) {
	fpos_t pos;	// Position in the input stream
	int res;
	res = fscanf(stream, "\n");
	fgetpos(stream, &pos);
	while (getc(stream) == '%') {
		while (getc(stream) != 10);
		res = fscanf(stream, "\n");
		fgetpos(stream, &pos);
	};
	fsetpos(stream, &pos);
};

inline PT_float_T RndPositiveValue(PT_float_T rndMax) {
	return (((PT_float_T)rand() / ((PT_float_T)RAND_MAX + 1)) * rndMax);
}

inline void RndPoint(PT_vector_T point, PT_float_T maxValue) {
	for (int i = 0; i < PD_n; i++)
		point[i] = RndPositiveValue(maxValue);
}

inline void WriteTrace(PT_vector_T x) {
	for (int j = 0; j < PD_n; j++)
		fprintf(PD_traceStream, "%f\t", x[j]);
	fprintf(PD_traceStream, "\n");
}

void PC_bsf_MainArguments(int argc, char* argv[]) {

	if (argc > 1) {
		PD_lppFile = argv[1];
		if (argc > 2) {
			PD_solutionFile = argv[2];
			if (argc > 2)
				PD_traceFile = argv[3];
		}
	}
	else {
		PD_lppFile = PP_PATH;
		PD_lppFile += PP_LPP_FILE;
		PD_solutionFile = PP_PATH;
		PD_solutionFile += PP_SOLUTION_FILE;
		PD_traceFile = PP_PATH;
		PD_traceFile += PP_TRACE_FILE;
	}
}