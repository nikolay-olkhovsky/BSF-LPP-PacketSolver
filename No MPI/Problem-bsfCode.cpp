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
		parameter->x[j] = PD_apex[j];
};

void PC_bsf_Init(bool* success) {
	PD_state = PP_STATE_START;
	

	// ------------- Load LPP data -------------------
	PD_lppFile = PP_PATH;
	PD_lppFile += PP_LPP_FILE;
	const char* lppFile = PD_lppFile.c_str();

	FILE* stream;
	float buf;

	stream = fopen(lppFile, "r");

	if (stream == NULL) {
		cout << "Failure of opening file '" << lppFile << "'.\n";
		*success = false;
		system("pause"); return;
	}

	if (fscanf(stream, "%d%d", &PD_m_init, &PD_n) == 0) { cout << "Unexpected end of file" << endl; *success = false; return; }
	
	if (PD_n > PP_MAX_N) {
		cout << "Invalid input data: Space dimension n = " << PD_n << " must be < " << PP_MAX_N + 1 << "\n";
		*success = false;
		system("pause");
		return;
	}

	if (PD_n >= PD_m_init) {
		cout << "Invalid input data: Space dimension n = " << PD_n 
			<< " must be less then number of inequalities m = " << PD_m_init << "\n";
		*success = false;
		system("pause");
		return;
	}

	for (int i = 0; i < PD_m_init; i++) {
		for (int j = 0; j < PD_n; j++) {
			if (fscanf(stream, "%f", &buf) == 0) { cout << "Unexpected end of file" << endl; *success = false; return; }
			PD_A[i][j] = buf;
		}
		if (fscanf(stream, "%f", &buf) == 0) { cout << "Unexpected end of file" << endl; *success = false; return; }
		PD_b[i] = buf;
	}

	for (int j = 0; j < PD_n; j++) {
		if (fscanf(stream, "%f", &buf) == 0) { cout << "Unexpected end of file" << endl; *success = false; return; }
		PD_c[j] = buf;
	}

	if (PD_c[0] != (PT_float_T)PD_n * PP_RHO) {
		cout << "Invalid input data!\n";
		*success = false;
		system("pause");
		return;
	}

	fclose(stream);

	PD_solutionFile = PP_PATH;
	PD_solutionFile += PP_SOLUTION_FILE;

	// Open trace file
	PD_traceFile = PP_PATH;
	PD_traceFile += PP_TRACE_FILE;
	const char* traceFile = PD_traceFile.c_str();
	PD_traceStream = fopen(traceFile, "w");
	if (stream == NULL) {
		cout << "Failure of opening file '" << PD_traceFile << "'.\n";
		*success = false;
		system("pause");
		return;
	}

	// Generating a point inside the polytope
	for (int j = 0; j < PD_n; j++)
		PD_basePoint[j] = (double)PP_SF / 2;

	// Generating Coordinates of Apex Point
	ObjectiveUnitVector(PD_objectiveUnitVector);
	Vector_MultiplyByNumber(PD_objectiveUnitVector, (PT_float_T)PP_DIST_TO_APEX, PD_apex);
	Vector_PlusEquals(PD_apex, PD_basePoint);

	/*if (PD_apex[PD_n - 1] < (PT_float_T)PP_SF * 1.1) {
		if (BSF_sv_mpiRank == 0)
			cout << "PC_bsf_Init.Error: PD_apex[PD_n - 1] = " << PD_apex[PD_n - 1] << " < PP_SF * 1.1 = " << PP_SF * 1.1 << endl;
		*success = false;
		return;
	}/**/
	PD_m = PP_MAX_M;
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
	const char* traceFile = PD_traceFile.c_str();
	static bool fistTime = true;

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

		if (fistTime) {
			StoreTracePoint(parameter->x);
			fistTime = false;
		}

#ifdef PP_DEBUG
		cout << "\t\t\tu = ";
		for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
			cout << setw(PP_SETW) << parameter->x[j];
		if (PP_OUTPUT_LIMIT < PD_n) cout << "	..." << setw(PP_SETW) << parameter->x[PD_n - 1];
		cout << "\tF(u) = " << ObjectiveF(parameter->x);
		cout << endl;
#endif // PP_DEBUG

		// Preparations for determining direction
		ObjectiveUnitVector(PD_objectiveUnitVector);
		Vector_MultiplyByNumber(PD_objectiveUnitVector, PP_OBJECTIVE_VECTOR_LENGTH, objectiveVector);
		Vector_Copy(parameter->x, PD_basePoint);
		Vector_PlusEquals(parameter->x, objectiveVector);
/*#ifdef PP_DEBUG
		cout << "\t\t\tv = ";
		for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
			cout << setw(PP_SETW) << parameter->x[j];
		if (PP_OUTPUT_LIMIT < PD_n) cout << "	..." << setw(PP_SETW) << parameter->x[PD_n - 1];
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
		fflush(PD_traceStream);
		cout << "\t\t\tw = ";
		for (int j = 0; j < PD_n; j++)
			cout << setw(PP_SETW) << parameter->x[j];
		cout << "\tF(w) = " << ObjectiveF(parameter->x);
		cout << endl;/**/
#endif // PP_DEBUG

		if (ObjectiveF(parameter->x) <= ObjectiveF(PD_basePoint) - PP_EPS_OBJECTIVE) {
#ifdef PP_DEBUG
			cout << setw(PP_SETW) << "F(u) = " << ObjectiveF(PD_basePoint) << " >= F(w) = " << ObjectiveF(parameter->x) << "\n";
			//system("pause");
#endif // PP_DEBUG
			* exit = true;
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
		for (int j = 0; j < PD_n; j++)					//
			cout << setw(PP_SETW) << PD_direction[j];	//
		cout << endl;									//
		//system("pause");								//
		cout << "-----------------------------------\n";//
#endif // PP_DEBUG -------------------------------------//

#ifdef PP_MAJOR_COORDINATES_CAN_NOT_DECREASE
		PD_newInequations = false;
		for (int j = 0; j < PD_n - 1; j++) {
			if (PD_direction[j] > 0)
				break;
			if (PD_direction[j] == 0)
				continue;
			PD_A[PD_m + 1][j] = -1;
			PD_b[PD_m + 1] = -PD_basePoint[j];
			PD_m++;
			PD_newInequations = true;
			break;
		}
		if (PD_newInequations) {
			PD_state = PP_STATE_FIND_BEGINNING_OF_PATH;
			*job = PP_JOB_PSEUDOPOJECTION;
			return;
		}
#endif

		WriteTrace(PD_tracePoint); // Trace!!! --------->>>

		PD_numShiftsSameLength = 0;
		Shift(PD_basePoint, PD_direction, PD_shiftLength, parameter->x);
		*job = PP_JOB_CHECK;
		PD_state = PP_STATE_MOVE_AND_CHECK;
		return;

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
			for (int j = 0; j < PD_n; j++)
				cout << setw(PP_SETW) << parameter->x[j];
			// cout << "\tF(t) = " << setw(PP_SETW) << ObjectiveF(parameter->x);
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
		if (!PointInPolytope_s(PD_basePoint))
			Vector_Copy(PD_basePoint, parameter->x);
		*job = PP_JOB_PSEUDOPOJECTION;
		PD_state = PP_STATE_LANDING;
		break;
	case PP_STATE_LANDING://-------------------------- Landing -----------------------------
		if (Vector_NormSquare(PD_relaxationVector) >= PP_EPS_RELAX * PP_EPS_RELAX)
			return;

		Vector_Copy(parameter->x, PD_basePoint);
		StoreTracePoint(PD_basePoint);

		// Preparations for determining direction
		Vector_MultiplyByNumber(PD_objectiveUnitVector, PP_OBJECTIVE_VECTOR_LENGTH, objectiveVector);
		Vector_PlusEquals(parameter->x, objectiveVector);

#ifdef PP_DEBUG
		cout << "\n\t\t\tu = ";
		for (int j = 0; j < PD_n; j++)
			cout << setw(PP_SETW) << PD_basePoint[j];
		cout << "\tF(u) = " << ObjectiveF(PD_basePoint);
		cout << endl;
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
	PD_m = PD_m_init;
#ifdef PP_MATRIX_OUTPUT
	cout << "------- Matrix PD_A & Column PD_b -------" << endl;
	for (int i = 0; i < PD_m; i++) {
		cout << i << ")";
		for (int j = 0; j < PD_n; j++)
			cout << setw(PP_SETW) << PD_A[i][j];
		cout << "\t<=" << setw(PP_SETW) << PD_b[i] << endl;
	}
#endif // PP_MATRIX_OUTPUT
	cout << "-------------------------------------------" << endl;
}

void PC_bsf_CopyParameter(PT_bsf_parameter_T parameterIn, PT_bsf_parameter_T* parameterOutP) {
	for (int j = 0; j < PD_n; j++) { 
		parameterOutP->x[j] = parameterIn.x[j];
	}
	if (parameterIn.i > 0) {
		for (int j = 0; j < PD_n; j++) {
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
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++) cout << setw(PP_SETW) << parameter.x[j];
	if (PP_OUTPUT_LIMIT < PD_n) cout << "	..." << setw(PP_SETW) << parameter.x[PD_n - 1];
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
	if(PP_OUTPUT_LIMIT < PD_n) cout << "	..." << setw(PP_SETW) << PD_basePoint[PD_n - 1];
	cout << endl;
	cout << "PD_direction:";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
		cout << setw(PP_SETW) << PD_direction[j];
	if(PP_OUTPUT_LIMIT < PD_n) cout << "	..." << setw(PP_SETW) << PD_direction[PD_n - 1];
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
	for (int j = 0; j < PD_n; j++) {
		sum += x[j] * y[j];
	}
	return sum;
}

inline void Vector_Relaxation(PT_vector_T sumOfProjections, int numberOfProjections, PT_vector_T relaxationVector) {
	for (int j = 0; j < PD_n; j++) {
		relaxationVector[j] = sumOfProjections[j] / (double)numberOfProjections;
	}
}

inline double Vector_NormSquare(PT_vector_T x) { 
	double sum = 0;

	for (int j = 0; j < PD_n; j++) 
		sum += x[j] * x[j];
	return sum;
}

inline bool PointInHalfspace // If the point belongs to the Halfspace with prescigion of PP_EPS_IN
(PT_vector_T point, PT_vector_T a, PT_float_T b) {
	return Vector_DotProductSquare(a, point) <= b + PP_EPS_IN;
}

inline bool PointInPolytope_s(PT_vector_T point) { // If the point belongs to the Polytope with prescigion of PP_EPS_ZERO
	for (int i = 0; i < PD_m; i++)
			if (Vector_DotProductSquare(point, PD_A[i]) > PD_b[i] + PP_EPS_ZERO)
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

inline void StoreTracePoint(PT_vector_T point) {
	for (int j = 0; j < PD_n; j++) {
		PD_tracePoint[j] = point[j];
	}
	Vector_Round(PD_tracePoint);
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

inline void WriteTrace(PT_vector_T x) {
	for (int j = 0; j < PD_n; j++)
		fprintf(PD_traceStream, "%f\t", x[j]);
	fprintf(PD_traceStream, "\n");
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

inline void ProblemOutput(double elapsedTime) {
	WriteTrace(PD_basePoint);
	for (int j = 0; j < PD_n; j++)
		cout << setw(PP_SETW) << PD_basePoint[j];
	cout << endl;
	fclose(PD_traceStream);
	cout << "=============================================" << endl;
	const char* solutionFile = PD_solutionFile.c_str();
	if (SaveSolution(PD_basePoint, solutionFile))
		cout << "Solution is saved into the file '" << PD_solutionFile << "'." << endl;
	cout << "Trace is saved into the file '" << PD_traceFile << "'." << endl;
	//system("pause");
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