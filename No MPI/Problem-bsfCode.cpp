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
	PD_ObjectiveVectorLength = PP_OBJECTIVE_VECTOR_LENGTH;	// Initiat length of objective vector
	PD_listSize = PP_MM;

	//
	*success = LoadMatrixFormat();
	if (*success == false)
		return;

	*success = OpenTraceFile();
	if (*success == false)
		return;

	PD_MTX_File_so = PP_PATH;
	PD_MTX_File_so += PP_MTX_PREFIX;
	PD_MTX_File_so += PD_problemName;
	PD_MTX_File_so += PP_MTX_POSTFIX_SO;

	ObjectiveUnitVector(PD_objectiveUnitVector);
	Vector_MultiplyByNumber(PD_objectiveUnitVector, PD_ObjectiveVectorLength, PD_objectiveVector);
}

void PC_bsf_SetListSize(int* listSize) {
	*listSize = PD_listSize;
}

void PC_bsf_SetMapListElem(PT_bsf_mapElem_T* elem, int i) {
	elem->a = PD_A[i];
	elem->b = &(PD_b[i]);
}

// 0. Pseudo-pojection
void PC_bsf_MapF(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T* reduceElem, int* success // 1 - reduceElem was produced successfully; 0 - otherwise
) {
	*success = Vector_ProjectOnHalfspace(BSF_sv_parameter.x, mapElem->a, *mapElem->b, reduceElem->projection);
	reduceElem->pointIn = PointInHalfspace_s(BSF_sv_parameter.x, mapElem->a, *mapElem->b);
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
	z->pointIn = x->pointIn && y->pointIn;
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
	int reduceCounter, // Number of successfully produced Elements of Reduce List
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

	PD_pointIn = reduceResult->pointIn;
	Vector_Relaxation(reduceResult->projection, reduceCounter, PD_relaxationVector);
	Vector_PlusEquals(parameter->x, PD_relaxationVector);
}

// 1. Movement on Polytope  ========================================================
void PC_bsf_ProcessResults_1(
	PT_bsf_reduceElem_T_1* reduceResult,
	int reduceCounter, // Number of successfully produced Elements of Reduce List
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
	int reduceCounter, // Number of successfully produced Elements of Reduce List
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
	int reduceCounter, // Number of successfully produced Elements of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* nextJob,
	bool* exit // "true" if Stopping Criterion is satisfied, and "false" otherwise
) {
	// not used
}

void PC_bsf_JobDispatcher(
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* job,
	bool* exit,
	double t
) {
	const char* traceFile = PD_MTX_File_tr.c_str();
	const char* x0_File = PD_MTX_File_x0.c_str();
	static int indexToBlock;
	static double prevLandingObjVal = INFINITY;
	int sign;

	parameter->i = 0;

	switch (PD_state) {
	case PP_STATE_START://-------------------------- Start -----------------------------
		PD_listSize = PD_m;
		if (PointInPolytope_s(PD_basePoint)) {
			Vector_Copy(PD_objectiveUnitVector, PD_direction);
			PD_shiftLength = PP_START_SHIFT_LENGTH;
			// Preparations for moving inside the polytope
			PD_numShiftsSameLength = 0;
			Shift(PD_basePoint, PD_direction, PD_shiftLength, parameter->x);
			if (PointInPolytope_s(parameter->x)) {
				*job = PP_JOB_CHECK_S;
				PD_state = PP_MOVE_INSIDE_POLYTOPE;
			}
			else {
#ifdef PP_DEBUG
				cout << "\t\t\tu0 =";
				for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
					cout << setw(PP_SETW) << PD_basePoint[PD_objI[j]];
				if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
				cout << "\tF(t) = " << setw(PP_SETW) << ObjectiveF(PD_basePoint);
				cout << endl;
#endif // PP_DEBUG
				// Preparations for determining direction
				Vector_Copy(PD_basePoint, parameter->x);
				Vector_PlusEquals(parameter->x, PD_objectiveVector);
				*job = PP_JOB_PSEUDOPOJECTION;
				PD_state = PP_STATE_DETERMINE_DIRECTION;
				PD_numDetDir = 0;
#ifdef PP_DEBUG
				cout << "--------- Determine Direction ------------\n";
#ifdef PP_PAUSE
				system("pause");
#endif 
#endif
			}
		}
		else {
			// Preparations for finding a start point
			PD_state = PP_STATE_FIND_START_POINT;
			*job = PP_JOB_PSEUDOPOJECTION;
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
			cout << "Point in! Sift = " << setw(PP_SETW) << PD_shiftLength << "\tt = ";
			for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
				cout << setw(PP_SETW) << parameter->x[PD_objI[j]];
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

		/*debug* WriteTrace(parameter->x);/**/

		// Preparations for determining direction
		Vector_Copy(parameter->x, PD_basePoint);
		Vector_PlusEquals(parameter->x, PD_objectiveVector);
		*job = PP_JOB_PSEUDOPOJECTION;
		PD_state = PP_STATE_DETERMINE_DIRECTION;
		PD_numDetDir = 0;
#ifdef PP_DEBUG
		cout << "--------- Determine Direction ------------\n";
#ifdef PP_PAUSE
		system("pause");
#endif 
#endif 
		break;
	case PP_STATE_DETERMINE_DIRECTION://------------------------- Determine Direction -----------------------------
		if (Vector_Norm(PD_relaxationVector) >= PP_EPS_RELAX)
			return;
		if (!PointInPolytope_s(parameter->x))
			return;

#ifdef PP_DEBUG
		//		static int counterIter;
		//		if (counterIter % PP_BSF_TRACE_COUNT == 0) {
		//			cout << "Iter # " << counterIter << "\tTime: " << round(t) << endl;
		cout << "\t\t\tu = ";
		for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
			cout << setw(PP_SETW) << parameter->x[PD_objI[j]];
		if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
		cout << "\tF(u) = " << setw(PP_SETW) << ObjectiveF(parameter->x);
		cout << endl;
		//		}
		//		counterIter++;
				/**cout << "\t\t\tw = ";
				for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
					cout << setw(PP_SETW) << parameter->x[j];
				if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
				cout << "\tF(w) = " << ObjectiveF(parameter->x);
				cout << endl;/**/
#endif // PP_DEBUG

		if (fabs(ObjectiveF(parameter->x) - ObjectiveF(PD_basePoint)) < PP_EPS_OBJECTIVE) {
			cout << setw(PP_SETW) << "F(u) = " << ObjectiveF(PD_basePoint) << " == F(w) = " << ObjectiveF(parameter->x) << "\n";
			cout << "Maybe, you should decreas PP_EPS_OBJECTIVE.\n";
			*exit = true;
			return;
		}

		PD_numDetDir++;
		if (PD_numDetDir > 2) {
			cout << "Number of Squential States 'Determine Direction' is greater than 2!\n";
			cout << "Maybe, you should increase PP_OBJECTIVE_VECTOR_LENGTH.\n";
			*exit = true;
			return;
		}

		if (ObjectiveF(parameter->x) <= ObjectiveF(PD_basePoint) - PP_EPS_OBJECTIVE) {
			cout << setw(PP_SETW) << "F(u) = " << ObjectiveF(PD_basePoint) << " >= F(w) = "
				<< ObjectiveF(parameter->x) << endl;
			cout << "Maybe, you should decreas PP_OBJECTIVE_VECTOR_LENGTH.\n";
			*exit = true;
			return;
		}

		// Preparations for motion
		PD_shiftLength = PP_START_SHIFT_LENGTH;
		if (!GetDirection(PD_basePoint, parameter->x, PD_direction)) {
			cout << "Direction is too small!\n";
			cout << "Maybe, you should increas PP_OBJECTIVE_VECTOR_LENGTH.\n";
			*exit = true;
			return;
		}

		Vector_EpsZero(PD_direction);

#ifdef PP_DEBUG //----------------------------------------------//
		cout << "\t\t\tD = ";									//
		for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)	//
			cout << setw(PP_SETW) << PD_direction[PD_objI[j]];	//
		if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";			//
		cout << endl;											//
#endif // PP_DEBUG ---------------------------------------------//

		if (PP_STRAIGHT_TRACE) {
			sign = (PD_c[PD_objI[indexToBlock]] >= 0 ? 1 : -1);
			if (sign * PD_direction[PD_objI[indexToBlock]] <= 0 && indexToBlock < PD_n - 1) {
				PD_A[PD_m][PD_objI[indexToBlock]] = -sign;
				PD_b[PD_m] = -sign * PD_basePoint[PD_objI[indexToBlock]];
				//
				//
				PD_m++;
				PD_listSize = PD_m;
#ifdef PP_DEBUG
				cout << "Variable " << indexToBlock << " is blocked.\n";
#endif // PP_DEBUG
				PD_numDetDir = 0;
				indexToBlock++;
				if (fabs(PD_direction[PD_objI[indexToBlock - 1]]) > PP_EPS_ZERO_COMPARE) {
					Vector_Copy(PD_basePoint, parameter->x);
					Vector_PlusEquals(parameter->x, PD_objectiveVector);
					return;
				}
			}
		}

		PD_numShiftsSameLength = 0;
		Shift(PD_basePoint, PD_direction, PD_shiftLength, parameter->x);
		*job = PP_JOB_CHECK;
		PD_state = PP_STATE_MOVE_AND_CHECK;
		/**#ifdef PP_DEBUG
				cout << "--------- Movement on surface ------------\n";
		#ifdef PP_PAUSE
				system("pause");
		#endif
		#endif // PP_DEBUG /**/
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
				cout << setw(PP_SETW) << parameter->x[PD_objI[j]];
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
#endif // PP_DEBUG/**/
		break;
	case PP_STATE_LANDING://-------------------------- Landing -----------------------------
		if (Vector_Norm(PD_relaxationVector) >= PP_EPS_RELAX)
			return;
		if (!PointInPolytope_s(parameter->x))
			return;

#ifdef PP_DEBUG
		cout << "Iter # " << BSF_sv_iterCounter << ". Elapsed time: " << round(t) << endl;
		cout << "\t\t\tu = ";
		for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
			cout << setw(PP_SETW) << parameter->x[PD_objI[j]];
		if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
		cout << "\tF(t) = " << setw(PP_SETW) << ObjectiveF(parameter->x);
		cout << endl;
#endif // PP_DEBUG /**/

		Vector_Copy(parameter->x, PD_basePoint);
		//WriteTrace(PD_basePoint);
		/*debug
		SavePoint(PD_basePoint, x0_File, t);
		/*end debug*/

		// Preparations for determining direction
		Vector_PlusEquals(parameter->x, PD_objectiveVector);

		*job = PP_JOB_PSEUDOPOJECTION;
		PD_state = PP_STATE_DETERMINE_DIRECTION;
#ifdef PP_DEBUG
		cout << "--------- Determining direction ------------\n";
#ifdef PP_PAUSE
		system("pause");
#endif 
#endif/**/
		break;
	case PP_STATE_FIND_START_POINT://-------------------------- Finding a start point -----------------------------
		/**if (Vector_Norm(PD_relaxationVector) >= PP_EPS_RELAX)
			return;/**/
		/**if (!PointInPolytope_s(parameter->x))
			return;/**/
		/**/if (!PD_pointIn)
			return;/**/
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
			//WriteTrace(PD_basePoint);
#ifdef PP_DEBUG
			cout << "\t\t\tu0 =";
			for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
				cout << setw(PP_SETW) << PD_basePoint[PD_objI[j]];
			if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
			cout << "\tF(t) = " << setw(PP_SETW) << ObjectiveF(PD_basePoint);
			cout << endl;
#endif // PP_DEBUG

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
#endif
#endif
		}
		break;
	default://------------------------------------- default -----------------------------------
		cout << "PC_bsf_JobDispatcher: Undefined state!" << endl;
		break;
	}
}

void PC_bsf_ParametersOutput(PT_bsf_parameter_T parameter) {
	cout << "=================================================== Target ====================================================" << endl;
	cout << "Problem name: " << PD_problemName << endl;
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
	cout << "Before conversion: m =\t" << PP_M << "\tn = " << PP_N << endl;
	cout << "After conversion:  m =\t" << PD_m << "\tn = " << PD_n << endl;
	cout << "Eps Relax:\t\t" << PP_EPS_RELAX << endl;
	cout << "Eps Min Dir Length:\t" << PP_EPS_DIR_LENGTH << endl;
	cout << "Eps Objective:\t\t" << PP_EPS_OBJECTIVE << endl;
	cout << "Eps Shift:\t\t" << PP_EPS_SHIFT << endl;
	cout << "Eps Zero Direction:\t" << PP_EPS_ZERO_DIR << endl;
	cout << "Eps Zero Compare:\t" << PP_EPS_ZERO_COMPARE << endl;
	cout << "Exact Obj Value:\t" << PP_EXACT_OBJ_VALUE << endl;
	cout << "Gap Max:\t\t" << PP_GAP << endl;
	cout << "Lambda:\t\t\t" << PP_LAMBDA << endl;
	cout << "Obj Vector Length:\t" << PP_OBJECTIVE_VECTOR_LENGTH << endl;
	cout << "Start Shift Lengt:\t" << PP_START_SHIFT_LENGTH << endl;
	cout << "Straight Trace:\t\t" << (PP_STRAIGHT_TRACE ? "true" : "false") << endl;

#ifdef PP_MATRIX_OUTPUT
	cout << "------- Matrix PD_A & Column PD_b -------" << endl;
	for (int i = 0; i < PD_m; i++) {
		cout << i << ")";
		for (int j = 0; j < PD_n; j++)
			cout << setw(PP_SETW) << PD_A[i][j];
		cout << "\t<=" << setw(PP_SETW) << PD_b[i] << endl;
	}
#endif // PP_MATRIX_OUTPUT

	cout << "Objective Function:\t";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++) cout << setw(PP_SETW) << PD_c[PD_objI[j]];
	cout << (PP_OUTPUT_LIMIT < PD_n ? "	..." : "") << endl;
	cout << "Starting point:\t\t";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++) cout << setw(PP_SETW) << PD_basePoint[PD_objI[j]];
	if (PP_OUTPUT_LIMIT < PD_n)
		cout << "	...";
	cout << "\tF(x) = " << ObjectiveF(PD_basePoint);
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
	cout << "Elapsed time: " << round(elapsedTime) << endl;
	cout << "Approximat. :";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++) cout << setw(PP_SETW) << parameter.x[PD_objI[j]];
	if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
	cout << "\tF(x) = " << setw(PP_SETW) << ObjectiveF(parameter.x) << endl;
	/* debug *SavePoint(parameter.x, PD_MTX_File_x0.c_str(), elapsedTime);/**/
}

// 1. Movement on Polytope
void PC_bsf_IterOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	/*cout << "------------------ 1. Movement on Polytope. Iter # " << BSF_sv_iterCounter << " ------------------" << endl;
	cout << "Elapsed time: " << round(elapsedTime) << endl;
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
	cout << "Sift Length = " << PD_shiftLength << endl;/**/
};

// 2. Movement iside Polytope
void PC_bsf_IterOutput_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	/* cout << "------------------ 2. Movement inside Polytope. Iter # " << BSF_sv_iterCounter << " ------------------" << endl;
	cout << "Elapsed time: " << round(elapsedTime) << endl;
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
	cout << "Sift Length = " << PD_shiftLength << endl;/**/
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
void PC_bsf_ProblemOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	ProblemOutput(t);
}

// 2. Movement inside Polytope
void PC_bsf_ProblemOutput_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	ProblemOutput(t);
}

void PC_bsf_ProblemOutput_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
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
inline PT_float_T Vector_DotProductSquare(PT_vector_T x, PT_vector_T y) {
	PT_float_T sum = 0;
	for (int j = 0; j < PD_n; j++)
		sum += x[j] * y[j];
	return sum;
}

inline void Vector_Relaxation(PT_vector_T sumOfProjections, int numberOfProjections, PT_vector_T relaxationVector) {
	for (int j = 0; j < PD_n; j++) {
		relaxationVector[j] = PP_LAMBDA * sumOfProjections[j] / (double)numberOfProjections;
	}
}

inline PT_float_T Vector_Norm(PT_vector_T x) {
	return sqrt(Vector_NormSquare(x));
}

inline PT_float_T Vector_NormSquare(PT_vector_T x) {
	PT_float_T sum = 0;

	for (int j = 0; j < PD_n; j++) {
		sum += x[j] * x[j];
	}
	return sum;
}

inline bool PointInHalfspace // If the point belongs to the Halfspace with prescigion of PD_Gap
(PT_vector_T point, PT_vector_T a, PT_float_T b) {
	return Vector_DotProductSquare(a, point) <= b + PP_GAP;
}

inline bool PointInPolytope(PT_vector_T x) { // If the point belongs to the polytope
	for (int i = 0; i < PD_m; i++)
		if (!PointInHalfspace(x, PD_A[i], PD_b[i]))
			return false;
	return true;
}

inline bool PointInHalfspace_s // If the point belongs to the Halfspace with prescigion of PP_EPS_ZERO_COMPARE
(PT_vector_T point, PT_vector_T a, PT_float_T b) {
	return Vector_DotProductSquare(a, point) <= b + PP_EPS_ZERO_COMPARE;
}

inline bool PointInPolytope_s(PT_vector_T x) { // If the point belongs to the polytope with prescigion of PP_EPS_ZERO_COMPARE
	for (int i = 0; i < PD_m; i++)
		if (!PointInHalfspace_s(x, PD_A[i], PD_b[i]))
			return false;
	return true;
}

inline void Shift(PT_vector_T basePoint, PT_vector_T direction, PT_float_T siftLength, PT_vector_T endPoint) {
	for (int j = 0; j < PD_n; j++)
		endPoint[j] = basePoint[j] + direction[j] * siftLength;
}

// Calculating unit vector of direction from startPoint to endPoint
inline bool GetDirection(PT_vector_T startPoint, PT_vector_T endPoint, PT_vector_T unitVector) {
	for (int j = 0; j < PD_n; j++) {
		unitVector[j] = endPoint[j] - startPoint[j];
	}
	double normOfUnitVector = Vector_Norm(unitVector);

	if (normOfUnitVector < PP_EPS_DIR_LENGTH)
		return false;

	for (int j = 0; j < PD_n; j++) {
		unitVector[j] /= normOfUnitVector;
	}
	return true;
}

inline void Vector_Copy(PT_vector_T fromPoint, PT_vector_T toPoint) { // toPoint = fromPoint
	for (int j = 0; j < PD_n; j++)
		toPoint[j] = fromPoint[j];
}

inline void Vector_PlusEquals(PT_vector_T equalVector, PT_vector_T plusVector) { // equalVector += plusVector
	for (int j = 0; j < PD_n; j++)
		equalVector[j] += plusVector[j];
}

inline void Vector_MinusEquals(PT_vector_T equalPoint, PT_vector_T minusVector) { // equalPoint += minusVector
	for (int j = 0; j < PD_n; j++)
		equalPoint[j] -= minusVector[j];
}

inline void Vector_Addition(PT_vector_T x, PT_vector_T y, PT_vector_T z) {  // z = x + y
	for (int j = 0; j < PD_n; j++)
		z[j] = x[j] + y[j];
}

inline void Vector_Subtraction(PT_vector_T x, PT_vector_T y, PT_vector_T z) {  // z = x - y
	for (int j = 0; j < PD_n; j++)
		z[j] = x[j] - y[j];
}

inline void Vector_MultiplyByNumber(PT_vector_T x, double r, PT_vector_T y) {  // y = r*x
	for (int j = 0; j < PD_n; j++)
		y[j] = x[j] * r;
}

inline void Vector_MultiplyEquals(PT_vector_T x, double r) {  // x = r*x
	for (int j = 0; j < PD_n; j++)
		x[j] *= r;
}

inline void Vector_DivideEquals(PT_vector_T x, double r) {  // x = x/r
	for (int j = 0; j < PD_n; j++)
		x[j] /= r;
}

inline void Vector_ResetToZero(PT_vector_T x) {  // x = 0
	for (int j = 0; j < PD_n; j++) x[j] = 0;
}

inline void Vector_DivideByNumber(PT_vector_T x, double r, PT_vector_T y) {  // x = x/r
	for (int j = 0; j < PD_n; j++)
		y[j] = x[j] / r;
}

inline void Vector_Round(PT_vector_T x) {
	double floorValue;
	double fractionalPart;
	double sign;
	double absValue;

	for (int j = 0; j < PD_n; j++) {
		if (fabs(x[j]) < PP_EPS_ZERO_DIR) {
			x[j] = 0;
			continue;
		}
		absValue = fabs(x[j]);
		sign = x[j] > 0 ? 1 : -1;
		floorValue = floor(absValue);
		fractionalPart = absValue - floorValue;
		if (1 - fractionalPart < PP_EPS_ZERO_DIR) {
			x[j] = sign * (floorValue + 1);
			continue;
		}
		if (fractionalPart < PP_EPS_ZERO_DIR)
			x[j] = sign * floorValue;
	}
}

inline void Vector_EpsZero(PT_vector_T x) { // If x[j] < PP_EPS_ZERO_DIR then x[j] = 0
	for (int j = 0; j < PD_n; j++)
		if (fabs(x[j]) < PP_EPS_ZERO_DIR)
			x[j] = 0;
}

inline void Vector_Unit(PT_vector_T vector) { // Calculating unit vector
	double normOfVector = Vector_Norm(vector);
	for (int j = 0; j < PD_n; j++) {
		vector[j] /= normOfVector;
	}
};

inline PT_float_T ObjectiveF(PT_vector_T x) {
	PT_float_T s = 0;
	for (int j = 0; j < PD_n; j++)
		s += PD_c[j] * x[j];
	return s;
}

static bool LoadMatrixFormat() {
	int nor,	// Number of matrix rows
		noc,	// Number of matrix columns
		non,	// Number of non-zero elements
		noe;	// Number of equations
	const char* mtxFile;
	FILE* stream;// Input stream
	float buf;

	//--------------- Reading A ------------------
	PD_MTX_File_A = PP_PATH;
	PD_MTX_File_A += PP_MTX_PREFIX;
	PD_MTX_File_A += PD_problemName;
	PD_MTX_File_A += PP_MTX_POSTFIX_A;
	mtxFile = PD_MTX_File_A.c_str();
	stream = fopen(mtxFile, "r+b");

	if (stream == NULL) {
		//
		cout
			<< "Failure of opening file '" << mtxFile << "'.\n";
		return false;
	}

	SkipComments(stream);
	if (fscanf(stream, "%d%d%d", &nor, &noc, &non) < 3) {
		//
		cout
			<< "Unexpected end of file " << mtxFile << endl;
		return false;
	}

	if (nor >= noc) {
		//
		cout
			<< "Number of rows m = " << nor << " must be < " << "Number of columns n = " << noc << "\n";
		return false;
	}

	if (noc != PP_N) {
		//
		cout
			<< "Invalid input data: PP_N must be = " << noc << "\n";
		return false;
	}

	if (nor != PP_M) {
		//
		cout
			<< "Invalid input data: PP_M must be = " << nor << "\n";
		return false;
	}

	PD_m = noe = nor;
	PD_n = noc;

	if (2 * nor + noc > PP_MM) {
		//
		cout
			<< "Invalid input data: number of inequalities m = " << 2 * nor + noc
			<< " must be < PP_MM + 1 =" << PP_MM + 1 << "\n";
		return false;
	}

	for (int k = 0; k < non; k++) {
		int i, j;

		if (fscanf(stream, "%d%d%f", &i, &j, &buf) < 3) {
			//	
			cout
				<< "Unexpected end of file'" << mtxFile << "'." << endl;
			return false;
		}

		i -= 1;
		j -= 1;
		if (i < 0) {
			//
			cout
				<< "Negative row index in'" << mtxFile << "'.\n" << endl;
			return false;
		}
		if (j < 0) {
			//
			cout
				<< "Negative column index in'" << mtxFile << "'.\n" << endl;
			return false;
		}
		PD_A[i][j] = buf;
	}

	fclose(stream);

	//--------------- Reading b ------------------
	PD_MTX_File_b = PP_PATH;
	PD_MTX_File_b += PP_MTX_PREFIX;
	PD_MTX_File_b += PD_problemName;
	PD_MTX_File_b += PP_MTX_POSTFIX_B;
	mtxFile = PD_MTX_File_b.c_str();
	stream = fopen(mtxFile, "r+b");

	if (stream == NULL) {
		//
		cout
			<< "Failure of opening file '" << mtxFile << "'.\n";
		return false;
	}

	SkipComments(stream);
	if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
		//
		cout
			<< "Unexpected end of file'" << mtxFile << "'." << endl;
		return false;
	}
	if (noe != nor) {
		//
		cout
			<< "Incorrect number of rows in'" << mtxFile << "'.\n";
		return false;
	}
	if (noc != 1) {
		//
		cout
			<< "Incorrect number of columnws in'" << mtxFile << "'.\n";
		return false;
	}

	for (int i = 0; i < noe; i++) {
		if (fscanf(stream, "%f", &buf) < 1) {
			//
			cout
				<< "Unexpected end of file '" << mtxFile << "'." << endl;
			return false;
		}
		PD_b[i] = buf;
	}
	fclose(stream);

	//--------------- Reading lo ------------------
	PD_MTX_File_lo = PP_PATH;
	PD_MTX_File_lo += PP_MTX_PREFIX;
	PD_MTX_File_lo += PD_problemName;
	PD_MTX_File_lo += PP_MTX_POSTFIX_LO;
	mtxFile = PD_MTX_File_lo.c_str();
	stream = fopen(mtxFile, "r+b");

	if (stream == NULL) {
		//
		cout
			<< "Failure of opening file '" << mtxFile << "'.\n";
		return false;
	}

	SkipComments(stream);
	if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
		//
		cout
			<< "Unexpected end of file'" << mtxFile << "'." << endl;
		return false;
	}
	if (nor != PD_n) {
		//
		cout
			<< "Incorrect number of rows in'" << mtxFile << "'.\n";
		return false;
	}
	if (noc != 1) {
		//
		cout
			<< "Incorrect number of columnws in'" << mtxFile << "'.\n";
		return false;
	}

	for (int j = 0; j < PD_n; j++) {
		if (fscanf(stream, "%f", &buf) < 1) {
			//
			cout
				<< "Unexpected end of file '" << mtxFile << "'." << endl;
			return false;
		}
		PD_lo[j] = buf;
	}

	fclose(stream);

	//--------------- Reading c ------------------
	PD_MTX_File_c = PP_PATH;
	PD_MTX_File_c += PP_MTX_PREFIX;
	PD_MTX_File_c += PD_problemName;
	PD_MTX_File_c += PP_MTX_POSTFIX_C;
	mtxFile = PD_MTX_File_c.c_str();
	stream = fopen(mtxFile, "r+b");

	if (stream == NULL) {
		//
		cout
			<< "Failure of opening file '" << mtxFile << "'.\n";
		return false;
	}

	SkipComments(stream);
	if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
		//
		cout
			<< "Unexpected end of file'" << mtxFile << "'." << endl;
		return false;
	}
	if (nor != PD_n) {
		//
		cout
			<< "Incorrect number of rows in'" << mtxFile << "'.\n";
		return false;
	}
	if (noc != 1) {
		//
		cout
			<< "Incorrect number of columnws in'" << mtxFile << "'.\n";
		return false;
	}

	for (int j = 0; j < PD_n; j++) {
		if (fscanf(stream, "%f", &buf) < 0) {
			//
			cout
				<< "Unexpected end of file" << endl;
			return false;
		}
		PD_c[j] = -buf;
	}
	fclose(stream);

	//--------------- Reading hi ------------------
	PD_MTX_File_hi = PP_PATH;
	PD_MTX_File_hi += PP_MTX_PREFIX;
	PD_MTX_File_hi += PD_problemName;
	PD_MTX_File_hi += PP_MTX_POSTFIX_HI;
	mtxFile = PD_MTX_File_hi.c_str();
	stream = fopen(mtxFile, "r+b");

	if (stream == NULL) {
		//
		cout
			<< "Failure of opening file '" << mtxFile << "'.\n";
		return false;
	}

	SkipComments(stream);
	if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
		//
		cout
			<< "Unexpected end of file'" << mtxFile << "'." << endl;
		return false;
	}
	if (nor != PD_n) {
		//
		cout
			<< "Incorrect number of rows in'" << mtxFile << "'.\n";
		return false;
	}
	if (noc != 1) {
		//
		cout
			<< "Incorrect number of columnws in'" << mtxFile << "'.\n";
		return false;
	}

	for (int j = 0; j < PD_n; j++) {
		if (fscanf(stream, "%f", &buf) < 1) {
			//
			cout
				<< "Unexpected end of file '" << mtxFile << "'." << endl;
			return false;
		}
		PD_hi[j] = buf;
	}
	fclose(stream);

	bool error = !Conversion();
	if (error) return false;

	SortObjVarI();

	//--------------- Reading x0 ------------------
	PD_MTX_File_x0 = PP_PATH;
	PD_MTX_File_x0 += PP_MTX_PREFIX;
	PD_MTX_File_x0 += PD_problemName;
	PD_MTX_File_x0 += PP_MTX_POSTFIX_X0;
	mtxFile = PD_MTX_File_x0.c_str();
	stream = fopen(mtxFile, "r+b");

	if (stream == NULL) {
		// Generating Coordinates of starting point
		for (int j = 0; j < PD_n; j++)
			PD_basePoint[j] = 0;
		return true;
	}

	SkipComments(stream);
	if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
		//
		cout
			<< "Unexpected end of file'" << mtxFile << "'." << endl;
		return false;
	}
	if (nor != PD_n) {
		//
		cout
			<< "Incorrect number of rows in'" << mtxFile << "'.\n";
		return false;
	}
	if (noc != 1) {
		//
		cout
			<< "Incorrect number of columnws in'" << mtxFile << "'.\n";
		return false;
	}

	for (int j = 0; j < PD_n; j++) {
		if (fscanf(stream, "%f", &buf) < 0) {
			//
			cout
				<< "Unexpected end of file" << endl;
			return false;
		}
		PD_basePoint[j] = buf;
	}
	fclose(stream);
	return true;
}

static bool Conversion() { // Transformation to inequalities & dimensionality reduction
	static PT_float_T iA[PP_M]; // Free variable coefficients
	static bool Flag[PP_N];		// Flags of free variables to delete
	bool single;

	for (int jc = 0; jc < PD_n; jc++) { // Detecting free variables
		if (PD_c[jc] != 0) continue;
		for (int i = 0; i < PD_m; i++) { // Find corresponding equation
			if (PD_A[i][jc] == 0) continue;
			single = true;
			for (int ii = i + 1; ii < PD_m; ii++)
				if (PD_A[ii][jc] != 0) {
					single = false;
					break;
				}
			if (!single) break;
			for (int j = jc + 1; j < PD_n; j++) {
				if (PD_A[i][j] != 0 && PD_c[jc] == 0) {
					single = false;
					break;
				}
			}
			if (!single) break;
			Flag[jc] = true;
			iA[i] = PD_A[i][jc];
			PD_A[i][jc] = 0;
		}
	}

	static bool PD_delete[PP_MM]; // Columns to delete
	PT_float_T s;

	for (int i = 0; i < PD_m; i++) { // Marking columns for deletion
		s = 0;
		for (int j = 0; j < PD_n; j++)
			s += fabs(PD_A[i][j]);
		if (s < PP_EPS_ZERO_COMPARE) {
			if (PD_b[i] != 0) {
				//
				cout
					<< "Inconsistent equation " << i << ": " << s << " = " << PD_b[i] << endl;
				return false;
				PD_delete[i] = true;
			}
		}
	}

	for (int i = 0; i < PD_m; i++) { // Removing null equations
		if (!PD_delete[i]) continue;
		for (int ii = i; ii < PD_m - 1; ii++) {  // Remove null equation
			for (int j = 0; j < PD_n; j++)
				PD_A[ii][j] = PD_A[ii + 1][j];
			iA[ii] = iA[ii + 1];
			PD_b[ii] = PD_b[ii + 1];
		}
		PD_m--;
	}

	for (int jc = 0; jc < PD_n; jc++) { // Delete free variables
		if (!Flag[jc]) continue;
		for (int j = jc; j < PD_n - 1; j++) { // Delete column
			PD_c[j] = PD_c[j + 1];
			Flag[j] = Flag[j + 1];
			for (int i = 0; i < PD_m; i++)
				PD_A[i][j] = PD_A[i][j + 1];
		}
		PD_n--;
		jc--;
		Flag[PD_n] = false;
		PD_c[PD_n] = 0;
		for (int i = 0; i < PD_m; i++)
			PD_A[i][PD_n] = 0;
	}

	int m = PD_m;
	for (int i = 0; i < m; i++) { // Conversion to inequalities

		if (iA[i] == 0) { // Equation without free variable => adding inequality.
			for (int j = 0; j < PD_n; j++)
				PD_A[PD_m][j] = -PD_A[i][j];
			PD_b[PD_m] = -PD_b[i];
			PD_m++;
		}
		else {
			if (iA[i] < 0) { // Free variable is negative => change sign to opposite.
				for (int j = 0; j < PD_n; j++)
					PD_A[i][j] = -PD_A[i][j];
				PD_b[i] = -PD_b[i];
			}
		}
	}

	for (int i = 0; i < PD_n; i++) { // Adding lower bound conditions
		for (int j = 0; j < PD_n; j++)
			PD_A[i + PD_m][j] = 0;
		PD_A[i + PD_m][i] = -1;
		PD_b[i + PD_m] = -PD_lo[i];
	}
	PD_m += PD_n;

	for (int i = 0; i < PD_n; i++) { // Adding higher bound conditions
		if (PD_hi[i] != INFINITY) {
			for (int j = 0; j < PD_n; j++)
				PD_A[PD_m][j] = 0;
			PD_A[PD_m][i] = 1;
			PD_b[PD_m] = PD_hi[i];
			PD_m++;
		}
	}

	return true;
}

static bool SavePoint(PT_vector_T x, const char* filename, double elapsedTime) {
	FILE* stream;

	stream = fopen(filename, "w");
	if (stream == NULL) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Failure of opening file '" << filename << "'.\n";
		return false;
	}

	fprintf(stream, "%c Elapsed time: %.0f\n%d %d\n", '%', elapsedTime, PD_n, 1);

	for (int j = 0; j < PD_n; j++)
		fprintf(stream, "%.14f\n", x[j]);

	fclose(stream);
	return true;
}

// Point projection onto Half-space <a,x> <= b
inline bool // true if the point does not belong to the half-space and false otherwise 
Vector_ProjectOnHalfspace(PT_vector_T point, PT_vector_T a, PT_float_T b, PT_vector_T projection) {
	double factor;
	double aNormSquare = Vector_NormSquare(a);

	if (sqrt(aNormSquare) < PP_EPS_ZERO_COMPARE)
		return false;

	factor = (b - Vector_DotProductSquare(point, a)) / aNormSquare;

	if (factor > PP_EPS_ZERO_COMPARE)
		return false;

	for (int j = 0; j < PD_n; j++) {
		projection[j] = factor * a[j];
	}

	return true;
}

inline PT_float_T Distance(PT_vector_T x, PT_vector_T y) {
	PT_vector_T z;
	Vector_Subtraction(x, y, z);
	return Vector_Norm(z);
}

inline void ObjectiveUnitVector(PT_vector_T objectiveUnitVector) { // Calculating Objective Unit Vector
	double c_norm = Vector_Norm(PD_c);
	Vector_DivideByNumber(PD_c, c_norm, objectiveUnitVector);
}

inline void ProblemOutput(double elapsedTime) {
	PT_float_T F_basePoint = ObjectiveF(PD_basePoint);
	cout << "=============================================" << endl;
	cout << "Elapsed time: " << elapsedTime << endl;
	cout << "Iterations: " << BSF_sv_iterCounter << endl;
	//cout << "Optimal objective value: " << setprecision(PP_SETW) << setw(PP_SETW) << ObjectiveF(PD_basePoint) << endl;
	cout << "Optimal objective value: " << setw(PP_SETW) << F_basePoint << endl;
	cout << "Exact objective value: " << setw(PP_SETW) << PP_EXACT_OBJ_VALUE << endl;
	cout << "Relative error = " << fabs(F_basePoint - PP_EXACT_OBJ_VALUE) / fabs(PP_EXACT_OBJ_VALUE) << endl;
	cout << "=============================================" << endl;
	fclose(PD_traceStream);
	const char* solutionFile = PD_MTX_File_so.c_str();
	Vector_EpsZero(PD_basePoint);
	if (SavePoint(PD_basePoint, solutionFile, elapsedTime))
		cout << "Solution is saved into the file '" << solutionFile << "'." << endl;
	cout << "Solution:\t";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++) cout << setw(PP_SETW) << PD_basePoint[j];
	if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
	cout << endl;
	cout << "Ordered:\t";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++) cout << setw(PP_SETW) << PD_basePoint[PD_objI[j]];
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

inline void SortObjVarI() { // Sorting objective variables in absolute descending order
	PT_float_T bigestAbsVal;
	int iBig;

	for (int j = 0; j < PD_n; j++)
		PD_objI[j] = j;

	for (int J = 0; J < PD_n - 1; J++) {
		bigestAbsVal = fabs(PD_c[PD_objI[J]]);
		iBig = J;
		for (int j = J + 1; j < PD_n; j++) {
			if (fabs(PD_c[PD_objI[j]]) > bigestAbsVal) {
				bigestAbsVal = fabs(PD_c[PD_objI[j]]);
				iBig = j;
			}
		}
		if (iBig != J) {
			int j;
			j = PD_objI[iBig];
			PD_objI[iBig] = PD_objI[J];
			PD_objI[J] = j;
		}
	}
}

inline void WriteTrace(PT_vector_T x) {
	for (int j = 0; j < PD_n; j++)
		fprintf(PD_traceStream, "%f\t", x[j]);
	fprintf(PD_traceStream, "\n");
}

void PC_bsf_MainArguments(int argc, char* argv[]) {

	if (argc > 1)
		PD_problemName = argv[1];
	else
		PD_problemName = PP_PROBLEM_NAME;
}

inline bool OpenTraceFile() {
	PD_MTX_File_tr = PP_PATH;
	PD_MTX_File_tr += PP_MTX_PREFIX;
	PD_MTX_File_tr += PD_problemName;
	PD_MTX_File_tr += PP_MTX_POSTFIX_TR;

	const char* traceFile = PD_MTX_File_tr.c_str();
	PD_traceStream = fopen(traceFile, "w");
	if (PD_traceStream == NULL) {
		cout << "Failure of opening file '" << PD_MTX_File_tr << "'.\n";
		return false;
	}
	return true;
}