/*==============================================================================
Project: LiFe
Theme: Apex Method
Module: Problem-Parameters.h (Problem Parameters)
Prefix: PP
Author: Leonid B. Sokolinsky
This source code has been produced with using BSF-skeleton
==============================================================================*/
//-------------------------- Compilation Modes -----------------------
//#define PP_DEBUG
//#define PP_PATH "Problems/"
#define PP_FILE_INI "config.ini"
static std::string PP_PATH;
//------------------------------------------------------------------/**/
static int PP_M;		// Number of equations (number of rows in *.mtx)
static int PP_N;		// Number of variables (number of cols in *.mtx)

#define PP_MAX_M 61
#define PP_MAX_N 81
#define PP_MAX_MM 365
//------------------------------------------------------------------
#define EPS							1E-8
#define PP_EPS_DIR_LENGTH			EPS		// Minimal Length of Direction Vector
#define PP_EPS_OBJ					1E-1	// |F(u) - PP_EXACT_OBJ_VALUE| <= PP_EPS_OBJ
#define PP_EPS_SHIFT				EPS		// Minimal shift to stop motion
#define PP_EPS_ZERO_COMPARE			EPS		// Zero comparison precision
#define PP_EPS_ZERO_DIR				1E-5	// Accuracy of zero value for determining direction vector coordinates
#define PP_EPS_DISTANCE				1.0		// Minimal distance between points in a trace
#define PP_EXACT_OBJ_VALUE			50000	// Exact maximum value of objective function
#define PP_INFINITY					1E+308	// Highest bound in *_hi.mtx
#define PP_SIGMA_TO_APEX			0	// Distance from apex base to apex point
#define PP_LOW_COST_PERCENTILE		0.01	// Percentile for low cost variable (must be in [0,1])
#define PP_GAP						1E-2	// Maximum gap from polytope surface (not too small!!!)
#define PP_OBJECTIVE_VECTOR_LENGTH	0.1		// Length of Objective Vector
#define PP_START_SHIFT_LENGTH		5		// Start length of shift vector
#define PP_MODE_BLOCK_HCV_VARIABLE	true	// block high cost variables
#define PP_MODE_USE_LCV_VARIABLE	false	// utilize low cost and zero variables
//------------------------------------------------------------------/**/
//================================ Common Paramrters ===========================
static int PP_MM;
#define PP_MAX_NUM_SHIFTS_SAME_LENGTH	5 // Maximal number of shifts with the same length
#define PP_MAX_ITER_COUNT				40000 // Maximal count of iterations
static int PP_ADD_FLAG;
//#define PP_RANDOM_X0
//-------------------------- Input/Outpoot Parameters ---------------------------
#define PP_OUTPUT_LIMIT	8	// Number of Elements to output
#define PP_TRACE_LIMIT 50	// Maximum points in a problem trace
#define PP_SETW 18
//
//------------------------- Matrix format ----------------
static std::string PP_PROBLEM_NAME;
static std::string PP_MTX_PREFIX;
static std::string PP_MTX_POSTFIX_A;
static std::string PP_MTX_POSTFIX_B;
static std::string PP_MTX_POSTFIX_LO;
static std::string PP_MTX_POSTFIX_HI;
static std::string PP_MTX_POSTFIX_C;
static std::string PP_MTX_POSTFIX_X0;
static std::string PP_MTX_POSTFIX_TR;
static std::string PP_MTX_POSTFIX_SO;
//
//-------------------------- Jobs  -----------------------
#define PP_JOB_PSEUDOPOJECTION	0 
#define PP_JOB_CHECK			1		
//-------------------------- Process states --------------------------
#define PP_STATE_START						0
#define PP_STATE_FIND_INTERIOR_POINT		1
#define PP_STATE_FIND_INITIAL_APPROXIMATION	2
#define PP_STATE_DETERMINE_DIRECTION		3
#define PP_STATE_MOVING_ALONG_SURFACE		4
#define PP_STATE_LANDING					5

//------------- Vector_ProjectOnHalfspace exit codes -------------
#define PP_EXITCODE_DEGENERATE_INEQUALITY		0
#define PP_EXITCODE_POINT_BELONGS_TO_HALFSPACE	1
#define PP_EXITCODE_SUCCESSFUL_PROJECTING		2