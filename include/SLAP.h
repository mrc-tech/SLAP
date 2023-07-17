/*
	Simple Linear Algebra Package (SLAP)
*/
#include <stdlib.h>
#include <string.h>	// for memory menagement functions
#include <stdio.h>	// for error messages
#include <math.h>


// INITIAL DEFs:
#include "core/defs.h"

// CORE:
#include "core/datatype.h" // matrix and other data type and menagement
#include "core/basicops.h" // basic operations
#include "core/strmod.h" // matrix structure modification

// DECOMPOSITIONS:
#include "solvers/LUp.h" // LU(P) factorialization with partial pivoting
#include "solvers/QR.h" // QR decomposition

// OPERATIONS:
#include "fun/FUN.h" // additional functions (det, inv, etc.)

// SOLVERS:
#include "solvers/Gauss.h"

// EIGEN:
#include "eigen/eigen_qr.h" // Eigen-analysis using QR decomposition



// UTILITIES:
#include "utils/utils.h"
#include "utils/random.h"
#include "utils/constant.h" // scientific constants

