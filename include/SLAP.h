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

// OPERATIONS:
#include "core/strmod.h" // matrix structure modification
#include "solvers/Gauss.h"
#include "solvers/LUp.h" // LU(P) factorialization with partial pivoting
#include "solvers/QR.h" // QR decomposition
#include "eigen/eigen_qr.h" // Eigen-analysis using QR decomposition



// UTILITIES:
#include "utils/utils.h"
#include "utils/random.h"
#include "utils/constant.h" // scientific constants

