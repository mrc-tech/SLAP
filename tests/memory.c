

#define SLAP_DEBUG 1 // debug analysis

#include "../include/SLAP.h"


void transpose();
void QR();



int main()
{
	
//	transpose();
	QR();
	
	
	return 0;
}

// -------------------------------------------------------------


void transpose()
{
	mat *A = mat_init(3,3, (double[]){1,2,3,4,5,6,7,8,9});
	mat_print(mat_transpose(A)); // MEMORY LEAKAGE!!!
	mat_transpose_r(A);
	mat_print(A);
	mat_free(A);
}


void QR()
{
	mat *A = mat_init(3,3, (double[]){1,3,3,4,5,6,7,8,9});
	mat_qr *qr;
	
	qr = mat_qr_solve(A);
	
	mat_qr_free(qr);
	mat_free(A);
}
