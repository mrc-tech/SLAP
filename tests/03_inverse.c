/*
	calculate the inverse of a square matrix
 */
#include <stdio.h>

#include "../include/SLAP.h"


int main()
{
	matd *A = matd_init(3,3, 1.,3.,3.,4.,5.,6.,7.,8.,9.);
	matd *b = matd_init(3,1, 1.,1.,1.);
	
	printf("A =\n"); print_mat(A);
	printf("b =\n"); print_mat(b);
	printf("Solve  A * x = b  for x\n\n");
	
	printf("LU(P) decomposition:\n");
	matd_lup *lu = matd_lup_solve(A);
	matd* x = lu_solve(lu, b);
	printf("x =\n"); print_mat(x);
	
	
	printf("\n\n\nGauss elimination of A:\n");
	print_mat(matd_GaussJordan(A));
	
	
	return 0;
}
