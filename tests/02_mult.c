/*
	matrix multiplication
*/
#include <stdio.h>

#include "../include/SLAP.h"


int main()
{
	matd *A, *B, *u, *v, *tmp;
	
	A = new_matd(3,3);
	A = matd_init(3,3, 1.,2.,3.,4.,5.,6.,7.,8.,9.);
	B = new_matd(3,3);
	B = matd_init(3,3, 1.,0.,0.,0.,1.,0.,0.,0.,1.); // identity matrix
	
	printf("A = \n");
	print_mat(A);
	printf("B = \n");
	print_mat(B);
	printf("2.5 * A =\n");
	print_mat(matd_smul(A, 2.5));
	
	printf("A * B =\n");
	print_mat(matd_mul(A,B));
	printf("A * A =\n");
	print_mat(matd_mul(A,A));
	
	u = new_matd(1,3); u = matd_init(1,3, 1.,2.,3.);
	v = new_matd(1,3); v = matd_init(1,3, 2.,2.,2.);
	printf("u = "); print_mat(u);
	printf("v = "); print_mat(v);
	printf("u * v^T = "); print_mat(matd_mul(u, matd_transpose(v)));
	printf("u^T * v=\n"); print_mat(matd_mul(matd_transpose(u), v));

	
	free_mat(A);
	free_mat(B);
	free_mat(tmp);
	
	return 0;
}
