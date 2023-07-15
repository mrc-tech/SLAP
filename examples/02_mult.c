/*
	matrix multiplication
*/
#include <stdio.h>

#include "../include/SLAP.h"


int main()
{
	mat *A, *B, *u, *v, *tmp;
	
	A = mat_new(3,3);
	A = mat_init2(3,3, (double[]){1,2,3,4,5,6,7,8,9});
	B = mat_new(3,3);
	B = mat_init2(3,3, (double[]){1,0,0,0,1,0,0,0,1}); // identity matrix
	
	printf("A = \n");
	mat_print(A);
	printf("B = \n");
	mat_print(B);
	printf("2.5 * A =\n");
	mat_print(mat_smul(A, 2.5));
	
	printf("A * B =\n");
	mat_print(mat_mul(A,B));
	printf("A * A =\n");
	mat_print(mat_mul(A,A));
	
	u = mat_new(1,3); u = mat_init2(1,3, (double[]){1,2,3});
	v = mat_new(1,3); v = mat_init2(1,3, (double[]){2,2,2});
	printf("u = "); mat_print(u);
	printf("v = "); mat_print(v);
	printf("u * v^T = "); mat_print(mat_mul(u, mat_transpose(v)));
	printf("u^T * v=\n"); mat_print(mat_mul(mat_transpose(u), v));

	
	mat_free(A);
	mat_free(B);
	mat_free(tmp);
	
	return 0;
}
