/*
	manipulations for matrices
 */
#include <stdio.h>

#include "../include/SLAP.h"

int main()
{
	mat *A, *B, *tmp;
	printf("A = \n");
	A = mat_fromfilename("matrix.txt");
	mat_print(A);
	printf("matd_getcol(A,1) = \n");
	tmp = mat_getcol(A,1);
	mat_print(tmp);
	printf("matd_getrow(A,1) = \n");
	tmp = mat_getrow(A,1);
	mat_print(tmp);
	printf("matd_getrow(A,2) = \n");
	tmp = mat_getrow(A,2);
	mat_print(tmp);
	
	printf("\n");
	A = mat_remcol(A,1);
	A = mat_remrow(A,1);
	mat_print(A);
	mat *mat_arr[2];
	mat_arr[0] = mat_copy(A);
	mat_arr[1] = mat_copy(A);
	B = mat_cathor(2,mat_arr); // ERROREEE!!!!!!!!!!!!!!!!!!!!
	mat_print(B);
	B = mat_catver(2,mat_arr);
	mat_print(B);
	
	
	return 0;
}
