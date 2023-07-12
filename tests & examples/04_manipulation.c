/*
	manipulations for matrices
 */
#include <stdio.h>

#include "../include/SLAP.h"

int main()
{
	matd *A, *B, *tmp;
	printf("A = \n");
	A = matd_fromfilename("matrix.txt");
	print_mat(A);
	printf("matd_getcol(A,1) = \n");
	tmp = matd_getcol(A,1);
	print_mat(tmp);
	printf("matd_getrow(A,1) = \n");
	tmp = matd_getrow(A,1);
	print_mat(tmp);
	printf("matd_getrow(A,2) = \n");
	tmp = matd_getrow(A,2);
	print_mat(tmp);
	
	printf("\n");
	A = matd_remcol(A,1);
	A = matd_remrow(A,1);
	print_mat(A);
	matd *mat_arr[2];
	mat_arr[0] = matd_copy(A);
	mat_arr[1] = matd_copy(A);
	B = matd_cathor(2,mat_arr); // ERROREEE!!!!!!!!!!!!!!!!!!!!
	print_mat(B);
	B = matd_catver(2,mat_arr);
	print_mat(B);
	
	
	return 0;
}
