/*
	calculate the inverse of a square matrix
 */
#include <stdio.h>

#include "../include/SLAP.h"


int main()
{
	matd *A = matd_init(3,3, 1.,2.,3.,4.,5.,6.,7.,8.,9.);
	matd *b = matd_init(3,1, 1.,1.,1.);
	
	matd_lup *lu = matd_lup_solve(A);
	matd* x = lu_solve(lu, b);
	
	printf("A =\n"); print_mat(A);
	printf("b =\n"); print_mat(b);
	printf("A * x = b\nx =\n"); print_mat(x); // ERROREEEEEE!!!!!!!!!!!!!!!!!! RISULTATO SBAGLIATO!!!!
	
	
	return 0;
}
