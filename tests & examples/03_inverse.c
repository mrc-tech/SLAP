/*
	calculate the inverse of a square matrix
 */
#include <stdio.h>

#include "../include/SLAP.h"


int main()
{
//	mat *A = mat_init(3,3, 1.,3.,3.,4.,5.,6.,7.,8.,9.);
//	mat *b = mat_init(3,1, 1.,1.,1.);
	mat *A = mat_init2(3,3, (double[]){1,3,3,4,5,6,7,8,9});
	mat *b = mat_init2(3,1, (double[]){1,1,1});
	
	printf("A =\n"); mat_print(A);
	printf("b =\n"); mat_print(b);
	printf("Solve  A * x = b  for x\n\n");
	
	printf("LU(P) decomposition:\n");
	mat_lup *lu = mat_lup_solve(A);
	mat* x = solve_lu(lu, b);
	printf("x =\n"); mat_print(x);
	
	printf("\nQR decomposition:\nA=Q*R\n");
	mat_qr *qr = mat_qr_solve(A); // QUALCHE VOLTA DA ERRORE DI MEMORIA!!!!! inoltre il risultato va controllato bene!@!!!!!
	mat_print(qr->Q);
	mat_print(qr->R);
	
	printf("\n\n\nGauss elimination of A:\n");
	mat_print(mat_GaussJordan(A));
	
	
	return 0;
}
