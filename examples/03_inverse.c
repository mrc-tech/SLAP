/*
	calculate the inverse of a square matrix
 */
#include <stdio.h>

#include "../include/SLAP.h"


int main()
{
	mat *A = mat_init(3,3, (double[]){1,3,3,4,5,6,7,8,9});
	mat *b = mat_init(3,1, (double[]){1,1,1});
	mat *x;
	
	printf("A =\n"); mat_print(A);
	printf("b =\n"); mat_print(b);
	printf("Solve  A * x = b  for x\n\n");
	
	printf("LU(P) decomposition:\n");
	mat_lup *lu = mat_lup_solve(A);
	x = solve_lu(lu, b);
	printf("x =\n"); mat_print(x);
	
//	A = mat_init(2,2, (TYPE[]){4,1,1,3});
//	b = mat_init(2,1, (TYPE[]){1,2});
	printf("\nConjugate gradient:\n");
	x = mat_conjgrad(A, b);
	printf("x =\n"); mat_print(x);
	
	return 0;

/*qr(A) =
   -8.1240   -9.7242  -11.0782
         0    1.8546    1.7647
         0         0   -0.3982
*/
	printf("\nQR decomposition:\nA=Q*R\n");
	mat_qr *qr = mat_qr_solve(A);
	mat_print(qr->Q);
	mat_print(qr->R);
	
	
	printf("Example (Wikipedia):\n");
	A = mat_init(3,3, (TYPE[]){12,-51,4,6,167,-68,-4,24,-41});
	qr = mat_qr_solve(A);
	printf("Q =\n"); mat_print(qr->Q);
	printf("R =\n"); mat_print(qr->R);
	
	
	printf("\n\n\nGauss elimination of A:\n");
	mat_print(mat_GaussJordan(A));
	
	
	return 0;
}
