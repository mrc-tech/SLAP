#include "../include/SLAP.h"


int main()
{
	mat *A, *x, *b;
	mat_lup *lu;
	mat *tmp;
	
	printf("solve_lu : ");
	A = mat_init(3,3, (TYPE[]){1,3,3,4,5,6,7,8,9});
	b = mat_init(3,1, (TYPE[]){1,1,1});
	lu = mat_lup_solve(A);
	x = solve_lu(lu, b);
	tmp = mat_init(3,1, (double[]){-0.5,0,0.5}); // solution
	if(mat_equal(x,tmp,1e-15)) printf("passed\n"); else printf("NOT PASSED!!!\n");
	A = mat_init(3,3, (TYPE[]){1,2,3,4,5,6,7,8,9}); // nearly singular (UNSTABLE MATRIX)
	lu = mat_lup_solve(A);
	x = solve_lu(lu, b);
	tmp = mat_init(3,1, (double[]){-2,5,1.5}); // solution (MATLAB with RCOND = 1.541976e-18)
	if(mat_equal(x,tmp,1e-15)) printf("unstable-> passed\n"); else printf("unstable-> NOT PASSED!!!\n");
//	mat_print(x);
	
	
	
	mat_free(A); mat_free(x); mat_free(b); mat_lup_free(lu);
	getch();
	
	return 0;
}


/* fortran code for LINPACKD benchmark

	subroutine matgen(n,a,lda)
	real a(lda,*)
	init = 1325
	do 10 j = 1,n
		do 20 i = 1,n
			init = mod(3125*init,65536)
			a(i,j) = (init - 32768.0)/16384.0
20		continue
10	continue
	end

The period of this pseudo–random number generator is: 214 = 16, 384
*/
