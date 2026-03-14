#include <stdio.h>

#include "../include/SLAP.h" // Simple Linear Algebra Package


int main()
{
	int i, j;
	mat *A, *B, *tmp;
	mat *u, *v;
	
	
	printf("mat_det       : ");
	A = mat_init(3,3, (TYPE[]){1,3,3,4,5,6,7,8,9});
	if(fabs(mat_det(A) - 6) < 1e-10) printf("passed\n"); else printf("NOT PASSED!!!\n");
	A = mat_init(3,3, (TYPE[]){1.,1/2.,1/3.,1/2.,1/3.,1/4.,1/3.,1/4.,1/5.}); // matrice di Hilbert, instabile (det = 1/2160)
	if(fabs(mat_det(A) - 1/2160.) < 1e-10) printf("              : passed\n"); else printf("              : NOT PASSED!!!\n");
	A = mat_init(3,3, (TYPE[]){1,1,1,1,1.0001,1,1,1,1.0001}); // matrice quasi instabile (det approx. 1e-8)
	if(fabs(mat_det(A) - 1e-8) < 1e-15) printf("              : passed\n"); else printf("              : NOT PASSED!!!\n");
	
	
	mat_free(A); mat_free(B); mat_free(tmp); mat_free(u); mat_free(v);
	printf("\npress a key to exit...");
	getch();
	
	return 0;
}
