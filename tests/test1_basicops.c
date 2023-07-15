#include <stdio.h>

#include "../include/SLAP.h" // Simple Linear Algebra Package


int main()
{
	/*
	    | 1 2 3 |       |  |
	A = | 4 5 6 |   B = |  |
	    | 7 8 9 |       |  |
	*/
	
	int i, j;
	mat *A = mat_new(3,3);
	mat *B = mat_new(3,3);
	mat *tmp = mat_new(3,3);
	mat *u, *v;
	
	printf("mat_transpose : ");
	A = mat_init2(3,3, (double[]){1,2,3,4,5,6,7,8,9});
	tmp = mat_init2(3,3, (double[]){1,4,7,2,5,8,3,6,9});
	if(mat_equal(mat_transpose(A),tmp,0))				printf("passed\n"); else printf("NOT PASSED!!!\n");
	
	printf("mat_smul      : ");
	tmp = mat_init2(3,3, (double[]){2,4,6,8,10,12,14,16,18});
	mat_smul_r(tmp,0.5); // divide by two
	if(mat_equal(A,tmp,0))								printf("passed\n"); else printf("NOT PASSED!!!\n");
	
	printf("mat_add       : ");
	tmp = mat_init2(3,3, (double[]){2,4,6,8,10,12,14,16,18});
	if(mat_equal(mat_add(A,A),tmp,0))					printf("passed\n"); else printf("NOT PASSED!!!\n");
	
	printf("mat_sub       : ");
	if(mat_equal(mat_sub(tmp,tmp),mat_new(3,3),0))		printf("passed\n"); else printf("NOT PASSED!!!\n");
	
	printf("mat_mul       : ");
	tmp = mat_init2(3,3, (double[]){30,36,42,66,81,96,102,126,150});
	if(mat_equal(mat_mul(A,A),tmp,0))					printf("passed\n"); else printf("NOT PASSED!!!\n");
	u = mat_new(1,3); u = mat_init2(1,3, (double[]){1,2,3});
	v = mat_new(1,3); v = mat_init2(1,3, (double[]){2,2,2});
	tmp = mat_new(1,1); tmp = mat_init2(1,1, (double[]){12}); // scalar product (vectors)
	if(mat_equal(mat_mul(u,mat_transpose(v)),tmp,0)) printf("              : passed\n"); else printf("                 NOT PASSED!!!\n");
	tmp = mat_new(3,3); tmp = mat_init2(3,3, (double[]){2,2,2,4,4,4,6,6,6}); // tensor product (vectors)
	if(mat_equal(mat_mul(mat_transpose(u),v),tmp,0)) printf("              : passed\n"); else printf("                 NOT PASSED!!!\n");
	
	
	
	
	mat_free(A);
	mat_free(B);
	mat_free(tmp);
	
	getch();
	
	return 0;
}
