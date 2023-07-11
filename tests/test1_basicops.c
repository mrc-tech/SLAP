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
	matd *A = new_matd(3,3);
	matd *B = new_matd(3,3);
	matd *tmp = new_matd(3,3);
	matd *u, *v;
	
	printf("matd_transpose : ");
	A = matd_init(3,3, 1.,2.,3.,4.,5.,6.,7.,8.,9.);
//	tmp = matd_init(3,3, 1.,4.,7.,2.,5.,8.,3.,6.,9.);
	tmp = matd_init2(3,3, (double[]){1,4,7,2,5,8,3,6,9});
	if(matd_equal(matd_transpose(A),tmp,0))				printf("passed\n"); else printf("NOT PASSED!!!\n");
	
	printf("matd_smul      : ");
	tmp = matd_init(3,3, 2.,4.,6.,8.,10.,12.,14.,16.,18.);
	matd_smul_r(tmp,0.5); // divide by two
	if(matd_equal(A,tmp,0))								printf("passed\n"); else printf("NOT PASSED!!!\n");
	
	printf("matd_add       : ");
	tmp = matd_init(3,3, 2.,4.,6.,8.,10.,12.,14.,16.,18.);
	if(matd_equal(matd_add(A,A),tmp,0))					printf("passed\n"); else printf("NOT PASSED!!!\n");
	
	printf("matd_sub       : ");
	if(matd_equal(matd_sub(tmp,tmp),new_matd(3,3),0))	printf("passed\n"); else printf("NOT PASSED!!!\n");
	
	printf("matd_mul       : ");
	tmp = matd_init(3,3, 30.,36.,42.,66.,81.,96.,102.,126.,150.);
	if(matd_equal(matd_mul(A,A),tmp,0))					printf("passed\n"); else printf("NOT PASSED!!!\n");
	u = new_matd(1,3); u = matd_init(1,3, 1.,2.,3.);
	v = new_matd(1,3); v = matd_init(1,3, 2.,2.,2.);
	tmp = new_matd(1,1); tmp = matd_init(1,1, 12.); // scalar product (vectors)
	if(matd_equal(matd_mul(u,matd_transpose(v)),tmp,0)) printf("               : passed\n"); else printf("                 NOT PASSED!!!\n");
	tmp = new_matd(3,3); tmp = matd_init(3,3, 2.,2.,2.,4.,4.,4.,6.,6.,6.); // tensor product (vectors)
	if(matd_equal(matd_mul(matd_transpose(u),v),tmp,0)) printf("               : passed\n"); else printf("                 NOT PASSED!!!\n");
	
	
	
	
	free_mat(A);
	free_mat(B);
	free_mat(tmp);
	
	getch();
	
	return 0;
}
