#include <stdio.h>

#include "../include/SLAP.h" // Simple Linear Algebra Package


int main()
{
	FILE *f;
	mat *A = mat_new(3,3); // double matrix (3x3)
	mat *B = mat_new(4,3); // 4 rows, 3 cols
	
	mat_set(A, 0,0, 2); // assign values on the main diagonal
	mat_set(A, 1,1, 4);
	mat_set(A, 2,2, 8);
	printf("A = \n");
	mat_print(A); // print the matrix
	mat_set(B, 0,0, 1); mat_set(B, 0,1, 2); mat_set(B, 0,2, 3); mat_set(B, 1,0, 4); mat_set(B, 1,1, 5); mat_set(B, 1,2, 6); mat_set(B, 2,0, 7); 
	printf("B = \n");
	mat_print(B);
	printf("\n");
	mat *test = mat_transpose(B); // QUESTO COMANDO DA PROBLEMI ALLA RIGA 24 !!!!!
	mat_transpose_r(B); mat_print(B); // questo funziona
	
	printf("matd_equal: %d\n", mat_equal(A,B, 0));
	B = mat_new(3,3); mat_set(B,0,0,2); mat_set(B,1,1,4); mat_set(B,2,2,8); // DA PROBLEMI SE SI FA LA TRASPOSTA DI B !!!!!!
	printf("matd_equal: %d\n", mat_equal(A,B, 0)); // exactly equal
	printf("matd_equal: %d\n", mat_equal(A,B, 1e-20)); // almost equal
	
	mat_free(A);
	mat_free(B);
	
	f = fopen("matrix.txt","r");
	A = mat_fromfile(f);
	fclose(f);
	
	printf("A(from file) = \n");
	mat_print(A);
	mat *C = mat_transpose(A);
	printf("A^(T) = \n");
	mat_print(C);
	printf("matd_equal: %d\n", mat_equal(A,mat_transpose(C), 0));
	A = mat_remcol(A, 1);
	mat_print(A);
	A = mat_remrow(A, 1);
	mat_print(A);
	
	
	
	printf("end program.\n");
	return 0;
}
