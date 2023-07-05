#include <stdio.h>

#include "../include/SLAP.h" // Simple Linear Algebra Package


int main()
{
	FILE *f;
	matd *A = new_matd(3,3); // double matrix (3x3)
	matd *B = new_matd(4,3); // 4 rows, 3 cols
	
	matd_set(A, 0,0, 2); // assign values on the main diagonal
	matd_set(A, 1,1, 4);
	matd_set(A, 2,2, 8);
	printf("A = \n");
	print_mat(A); // print the matrix
	matd_set(B, 0,0, 1); matd_set(B, 0,1, 2); matd_set(B, 0,2, 3); matd_set(B, 1,0, 4); matd_set(B, 1,1, 5); matd_set(B, 1,2, 6); matd_set(B, 2,0, 7); 
	printf("B = \n");
	print_mat(B);
	printf("\n");
//	matd *test = matd_transpose(B); // QUESTO COMANDO DA PROBLEMI ALLA RIGA 24 !!!!!
	
	printf("is_equal: %d\n", is_equal(A,B, 0));
	B = new_matd(3,3); matd_set(B,0,0,2); matd_set(B,1,1,4); matd_set(B,2,2,8); // DA PROBLEMI SE SI FA LA TRASPOSTA DI B !!!!!!
	printf("is_equal: %d\n", is_equal(A,B, 0)); // exactly equal
	printf("is_equal: %d\n", is_equal(A,B, 1e-20)); // almost equal
	
	free_mat(A);
	free_mat(B);
	
	f = fopen("matrix.txt","r");
	A = matd_fromfile(f);
	fclose(f);
	
	printf("A(from file) = \n");
	print_mat(A);
	matd *C = matd_transpose(A);
	printf("A^(T) = \n");
	print_mat(C);
	printf("is_equal: %d\n", is_equal(A,matd_transpose(C), 0));
	
	
	
	
	return 0;
}
