#include <stdio.h>

#include "../include/SLAP.h" // Simple Linear Algebra Package


int main()
{
	
	matd *A = new_matd(3,3); // double matrix (3x3)
	matd *B = new_matd(4,3); // 4 rows, 3 cols
	
	matd_set(*A, 0,0, 2); // assign values on the main diagonal
	matd_set(*A, 1,1, 4);
	matd_set(*A, 2,2, 8);
	print_mat(*A); // print the matrix
	printf("\n");
	print_mat(*B);
	
	free_mat(A);
	free_mat(B);
	
	return 0;
}
