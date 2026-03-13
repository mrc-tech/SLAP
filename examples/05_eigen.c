//#define SLAP_DEBUG 1
#include "../include/SLAP.h"

int main()
{
	mat *A = mat_init(3,3, (TYPE[]){1,3,3,4,5,6,7,8,9});
	
	printf("A = \n"); mat_print(A);
	printf("benchmark: eig(A) = [16.366, -1.000, -0.366]\n");
	
	printf("\neigen_qr(A) = \n");
	mat *eig = eigen_qr(A);
	mat_print(eig);
	
	mat_free(A); mat_free(eig);
	
	printf("\nend program.\n");
	return 0;
}
