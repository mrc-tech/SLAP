//#define SLAP_DEBUG 1
#include "../include/SLAP.h"

int main()
{
	mat *A = mat_init(3,3, (TYPE[]){1,3,3,4,5,6,7,8,9});
	
	printf("A = \n"); mat_print(A);
	printf("benchmark: eig(A) = [16.366, -1.000, -0.366]\n");
	
	printf("\neigen_qr(A) = \n");
	mat *eig = eigen_qr(A, 1e-6, 1000);
	mat_print(eig);
	
	printf("\neigen_power_method(A) = \n");
	mat *autovettore = mat_new(A->n_rows, 1);
	TYPE max_eig = eigen_power_method(A, autovettore, 1e-6, 1000);
	printf("massimo autovalore: %f\n", max_eig);
	printf("equivalente autovalore:\n");
	mat_print(autovettore);
	
	mat_free(A); mat_free(eig);
	printf("\nend program.\n");
	return 0;
}
