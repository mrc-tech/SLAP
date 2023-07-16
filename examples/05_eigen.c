//#define SLAP_DEBUG 1
#include "../include/SLAP.h"

int main()
{
	mat *A = mat_init2(3,3, (double[]){1,3,3,4,5,6,7,8,9});
	
	// eig(A) = [16.3666, -1.0000, -0.3666]
	mat *eig = eigen_qr(A);
	mat_print(eig);
	
	mat_free(A); mat_free(eig);
	return 0;
}
