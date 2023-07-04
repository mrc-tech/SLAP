/*
	Simple Linear Algebra Package (SLAP)
	basic operations
*/
#ifndef SLAP_BASICOPS
#define SLAP_BASICOPS

void print_mat(matd matrix)
{
	// FARE IN MODO CHE STAMPA COME UNA TABELLA PREORDINATA DAL NUMERO DELLE CIFRE (tutto compatto)
	int r,c;
	for(r=0; r<matrix.n_rows; r++){
		for(c=0; c<matrix.n_cols; c++){
			printf("%lf\t", matd_get(matrix, r,c));
		}
		printf("\n");
	}
}



/*
	TODO:
		- transpose
		- equality
		- sum
		- multiply
		- extract row/column
		
*/


#endif // SLAP_BASICOPS
