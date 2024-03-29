/*
	LU factorialization (or decomposition) with partial pivoting
	
	P * A = L * U
		A : square matrix to be decomposed
		P : represents any valid (row) permutation of the Identity I matrix, and it�s computed during the process
		L : is a lower diagonal matrix, with all the elements of the first diagonal = 1
		U : is an upper diagonal matrix
*/
#ifndef SLAP_LUP
#define SLAP_LUP

#include "Gauss.h" // for the Gaussian elimination routines


typedef struct _mat_lup {
	mat *L;
	mat *U;
	mat *P;
	unsigned int num_permutations; // useful when computing the determinant
} mat_lup;


mat_lup* mat_lup_new(mat *L, mat *U, mat *P, unsigned int num_permutations)
{
	mat_lup *m = malloc(sizeof(*m));
	MEM_CHECK(m);
	m->L = L;
	m->U = U;
	m->P = P;
	m->num_permutations = num_permutations;
	return m;
}

void mat_lup_free(mat_lup* lu)
{
	if(lu){
		mat_free(lu->L);
		mat_free(lu->U);
		mat_free(lu->P);
		free(lu);
	}
}


int mat_setdiag(mat *m, TYPE value)
{
	// Sets all elements of the matrix to given value
	int i;
	if(m->n_rows != m->n_cols){ /* SLAP_ERROR(CANNOT_SET_DIAG, value); */ return 0; }
	for(i=0; i<m->n_rows; i++) m->data[i*m->n_cols+i] = value;
	return 1;
}


int mat_absmaxr(mat *m, unsigned int k)
{
	// Finds the id of the max on the column (starting from k -> num_rows)
	int i;
	TYPE max = m->data[k*m->n_cols+k];
	int maxIdx = k;
	for(i=k+1; i<m->n_rows; i++){
		if(fabs(m->data[i*m->n_cols+k]) > max){
			max = fabs(m->data[i*m->n_cols+k]);
			maxIdx = i;
		}
	}
	return maxIdx;
}


mat_lup* mat_lup_solve(mat *m)
{
	// perform the LU(P) factorization
	mat *L, *U, *P;
	int j,i, pivot;
	unsigned int num_permutations = 0;
	TYPE mul;
	if(m->n_rows != m->n_cols){
//		SLAP_ERROR(CANNOT_LU_MATRIX_SQUARE, m->num_rows, m-> num_cols);
		return NULL;
	}
	L = mat_new(m->n_rows, m->n_rows);
	U = mat_copy(m);
	P = mat_eye(m->n_rows);
	
	for(j=0; j<U->n_cols; j++){
		// Retrieves the row with the biggest element for column (j)
		pivot = mat_absmaxr(U, j);
//		if(fabs(U->data[pivot*U->n_cols+j]) < SLAP_MIN_COEF){ // DA PROBLEMI DI MEMORIA RUNTIME!!!!!!!
////			SLAP_ERROR(CANNOT_LU_MATRIX_DEGENERATE);
//			return NULL;
//		}
		if(pivot!=j){
			// Pivots LU and P accordingly to the rule
			mat_row_swap_r(U, j, pivot);
			mat_row_swap_r(L, j, pivot);
			mat_row_swap_r(P, j, pivot);
			num_permutations++; // Keep the number of permutations to easily calculate the determinant sign afterwards
		}
		for(i=j+1; i<U->n_rows; i++){
			mul = U->data[i*U->n_cols+j] / U->data[j*U->n_cols+j];
			mat_row_addrow_r(U, i, j, -mul); // Building the U upper rows
			L->data[i*L->n_cols+j] = mul; // Store the multiplier in L
		}
	}
	mat_setdiag(L, 1.0); // set the diagonal to 1.0
	
	return mat_lup_new(L, U, P, num_permutations);
}



// Forward substitution algorithm
// Solves the linear system L * x = b
//
// L is lower triangular matrix of size NxN
// B is column matrix of size Nx1
// x is the solution column matrix of size Nx1
//
// Note: In case L is not a lower triangular matrix, the algorithm will try to
// select only the lower triangular part of the matrix L and solve the system
// with it.
//
// Note: In case any of the diagonal elements (L[i][i]) are 0 the system cannot
// be solved
//
// Note: This function is usually used with an L matrix from a LU decomposition
mat *solvefwd_lu(mat *L, mat *b)
{
	mat *x = mat_new(L->n_cols, 1);
	int i,j;
	TYPE tmp;
	for(i=0; i<L->n_cols; i++){
		tmp = b->data[i*b->n_cols];
		for(j=0; j<i; j++){
			tmp -= L->data[i*L->n_cols+j] * x->data[j*x->n_cols];
		}
		x->data[i*x->n_cols] = tmp / L->data[i*L->n_cols+i];
	}
	return x;
}

// Back substition algorithm
// Solves the linear system U *x = b
//
// U is an upper triangular matrix of size NxN
// B is a column matrix of size Nx1
// x is the solution column matrix of size Nx1
//
// Note in case U is not an upper triangular matrix, the algorithm will try to
// select only the upper triangular part of the matrix U and solve the system
// with it
//
// Note: In case any of the diagonal elements (U[i][i]) are 0 the system cannot
// be solved
mat *solvebck_lu(mat *U, mat *b)
{
	mat *x = mat_new(U->n_cols, 1);
	int i = U->n_cols, j;
	TYPE tmp;
	while(i-- > 0){
		tmp = b->data[i*b->n_cols];
		for(j=i; j<U->n_cols; j++) tmp -= U->data[i*U->n_cols+j] * x->data[j*x->n_cols];
		x->data[i*x->n_cols] = tmp / U->data[i*U->n_cols+i];
	}
	return x;
}



mat *solve_lu(mat_lup *lu, mat* b)
{
	mat *Pb, *x, *y;
	if(lu->U->n_rows != b->n_rows || b->n_cols != 1){
//		SLAP_ERROR(CANNOT_SOLVE_LIN_SYS_INVALID_B,b->n_rows,b->n_cols,lu->U->n_rows,1);
		return NULL;
	}
	Pb = mat_mul(lu->P, b);
	
	y = solvefwd_lu(lu->L, Pb); // We solve L*y = P*b using forward substition
	x = solvebck_lu(lu->U, y); // We solve U*x=y
	
	mat_free(y);
	mat_free(Pb);
	return x;
}





#endif // SLAP_LUP

