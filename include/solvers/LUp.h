/*
	LU factorialization (or decomposition) with partial pivoting
	
	P * A = L * U
		A : square matrix to be decomposed
		P : represents any valid (row) permutation of the Identity I matrix, and it’s computed during the process
		L : is a lower diagonal matrix, with all the elements of the first diagonal = 1
		U : is an upper diagonal matrix
*/
#ifndef SLAP_LUP
#define SLAP_LUP


typedef struct _matd_lup {
	matd *L;
	matd *U;
	matd *P;
	unsigned int num_permutations; // useful when computing the determinant
} matd_lup;


matd_lup* matd_lup_new(matd *L, matd *U, matd *P, unsigned int num_permutations)
{
	matd_lup *m = malloc(sizeof(*m));
	MEM_CHECK(m);
	m->L = L;
	m->U = U;
	m->P = P;
	m->num_permutations = num_permutations;
	return m;
}

void matd_lup_free(matd_lup* lu)
{
	if(lu){
		free_mat(lu->L);
		free_mat(lu->U);
		free_mat(lu->P);
		free(lu);
	}
}


int matd_setdiag(matd *m, double value)
{
	// Sets all elements of the matrix to given value
	int i;
	if(m->n_rows != m->n_cols){ /* SLAP_ERROR(CANNOT_SET_DIAG, value); */ return 0; }
	for(i=0; i<m->n_rows; i++) m->data[i*m->n_cols+i] = value;
	return 1;
}

int matd_row_addrow_r(matd *m, unsigned int where, unsigned int row, double multiplier)
{
	int i = 0;
	if(where >= m->n_rows || row >= m->n_rows){
//		SLAP_ERROR(CANNOT_ADD_TO_ROW, multiplier, row, where, m->num_rows);
		return 0;
	}
	for(i=0; i<m->n_cols; i++) m->data[where*m->n_cols+i] += multiplier * m->data[row*m->n_cols+i];
	return 1;
}

int matd_row_swap_r(matd *m, unsigned int row1, unsigned int row2)
{
	// swap two rows of matrix m
	int i;
	double tmp;
	if(row1 >= m->n_rows || row2 >= m->n_rows){
//		SLAP_ERROR(CANNOT_SWAP_ROWS, row1, row2, m->num_rows);
		return 0;
	}
	for(i=0; i<m->n_cols; i++){
		tmp = m->data[row2*m->n_cols+i];
		m->data[row2*m->n_cols+i] = m->data[row1*m->n_cols+i];
		m->data[row1*m->n_cols+i] = tmp;
	}
	return 1;
}


int matd_absmaxr(matd *m, unsigned int k)
{
	// Finds the id of the max on the column (starting from k -> num_rows)
	int i;
	double max = m->data[k*m->n_cols+k];
	int maxIdx = k;
	for(i=k+1; i<m->n_rows; i++){
		if(fabs(m->data[i*m->n_cols+k]) > max){
			max = fabs(m->data[i*m->n_cols+k]);
			maxIdx = i;
		}
	}
	return maxIdx;
}


matd_lup* matd_lup_solve(matd *m)
{
	// perform the LU(P) factorization
	if(m->n_rows != m->n_cols){
//		SLAP_ERROR(CANNOT_LU_MATRIX_SQUARE, m->num_rows, m-> num_cols);
		return NULL;
	}
	matd *L = new_matd(m->n_rows, m->n_rows);
	matd *U = matd_copy(m);
	matd *P = matd_eye(m->n_rows);
	
	int j,i, pivot;
	unsigned int num_permutations = 0;
	double mul;
	
	for(j=0; j<U->n_cols; j++){
		// Retrieves the row with the biggest element for column (j)
		pivot = matd_absmaxr(U, j);
//		if(fabs(U->data[pivot][j]) < SLAP_MIN_COEF){
////			SLAP_ERROR(CANNOT_LU_MATRIX_DEGENERATE);
//			return NULL;
//		}
		if(pivot!=j){
			// Pivots LU and P accordingly to the rule
			matd_row_swap_r(U, j, pivot);
			matd_row_swap_r(L, j, pivot);
			matd_row_swap_r(P, j, pivot);
			num_permutations++; // Keep the number of permutations to easily calculate the determinant sign afterwards
		}
		for(i=j+1; i<U->n_rows; i++){
			mul = U->data[i*U->n_cols+j] / U->data[j*U->n_cols+j];
			matd_row_addrow_r(U, i, j, -mul); // Building the U upper rows
			L->data[i*L->n_cols+j] = mul; // Store the multiplier in L
		}
	}
	matd_setdiag(L, 1.0); // set the diagonal to 1.0
	
	return matd_lup_new(L, U, P, num_permutations);
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
matd *lu_solvefwd(matd *L, matd *b)
{
	matd *x = new_matd(L->n_cols, 1);
	int i,j;
	double tmp;
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
matd *lu_solvebck(matd *U, matd *b)
{
	matd *x = new_matd(U->n_cols, 1);
	int i = U->n_cols, j;
	double tmp;
	while(i-- > 0){
		tmp = b->data[i*b->n_cols];
		for(j=i; j<U->n_cols; j++) tmp -= U->data[i*U->n_cols+j] * x->data[j*x->n_cols];
		x->data[i*x->n_cols] = tmp / U->data[i*U->n_cols+i];
	}
	return x;
} 



matd *lu_solve(matd_lup *lu, matd* b)
{
	if(lu->U->n_rows != b->n_rows || b->n_cols != 1){
//		SLAP_ERROR(CANNOT_SOLVE_LIN_SYS_INVALID_B,b->n_rows,b->n_cols,lu->U->n_rows,1);
		return NULL;
	}
	matd *Pb = matd_mul(lu->P, b);
	
	matd *y = lu_solvefwd(lu->L, Pb); // We solve L*y = P*b using forward substition
	matd *x = lu_solvebck(lu->U, y); // We solve U*x=y
	
	free_mat(y);
	free_mat(Pb);
	return x;
}





#endif // SLAP_LUP

