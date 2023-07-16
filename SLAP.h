/*
	Simple Linear Algebra Package (SLAP)
*/
#include <stdlib.h>
#include <string.h>	// for memory menagement functions
#include <stdio.h>	// for error messages
#include <math.h>


// INITIAL DEFs:
#ifndef SLAP_DEFS
#define SLAP_DEFS



#ifndef TYPE // data type of the matrix
#ifdef __MSDOS__ // MS-DOS system
#pragma message You are compiling using Borland C++ version __BORLANDC__.
#define SLAP_DOS
#define TYPE long double
#else // other non-DOS systems
#define TYPE double
#endif
#endif // TYPE

#ifndef SLAP_DEBUG
#define SLAP_DEBUG 0 // no debug
#endif


//#define NULL 0
#define SLAP_MIN_COEF 0.000000000000001 // DIPENDE DAL SISTEMA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#define MEM_CHECK(ptr) \
	if (!(ptr)) { \
		fprintf(stderr, "%s:%d NULL POINTER: %s n", __FILE__, __LINE__, (#ptr)); \
		exit(-1); \
	}


#define SWAP(a,b) \
	TYPE _temp = a; \
	a = b; \
	b = _temp;


#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
	



#endif // SLAP_DEFS

// CORE:
/*
	Simple Linear Algebra Package (SLAP)
	data type definition
	basic type is defined through TYPE preprocessor definition
*/
#ifndef SLAP_DATATYPE
#define SLAP_DATATYPE

#include <stdarg.h> // for variable argument list ("va_list")




typedef struct _mat{
	unsigned int n_rows;
	unsigned int n_cols;
	TYPE *data; // row-major matrix data array
} mat;



// --------------------------------------------------------------------


//#define MAT(m, row,col) ( &m + row*m.n_cols + col ) // row-major accessing member NOT WORKING!!


mat* mat_new(unsigned int num_rows, unsigned int num_cols)
{
	// create a new double matrix
	mat * m; // return matrix
	int i;
	
	if(num_rows == 0) { /*SLAP_ERROR(INVALID_ROWS);*/ return NULL; }
	if(num_cols == 0) { /*SLAP_ERROR(INVALID_COLS);*/ return NULL; }
	
	m = calloc(1, sizeof(*m)); // allocate space for the struct
	MEM_CHECK(m);
	m->n_rows = num_rows;
	m->n_cols = num_cols;
	m->data = calloc(num_rows*num_cols, sizeof(*m->data));
	MEM_CHECK(m->data);
	for(i=0; i<num_rows*num_cols; i++) m->data[i] = 0; // set to zero
	if(SLAP_DEBUG) printf("mat_new\n"); // for memory allocation check
	
	return m;
}

void mat_free(mat* matrix)
{
	if(matrix){
		if(matrix->data) free(matrix->data); // delete the data
		free(matrix); // delete the data structure
		if(SLAP_DEBUG) printf("mat_free\n"); // for memory allocation check
	}
}


TYPE mat_get(mat* M, unsigned int row, unsigned int col) { return M->data[row*M->n_cols + col]; } // row-major CONTROLLARE LA VALIDITA` DEGLI INDICI!!!!!
void mat_set(mat* M, unsigned int row, unsigned int col, TYPE val) { M->data[row*M->n_cols + col] = val; } // row-major CONTROLLARE LA VALIDITA` DEGLI INDICI!!!!!

unsigned int mat_size(mat* M) { return M->n_rows * M->n_cols; }


mat* mat_init(unsigned int num_rows, unsigned int num_cols, TYPE data[])
{
	int i;
	mat *m = mat_new(num_rows,num_cols);
	for (i=0; i<num_rows*num_cols; i++) m->data[i] = data[i];
	return m;
}

#ifdef SLAP_DOS
mat* mat_init_DOS(unsigned int num_rows, unsigned int num_cols, ...)
{
	// non funziona con i vecchi compilatori perche` il secondo parametro di va_start(,) deve essere l'ultimo parametro passato alla funzione
	// RISOLVERER LEGGENDO HELP DI BORLAND TURBO C!!!!!!!
	va_list valist;
	int i = num_rows*num_cols;
	mat *m = mat_new(num_rows,num_cols);
	va_start(valist, i); // initialize valist for num number of arguments
	for (i=0; i<num_rows*num_cols; i++) m->data[i] = va_arg(valist, double);
	va_end(valist); // clean memory reserved for valist
	return m;
}
//void matd_init(matd *m, unsigned num, ...)
//{
//	va_list valist;
//	int i;
//	double tmp;
//	va_start(valist, num); // initialize valist for num number of arguments
//	for (i=0; i<m->n_rows*m->n_cols; i++){
//		tmp = va_arg(valist, double);
//		m->data[i] = tmp;
//		printf("%lf\n", tmp);
//	}
//	va_end(valist); // clean memory reserved for valist
//}
#endif




#endif // SLAP_DATATYPE
/*
	Simple Linear Algebra Package (SLAP)
	basic operations
	
	TODO:
		- sum
		- multiply
		- extract row/column
*/
#ifndef SLAP_BASICOPS
#define SLAP_BASICOPS




void mat_print(mat* matrix)
{
	// FARE IN MODO CHE STAMPA COME UNA TABELLA PREORDINATA DAL NUMERO DELLE CIFRE (tutto compatto)
	int r,c;
	for(r=0; r<matrix->n_rows; r++){
		for(c=0; c<matrix->n_cols; c++){
			printf("%lf\t", mat_get(matrix, r,c));
		}
		printf("\n");
	}
}


short mat_equal(mat* m1, mat* m2, TYPE tolerance)
{
	// return 1 if m1 = m2, else returns 0
	int i;
	if((m1->n_rows != m2->n_rows) || (m1->n_cols != m2->n_cols)) return 0; // check dimensions
	for(i=0; i<m1->n_rows*m1->n_cols; i++){
		if(fabs(m1->data[i] - m2->data[i]) > tolerance) return 0;
	}
	return 1;
}


mat* mat_copy(mat *m)
{
	// Dynamically allocates a new Matrix
	// Initialise the matrix by copying another one
	mat *res  = mat_new(m->n_rows, m->n_cols);
	int i;
	for(i=0; i<res->n_rows*res->n_cols; i++) res->data[i] = m->data[i];
	return res;
}


mat* mat_transpose(mat* matrix)
{
	// QUALCOSA NON MI CONVINCE CON LA GESTIONE DELLA MEMORIA!!!!
	mat *m = mat_new(matrix->n_cols, matrix->n_rows); // return matrix
	int r,c;
	for(r=0; r<m->n_rows; r++){
		for(c=0; c<m->n_cols; c++){
			m->data[r * m->n_cols + c] = matrix->data[c * m->n_rows + r];
		}
	}
	return m;
}

int mat_transpose_r(mat* m)
{
	// change the matrix by reference ("_r")
	// without swap dimensions this is a conversion between row-major and column-major
	int i, j;
	TYPE temp;
	TYPE *tmp = (TYPE*) malloc(m->n_rows*m->n_cols * sizeof(TYPE)); // allocate temporary array
	for(i=0; i<m->n_rows*m->n_cols; i++){
		j = m->n_cols * (i % m->n_rows) + (i / m->n_rows);
		tmp[i] = m->data[j];
	}
	free(m->data); // free the previous matrix
	m->data = tmp; // set the new matrix data
	temp = m->n_cols; m->n_cols = m->n_rows; m->n_rows = temp; // swap dimensions
	return 1;
}



int mat_smul_r(mat *m, TYPE num)
{
	// multiply matrix by a scalar (by reference)
	int i;
	for(i=0; i<m->n_rows*m->n_cols; i++) m->data[i] *= num;
	// CONTROLLO DELLA MEMORIA E RITORNARE VALORI DI CONTROLLO
	return 1;
}

mat* mat_smul(mat *m, TYPE num)
{
	// multiply matrix by a scalar
	mat* res = mat_copy(m);
	mat_smul_r(res,num);
	return res;
}


int mat_add_r(mat *m1, mat *m2)
{
	// reference version (return value in matrix m1)
	int i;
	if((m1->n_rows != m2->n_rows) && (m1->n_cols != m2->n_cols)){
//		SLAP_ERROR(CANNOT_ADD);
		return 0;
	}
	for(i=0; i<m1->n_rows*m1->n_cols; i++) m1->data[i] += m2->data[i];
	
	return 1;
}

mat* mat_add(mat *m1, mat *m2)
{
	mat *m = mat_copy(m1);
	if(!mat_add_r(m, m2)) { mat_free(m); return NULL; }
	return m;
}

int mat_sub_r(mat *m1, mat *m2)
{
	// reference version (return value in matrix m1)
	int i;
	if((m1->n_rows != m2->n_rows) && (m1->n_cols != m2->n_cols)){
//		SLAP_ERROR(CANNOT_SUBTRACT);
		return 0;
	}
	for(i=0; i<m1->n_rows*m1->n_cols; i++) m1->data[i] -= m2->data[i];
	
	return 1;
}

mat* mat_sub(mat *m1, mat *m2)
{
	mat *m = mat_copy(m1);
	if(!mat_sub_r(m, m2)) { mat_free(m); return NULL; }
	return m;
}



mat* mat_mul(mat* m1, mat* m2)
{
	// multiply two matrices
	mat *m;
	int r, c, i;
	if(!(m1->n_cols == m2->n_rows)){
//		SLAP_ERROR(CANNOT_MULTIPLY);
		return NULL;
	}
	m = mat_new(m1->n_rows, m2->n_cols);
	for(r=0; r<m->n_rows; r++){
		for(c=0; c<m->n_cols; c++){
			for(i=0; i<m1->n_cols; i++){
				m->data[r*m->n_cols+c] += m1->data[r*m1->n_cols+i] * m2->data[i*m2->n_cols+c];
			}
		}
	}
	
	return m;
}




mat* mat_eye(unsigned int size)
{
	// identity square matrix
	int i;
	mat *m = mat_new(size, size);
	for(i=0; i<m->n_rows; i++) m->data[i*m->n_cols+i] = 1.0;
	return m;
}




#endif // SLAP_BASICOPS

// OPERATIONS:
/*
	matrix structure modification
*/
#ifndef SLAP_STRMOD
#define SLAP_STRMOD


mat* mat_remcol(mat *m, unsigned int column)
{
	// remove the i-th column (start counting from zero)
	mat *ret;
	int i, j, k;
	if(column >= m->n_cols){
//    	SLAP_FERROR(CANNOT_REMOVE_COLUMN, column, m->num_cols);
		return NULL;
	}
	ret = mat_new(m->n_rows, m->n_cols-1);
	for(i=0; i<m->n_rows; i++){
		for(j=0,k=0; j<m->n_cols; j++){
			if(column != j) ret->data[i*ret->n_cols + k++] = m->data[i*m->n_cols + j];
		}
	}
	return ret;
}

mat* mat_remrow(mat *m, unsigned int row)
{
	// remove the i-th row (start counting from zero)
	mat *ret;
	int i, j, k;
	if(row >= m->n_rows){
//    	SLAP_FERROR(CANNOT_REMOVE_ROW, row, m->num_rows);
		return NULL;
	}
	ret = mat_new(m->n_rows-1, m->n_cols);
	for(i=0,k=0; i<m->n_rows; i++){
		if(row != i){
			for(j=0; j<m->n_cols; j++){
				ret->data[k*ret->n_cols + j] = m->data[i*ret->n_cols + j];
			}
			k++;
		}
	}
	return ret;
}



mat *mat_getcol(mat *m, unsigned int col)
{
	// return matrix column
	int j;
	mat *res;
	if(col >= m->n_cols){
//		SLAP_FERROR(CANNOT_GET_COLUMN, col, m->num_cols);
		return NULL;
	}
	res = mat_new(m->n_rows, 1);
	for(j=0; j<res->n_rows; j++) res->data[j*res->n_cols] = m->data[j*m->n_cols+col];
	return res;
}

double *mat_getcol_array(mat *m, unsigned int col)
{
	// return column via array
	int i;
	TYPE *res;
	if(col >= m->n_cols) return NULL;
	res = (TYPE*)malloc(m->n_rows * sizeof(TYPE));
	for(i=0; i<m->n_rows; i++) res[i] = m->data[i*m->n_cols+col];
	return res;
}

mat *mat_getrow(mat *m, unsigned int row)
{
	// return matrix row
	mat *res;
	if(row >= m->n_rows){
//		SLAP_FERROR(CANNOT_GET_ROW, row, m->num_rows);
		return NULL;
	}
	res = mat_new(1, m->n_cols);
	memcpy(&res->data[0], &m->data[row*m->n_cols], m->n_cols * sizeof(res->data[0]));
	return res;
}

double *mat_getrow_array(mat *m, unsigned int row)
{
	// return a row via array
	TYPE *res;
	if(row >= m->n_rows) return NULL;
	res = (TYPE*)malloc(m->n_cols * sizeof(TYPE));
	memcpy(&res, &m->data[row*m->n_cols], m->n_cols * sizeof(TYPE));
	return res;
}








mat* mat_cathor(int N, mat **marr) // NON VERIFICATO!!!!!!
{
	// concatenate matrices horizontally (same number of rows, aumented number of columns)
	mat *m;
	int i, j, k, offset;
	unsigned int lrow, ncols;
	if (N == 0) return NULL; // No matrices, nothing to return
	if (N == 1) return mat_copy(marr[0]); // no need for additional computations
    
	// We calculate the total number of columns to know how to allocate memory for the resulting matrix:
	lrow  = marr[0]->n_rows;
	ncols = marr[0]->n_cols;
	for(k=1; k<N; k++){
		if (NULL == marr[k]){
//			SLAP_ERROR(INCONSITENT_ARRAY, k, mnum);
			return 0;
		}
		if (lrow != marr[k]->n_rows){
//			SLAP_ERROR(CANNOT_CONCATENATE_H, lrow, marr[k]->num_rows);
			return 0;
		}
		ncols += marr[k]->n_cols;
	}
	// allocate memory for the resulting matrix
	m = mat_new(lrow, ncols);
	for(i=0; i<m->n_rows; i++){
		k = 0;
		offset = 0;
		for(j=0; j<m->n_cols; j++){
			// If the column index of marr[k] overflows
			if(j-offset == marr[k]->n_cols){
				offset += marr[k]->n_cols;
				k++; // jump to the next matrix in the array
			}
		m->data[i*m->n_cols+j] = marr[k]->data[i*marr[k]->n_cols + j - offset];
		}
	}
	return m;
}


mat* mat_catver(unsigned int N, mat **marr)
{
	// concatenate vertically N matrices
	mat *res;
	unsigned int numrows = 0;
	int lcol, i, j, k, offset;
	if(N == 0) return NULL;
	if(N == 1) return mat_copy(marr[0]);
	
	lcol = marr[0]->n_cols;
	for(i=0; i<N; i++){
		if(marr[i] == 0){
//			SLAP_FERROR(INCONSITENT_ARRAY, i, mnum);
			return NULL;
		}
		if(lcol != marr[i]->n_cols){
//			SLAP_FERROR(CANNOT_CONCATENATE_V,lcol,marr[i]->num_cols);
			return NULL;
		}
		numrows += marr[i]->n_rows;
	}
	res = mat_new(numrows, lcol);
	for(j=0; j<res->n_cols; j++){
		offset = 0;
		k = 0;
		for(i=0; i<res->n_rows; i++){
			if(i - offset == marr[k]->n_rows){
				offset += marr[k]->n_rows;
				k++;
			}
			res->data[i * res->n_cols + j] = marr[k]->data[(i-offset) * res->n_cols + j];
		}
	}
	return res;
}



#endif // SLAP_STRMOD
/*
	Gauss Elimination
*/
#ifndef SLAP_GAUSS
#define SLAP_GAUSS



int mat_row_smul_r(mat *m, unsigned int row, TYPE num)
{
	int i;
	if(row >= m->n_rows){
//		SLAP_FERROR(CANNOT_ROW_MULTIPLY, row, m->num_rows);
		return 0;
	}
	for(i=0; i<m->n_cols; i++) m->data[row*m->n_cols+i] *= num;
	return 1;
}

int mat_col_smul_r(mat *m, unsigned int col, TYPE num)
{
	int i;
	if(col >= m->n_cols){
//		SLAP_FERROR(CANNOT_COL_MULTIPLY, row, m->num_rows);
		return 0;
	}
	for(i=0; i<m->n_rows; i++) m->data[i*m->n_cols+col] *= num;
	return 1;
}


int mat_row_addrow_r(mat *m, unsigned int where, unsigned int row, TYPE multiplier)
{
	int i = 0;
	if(where >= m->n_rows || row >= m->n_rows){
//		SLAP_ERROR(CANNOT_ADD_TO_ROW, multiplier, row, where, m->num_rows);
		return 0;
	}
	for(i=0; i<m->n_cols; i++) m->data[where*m->n_cols+i] += multiplier * m->data[row*m->n_cols+i];
	return 1;
}


int mat_row_swap_r(mat *m, unsigned int row1, unsigned int row2)
{
	// swap two rows of matrix m
	int i;
	TYPE tmp;
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


// Finds the first non-zero element on the col column, under the row row.
// Used to determine the pivot
// If not pivot is found, returns -1
int mat_pivot_id(mat *m, unsigned int col, unsigned int row)
{
	int i;
	for(i=row; i<m->n_rows; i++) if(fabs(m->data[i*m->n_cols+col]) > SLAP_MIN_COEF) return i;
	return -1;
}

// Find the max element from the column "col" under the row "row"
// This is needed to pivot in Gauss-Jordan elimination
// Return the maximum pivot for numerical stability. If pivot is not found, return -1
int mat_pivot_maxid(mat *m, unsigned int col, unsigned int row)
{
	int i, maxi;
	TYPE micol;
	TYPE max = fabs(m->data[row*m->n_cols+col]);
	maxi = row;
	for(i=row; i<m->n_rows; i++){
		micol = fabs(m->data[i*m->n_cols+col]);
		if(micol > max){
			max = micol;
			maxi = i;
		}
	}
	return(max < SLAP_MIN_COEF) ? -1 : maxi;
}


// Retrieves the matrix in Row Echelon form using Gauss Elimination
mat *mat_GaussJordan(mat *m)
{
	mat *r = mat_copy(m);
	int i=0, j=0, k, pivot;
	while(j < r->n_cols && i < r->n_cols){
		// Find the pivot - the first non-zero entry in the first column of the matrix
		pivot = mat_pivot_maxid(r, j, i);
		if(pivot<0){ // All elements on the column are zeros
			j++; // Move to the next column without doing anything
			continue;
		}
		if(pivot != i) mat_row_swap_r(r, i, pivot); // We interchange rows moving the pivot to the first row that doesn't have already a pivot in place
		mat_row_smul_r(r, i, 1/r->data[i*r->n_cols+j]); // Multiply each element in the pivot row by the inverse of the pivot
		for(k=i+1; k<r->n_rows; k++){
			if(fabs(r->data[k*r->n_cols+j]) > SLAP_MIN_COEF){
				mat_row_addrow_r(r, k, i, -(r->data[k*r->n_cols+j])); // We add multiplies of the pivot so every element on the column equals 0
			}
		}
		i++; j++;
	}
	return r;
}







#endif // SLAP_GAUSS
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

/*
	QR decomposition
*/
#ifndef SLAP_QR
#define SLAP_QR



typedef struct _mat_qr {
	mat *Q;
	mat *R;
} mat_qr;


mat_qr* mat_qr_new()
{
	mat_qr *qr = (mat_qr*)malloc(sizeof(*qr));
	MEM_CHECK(qr);
	return qr;
}

void mat_qr_free(mat_qr *qr)
{
	if(qr){
		if(qr->Q) mat_free(qr->Q);
		if(qr->R) mat_free(qr->R);
		free(qr);
	}
}


double mat_l2norm(mat* m)
{
	int i;
	TYPE sum = 0.0;
	if(m->n_cols != 1 && m->n_rows != 1) return -1; // only vectors
	for(i=0; i<MAX(m->n_rows,m->n_cols); i++) sum += m->data[i] * m->data[i];
	return sqrt(sum);
}



//mat_qr* mat_qr_solve(mat *m)
//{
//	// find the QR decomposition of the matrix m
//	mat_qr *qr = mat_qr_new();
//	mat *Q = mat_copy(m);
//	mat *R = mat_new(m->n_rows, m->n_cols); // n_cols and n_rows have to be equal
//	
//	int j, k;
//	TYPE l2norm;
//	mat *rkj = mat_new(1,1); // scalar MEMORY ALLOCATED!!!!
//	mat *aj, *qk;
//	for(j=0; j<m->n_cols; j++){
//		aj = mat_getcol(m, j); // j-th column of the matrix m
//		for(k=0; k<j; k++){
//			rkj = mat_mul(mat_transpose(mat_getcol(m,j)), mat_getcol(Q,k)); // scalar product MEMORY LEAKAGE
//			R->data[k*R->n_cols+j] = rkj->data[0];
//			qk = mat_getcol(Q, k);
//			mat_col_smul_r(qk, 0, rkj->data[0]);
//			mat_sub_r(aj, qk);
//			mat_free(rkj); mat_free(qk); // free rjk and qk each iteration
//		}
//		for(k=0; k<Q->n_rows; k++) Q->data[k*Q->n_cols+j] = aj->data[k]; // set the j-th column of Q
//		l2norm = mat_l2norm(mat_getcol(Q, j)); // L2-norm (Euclidean) o the j-th column of Q MEMORY LEAKAGE!!!!!
//		mat_col_smul_r(Q, j, 1/l2norm); // divide by the norm
//		R->data[j*R->n_cols+j] = l2norm;
//		mat_free(aj);
//	}
//	qr->Q = Q;
//	qr->R = R;
//	return qr;
//}
mat_qr* mat_qr_solve(mat *m) // without memory leakage
{
	// find the QR decomposition of the matrix m
	mat_qr *qr = mat_qr_new();
	mat *Q = mat_copy(m);
	mat *R = mat_new(m->n_rows, m->n_cols); // n_cols and n_rows have to be equal
	
	int j, k;
	TYPE l2norm;
	mat *rkj; // scalar
	mat *aj, *qk; // column vectors of A and Q
	mat *tmp1, *tmp2; // temporary matrices for correct memory menagement
	for(j=0; j<m->n_cols; j++){
		aj = mat_getcol(m, j); // j-th column of the matrix m
		for(k=0; k<j; k++){
			tmp1 = mat_getcol(m,j); mat_transpose_r(tmp1); // transpose of the j-th column of matrix m (row-vector)
			tmp2 = mat_getcol(Q,k); // k-th column of the matrix Q (column-vector)
			rkj = mat_mul(tmp1, tmp2); // scalar product
			mat_free(tmp1); mat_free(tmp2); // free temp mem
			R->data[k*R->n_cols+j] = rkj->data[0];
			qk = mat_getcol(Q, k);
			mat_col_smul_r(qk, 0, rkj->data[0]);
			mat_sub_r(aj, qk);
			mat_free(rkj); mat_free(qk); // free rjk and qk each iteration
		}
		for(k=0; k<Q->n_rows; k++) Q->data[k*Q->n_cols+j] = aj->data[k]; // set the j-th column of Q
		tmp1 = mat_getcol(Q, j); // j-th column of Q
		l2norm = mat_l2norm(tmp1); // L2-norm (Euclidean)
		mat_free(tmp1); // free temp mem
		mat_col_smul_r(Q, j, 1/l2norm); // divide by the norm
		R->data[j*R->n_cols+j] = l2norm;
		mat_free(aj);
	}
	qr->Q = Q;
	qr->R = R;
	return qr;
}


#endif // SLAP_QR
/*
	Eigen-analysis using QR decomposition
*/
#ifndef SLAP_EIGEN_QR
#define SLAP_EIGEN_QR

mat* mat_get_diag(mat *m) // DA METTERE IN STRMOD!!!!!!!!!!!!!
{
	// returns the diagonal of the matrix m as a column vector
	int N = MIN(m->n_rows,m->n_cols); // for non-square matrices
	mat *d = mat_new(N,1); // column vector
	int i;
	for(i=0; i<N; i++) d->data[i] = m->data[i*m->n_cols+i];
	return d;
}

mat* eigen_qr(mat *m)
{
	mat *A = mat_copy(m);
	mat_qr *qr;
	int i;
	if(m->n_rows != m->n_cols) return NULL; // ERROR!! nxn square matrix
	for(i=0; i<100; i++){ // IL NUMERO DI ITERAZIONI DIPENDE DALLA CONVERGENZA CERCATA!!!!!!
		qr = mat_qr_solve(A);
		mat_free(A); // avoid memory leakage
		A = mat_mul(qr->R, qr->Q); // update A matrix for i-th step
		mat_qr_free(qr); // free used memory for QR decomposition
	}
//	mat_print(A);
	return mat_get_diag(A);
}


#endif // SLAP_EIGEN_QR



// UTILITIES:
#ifndef SLAP_UTILS
#define SLAP_UTILS

#include <stdio.h> // per "FILE"


mat* mat_fromfile(FILE *f)
{
	mat *m;
	int r,c;
	unsigned int num_rows = 0, num_cols = 0;
	fscanf(f, "%d", &num_rows);
	fscanf(f, "%d", &num_cols);
	m = mat_new(num_rows, num_cols);
	for(r=0; r<m->n_rows; r++){
		for(c=0; c<m->n_cols; c++){
			fscanf(f, "%lf\t", &m->data[r * m->n_cols + c]);
		}
	}
	return m;
}

mat* mat_fromfilename(const char *file)
{
	mat *m;
	FILE *m_file = fopen(file, "r");
	if (!m_file) {
//		SLAP_FERROR(CANNOT_OPEN_FILE, file);
		return NULL;
	}
	m = mat_fromfile(m_file);
	fclose(m_file);
	return m;
}

#endif // SLAP_UTILS
/*
	random number generator

	TODO:
		- Random number generator (mod-type)
		- Uniform distribution sampler
		- Normal distribution sampler (Box-Muller based)
		- Reliable random generator (research)
*/
#ifndef SLAP_RANDOM
#define SLAP_RANDOM





#endif // SLAP_RANDOM
/*
	scientific constants
*/
#ifndef SLAP_CONSTANTS
#define SLAP_CONSTANTS

#define PI 3.14159265359
#define E  2.718281828459045
#define FI 1.61803398875 // golden ratio

#endif // SLAP_CONSTANTS

