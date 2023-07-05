/*
	Simple Linear Algebra Package (SLAP)
*/

// INITIAL DEFs:
#ifndef SLAP_DEFS
#define SLAP_DEFS


#define MEM_CHECK(ptr) \
	if (!(ptr)) { \
		fprintf(stderr, "%s:%d NULL POINTER: %s n", __FILE__, __LINE__, (#ptr)); \
		exit(-1); \
	}


#endif // SLAP_DEFS


// CORE:
/*
	Simple Linear Algebra Package (SLAP)
	data type definition
		d: double
		f: float
		i: int
		b: unsigned char
*/
#ifndef SLAP_DATATYPE
#define SLAP_DATATYPE

#include <stdlib.h> // for "calloc"
#include <stdarg.h> // for variable argument list ("va_list")


// ####################################################################
// ######################### MATD (DOUBLE)  ###########################
// ####################################################################


typedef struct _matd{
	unsigned int n_rows;
	unsigned int n_cols;
	double *data; // row-major matrix data array
} matd;



// --------------------------------------------------------------------


//#define MAT(m, row,col) ( &m + row*m.n_cols + col ) // row-major accessing member NOT WORKING!!


matd* new_matd(unsigned int num_rows, unsigned int num_cols)
{
	// create a new double matrix
	matd * m; // return matrix
	int i;
	
	if(num_rows == 0) { /*SLAP_ERROR(INVALID_ROWS);*/ return 0; } // dovrebbe ritornare NULL
	if(num_cols == 0) { /*SLAP_ERROR(INVALID_COLS);*/ return 0; } // dovrebbe ritornare NULL
	
	m = calloc(1, sizeof(*m)); // allocate space for the struct
	MEM_CHECK(m);
	m->n_rows = num_rows;
	m->n_cols = num_cols;
	m->data = calloc(num_rows*num_cols, sizeof(*m->data));
	MEM_CHECK(m->data);
	for(i=0; i<num_rows*num_cols; i++) m->data[i] = 0; // set to zero
	
	return m;
}

void free_mat(matd* matrix)
{
	free(matrix->data); // delete the data
	free(matrix); // delete the data structure
}


double matd_get(matd* M, unsigned int row, unsigned int col) { return M->data[row*M->n_cols + col]; } // row-major CONTROLLARE LA VALIDITA` DEGLI INDICI!!!!!
void   matd_set(matd* M, unsigned int row, unsigned int col, double val) { M->data[row*M->n_cols + col] = val; } // row-major CONTROLLARE LA VALIDITA` DEGLI INDICI!!!!!

unsigned int matd_size(matd* M) { return M->n_rows * M->n_cols; }


matd* matd_init(unsigned int num_rows, unsigned int num_cols, ...)
{
	// non funziona con i vecchi compilatori perche` il secondo parametro di va_start(,) deve essere l'ultimo parametro passato alla funzione
	va_list valist;
	int i = num_rows*num_cols;
	matd *m = new_matd(num_rows,num_cols);
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


// ####################################################################
// ######################### MATF (float)  ############################
// ####################################################################


typedef struct _matf{
	unsigned int n_rows;
	unsigned int n_cols;
	float *data;
} matf;





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

#include <math.h>


void print_mat(matd* matrix)
{
	// FARE IN MODO CHE STAMPA COME UNA TABELLA PREORDINATA DAL NUMERO DELLE CIFRE (tutto compatto)
	int r,c;
	for(r=0; r<matrix->n_rows; r++){
		for(c=0; c<matrix->n_cols; c++){
			printf("%lf\t", matd_get(matrix, r,c));
		}
		printf("\n");
	}
}


short matd_equal(matd* m1, matd* m2, double tolerance)
{
	// return 1 if m1 = m2, else returns 0
	int i;
	if((m1->n_rows != m2->n_rows) || (m1->n_cols != m2->n_cols)) return 0; // check dimensions
	for(i=0; i<m1->n_rows*m1->n_cols; i++){
		if(fabs(m1->data[i] - m2->data[i]) > tolerance) return 0;
	}
	return 1;
}


matd* matd_copy(matd *m)
{
	// Dynamically allocates a new Matrix
	// Initialise the matrix by copying another one
	matd *res  = new_matd(m->n_rows, m->n_cols);
	int i;
	for(i=0; i<res->n_rows*res->n_cols; i++) res->data[i] = m->data[i];
	return res;
}


matd* matd_transpose(matd* matrix)
{
	// QUALCOSA NON MI CONVINCE CON LA GESTIONE DELLA MEMORIA!!!!
	matd * m = new_matd(matrix->n_cols, matrix->n_rows); // return matrix
	int r,c;
	for(r=0; r<m->n_rows; r++){
		for(c=0; c<m->n_cols; c++){
			m->data[c * m->n_cols + r] = matrix->data[r * matrix->n_cols + c];
		}
	}
	
	return m;
}



int matd_smul(matd *m, double num)
{
	// nultiply matrix by a scalar
	int i;
	for(i=0; i<m->n_rows*m->n_cols; i++) m->data[i] *= num;
	// CONTROLLO DELLA MEMORIA E RITORNARE VALORI DI CONTROLLO
	return 1;
}


matd* matd_add(matd *m1, matd *m2)
{
	matd *m = matd_copy(m1);
	if(!matd_add_r(m, m2)) { free_mat(m); return 0; }
	return m;
}

int matd_add_r(matd *m1, matd *m2)
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

matd* matd_sub(matd *m1, matd *m2)
{
	matd *m = matd_copy(m1);
	if(!matd_sub_r(m, m2)) { free_mat(m); return 0; }
	return m;
}

int matd_sub_r(matd *m1, matd *m2)
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





#endif // SLAP_BASICOPS



// OPERATIONS:
/*
	matrix structure modification
*/
#ifndef SLAP_STRMOD
#define SLAP_STRMOD


matd* matd_remcol(matd *m, unsigned int column)
{
	// remove the i-th column (start counting from zero)
	matd *ret;
	int i, j, k;
	if(column >= m->n_cols){
//    	SLAP_FERROR(CANNOT_REMOVE_COLUMN, column, m->num_cols);
		return 0; // return NULL
	}
	ret = new_matd(m->n_rows, m->n_cols-1);
	for(i=0; i<m->n_rows; i++){
		for(j=0,k=0; j<m->n_cols; j++){
			if(column != j) ret->data[i*ret->n_cols + k++] = m->data[i*m->n_cols + j];
		}
	}
	return ret;
}

matd* matd_remrow(matd *m, unsigned int row)
{
	// remove the i-th row (start counting from zero)
	matd* ret;
	int i, j, k;
	if(row >= m->n_rows){
//    	SLAP_FERROR(CANNOT_REMOVE_ROW, row, m->num_rows);
		return 0; // return NULL
	}
	ret = new_matd(m->n_rows-1, m->n_cols);
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



matd* matd_catver(int N, matd **marr) // NON VERIFICATO!!!!!!
{
	// concatenate vertically N matrices
	matd *r;
	unsigned int lrow, ncols;
	int i,j,k,offset;
	if (N == 0) return 0; // No matrices, nothing to return
	if (N == 1) return matd_copy(marr[0]); // no need for additional computations

	// We calculate the total number of columns to know how to allocate memory for the resulting matrix:
	lrow  = marr[0]->n_rows;
	ncols = marr[0]->n_cols;
	for(k=1; k<N; k++){
		if (NULL == marr[k]){
//			SLAP_FERROR(INCONSITENT_ARRAY, k, mnum);
			return 0;
		}
		if (lrow != marr[k]->n_rows){
//			SLAP_FERROR(CANNOT_CONCATENATE_H, lrow, marr[k]->num_rows);
			return 0;
		}
		ncols += marr[k]->n_cols;
	}
	// allocate memory for the resulting matrix
	r = new_matd(lrow, ncols);
	for(i = 0; i<r->n_rows; i++){
		k = 0;
		offset = 0;
		for(j = 0; j<r->n_cols; j++){
			// If the column index of marr[k] overflows
			if(j-offset == marr[k]->n_cols){
				offset += marr[k]->n_cols;
				k++; // jump to the next matrix in the array
			}
		r->data[i*r->n_cols+j] = marr[k]->data[i*r->n_cols + j - offset];
		}
	}
	return r;
}


matd* matd_cathor(unsigned int N, matd **marr) // NON VERIFICATO!!!!!
{
	// concatenate matrices horizontally (same number of rows, aumented number of columns)
	matd *res;
	int lcol, i, j, k, offset;
	unsigned int numrows = 0;
	if(N == 0) return 0;
	if(N == 1) return matd_copy(marr[0]);

	lcol = marr[0]->n_cols;
	for(i=0; i<N; i++){
		if(marr[i] == 0){
//			SLAP_FERROR(INCONSITENT_ARRAY, i, mnum);
			return 0;
		}
		if(lcol != marr[i]->n_cols){
//			SLAP_FERROR(CANNOT_CONCATENATE_V,lcol,marr[i]->num_cols);
			return 0;
		}
		numrows += marr[i]->n_rows;
	}
	res = new_matd(numrows, lcol);
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


// UTILITIES:
#ifndef SLAP_UTILS
#define SLAP_UTILS

#include <stdio.h> // per "FILE"


matd* matd_fromfile(FILE *f)
{
	matd *m;
	int r,c;
	unsigned int num_rows = 0, num_cols = 0;
	fscanf(f, "%d", &num_rows);
	fscanf(f, "%d", &num_cols);
	m = new_matd(num_rows, num_cols);
	for(r=0; r<m->n_rows; r++){
		for(c=0; c<m->n_cols; c++){
			fscanf(f, "%lf\t", &m->data[r * m->n_cols + c]);
		}
	}
	return m;
}

matd* matd_fromfilename(const char *file)
{
	matd *res;
	FILE *m_file = fopen(file, "r");
	if (!m_file) {
//		SLAP_FERROR(CANNOT_OPEN_FILE, file);
		return 0;
	}
	res = matd_fromfile(m_file);
	fclose(m_file);
	return res;
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

