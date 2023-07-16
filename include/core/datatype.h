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
