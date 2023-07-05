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
	// CONTROLLARE LA MATRICE CREATA ( SLAP_CHECK(m) )
	m->n_rows = num_rows;
	m->n_cols = num_cols;
	m->data = calloc(m->n_rows*m->n_cols, sizeof(*m->data));
	// CONTROLLARE I DATI CREATI ( SLAP_CHECK(m->data) )
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

unsigned int size(matd* M) { return M->n_rows * M->n_cols; }


// ####################################################################
// ######################### MATF (float)  ############################
// ####################################################################


typedef struct _matf{
	unsigned int n_rows;
	unsigned int n_cols;
	float *data;
} matf;





#endif // SLAP_DATATYPE
