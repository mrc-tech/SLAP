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
