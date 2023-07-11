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
	matd *m = new_matd(matrix->n_cols, matrix->n_rows); // return matrix
	int r,c;
	for(r=0; r<m->n_rows; r++){
		for(c=0; c<m->n_cols; c++){
			m->data[r * m->n_cols + c] = matrix->data[c * m->n_rows + r];
		}
	}
	return m;
}

int matd_transpose_r(matd* m)
{
	// change the matrix by reference ("_r")
	// without swap dimensions this is a conversion between row-major and column-major
	int i, j;
	double temp;
	double *tmp = (double*) malloc(m->n_rows*m->n_cols * sizeof(double)); // allocate temporary array
	for(i=0; i<m->n_rows*m->n_cols; i++){
		j = m->n_cols * (i % m->n_rows) + (i / m->n_rows);
		tmp[i] = m->data[j];
	}
	free(m->data); // free the previous matrix
	m->data = tmp; // set the new matrix data
	temp = m->n_cols; m->n_cols = m->n_rows; m->n_rows = temp; // swap dimensions
	return 1;
}



int matd_smul_r(matd *m, double num)
{
	// multiply matrix by a scalar (by reference)
	int i;
	for(i=0; i<m->n_rows*m->n_cols; i++) m->data[i] *= num;
	// CONTROLLO DELLA MEMORIA E RITORNARE VALORI DI CONTROLLO
	return 1;
}

matd* matd_smul(matd *m, double num)
{
	// multiply matrix by a scalar
	matd* res = matd_copy(m);
	matd_smul_r(res,num);
	return res;
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

matd* matd_add(matd *m1, matd *m2)
{
	matd *m = matd_copy(m1);
	if(!matd_add_r(m, m2)) { free_mat(m); return NULL; }
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

matd* matd_sub(matd *m1, matd *m2)
{
	matd *m = matd_copy(m1);
	if(!matd_sub_r(m, m2)) { free_mat(m); return NULL; }
	return m;
}



matd* matd_mul(matd* m1, matd* m2)
{
	// multiply two matrices
	matd *m;
	int r, c, i;
	if(!(m1->n_cols == m2->n_rows)){
//		SLAP_ERROR(CANNOT_MULTIPLY);
		return NULL;
	}
	m = new_matd(m1->n_rows, m2->n_cols);
	for(r=0; r<m->n_rows; r++){
		for(c=0; c<m->n_cols; c++){
			for(i=0; i<m1->n_cols; i++){
				m->data[r*m->n_cols+c] += m1->data[r*m1->n_cols+i] * m2->data[i*m2->n_cols+c];
			}
		}
	}
	
	return m;
}




matd* matd_eye(unsigned int size)
{
	// identity square matrix
	int i;
	matd *m = new_matd(size, size);
	for(i=0; i<m->n_rows; i++) m->data[i*m->n_cols+i] = 1.0;
	return m;
}




#endif // SLAP_BASICOPS
