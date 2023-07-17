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
	m = mat_new(m1->n_rows, m2->n_cols); // also set all values to zero
	for(r=0; r<m->n_rows; r++){
		for(c=0; c<m->n_cols; c++){
			for(i=0; i<m1->n_cols; i++){
				m->data[r*m->n_cols+c] += m1->data[r*m1->n_cols+i] * m2->data[i*m2->n_cols+c];
			}
		}
	}
	
	return m;
}


TYPE mat_trace(mat *m)
{
	// trace of the matrix m
	TYPE tr = 0;
	int i;
	for(i=0; i<MIN(m->n_rows,m->n_cols); i++) tr += m->data[i * m->n_cols + i];
	return tr;
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
