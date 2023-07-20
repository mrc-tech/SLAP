// --- efficient matrix multiplication ---

mat* mat_mul_cache1(mat* m1, mat* m2)
{
	// multiply two matrices
	// cache optimization with loop interchanging
	mat *m;
	int r, c, i;
	if(!(m1->n_cols == m2->n_rows)){
//		SLAP_ERROR(CANNOT_MULTIPLY);
		return NULL;
	}
	m = mat_new(m1->n_rows, m2->n_cols); // also set all values to zero
	for(r=0; r<m->n_rows; r++){
		for(i=0; i<m1->n_cols; i++){ // swapped with the for below
			for(c=0; c<m->n_cols; c++){
				m->data[r*m->n_cols+c] += m1->data[r*m1->n_cols+i] * m2->data[i*m2->n_cols+c];
			}
		}
	}
	
	return m;
}


int mat_mul_Hadamard_r{mat* A, mat* B}
{
	// tensor inner product (element-wise)
	// put the result into A (reference)
	int i;
	if((A->n_rows != B->n_rows) || (A->n_cols != B->n_cols)) return 0; // dimension of A and B must be the same
	for(i=0; i<A->n_rows*A->n_cols; i++){
		A->data[i] *= B->data[i];
	}
	return 1;
}
