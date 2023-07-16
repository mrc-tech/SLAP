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
	for(i=0; i<22; i++){ // IL NUMERO DI ITERAZIONI DIPENDE DALLA CONVERGENZA CERCATA!!!!!!
		qr = mat_qr_solve(m);
		A = mat_mul(qr->R, qr->Q); // update A matrix for i-th step
	}
	mat_qr_free(qr); // free used memory for QR decomposition
	mat_print(A);
	return mat_get_diag(A);
}


#endif // SLAP_EIGEN_QR
