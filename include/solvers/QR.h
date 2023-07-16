/*
	QR decomposition
*/
#ifndef SLAP_QR
#define SLAP_QR

#include "Gauss.h"


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
