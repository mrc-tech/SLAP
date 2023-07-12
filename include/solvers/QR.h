/*
	QR decomposition
*/
#ifndef SLAP_QR
#define SLAP_QR

#include "Gauss.h"


typedef struct _matd_qr {
	matd *Q;
	matd *R;
} matd_qr;


matd_qr* matd_qr_new()
{
	matd_qr *qr = malloc(sizeof(*qr));
	MEM_CHECK(qr);
	return qr;
}

void matd_qr_free(matd_qr *qr)
{
	if(qr){
		if(qr->Q) free_mat(qr->Q);
		if(qr->R) free_mat(qr->R);
		free(qr);
	}
}


double matd_l2norm(matd* m)
{
	int i;
	double doublesum = 0.0;
	if(m->n_cols != 1 || m->n_rows != 1) return -1; // only vectors
	for(i=0; i<MAX(m->n_rows,m->n_cols); i++){
		doublesum += (m->data[i]*m->data[i]);
	}
	return sqrt(doublesum);
}



matd_qr* matd_qr_solve(matd *m)
{
	matd_qr *qr = matd_qr_new();
	matd *Q = matd_copy(m);
	matd *R = new_matd(m->n_rows, m->n_cols); // n_cols and n_rows have to be equal
	
	int j, k;
	double l2norm;
	matd *rkj = new_matd(1,1); // scalar
	matd *aj, *qk;
	for(j=0; j<m->n_cols; j++){
		rkj->data[0] = 0.0;
		aj = matd_getcol(m, j);
		for(k=0; k<j; k++){
//			rkj = nml_vect_dot(m, j, Q, k);
			rkj = matd_mul(matd_transpose(matd_getcol(m,j)), matd_getcol(Q,k)); // scalar product
			R->data[k*R->n_cols+j] = rkj->data[0];
			qk = matd_getcol(Q, k);
			matd_col_smul_r(qk, 0, rkj->data[0]);
			matd_sub_r(aj, qk);
			free_mat(qk);
		}
		for(k=0; k<Q->n_rows; k++) Q->data[k*Q->n_cols+j] = aj->data[k*aj->n_cols];
		l2norm = matd_l2norm(matd_getcol(Q, j));
		matd_col_smul_r(Q, j, 1/l2norm);
		R->data[j*R->n_cols+j] = l2norm;
		free_mat(aj);
	}
	qr->Q = Q;
	qr->R = R;
	free_mat(rkj);
	return qr;
}


#endif // SLAP_QR
