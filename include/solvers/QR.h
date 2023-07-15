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
	mat_qr *qr = malloc(sizeof(*qr));
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
	TYPE doublesum = 0.0;
	if(m->n_cols != 1 || m->n_rows != 1) return -1; // only vectors
	for(i=0; i<MAX(m->n_rows,m->n_cols); i++){
		doublesum += (m->data[i]*m->data[i]);
	}
	return sqrt(doublesum);
}



mat_qr* mat_qr_solve(mat *m)
{
	mat_qr *qr = mat_qr_new();
	mat *Q = mat_copy(m);
	mat *R = mat_new(m->n_rows, m->n_cols); // n_cols and n_rows have to be equal
	
	int j, k;
	TYPE l2norm;
	mat *rkj = mat_new(1,1); // scalar
	mat *aj, *qk;
	for(j=0; j<m->n_cols; j++){
		rkj->data[0] = 0.0;
		aj = mat_getcol(m, j);
		for(k=0; k<j; k++){
//			rkj = nml_vect_dot(m, j, Q, k);
			rkj = mat_mul(mat_transpose(mat_getcol(m,j)), mat_getcol(Q,k)); // scalar product
			R->data[k*R->n_cols+j] = rkj->data[0];
			qk = mat_getcol(Q, k);
			mat_col_smul_r(qk, 0, rkj->data[0]);
			mat_sub_r(aj, qk);
			mat_free(qk);
		}
		for(k=0; k<Q->n_rows; k++) Q->data[k*Q->n_cols+j] = aj->data[k*aj->n_cols];
		l2norm = mat_l2norm(mat_getcol(Q, j));
		mat_col_smul_r(Q, j, 1/l2norm);
		R->data[j*R->n_cols+j] = l2norm;
		mat_free(aj);
	}
	qr->Q = Q;
	qr->R = R;
	mat_free(rkj);
	return qr;
}


#endif // SLAP_QR
