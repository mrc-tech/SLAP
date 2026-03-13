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
mat_qr* mat_qr_solve(const mat *m) // without memory leakage
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

mat_qr* mat_qr_solve_Gemini(const mat *m) 
{
	if (!m || m->n_cols == 0 || m->n_rows == 0) return NULL;// Controlli di sicurezza di base
	mat_qr *qr = mat_qr_new();
	mat *Q = mat_copy(m); // Q inizia come una copia speculare di A (m) e verrà modificata "in-place"
	
	// R è una matrice triangolare superiore.
	// Nota: idealmente R per una decomposizione "thin" è di dimensioni (n_cols x n_cols).
	// Mantengo le dimensioni originali per non rompere il resto del tuo codice.
	mat *R = mat_new(m->n_rows, m->n_cols); 

	int i, j, k;
	unsigned int rows = Q->n_rows;
	unsigned int cols = Q->n_cols;

	// Modified Gram-Schmidt (MGS) Algorithm
	/* Il passaggio dalla versione Gram-Schmidt Classica (CGS) alla Gram-Schmidt Modificata (MGS)
	consiste in una sottile differenza matematica . Nel metodo Classico, calcoli tutte le proiezioni
	rispetto ai vettori originali. Nel metodo Modificato, aggiorni i vettori in place mano a mano che
	trovi le basi ortogonali, riducendo immensamente l'errore di arrotondamento (cancellazione numerica). */
	for(k=0; k<cols; k++) {
		
		// 1. Calcolo della norma L2 della colonna k-esima di Q
		TYPE norm_sq = 0.0;
		for(i=0; i<rows; i++) {
			TYPE val = Q->data[i*cols+k];
			norm_sq += val * val;
		}
		TYPE l2norm = sqrt(norm_sq);
		R->data[k*R->n_cols+k] = l2norm; // R_{k,k} = ||q_k||

		// Prevenzione della divisione per zero in caso di matrici singolari/dipendenti
		if(l2norm < SLAP_MIN_COEF) {
			// Se la colonna è un vettore nullo, evitiamo NaN e saltiamo la normalizzazione
			l2norm = 1.0; 
		}else{
			// 2. Normalizzazione della colonna k-esima di Q: q_k = v_k / ||v_k||
			for(i=0; i<rows; i++) {
				Q->data[i*cols+k] /= l2norm;
			}
		}

		// 3. Ortogonalizzazione di tutte le colonne successive j contro la colonna ortogonale k
		for(j=k+1; j<cols; j++) {
			TYPE dot_product = 0.0;
			
			// Prodotto scalare: R_{k,j} = q_k^T * v_j
			for(i=0; i<rows; i++) {
				dot_product += Q->data[i*cols+k] * Q->data[i*cols+j];
			}
			R->data[k*R->n_cols+j] = dot_product;

			// Sottrazione della proiezione: v_j = v_j - R_{k,j} * q_k
			for(i=0; i<rows; i++) {
				Q->data[i*cols+j] -= dot_product * Q->data[i*cols+k];
			}
		}
	}

	qr->Q = Q;
	qr->R = R;
	
	return qr;
}
/*
- Zero Allocazioni all'interno del ciclo: La tua versione originale chiamava mat_getcol, mat_transpose, e mat_mul decine o centinaia di volte. Sotto al cofano, il C invocava malloc e free a ripetizione, operazione lentissima e prona alla frammentazione della memoria. Qui allochiamo solo Q e R in partenza: zero malloc nascoste nel core loop.
- Accesso diretto in Memoria (Pointer Math): Accedendo ai valori con Q->data[i * cols + j], facciamo calcoli scalari purissimi. Il compilatore (anche quelli vecchi per MS-DOS) riesce a ottimizzare o vettorizzare questi cicli nativamente (se compili con -O2 o -O3 in GCC).
- Stabilità (MGS): Noterai che l'indice j parte da k + 1 (guarda avanti) e modifichiamo direttamente Q per sottrarre le componenti. Questo è l'algoritmo Modified Gram-Schmidt, che preserva l'ortogonalità della matrice $Q$ molto meglio rispetto alla versione in cui k va da 0 a j.
- Gestione dello zero: Ho aggiunto un controllo if (l2norm < SLAP_MIN_COEF). Senza questo, se calcoli il QR di una matrice non a rango massimo (con vettori linearmente dipendenti), divideresti per zero ottenendo NaN in tutta la matrice.
*/


#endif // SLAP_QR
