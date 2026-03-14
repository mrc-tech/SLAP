/*
	Eigen-analysis using QR decomposition
*/
#ifndef SLAP_EIGEN_QR
#define SLAP_EIGEN_QR


mat* eigen_qr(const mat *m, TYPE tolerance, int max_iter)
{
	mat *A = mat_copy(m);
	mat_qr *qr;
	int i, r, c;
	short converged;
	
	if(m->n_rows != m->n_cols) {
		mat_free(A);
		return NULL; // ERROR!! not square matrix
	}
	
	for(i=0; i<max_iter; i++){
		qr = mat_qr_solve_Gemini(A);
		mat_free(A); // avoid memory leakage
		A = mat_mul(qr->R, qr->Q); // update A matrix for i-th step
		mat_qr_free(qr); // free used memory for QR decomposition
		
		converged = 1; // Assumiamo sia convergente
		
		// Esaminiamo solo la parte SOTTO la diagonale principale
		for(c=0; c<A->n_cols-1; c++) {
			for(r=c+1; r<A->n_rows; r++) {
				if(fabs(A->data[r*A->n_cols+c]) > tolerance) {
					converged = 0; // Trovato un elemento troppo grande, non ha ancora finito!
					break; 
				}
			}
			if(!converged) break;
		}

		if(converged) {
			#if SLAP_DEBUG
			printf("QR Autovalori convergenza raggiunta in %d iterazioni.\n", iter);
			#endif
			break;
		}
	}
//	mat_print(A);
	mat *diag = mat_get_diag(A);
	mat_free(A); // libera la memoria prima di uscire
	return diag;
}


#endif // SLAP_EIGEN_QR
