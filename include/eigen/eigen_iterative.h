/*
	Eigen-analysis using iterative methods
*/
#ifndef SLAP_EIGEN_ITERATIVE
#define SLAP_EIGEN_ITERATIVE




TYPE eigen_power_method(const mat *A, mat *eigenvec_out, TYPE tolerance, int max_iter) 
{
	/* * Trova l'autovalore dominante e il suo autovettore.
	 * A: Matrice quadrata
	 * eigenvec_out: (Opzionale) Vettore colonna pre-allocato nx1 dove salvare l'autovettore
	 * tolerance: Tolleranza per la convergenza (es. 1e-6)
	 * max_iter: Limite di sicurezza (es. 1000)
	 */
	int n = A->n_rows;
	mat *v = mat_new(n, 1);
	mat *w = mat_new(n, 1);
	int i, r, c, iter;
	TYPE lambda = 0.0, lambda_old = 0.0, norm;

	if(A->n_rows != A->n_cols) return 0.0; // Solo matrici quadrate

	// 1. Inizializziamo il vettore v con tutti 1.0 e lo normalizziamo
	for(i=0; i<n; i++) v->data[i] = 1.0;
	norm = mat_l2norm(v);
	for(i=0; i<n; i++) v->data[i] /= norm;

	for(iter = 0; iter < max_iter; iter++) {
		
		// 2. w = A * v (Fatto "inline" per evitare malloc/free lente in DOS)
		for(r = 0; r < n; r++) {
			w->data[r] = 0.0;
			for(c = 0; c < n; c++) {
				w->data[r] += A->data[r * A->n_cols + c] * v->data[c];
			}
		}

		// 3. Calcolo del Quoziente di Rayleigh (lambda = v^T * w)
		lambda = 0.0;
		for(i = 0; i < n; i++) {
			lambda += v->data[i] * w->data[i];
		}

		// 4. Controllo di convergenza
		if (fabs(lambda - lambda_old) < tolerance) break;
		lambda_old = lambda;

		// 5. Normalizziamo w per ottenere il nuovo v
		norm = mat_l2norm(w);
		if(norm < SLAP_ALMOST_ZERO) break; // Evita divisioni per zero
		for(i = 0; i < n; i++) {
			v->data[i] = w->data[i] / norm;
		}
	}

	// Se l'utente ha passato un vettore output, salviamo l'autovettore
	if(eigenvec_out && eigenvec_out->n_rows == n && eigenvec_out->n_cols == 1) {
		for(i = 0; i < n; i++) eigenvec_out->data[i] = v->data[i];
	}

	mat_free(v);
	mat_free(w);
	return lambda;
}


#endif // SLAP_EIGEN_ITERATIVE
