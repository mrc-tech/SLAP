/*
	Conjugate gradient solver
*/
#ifndef SLAP_CONJ_GRAD
#define SLAP_CONJ_GRAD

inline TYPE first_member(mat *m){ return m->data[0]; }


mat* mat_conjgrad(const mat *A, const mat *b)
{
	/* * Risolve il sistema lineare Ax = b usando il Gradiente Coniugato.
	 * A deve essere Simmetrica e Definita Positiva.
	 * x è il vettore "initial guess" (di solito inizializzato a zero) e 
	 * conterrà il risultato finale.
	 */
	int n = A->n_rows;
	int i, k, c, r;
	int max_iter = n * 2; // Limite di sicurezza classico per CG
	mat *x = mat_new(b->n_rows,1); // column vector (all zero)
	
	TYPE alpha, beta, r_dot_r, r_dot_r_new, p_Ap_dot, sum;
	TYPE norm_b = sqrt(mat_dot(b, b)); // Prima del ciclo calcola la norma di b
	
	mat *rv = mat_new(1, n); // Allocazioni iniziali (vengono fatte UNA SOLA VOLTA)
	mat *p  = mat_new(1, n); // vettori riga sono piu' veloci in row-major
	mat *Ap = mat_new(1, n); 
	
	// 1. Inizializzazione: r = b - A*x
	mat *Ax = mat_mul(A, x);
	for(i = 0; i < n; i++) {
		rv->data[i] = b->data[i] - Ax->data[i];
		p->data[i] = rv->data[i]; // All'inizio, la direzione p_0 è uguale al residuo r_0
	}
	mat_free(Ax); // Liberiamo subito Ax, non ci serve più
	
	r_dot_r = mat_dot(rv, rv); // Calcoliamo il prodotto scalare r^T * r iniziale
	
	if (norm_b == 0.0) norm_b = 1.0; // Evita divisioni per zero
	
	// CICLO PRINCIPALE
	for(k=0; k<max_iter; k++){
		
		if (sqrt(r_dot_r) / norm_b < SLAP_MIN_COEF) break; // Condizione di uscita: il residuo è abbastanza vicino a zero?
		// in questo caso SLAP_MIN_COEF e' usata come una tolleranza normalizzata
		
		// -----------------------------------------------------------
		// CALCOLO DI Ap = A * p (INLINE E IN-PLACE)
		// -----------------------------------------------------------
		// Niente chiamate a funzione, niente malloc. Calcolo puro.
		for(r = 0; r < n; r++) {
			sum = 0.0;
			for(c=0; c<n; c++) sum += A->data[r*A->n_cols+c] * p->data[c];
			Ap->data[r] = sum; // Scriviamo nel vettore pre-allocato!
		}
		// -----------------------------------------------------------
		
		p_Ap_dot = mat_dot(p, Ap); // Prodotto scalare p^T * A * p
		
		// Sicurezza per evitare divisioni per zero se la matrice è malcondizionata
		if (fabs(p_Ap_dot) < SLAP_ALMOST_ZERO) {
			mat_free(Ap);
			break;
		}
		alpha = r_dot_r / p_Ap_dot; // alpha = (r^T * r) / (p^T * A * p)
		
		// AGGIORNAMENTO VETTORIALE INLINE
		// Sostituisce mat_add, mat_sub, mat_scale risparmiando 4 allocazioni a ciclo!
		for(i = 0; i < n; i++) {
			x ->data[i] += alpha * p->data[i];  // x_{k+1} = x_k + alpha * p_k
			rv->data[i] -= alpha * Ap->data[i]; // r_{k+1} = r_k - alpha * A * p_k
		}
		r_dot_r_new = mat_dot(rv, rv); // r_{k+1}^T * r_{k+1}
		beta = r_dot_r_new / r_dot_r; // beta = (r_{k+1}^T * r_{k+1}) / (r_k^T * r_k)
		
		// AGGIORNAMENTO DIREZIONE INLINE
		for(i = 0; i < n; i++) {
			p->data[i] = rv->data[i] + beta * p->data[i]; // p_{k+1} = r_{k+1} + beta * p_k
		}
		r_dot_r = r_dot_r_new; // Prepariamo r_dot_r per il prossimo ciclo
	}
	
	// Pulizia finale della memoria
	mat_free(rv);
	mat_free(p);
	mat_free(Ap);
	
	return x;
}


//mat* mat_conjgrad(mat *A, mat *b) // forse dovrei usare const mat* ...
//{
//	// solve linear system A*x=b with conjugate gradient method
//	mat *x = mat_new(b->n_rows,1); // column vector (all zero)
//	mat *r, *p, *Ap;
//	TYPE rold, rnew, alpha;
//	int i;
//	
//	r = mat_mul(A, x); r = mat_sub(b, r);
//	p = mat_copy(r);
//	rold = first_member(mat_mul(mat_transpose(r), r)); // scalar product MEMORY LEAK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//	
//	for(i=0; i<x->n_rows; i++){
//		Ap = mat_mul(A, p);
//		alpha = rold / first_member(mat_mul(mat_transpose(p), Ap)); // MEMORY LEAK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//		mat_add_r(x, mat_scale(p,  alpha)); // update x
//		mat_sub_r(r, mat_scale(Ap, alpha)); // update r
//		rnew = first_member(mat_mul(mat_transpose(r),r)); // MEMORY LEAK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//		if (sqrt(rnew) < 1e-10) break; // convergence on desired precision
//		p = mat_add(r, mat_scale(p, rnew/rold)); // MEMORY LEAK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//		rold = rnew; // update (r^T * r)
//	}
//	
//	mat_free(r); mat_free(p); mat_free(Ap);
//	return x;
//}


#endif // SLAP_CONJ_GRAD
