/*
	Conjugate gradient solver
*/
#ifndef SLAP_CONJ_GRAD
#define SLAP_CONJ_GRAD

inline TYPE first_member(mat *m){ return m->data[0]; }


mat* mat_conjgrad(mat *A, mat *b)
{
	// solve linear system A*x=b with conjugate gradient method
	mat *x = mat_new(b->n_rows,1); // column vector (all zero)
	mat *r, *p, *Ap;
	TYPE rold, rnew, alpha;
	int i;
	
	r = mat_mul(A, x); r = mat_sub(b, r);
	p = mat_copy(r);
	rold = first_member(mat_mul(mat_transpose(r), r)); // scalar product MEMORY LEAK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	for(i=0; i<x->n_rows; i++){
		Ap = mat_mul(A, p);
		alpha = rold / first_member(mat_mul(mat_transpose(p), Ap)); // MEMORY LEAK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		mat_add_r(x, mat_smul(p,  alpha)); // update x
		mat_sub_r(r, mat_smul(Ap, alpha)); // update r
		rnew = first_member(mat_mul(mat_transpose(r),r)); // MEMORY LEAK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if (sqrt(rnew) < 1e-10) break; // convergence on desired precision
		p = mat_add(r, mat_smul(p, rnew/rold)); // MEMORY LEAK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		rold = rnew; // update (r^T * r)
	}
	
	mat_free(r); mat_free(p); mat_free(Ap);
	return x;
}


#endif // SLAP_CONJ_GRAD
