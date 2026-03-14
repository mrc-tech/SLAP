
mat* mat_inv(const mat *m) {
	mat_lup *lu;
	mat *inv, *e, *col_inv;
	unsigned int i, j;
	
	if(m->n_rows != m->n_cols) return NULL;
	
	lu = mat_lup_solve(m);
	if(!lu) return NULL; // Matrice singolare
	
	inv = mat_new(m->n_rows, m->n_cols);
	e = mat_new(m->n_rows, 1);
	
	for(i=0; i<m->n_cols; i++){
		// Genera il vettore colonna base e_i
		for(j=0; j<m->n_rows; j++){
			e->data[j] = (i == j) ? 1.0 : 0.0;
		}
		
		// Risolve il sistema per trovare la colonna i-esima dell'inversa
		col_inv = solve_lu(lu, e);
		
		if(col_inv){
			// Copia la soluzione nel posto giusto all'interno di inv
			for(j=0; j<m->n_rows; j++){
				inv->data[j * inv->n_cols + i] = col_inv->data[j];
			}
			mat_free(col_inv);
		}
	}
	
	mat_free(e);
	mat_lup_free(lu);
	return inv;
}


