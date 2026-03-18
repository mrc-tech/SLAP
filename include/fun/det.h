
TYPE mat_det(const mat *m)
{
	mat_lup *lu;
	TYPE det = 1.0;
	unsigned int i;
	
	if (m->n_rows != m->n_cols) return 0.0; // Solo matrici quadrate
	
	lu = mat_lup_solve(m);
	if(!lu) return 0.0; // Matrice singolare (non invertibile/degenere)
	
	for(i=0; i<lu->U->n_rows; i++){
		det *= lu->U->data[i * lu->U->n_cols + i];
	}
	
	// Se il numero di permutazioni delle righe e' dispari, il segno si inverte
	if(lu->num_permutations % 2 != 0){
		det = -det;
	}
	
	mat_lup_free(lu);
	return det;
}


