/*
	vector product
	multiply two vector (row or column vector) with dimension 3 (1x3 or 3x1)
*/
mat* vec_prod(mat *u, mat *v) // CHIAMARE LA FUNZIONE "cross" ?????????????!!!!!!!!!!!!!
{
	mat* w;
	if((u->n_rows != v->n_rows) || (u->n_cols != v->n_cols)) return NULL; // ERROR!! must be same size
	if((u->n_rows != 1) && (u->n_cols != 1)) return NULL; // u must be a vector (at least one dimension have to be 1)
	if((v->n_rows != 1) && (v->n_cols != 1)) return NULL; // v must be a vector (at least one dimension have to be 1)
	if((mat_size(u) < 3) || (mat_size(u) < 3)) return NULL; // not big enough
	if(u->n_rows > u->n_cols) w = mat_new(3,1); // row vector
	else w = mat_new(1,3); // column vector
	
	w->data[0] = u->data[1]*v->data[2] - u->data[2]*v->data[1];
	w->data[1] = u->data[2]*v->data[0] - u->data[0]*v->data[2];
	w->data[2] = u->data[0]*v->data[1] - u->data[1]*v->data[0];
	
	return w;
}
