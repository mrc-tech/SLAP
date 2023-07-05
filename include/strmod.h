/*
	matrix structure modification
*/
#ifndef SLAP_STRMOD
#define SLAP_STRMOD


matd* matd_remcol(matd *m, unsigned int column)
{
	// remove the i-th column (start counting from zero)
	if(column >= m->n_cols){
//    	SLAP_FERROR(CANNOT_REMOVE_COLUMN, column, m->num_cols);
		return 0; // return NULL
	}
	matd *ret = new_matd(m->n_rows, m->n_cols-1);
	int i, j, k;
	for(i=0; i<m->n_rows; i++){
		for(j=0,k=0; j<m->n_cols; j++){
			if(column != j) ret->data[i*ret->n_cols + k++] = m->data[i*m->n_cols + j];
		}
	}
	return ret;
}

matd* matd_remrow(matd *m, unsigned int row)
{
	// remove the i-th row (start counting from zero)
	if(row >= m->n_rows){
//    	SLAP_FERROR(CANNOT_REMOVE_ROW, row, m->num_rows);
		return 0; // return NULL
	}
	matd *ret = new_matd(m->n_rows-1, m->n_cols);
	int i, j, k;
	for(i=0,k=0; i<m->n_rows; i++){
		if(row != i){
			for(j=0; j<m->n_cols; j++){
				ret->data[k*ret->n_cols + j] = m->data[i*ret->n_cols + j];
			}
			k++;
		}
	}
	return ret;
}



matd* matd_catver(int N, matd **marr) // NON VERIFICATO!!!!!!
{
	// concatenate vertically N matrices
	if (N == 0) return 0; // No matrices, nothing to return
	if (N == 1) return matd_copy(marr[0]); // no need for additional computations
    
	// We calculate the total number of columns to know how to allocate memory for the resulting matrix:
	int i,j,k,offset;
	unsigned int lrow, ncols;
	lrow  = marr[0]->n_rows;
	ncols = marr[0]->n_cols;
	for(k=1; k<N; k++){
		if (NULL == marr[k]){
//			SLAP_FERROR(INCONSITENT_ARRAY, k, mnum);
			return 0;
		}
		if (lrow != marr[k]->n_rows){
//			SLAP_FERROR(CANNOT_CONCATENATE_H, lrow, marr[k]->num_rows);
			return 0;
		}
		ncols += marr[k]->n_cols;
	}
	// allocate memory for the resulting matrix
	matd *r = new_matd(lrow, ncols);
	for(i = 0; i<r->n_rows; i++){
		k = 0;
		offset = 0;
		for(j = 0; j<r->n_cols; j++){
			// If the column index of marr[k] overflows
			if(j-offset == marr[k]->n_cols){
				offset += marr[k]->n_cols;
				k++; // jump to the next matrix in the array
			}
		r->data[i*r->n_cols+j] = marr[k]->data[i*r->n_cols + j - offset];
		}
	}
	return r;
}


matd* matd_cathor(unsigned int N, matd **marr) // NON VERIFICATO!!!!!
{
	// concatenate matrices horizontally (same number of rows, aumented number of columns)
	if(N == 0) return 0;
	if(N == 1) return matd_copy(marr[0]);
	
	int lcol, i, j, k, offset;
	unsigned int numrows = 0;
	matd *res;
	lcol = marr[0]->n_cols;
	for(i=0; i<N; i++){
		if(marr[i] == 0){
//			SLAP_FERROR(INCONSITENT_ARRAY, i, mnum);
			return 0;
		}
		if(lcol != marr[i]->n_cols){
//			SLAP_FERROR(CANNOT_CONCATENATE_V,lcol,marr[i]->num_cols);
			return 0;
		}
		numrows += marr[i]->n_rows;
	}
	res = new_matd(numrows, lcol);
	for(j=0; j<res->n_cols; j++){
		offset = 0;
		k = 0;
		for(i=0; i<res->n_rows; i++){
			if(i - offset == marr[k]->n_rows){
				offset += marr[k]->n_rows;
				k++;
			}
			res->data[i * res->n_cols + j] = marr[k]->data[(i-offset) * res->n_cols + j];
		}
	}
	return res;
}



#endif // SLAP_STRMOD
