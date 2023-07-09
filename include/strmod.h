/*
	matrix structure modification
*/
#ifndef SLAP_STRMOD
#define SLAP_STRMOD


matd* matd_remcol(matd *m, unsigned int column)
{
	// remove the i-th column (start counting from zero)
	matd *ret;
	int i, j, k;
	if(column >= m->n_cols){
//    	SLAP_FERROR(CANNOT_REMOVE_COLUMN, column, m->num_cols);
		return 0; // return NULL
	}
	ret = new_matd(m->n_rows, m->n_cols-1);
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
	matd *ret;
	int i, j, k;
	if(row >= m->n_rows){
//    	SLAP_FERROR(CANNOT_REMOVE_ROW, row, m->num_rows);
		return 0; // return NULL
	}
	ret = new_matd(m->n_rows-1, m->n_cols);
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



matd *matd_getcol(matd *m, unsigned int col)
{
	// return matrix column
	int j;
	matd *res;
	if(col >= m->n_cols){
//		SLAP_FERROR(CANNOT_GET_COLUMN, col, m->num_cols);
		return NULL;
	}
	res = new_matd(m->n_rows, 1);
	for(j=0; j<res->n_rows; j++) res->data[j*res->n_cols] = m->data[j*m->n_cols+col];
	return res;
}

double *matd_getcol_array(matd *m, unsigned int col)
{
	// return column via array
	int i;
	double *res;
	if(col >= m->n_cols) return NULL;
	res = (double*)malloc(m->n_rows * sizeof(double));
	for(i=0; i<m->n_rows; i++) res[i] = m->data[i*m->n_cols+col];
	return res;
}

matd *matd_getrow(matd *m, unsigned int row)
{
	// return matrix row
	matd *res;
	if(row >= m->n_rows){
//		SLAP_FERROR(CANNOT_GET_ROW, row, m->num_rows);
		return NULL;
	}
	res = new_matd(1, m->n_cols);
	memcpy(&res->data[0], &m->data[row*m->n_cols], m->n_cols * sizeof(res->data[0]));
	return res;
}

double *matd_getrow_array(matd *m, unsigned int row)
{
	// return a row via array
	double *res;
	if(row >= m->n_rows) return NULL;
	res = (double*)malloc(m->n_cols * sizeof(double));
	memcpy(&res, &m->data[row*m->n_cols], m->n_cols * sizeof(double));
	return res;
} 



matd* matd_cathor(int N, matd **marr) // NON VERIFICATO!!!!!!
{
	// concatenate matrices horizontally (same number of rows, aumented number of columns)
	matd *m;
	int i, j, k, offset;
	unsigned int lrow, ncols;
	if (N == 0) return 0; // No matrices, nothing to return
	if (N == 1) return matd_copy(marr[0]); // no need for additional computations
    
	// We calculate the total number of columns to know how to allocate memory for the resulting matrix:
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
	m = new_matd(lrow, ncols);
	for(i = 0; i<m->n_rows; i++){
		k = 0;
		offset = 0;
		for(j = 0; j<m->n_cols; j++){
			// If the column index of marr[k] overflows
			if(j-offset == marr[k]->n_cols){
				offset += marr[k]->n_cols;
				k++; // jump to the next matrix in the array
			}
		m->data[i*m->n_cols+j] = marr[k]->data[i*m->n_cols + j - offset];
		}
	}
	return m;
}


matd* matd_catver(unsigned int N, matd **marr)
{
	// concatenate vertically N matrices
	matd *res;
	unsigned int numrows = 0;
	int lcol, i, j, k, offset;
	if(N == 0) return 0;
	if(N == 1) return matd_copy(marr[0]);
	
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
