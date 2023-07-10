/*
	Gauss Elimination
*/
#ifndef SLAP_GAUSS
#define SLAP_GAUSS



int matd_row_smul_r(matd *m, unsigned int row, double num)
{
	int i;
	if(row >= m->n_rows){
//		SLAP_FERROR(CANNOT_ROW_MULTIPLY, row, m->num_rows);
		return 0;
	}
	for(i=0; i<m->n_cols; i++) m->data[row*m->n_cols+i] *= num;
	return 1;
}

int matd_col_smul_r(matd *m, unsigned int col, double num)
{
	int i;
	if(col >= m->n_cols){
//		SLAP_FERROR(CANNOT_COL_MULTIPLY, row, m->num_rows);
		return 0;
	}
	for(i=0; i<m->n_rows; i++) m->data[i*m->n_cols+col] *= num;
	return 1;
}


int matd_row_addrow_r(matd *m, unsigned int where, unsigned int row, double multiplier)
{
	int i = 0;
	if(where >= m->n_rows || row >= m->n_rows){
//		SLAP_ERROR(CANNOT_ADD_TO_ROW, multiplier, row, where, m->num_rows);
		return 0;
	}
	for(i=0; i<m->n_cols; i++) m->data[where*m->n_cols+i] += multiplier * m->data[row*m->n_cols+i];
	return 1;
}


int matd_row_swap_r(matd *m, unsigned int row1, unsigned int row2)
{
	// swap two rows of matrix m
	int i;
	double tmp;
	if(row1 >= m->n_rows || row2 >= m->n_rows){
//		SLAP_ERROR(CANNOT_SWAP_ROWS, row1, row2, m->num_rows);
		return 0;
	}
	for(i=0; i<m->n_cols; i++){
		tmp = m->data[row2*m->n_cols+i];
		m->data[row2*m->n_cols+i] = m->data[row1*m->n_cols+i];
		m->data[row1*m->n_cols+i] = tmp;
	}
	return 1;
}


// Finds the first non-zero element on the col column, under the row row.
// Used to determine the pivot
// If not pivot is found, returns -1
int matd_pivot_id(matd *m, unsigned int col, unsigned int row)
{
	int i;
	for(i=row; i<m->n_rows; i++) if(fabs(m->data[i*m->n_cols+col]) > SLAP_MIN_COEF) return i;
	return -1;
}

// Find the max element from the column "col" under the row "row"
// This is needed to pivot in Gauss-Jordan elimination
// Return the maximum pivot for numerical stability. If pivot is not found, return -1
int matd_pivot_maxid(matd *m, unsigned int col, unsigned int row)
{
	int i, maxi;
	double micol;
	double max = fabs(m->data[row*m->n_cols+col]);
	maxi = row;
	for(i=row; i<m->n_rows; i++){
		micol = fabs(m->data[i*m->n_cols+col]);
		if(micol > max){
			max = micol;
			maxi = i;
		}
	}
	return(max < SLAP_MIN_COEF) ? -1 : maxi;
}


// Retrieves the matrix in Row Echelon form using Gauss Elimination
matd *matd_GaussJordan(matd *m)
{
	matd *r = matd_copy(m);
	int i=0, j=0, k, pivot;
	while(j < r->n_cols && i < r->n_cols){
		// Find the pivot - the first non-zero entry in the first column of the matrix
		pivot = matd_pivot_maxid(r, j, i);
		if(pivot<0){ // All elements on the column are zeros
			j++; // Move to the next column without doing anything
			continue;
		}
		if(pivot != i) matd_row_swap_r(r, i, pivot); // We interchange rows moving the pivot to the first row that doesn't have already a pivot in place
		matd_row_smul_r(r, i, 1/r->data[i*r->n_cols+j]); // Multiply each element in the pivot row by the inverse of the pivot
		for(k=i+1; k<r->n_rows; k++){
			if(fabs(r->data[k*r->n_cols+j]) > SLAP_MIN_COEF){
				matd_row_addrow_r(r, k, i, -(r->data[k*r->n_cols+j])); // We add multiplies of the pivot so every element on the column equals 0
			}
		}
		i++; j++;
	}
	return r;
}







#endif // SLAP_GAUSS
