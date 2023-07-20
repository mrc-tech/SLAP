// matrix norm

double mat_norm_col(mat *m)
{
	// matrix norm: max column sum
	double sum;
	double max = 0; // max sum
	int r,c;
	
	for(c=0; c<m->n_cols; c++){
		sum = 0;
		for(r=0; r<m->n_rows; r++){
			sum += fabs(m->data[r*m->n_cols+c]);
		}
		if(sum > max) max = sum; // update max
	}
	
	return max;
}


double mat_norm_row(mat *m)
{
	// matrix norm: max row sum
	double sum;
	double max = 0; // max sum
	int r,c;
	
	for(r=0; r<m->n_rows; r++){
		sum = 0;
		for(c=0; c<m->n_cols; c++){
			sum += fabs(m->data[r*m->n_cols+c]);
		}
		if(sum > max) max = sum; // update max
	}
	
	return max;
}


double mat_norm_Frobenius(mat *m)
{
	// matrix norm: Frobenius
	mat *tmp = mat_transpose(m);
	tmp = mat_mul(m,tmp);
	double res = mat_trace(tmp);
	mat_free(tmp);
	return res;
}

