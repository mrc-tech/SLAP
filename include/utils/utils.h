#ifndef SLAP_UTILS
#define SLAP_UTILS

#include <stdio.h> // per "FILE"


mat* mat_fromfile(FILE *f)
{
	mat *m;
	int r,c;
	unsigned int num_rows = 0, num_cols = 0;
	fscanf(f, "%d", &num_rows);
	fscanf(f, "%d", &num_cols);
	m = mat_new(num_rows, num_cols);
	for(r=0; r<m->n_rows; r++){
		for(c=0; c<m->n_cols; c++){
			fscanf(f, "%lf\t", &m->data[r * m->n_cols + c]);
		}
	}
	return m;
}

mat* mat_fromfilename(const char *file)
{
	mat *m;
	FILE *m_file = fopen(file, "r");
	if (!m_file) {
//		SLAP_FERROR(CANNOT_OPEN_FILE, file);
		return NULL;
	}
	m = mat_fromfile(m_file);
	fclose(m_file);
	return m;
}

#endif // SLAP_UTILS
