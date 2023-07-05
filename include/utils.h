#ifndef SLAP_UTILS
#define SLAP_UTILS

#include <stdio.h> // per "FILE"


matd* matd_fromfile(FILE *f)
{
	int r,c;
	unsigned int num_rows = 0, num_cols = 0;
	fscanf(f, "%d", &num_rows);
	fscanf(f, "%d", &num_cols);
	matd *m = new_matd(num_rows, num_cols);
	for(r=0; r<m->n_rows; r++){
		for(c=0; c<m->n_cols; c++){
			fscanf(f, "%lf\t", &m->data[r * m->n_cols + c]);
		}
	}
	return m;
}

matd* matd_fromfilename(const char *file)
{
	FILE *m_file = fopen(file, "r");
	if (!m_file) {
//		SLAP_FERROR(CANNOT_OPEN_FILE, file);
		return 0;
	}
	matd *r = matd_fromfile(m_file);
	fclose(m_file);
	return r;
}

#endif // SLAP_UTILS
