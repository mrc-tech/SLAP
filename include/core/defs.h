#ifndef SLAP_DEFS
#define SLAP_DEFS


//#define NULL 0
#define SLAP_MIN_COEF 0.000000000000001 // DIPENDE DAL SISTEMA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#define MEM_CHECK(ptr) \
	if (!(ptr)) { \
		fprintf(stderr, "%s:%d NULL POINTER: %s n", __FILE__, __LINE__, (#ptr)); \
		exit(-1); \
	}


#define SWAP(a,b) \
	double _temp = a; \
	a = b; \
	b = _temp;


#endif // SLAP_DEFS
