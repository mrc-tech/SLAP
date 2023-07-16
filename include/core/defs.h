#ifndef SLAP_DEFS
#define SLAP_DEFS



#ifndef TYPE // data type of the matrix
#ifdef __MSDOS__ // MS-DOS system
#pragma message You are compiling using Borland C++ version __BORLANDC__.
#define SLAP_DOS
#define TYPE long double
#else // other non-DOS systems
#define TYPE double
#endif
#endif // TYPE


#ifndef SLAP_DEBUG
#define SLAP_DEBUG 0 // no debug
#endif


//#define NULL 0
#define SLAP_MIN_COEF 0.000000000000001 // DIPENDE DAL SISTEMA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#define MEM_CHECK(ptr) \
	if (!(ptr)) { \
		fprintf(stderr, "%s:%d NULL POINTER: %s n", __FILE__, __LINE__, (#ptr)); \
		exit(-1); \
	}


#define SWAP(a,b) \
	TYPE _temp = a; \
	a = b; \
	b = _temp;


#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
	



#endif // SLAP_DEFS
