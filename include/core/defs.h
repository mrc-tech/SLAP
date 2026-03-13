#ifndef SLAP_DEFS
#define SLAP_DEFS



#include <float.h> // <-- FONDAMENTALE per i limiti numerici

#ifndef TYPE // data type of the matrix
    #ifdef __MSDOS__ // MS-DOS system
        #pragma message You are compiling using Borland C++ version __BORLANDC__.
        #define SLAP_DOS
        #define TYPE long double
        #define SLAP_FORMAT "%Lf" // L maiuscola per i long double
        
        // In MS-DOS long double è spesso 80-bit.
        // LDBL_EPSILON è circa 1.08e-19. 
        // Moltiplichiamo per 10 o 100 per avere una tolleranza sicura per gli arrotondamenti.
        #define SLAP_MIN_COEF (LDBL_EPSILON * 100.0L) 
        /* Se fai operazioni su una matrice (es. Gauss-Jordan o decomposizione QR), gli errori in 
		virgola mobile si accumulano a ogni moltiplicazione e addizione. Se la tua matrice è, per 
		esempio, 10x10 o 50x50, il rumore di fondo crescerà. Usare Epsilon * 100 (che per un double 
		moderno equivale a dire "considera zero tutto ciò che è più piccolo di ~2e-14") garantisce 
		che l'algoritmo non scambi del rumore di arrotondamento per un pivot valido, 
		evitando divisioni disastrose.*/
        
    #else // other non-DOS systems
        #define TYPE double
        #define SLAP_FORMAT "%lf" // l minuscola per i double
        
        // Nei sistemi moderni double è 64-bit.
        // DBL_EPSILON è circa 2.22e-16.
        #define SLAP_MIN_COEF (DBL_EPSILON * 100.0)
        
    #endif
#endif // TYPE

#ifndef SLAP_DEBUG
#define SLAP_DEBUG 0 // no debug
#endif



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
