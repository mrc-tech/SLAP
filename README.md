# SLAP
Simple Linear Algebra Package (SLAP)

## Features
- Written in `C` for maximum compatibility among various systems
- No external libraries (all-in-one header file)
- Small enough to fit inside **MS-DOS** (and eventually embedded systems)
- Tailored to be used in _Finite Element_ software

# Examples

### Matrix creation
```C++
#include <stdio.h>
#include "SLAP.h"

void main()
{
	matd *A = matd_init2(3,3, (double[]){1,2,3,4,5,6,7,8,9}); // allocate the matrix structure and set the values
  
	printf("A = \n"); print_mat(A); // print the matrix
	
	free_mat(A); // free the allocated memory
}
```

### Multiplication
```C++
#include <stdio.h>
#include "SLAP.h"

void main()
{
	matd *A, *I, *u, *v;
	
	A = matd_init2(3,3, (double[]){1,2,3,4,5,6,7,8,9});
	I = matd_init2(3,3, (double[]){1,0,0,0,1,0,0,0,1}); // identity matrix
	
	// scalar-matrix multiplication:
	print_mat(matd_smul(A, 2.5)); // print 2.5 * A
	
	// vector-vector multiplication:
	u = new_matd(1,3); u = matd_init2(1,3, (double[]){1,2,3}); // row vector
	v = new_matd(1,3); v = matd_init2(1,3, (double[]){2,2,2}); // row vector
	printf("u * v^T = ");  print_mat(matd_mul(u, matd_transpose(v))); // scalar product
	printf("u^T * v =\n"); print_mat(matd_mul(matd_transpose(u), v)); // tensor product
	
	// matrix-vector multiplication:
	print_mat(matd_mul(A,matd_transpose(u)));
	
	// matrix-matrix multiplication:
	print_mat(matd_mul(A,I)); // A * I
	print_mat(matd_mul(A,A)); // A * A
	
	free_mat(A);
	free_mat(I);
	free_mat(u);
	free_mat(v);
}
```

### Solve linear system
```C++
#include <stdio.h>
#include "SLAP.h"


int main()
{
	matd *A = matd_init2(3,3, (double[]){1,3,3,4,5,6,7,8,9});
	matd *b = matd_init2(3,1, (double[]){1,1,1});
	
	printf("A =\n"); print_mat(A);
	printf("b =\n"); print_mat(b);
	printf("Solve  A * x = b  for x\n\n");
	
	printf("LU(P) decomposition:\n");
	matd_lup *lu = matd_lup_solve(A);
	matd* x = lu_solve(lu, b);
	printf("x =\n"); print_mat(x);
	
	printf("\n\n\nGauss elimination of A:\n");
	print_mat(matd_GaussJordan(A));
  
	free_mat(A); free_mat(b);
}
```



# ToDo
- simple script that buld an ASCII file assessing the internal header-files structure (do inside the `header_merger`?)
- rendere coerenti i nomi delle funzioni (ad esempio mettendo `matd_` prima di ogni operazione sulle matrici double e poi il nome della funzione. Come ad esempio `matd_new`, `matd_free`, `matd_equal`, etc.)
- dividere la cartella `tests` da quella `examples`
- immagine di presentazione 1280×640px fatta carina con `SIL` che possa rendere accattivante cliccare su SLAP
- esempi nel README.md per come fare delle operazioni base con la libreria (come ad esempio risolvere un piccolo sistema lineare o qualche operazione base sui vettori)
