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
	// allocate the matrix structure and set the values:
	mat *A = mat_init2(3,3, (double[]){1,2,3,4,5,6,7,8,9});
	printf("A = \n"); mat_print(A); // print the matrix
	
	mat_free(A); // free the allocated memory
}
```

### Multiplication
```C++
#include <stdio.h>
#include "SLAP.h"

void main()
{
	mat *A, *I, *u, *v;
	
	A = mat_init2(3,3, (double[]){1,2,3,4,5,6,7,8,9});
	I = mat_init2(3,3, (double[]){1,0,0,0,1,0,0,0,1}); // identity matrix
	
	// scalar-matrix multiplication:
	mat_print(mat_smul(A, 2.5)); // print 2.5 * A
	
	// vector-vector multiplication:
	u = mat_init2(1,3, (double[]){1,2,3}); // row vector
	v = mat_init2(1,3, (double[]){2,2,2}); // row vector
	printf("u * v^T = ");  mat_print(mat_mul(u, mat_transpose(v))); // scalar product
	printf("u^T * v =\n"); mat_print(mat_mul(mat_transpose(u), v)); // tensor product
	
	// matrix-vector multiplication:
	mat_print(mat_mul(A,mat_transpose(u)));
	
	// matrix-matrix multiplication:
	mat_print(mat_mul(A,I)); // A * I
	mat_print(mat_mul(A,A)); // A * A
	
	mat_free(A); mat_free(I); mat_free(u); mat_free(v);
}
```

### Solve linear system
```C++
#include <stdio.h>
#include "SLAP.h"

void main()
{
	mat *A = mat_init2(3,3, (double[]){1,3,3,4,5,6,7,8,9});
	mat *b = mat_init2(3,1, (double[]){1,1,1});
	
	printf("A =\n"); mat_print(A);
	printf("b =\n"); mat_print(b);
	printf("Solve  A * x = b  for x\n\n");
	
	printf("LU(P) decomposition:\n");
	mat_lup *lu = mat_lup_solve(A); // LU(P) decomposition of matrix A
	mat* x = solve_lu(lu, b); // solve linear system using LU decomposition
	printf("x =\n"); mat_print(x);
	
	printf("\n\n\nGauss elimination of A:\n");
	print_mat(mat_GaussJordan(A));
	
	mat_free(A); mat_free(b);
}
```



# ToDo
- dividere la cartella `tests` da quella `examples`
