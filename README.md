<!--# SLAP
Simple Linear Algebra Package (SLAP) -->
![cover](doc/cover.png)

## Features
- Written in `C` for maximum compatibility among various systems
- No external libraries (all-in-one header file)
- Small enough to fit inside **MS-DOS** (and eventually embedded systems)
- Tailored to be used in _Finite Element_ software


# Functions

All the functions that return a matrix allocate new memory for the returned matrix.
| function | operation | math |
| --- | --- | --- |
| `mat_new(n,m)` | allocate memory for $n\times m$ matrix ($n$ rows and $m$ cols); set all values to zero | |
| `mat_free(A)` | free the allocated memory for matrix $\mathbf{A}$ | |
| `B = mat_copy(A)` | copy the matrix $\mathbf{A}$ into $\mathbf{B}$ allocating memory | |
| `mat_add_r(A,B)` | add matrix $\mathbf{A}$ and $\mathbf{B}$ and put the result in $\mathbf{A}$ (_reference_) | $\mathbf{A}=\mathbf{A}+\mathbf{B}$ |
| `mat_sub_r(A,B)` | subtract matrix $\mathbf{A}$ and $\mathbf{B}$ and put the result in $\mathbf{A}$ (_reference_) | $\mathbf{A}=\mathbf{A}-\mathbf{B}$ |
| `C = mat_add(A,B)` | add matrix $\mathbf{A}$ and $\mathbf{B}$ and put the result in $\mathbf{C}$ | $\mathbf{C}=\mathbf{A}+\mathbf{B}$ |
| `C = mat_sub(A,B)` | subtract matrix $\mathbf{A}$ and $\mathbf{B}$ and put the result in $\mathbf{C}$ | $\mathbf{C}=\mathbf{A}-\mathbf{B}$ |
| `B = mat_smul(A,c)` | scale the matrix $\mathbf{A}$ with the scalar $c$ and put into $\mathbf{C}$ | $\mathbf{B}=c\mathbf{A}$ |
| `mat_smul_r(A,c)` | scale the matrix $\mathbf{A}$ with the scalar $c$ and put into $\mathbf{A}$ (_reference_) | $\mathbf{A}=c\mathbf{A}$ |
| `C = mat_mul(A,B)` | multiply $\mathbf{A}$ and $\mathbf{B}$ and put the result in $\mathbf{C}$; also works with scalar (row times column vector) and tensor (column times row vector) product | $\mathbf{A}=\mathbf{A}\mathbf{B}\quad$  $\mathbf{u}\cdot\mathbf{v}=\mathbf{u}^T\mathbf{v}\quad$   $\mathbf{u}\otimes\mathbf{v}=\mathbf{u}\mathbf{v}^T$ |
| `mat_equal(A,B,tol)` | check if matrix $\mathbf{A}$ and $\mathbf{B}$ are equal | true if $abs(A_{ij}-B_{ij}) < tol$ |
| `B = mat_transpose(A)` | transpose the matrix $\mathbf{A}$ and put into matrix $\mathbf{B}$ | $\mathbf{B}=\mathbf{A}^T$ |
| `mat_transpose_r(A)` | transpose the matrix $\mathbf{A}$ without allocating memory (_reference_) | $\mathbf{A}=\mathbf{A}^T$ |
| `I = mat_eye(n)` | create an $n\times n$ identity matrix | |
| `U = mat_GaussJordan(A)` | transform matrix $\mathbf{A}$ in row echelon form $\mathbf{U}$ (upper triangular) via Gauss elimination | |
| `a = eigen_qr(A)` | calculate eigenvalues of $\mathbf{A}$ with QR decomposition and put the result in $\mathbf{a}$ | |
| `mat_print(A)` | print matrix $\mathbf{A}$ | |



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
- [ ] **CORE**
	- [ ] BTC++ `mat_init_DOS`
- [ ] **BASIC OPERATIONS**
	- [ ] trace
	- [ ] diagonal square matrix from vector
	- [ ] rename `smul` into `scale`?
	- [ ] ...
- [ ] **DECOMPOSITION**
	- [ ] QR
		- [ ] Householder method
		- [ ] Gibs rotations?
	- [ ] Singular Value Decomposition (SVD)
	- [ ] Principal Component Analysis (PCA)
	- [ ] ...
- [ ] **ADVANCED OPERATIONS**
	- [ ] determinant
	- [ ] inverse
	- [ ] Hessenberg form
	- [ ] FFT
	- [ ] ...
- [ ] **SOLVER**
	- [ ] controllare la stabilità (e la velocità)
	- [ ] QR decomposition
		- [ ] devo solamente implementare la routine che risolve il sistema
	- [ ] iterative?
	- [ ] ...
- [ ] **EIGEN**
	- [ ] QR decomposition
		- [ ] definire meglio quando finire la procedura iterativa
		- [ ] forma di Hessenberg per aumentare l'efficienza
	- [ ] Iterative methods
	- [ ] ...
- [ ] **UTILS**
	- [ ] random number generator
	- [ ] ...