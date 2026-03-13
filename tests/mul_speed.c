#include <stdio.h>
#include <time.h>

#include "../include/SLAP.h" // Simple Linear Algebra Package

int main()
{
	int i;
	mat *A, *B;
	clock_t start, end;
    double tempo;
    unsigned N = 7000; // dimensione della matrice da moltiplicare
    // Fare il test della moltiplicazione per una matrice di dimensione maggiore di 50000 crasha!
    // con N=20000 il programma occupa 6 Gb di RAM... Non mi stupisce che crasha!
    
//	A = mat_init(3,3, (double[]){1,2,3,4,5,6,7,8,9});
	A = mat_new(N,N);
	for(i=0; i<N*N; i++) A->data[i] = 1; // setta tutto a 1 per evitare lo skip a zero
	
    start = clock();
    // ------------------------------
	B = mat_mul(A,A);
//	for(i=0; i<100; i++){ B = mat_mul(A,A); }
    // ------------------------------
	end = clock();
	
    tempo = (double)(end - start) / CLOCKS_PER_SEC;
    printf("N = %d\n", N);
    printf("Tempo impiegato: %f millisecondi\n", tempo*1000.0);
//	getch();
	
	return 0;
}
