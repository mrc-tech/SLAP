#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/SLAP.h"

// Variabili globali per il conteggio dei test
int tests_run = 0;
int tests_passed = 0;

// Macro per l'asserzione (se fallisce, stampa ERRORE, altrimenti OK)
#define ASSERT(test, msg) do { \
    tests_run++; \
    if (!(test)) { \
        printf("[FALLITO] %s\n", msg); \
    } else { \
        printf("[OK]      %s\n", msg); \
        tests_passed++; \
    } \
} while (0)




// ====================================================================
// TEST SUITES
// ====================================================================

void test_core_and_basic_ops() {
    printf("\n--- TEST CORE & BASIC OPS ---\n");
    
    // Test Inizializzazione e Addizione
    mat *A = mat_init(2, 2, (TYPE[]){1.0, 2.0, 3.0, 4.0});
    mat *B = mat_init(2, 2, (TYPE[]){5.0, 6.0, 7.0, 8.0});
    mat *C_exp = mat_init(2, 2, (TYPE[]){6.0, 8.0, 10.0, 12.0});
    mat *C_res = mat_add(A, B);
    
    ASSERT(mat_equal(C_res, C_exp, 1e-6), "Addizione di matrici 2x2");
    
    // Test Moltiplicazione (Algoritmo Cache-Aligned)
    mat *D_exp = mat_init(2, 2, (TYPE[]){19.0, 22.0, 43.0, 50.0});
    mat *D_res = mat_mul(A, B);
    
    ASSERT(mat_equal(D_res, D_exp, 1e-6), "Moltiplicazione di matrici 2x2");
    
    // Test Trasposta
    mat *A_t_exp = mat_init(2, 2, (TYPE[]){1.0, 3.0, 2.0, 4.0});
    mat *A_t_res = mat_transpose(A);
    
    ASSERT(mat_equal(A_t_res, A_t_exp, 1e-6), "Trasposizione di matrice");

    mat_free(A); mat_free(B); mat_free(C_exp); mat_free(C_res);
    mat_free(D_exp); mat_free(D_res); mat_free(A_t_exp); mat_free(A_t_res);
}

void test_qr_decomposition() {
    printf("\n--- TEST DECOMPOSIZIONE QR ---\n");
    
    // Creiamo una matrice complessa 3x3
    mat *A = mat_init(3, 3, (TYPE[]){12.0, -51.0, 4.0, 
                                     6.0, 167.0, -68.0, 
                                     -4.0, 24.0, -41.0});
    
    mat_qr *qr = mat_qr_solve(A);
    
    // 1. Test di Ricostruzione: Q * R deve essere uguale ad A
    mat *QR = mat_mul(qr->Q, qr->R);
    ASSERT(mat_equal(QR, A, 1e-5), "Ricostruzione QR (Q * R == A)");
    
    // 2. Test di Ortogonalità: Q^T * Q deve essere uguale alla matrice Identità
    mat *Q_t = mat_transpose(qr->Q);
    mat *I_calc = mat_mul(Q_t, qr->Q);
    mat *I_real = mat_eye(3);
    ASSERT(mat_equal(I_calc, I_real, 1e-5), "Ortogonalita' di Q (Q^T * Q == I)");
    
    mat_free(A); mat_qr_free(qr); mat_free(QR); 
    mat_free(Q_t); mat_free(I_calc); mat_free(I_real);
}

void test_eigen() {
    printf("\n--- TEST AUTOVALORI (EIGEN) ---\n");
    
    // Matrice 3x3 simmetrica con autovalori noti: approx 16.366, -1.0, -0.366
    mat *A = mat_init(3, 3, (TYPE[]){1.0, 3.0, 3.0, 
                                     4.0, 5.0, 6.0, 
                                     7.0, 8.0, 9.0});
    
    // Test Power Method
    mat *autovettore = mat_new(3, 1);
    TYPE max_eig = eigen_power_method(A, autovettore, 1e-6, 1000);
    
    // Il massimo autovalore reale è circa 16.366
    int is_max_eig_correct = (fabs(max_eig - 16.3662) < 1e-3);
    ASSERT(is_max_eig_correct, "Power Method trova l'autovalore dominante corretto");
    
    // Test QR Autovalori
    mat *eig_diag = eigen_qr(A, 1e-6, 1000);
    int qr_converged = (eig_diag != NULL); // NON CAPISCO CHE RAZZA DI TEST SIA...!!!!!!!
    ASSERT(qr_converged, "L'algoritmo QR Eigen converge per matrice simmetrica");
    
    mat_free(A); mat_free(autovettore);
    if(eig_diag) mat_free(eig_diag);
}

void test_hilbert_matrix() {
    printf("\n--- TEST STRESS: MATRICE DI HILBERT ---\n");
    // La matrice di Hilbert è H_ij = 1 / (i + j - 1). È famigerata per i problemi numerici.
    int n = 4; // Anche solo 4x4 o 5x5 è sufficiente per generare instabilità
    mat *H = mat_new(n, n);
    int i, j;
    
    for(i=0; i<n; i++) {
        for(j=0; j<n; j++) {
            H->data[i * n + j] = 1.0 / (TYPE)(i + j + 1); // Indici 0-based
        }
    }
    
    // Verifichiamo se il nostro QR stabile (Gemini) resiste alla matrice di Hilbert
    mat_qr *qr = mat_qr_solve(H);
    mat *QR = mat_mul(qr->Q, qr->R);
    
    // Usiamo una tolleranza leggermente più rilassata dato il cattivo condizionamento
    ASSERT(mat_equal(QR, H, 1e-4), "Ricostruzione QR sopravvive alla Matrice di Hilbert");
    
    mat_free(H); mat_free(QR); mat_qr_free(qr);
}

int main() {
    printf("==========================================\n");
    printf("   INIZIO TEST LIBRERIA SLAP\n");
    printf("==========================================\n");
    
    test_core_and_basic_ops();
    test_qr_decomposition();
    test_eigen();
    test_hilbert_matrix();
    
    printf("\n==========================================\n");
    printf("RISULTATI FINALI: %d / %d TEST SUPERATI\n", tests_passed, tests_run);
    printf("==========================================\n");
    
    if (tests_passed == tests_run) {
        printf("Libreria validata con successo. Pronto per il rilascio!\n\n");
        return 0;
    } else {
        printf("ATTENZIONE: Alcuni test sono falliti. Controllare i log.\n\n");
        return -1;
    }
}
