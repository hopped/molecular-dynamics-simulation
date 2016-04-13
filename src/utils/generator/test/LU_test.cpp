#include "LU.h"

#include "LA.h"

int main() {
    double **A = new double*[3];
    double **L = new double*[3];
    double **R = new double*[3];
    for(int i = 0; i < 3; i++) {
        A[i] = new double[3];
        L[i] = new double[3];
        R[i] = new double[3];
    }
    /* A = [
     *      [1.000000,-1.000000,-1.000000],
     *      [-2.000000,6.000000,3.000000],
     *      [-1.000000,13.000000,6.000000]
     * ]
     */
    A[0][0] = 1;
    A[0][1] = -1;
    A[0][2] = -1;
    A[1][0] = -2;
    A[1][1] = 6;
    A[1][2] = 3;
    A[2][0] = -1;
    A[2][1] = 13;
    A[2][2] = 6;
    int n = 3;
    print_matrix(A,n);

    LU(A, n, L, R);
    print_matrix(L,n);
    print_matrix(R,n);
    return 0;
}
