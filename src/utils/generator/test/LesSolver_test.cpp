#include "LesSolver.h"

#include "LA.h"
#include "LU.h"


int main() {
    int n = 3;
    double **A = new double*[3];
    double *x = new double[3];
    for(int i = 0; i < 3; i++) {
        A[i] = new double[3];
    }
    A[0][0] = 4;
    A[0][1] = -1;
    A[0][2] = -1;
    A[1][0] = 8;
    A[1][1] = 0;
    A[1][2] = -1;
    A[2][0] = 4;
    A[2][1] = 1;
    A[2][2] = 4;
    double *a = new double[3];
    a[0] = 2;
    a[1] = 7;
    a[2] = 9;

    print_vector(a,n);
    print_matrix(A,n);
    LesSolve(A, a, n, x);
    print_vector(x,n);
    return 0;
}

