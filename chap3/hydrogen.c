#include <stdio.h>
#include <math.h>
#include "../headers/cmatrix.h"

const double pi = 3.141592653589793238462;
const double alpha[4] = {13.00773, 1.962079, 0.444529, 0.1219492};

extern void dsygv_(int *itype, char *jobz, char *uplo, int *n, 
    double *a, int *lda, double *b, int *ldb, double *w, double *work, int *lwork, int *info);

void call_dsygv(int n, double **a, double **b, double *eigenvals){
    int itype = 1;      // Az = λBz
    char jobz = 'N';    // only eigenvalues are computed
    char uplo = 'U';    // stored in upper triangles
    int lwork = 3*n;
    double *work = alloc_dvector(lwork);
    int info;

    dsygv_(&itype, &jobz, &uplo, &n, mat_ptr(a), &n, mat_ptr(b), &n, vec_ptr(eigenvals), vec_ptr(work), &lwork, &info);
    free_dvector(work);

    if(info < 0){
        fprintf(stderr, "argument %d had an illegal value.\n", -info);
        exit(1);
    }else if(1 <= info && info <= n){
        fprintf(stderr, "DSYGV failed to converge.\n");
        exit(1);
    }else if(info > n){
        fprintf(stderr, "returned error code %d.\n", info);
        exit(1);
    }
}

double calc_H(int m, int n){
    double T = 3.0*alpha[m]*alpha[n]*pow(pi, 1.5) / pow(alpha[m]+alpha[n], 2.5);
    double A = -2*pi / (alpha[m] + alpha[n]);
    return T + A;
}

double calc_S(int m, int n){
    return pow(pi/(alpha[m]+alpha[n]), 1.5);
}

int main(){
    int N = sizeof(alpha) / sizeof(alpha[0]);

    double **H = alloc_dmatrix(N, N);
    double **S = alloc_dmatrix(N, N);

    for(int i=0;i<N;i++){
        for(int j=i;j<N;j++){
            mat_elem(H, i, j) = calc_H(i, j);
            mat_elem(S, i, j) = calc_S(i, j);
        }
    }

    double *eigenvals = alloc_dvector(N);
    call_dsygv(N, H, S, eigenvals);

    printf("#計算結果\t解析解    \t差\n");
    printf("%.6f\t%.6f\t% .6f\n", eigenvals[0], -0.5, eigenvals[0]+0.5);

    free_dmatrix(H);
    free_dmatrix(S);
    free_dvector(eigenvals);

    return 0;
}
