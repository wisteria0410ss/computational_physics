#include <stdio.h>
#include <math.h>
#include "../headers/cmatrix.h"

const double pi = 3.141592653589793238462;

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
    if((n + m) % 2 == 1) return 0.0;
    return -8.0 * (1.0-m-n-2.0*m*n) / ((m+n+3.0)*(m+n+1.0)*(m+n-1.0));
}

double calc_S(int m, int n){
    if((n + m) % 2 == 1) return 0.0;
    return 2.0/(n+m+5.0) - 4.0/(n+m+3.0) + 2.0/(n+m+1.0);
}

int main(int argc, char** argv){
    int N = 5;
    if(argc >= 2) N = atoi(argv[1]);

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

    printf("#\t  計算結果  \t   解析解   \t     差\n");
    for(int i=0;i<N;i++){
        printf("%d\t%12.6f\t%12.6f\t% 12.6f\n", i+1, eigenvals[i], pow((i+1)*pi, 2)/4.0, eigenvals[i] - pow((i+1)*pi, 2)/4.0);
    }

    free_dmatrix(H);
    free_dmatrix(S);
    free_dvector(eigenvals);
    return 0;
}