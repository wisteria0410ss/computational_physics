#include <stdio.h>
#include <math.h>
#include <time.h>
#include "../headers/cmatrix.h"

extern void dgemv_(char *trans, int *m, int *n, 
    double *alph, double *A, int *ldA, double *x, int *incx,
    double *beta , double *y, int *incy);
extern double ddot_(int *n, double *x, int *incx, double *y, int *incy);
extern void dsygv_(int *itype, char *jobz, char *uplo, int *n, 
    double *a, int *lda, double *b, int *ldb, double *w, double *work, int *lwork, int *info);

void dgemv(int n, double **a, double *x, double *y){
    char trans = 'N';
    double alph = 1.0;
    int ld = n;
    int inc = 1;
    double beta = 0.0;
    dgemv_(&trans, &n, &n, &alph, mat_ptr(a), &ld, vec_ptr(x), &inc, &beta, vec_ptr(y), &inc);
}

double ddot(int n, double *x, double *y){
    int inc = 1;
    return ddot_(&n ,vec_ptr(x), &inc, vec_ptr(y), &inc);
}

void dsygv(int n, double **a, double **b, double *eigenvals){
    int itype = 1;      // Az = λBz
    char jobz = 'V';    // eigenvalues and eigenvectors are computed
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

const double pi = 3.14159265358979324;
const double alpha[] = {0.297104, 1.236745, 5.749982, 38.216677};
const int n_basis = sizeof(alpha) / sizeof(alpha[0]);

double get_h(int p, int q){
    double T = 3.0*alpha[p]*alpha[q]*pow(pi, 1.5) / pow(alpha[p]+alpha[q], 2.5);
    double A = -2*pi / (alpha[p] + alpha[q]);
    return T + 2.0*A;
}

double get_S(int p, int q){
    return pow(pi/(alpha[p]+alpha[q]), 1.5);
}

double get_Q(int p, int r, int q, int s){
    return 2.0*pow(pi, 2.5) / ((alpha[p]+alpha[q]) * (alpha[r]+alpha[s]) * sqrt(alpha[p]+alpha[q]+alpha[r]+alpha[s]));
}

int main(){
    double **h, **S, **Q[n_basis][n_basis], *c, **F;
    h = alloc_dmatrix(n_basis, n_basis);
    S = alloc_dmatrix(n_basis, n_basis);
    for(int i=0;i<n_basis;i++){
        for(int j=0;j<n_basis;j++){
            Q[i][j] = alloc_dmatrix(4, 4);
        }
    }
    c = alloc_dvector(n_basis);
    F = alloc_dmatrix(n_basis, n_basis);

    double *w, sc;
    w = alloc_dvector(n_basis);

    for(int i=0;i<n_basis;i++) c[i] = 1.0;

    for(int p=0;p<n_basis;p++){
        for(int q=0;q<n_basis;q++){
            mat_elem(h, p, q) = get_h(p, q);
            mat_elem(S, p, q) = get_S(p, q);
            for(int r=0;r<n_basis;r++)
                for(int s=0;s<n_basis;s++)
                    mat_elem(Q[p][q], r, s) = get_Q(p, r, q, s);
        }
    }

    dgemv(n_basis, S, c, w);
    sc = ddot(n_basis, c, w);
    for(int i=0;i<n_basis;i++) c[i] /= sqrt(sc);

    double E = 0.0, prev = 10.0;
    while(fabs(E-prev) >= 1e-20){
        for(int p=0;p<n_basis;p++){
            for(int q=0;q<n_basis;q++){
                dgemv(n_basis, Q[p][q], c, w);
                mat_elem(F, p, q) = mat_elem(h, p, q) + ddot(n_basis, c, w);
            }
        }

        dsygv(n_basis, F, S, w);

        prev = E;
        dgemv(n_basis, h, c, w);
        E = 2.0 * ddot(n_basis, c, w);
        for(int p=0;p<n_basis;p++){
            for(int q=0;q<n_basis;q++){
                dgemv(n_basis, Q[p][q], c, w);
                E += ddot(n_basis, c, w) * c[p] * c[q];
            }
        }

        for(int p=0;p<n_basis;p++){
            c[p] = mat_elem(F, p, 0);
            for(int q=0;q<n_basis;q++){
                mat_elem(S, p, q) = get_S(p, q);
            }
        }
    }

    printf("%.12f\n", E);

    FILE *fp = fopen("helium_wf.dat", "w");
    for(double r=0.0;r<=5.0;r+=0.01){
        double wf = 0.0;
        for(int i=0;i<n_basis;i++) wf += c[i]*exp(-alpha[i]*r*r);
        fprintf(fp, "%.12e\t%.12e\n", r, wf);
    }
    fclose(fp);

    free_dmatrix(h);
    free_dmatrix(S);
    for(int i=0;i<n_basis;i++){
        for(int j=0;j<n_basis;j++){
            free_dmatrix(Q[i][j]);
        }
    }
    free_dvector(c);
    free_dmatrix(F);
    free_dvector(w);
    return 0;
}