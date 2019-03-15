/* Utility for vector and matrix allocation/deallocation */

/* List of functions

  alloc_dvector : allocate vector of double
  alloc_fvector : allocate vector of float
  alloc_ivector : allocate vector of int

  alloc_dmatrix : allocate matrix of double
  alloc_fmatrix : allocate matrix of float
  alloc_imatrix : allocate matrix of int

  free_dvector : deallocate vector of double
  free_fvector : deallocate vector of float
  free_ivector : deallocate vector of int

  free_dmatrix : deallocate matrix of double
  free_fmatrix : deallocate matrix of float
  free_imatrix : deallocate matrix of int

  fprint_dvector: print out vector of double
  fprint_fvector: print out vector of float
  fprint_ivector: print out vector of int

  fprint_dmatrix: print out matrix of double
  fprint_fmatrix: print out matrix of float

  read_dvector: read vector of double from file
  read_fvector: read vector of float from file
  read_ivector: read vector of int from file

  read_dmatrix: read matrix of double from file
  read_fmatrix: read matrix of float from file
  read_imatrix: read matrix of int from file
*/

#ifndef CMATRIX_H
#define CMATRIX_H

#include <stdlib.h>
#include <stdio.h>

/* useful macros */
#define vec_ptr(vec) &(vec)[0]
#define mat_ptr(mat) &(mat)[0][0]
#define mat_elem(mat, i, j) (mat)[j][i]

/* allocate vector of double */
double *alloc_dvector(int n);

/* allocate vector of float */
float *alloc_fvector(int n);

/* allocate vector of int */
int *alloc_ivector(int n);

/* allocate m x n column-major matrix of double */
double **alloc_dmatrix(int m, int n);

/* allocate m x n column-major matrix of float */
float **alloc_fmatrix(int m, int n);

/* allocate m x n column-major matrix of int */
int **alloc_imatrix(int m, int n);

/* deallocate vector of double */
void free_dvector(double *vec);

/* deallocate vector of float */
void free_fvector(float *vec);

/* deallocate vector of int */
void free_ivector(int *vec);

/* deallocate matrix of double */
void free_dmatrix(double **mat);

/* deallocate float matrix of float */
void free_fmatrix(float **mat);

/* deallocate float matrix of int */
void free_imatrix(int **mat);

/* print out vector of double */
void fprint_dvector(FILE *fp, int n, double *vec);

/* print out vector of float */
void fprint_fvector(FILE *fp, int n, float *vec);

/* print out vector of int */
void fprint_ivector(FILE *fp, int n, int *vec);

/* print out matrix of double */
void fprint_dmatrix(FILE *fp, int m, int n, double **mat);

/* print out matrix of float */
void fprint_fmatrix(FILE *fp, int m, int n, float **mat);

/* print out matrix of int */
void fprint_imatrix(FILE *fp, int m, int n, int **mat);

/* read vector of double from file */
void read_dvector(FILE *fp, int *n, double **vec);

/* read vector of float from file */
void read_fvector(FILE *fp, int *n, float **vec);

/* read matrix of double from file */
void read_dmatrix(FILE *fp, int *m, int *n, double ***mat);

/* read matrix of float from file */
void read_fmatrix(FILE *fp, int *m, int *n, float ***mat);

#endif
