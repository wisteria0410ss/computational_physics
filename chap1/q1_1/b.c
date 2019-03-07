#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "duffing.h"

#define mesh    300
#define points  100000
#define steps   (mesh*points)

int main(int argc, char **argv){
    if(argc < 3){
        fprintf(stderr, "usage: ./b.out x0 v0");
        return 1;
    }
    double T = 2*PI/omega;
    double h = T/mesh;
    double *x, *v;
    x = (double *)malloc(sizeof(double)*(steps+1));
    v = (double *)malloc(sizeof(double)*(steps+1));
    x[0] = atof(argv[1]);
    v[0] = atof(argv[2]);

    rk4(F, x, v, h, steps);

    printf("# t\tx\tp\n");
    for(int i=0;i<=steps;i+=mesh){
        printf("%.8f\t%.8f\t%.8f\n", h*i, x[i], v[i]);
    }

    free(x);
    free(v);
    
    return 0;
}