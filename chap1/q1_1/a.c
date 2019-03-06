#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "duffing.h"

#define h       0.01
#define steps   10000

int main(int argc, char **argv){
    if(argc < 3){
        fprintf(stderr, "usage: ./a.out x0 v0");
        return 1;
    }
    double x[steps+1], v[steps+1];
    x[0] = atof(argv[1]);
    v[0] = atof(argv[2]);

    rk4(F, x, v, h, steps);

    printf("# t\tx\tv\n");
    for(int i=0;i<=steps;i++){
        printf("%.8f\t%.8f\t%.8f\n", h*i, x[i], v[i]);
    }

    return 0;
}