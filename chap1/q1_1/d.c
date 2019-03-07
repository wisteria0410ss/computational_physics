#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "duffing.h"

#define mesh    1000
#define points  25000
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

    bool in[128][128];
    printf("# b\tN\n");
    for(int l=2;l<=128;l*=2){
        int cnt = 0;
        for(int i=0;i<l;i++){
            for(int j=0;j<l;j++){
                in[i][j] = false;
            }
        }
        for(int i=0;i<=steps;i+=mesh){
            int xx = (x[i]+3.0)*l/6.0;
            int yy = (v[i]+3.0)*l/6.0;
            if(!in[xx][yy]){
                cnt++;
                in[xx][yy] = true;
            }
        }
        printf("%f\t%d\n", 6.0/l, cnt);
    }

    free(x);
    free(v);
    
    return 0;
}