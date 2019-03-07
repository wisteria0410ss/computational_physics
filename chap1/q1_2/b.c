#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include "../../headers/mt19937ar.h"

#define L 500
#define N 10000

int main(){
    init_genrand(time(NULL));

    bool site[L][L];
    for(int i=0;i<L;i++){
        for(int j=0;j<L;j++){
            site[i][j] = false;
        }
    }
    site[L/2][L/2] = true;

    int dx[4] = {1, -1, 0, 0}, dy[4] = {0, 0, 1, -1};
    double xsum = 0.0, ysum = 0.0, rsum = 0.0;

    for(int n=1;n<=N;n++){
        int x = 0, y = 0;

        // 境界上の初期位置の決定
        int r = genrand_int31() % ((L-1)*4);
        switch(r/(L-1)){
            case 0:
                x = r % (L-1);
                y = 0;
                break;
            case 1:
                x = (L-1);
                y = r % (L-1);
                break;
            case 2:
                x = (L-1) - r % (L-1);
                y = (L-1);
                break;
            case 3:
                x = 0;
                y = (L-1) - r % (L-1);
                break;
        }

        // ランダムウォーク
        while(true){
            for(int i=0;i<4;i++){
                if(x+dx[i] < 0 || x+dx[i] >= L || y+dy[i] < 0 || y+dy[i] >= L) continue;
                if(site[x+dx[i]][y+dy[i]]){
                    site[x][y] = true;
                    xsum += x;
                    ysum += y;
                    rsum += x*x + y*y;
                    goto next_particle;
                }
            }

            // 進む方向
            int d;
            for(d=genrand_int31()%4; x+dx[d]<0 || x+dx[d]>=L || y+dy[d]<0 || y+dy[d]>=L; d=genrand_int31()%4) ;
            x += dx[d];
            y += dy[d];
        }
        next_particle:
        if(n%200 == 0){
            double x0 = xsum/n, y0 = ysum/n;
            double rg = sqrt(rsum/n - (x0*x0 + y0*y0));
            printf("%f\t%d\n", rg, n);
        }
    }

    return 0;
}