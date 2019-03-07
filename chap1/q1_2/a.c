#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include "../../headers/mt19937ar.h"

#define L 150
#define N 9000

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

    for(int n=0;n<N;n++){
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
                    goto next_particle;
                }
            }

            // 進む方向
            int d;
            for(d=genrand_int31()%4; x+dx[d]<0 || x+dx[d]>=L || y+dy[d]<0 || y+dy[d]>=L; d=genrand_int31()%4) ;
            x += dx[d];
            y += dy[d];
        }
        next_particle: ;
    }
    
    for(int y=0;y<L;y++){
        for(int x=0;x<L;x++){
            printf("%d%c", site[x][y], " \n"[x==L-1]);
        }
    }

    return 0;
}