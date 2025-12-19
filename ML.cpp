/*
g++ -O2 -std=c++11 -fopenmp ML.cpp -o ML
*/
#pragma warning(disable:4710)
#pragma warning(disable:4711)
#pragma warning(disable:4820)
#include <vector>
#include <cfloat>
#include <cstdio>
#include <cstring>
#include <climits>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <string>
#include <iostream>
#include <cstdint>
#include <algorithm>
#include <cassert>
#include <random>
#include <queue>
#include <deque>
#include <list>
#include <map>
#include <array>
#include <chrono>
#include <fstream>
#include <functional>
#include <unordered_map>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;

// === 定数・マクロ定義 ===
#define DLT(ST,ED) ((double)((ED)-(ST))/CLOCKS_PER_SEC)
#define XX(PT)  ((PT)&15)
#define YY(PT)  XX((PT)>>4)
#define YX(Y,X) ((Y)<<4|(X))
#define DIR 4
#define ROW 5
#define COL 6
#define DROP 8
#define TRN  150
#define MAX_TURN 150
#define BEAM_WIDTH 10000
#define PROBLEM 10000
#define BONUS 10
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define NODE_SIZE MAX(500,4*BEAM_WIDTH)
#define TRAIN 100

// === ★追加・変更: 相対座標計算用マクロ ===
// ROW=5, COL=6 の場合、dyは-4~4, dxは-5~5 の範囲を取り得る
// これらを正のインデックスに変換するための定数
#define D_OFFSET_Y (ROW)       // Y方向オフセット
#define D_OFFSET_X (COL)       // X方向オフセット
#define D_stride_X (2*COL + 1) // 横幅のストライド（安全マージン込み）
// 配列サイズ： (縦幅 * 横幅) 以上の十分なサイズを確保
// (2*ROW+1) * (2*COL+1) = 11 * 13 = 143 なので 200あれば十分
#define REL_SIZE 256

typedef char F_T;
typedef char T_T;
typedef signed char sc;
typedef unsigned char uc;
typedef unsigned long long ll;
enum { EVAL_NONE = 0, EVAL_FALL, EVAL_SET, EVAL_FS, EVAL_COMBO };

// --- 関数プロトタイプ ---
void init(F_T field[ROW][COL]);
void fall(int x,int h,F_T field[ROW][COL]);
void set(F_T field[ROW][COL], int force);
void show_field(F_T field[ROW][COL]);
int rnd(int mini, int maxi);
int chain(int nrw, int ncl, F_T d, F_T field[ROW][COL], F_T chkflag[ROW][COL], F_T delflag[ROW][COL]);
int evaluate(F_T field[ROW][COL], int flag);
int sum_e(F_T field[ROW][COL]);
int sum_evaluate(F_T field[ROW][COL]);
void operation(F_T field[ROW][COL], T_T first_te,ll route[(TRN/21)+1]);
int evaluate2(F_T field[ROW][COL], int flag, sc* combo, ll* hash,int p_maxcombo[DROP+1]);
int sum_e2(F_T field[ROW][COL], sc* combo, ll* hash,int p_maxcombo[DROP+1]);
ll xor128();
int sum_e3(F_T field[ROW][COL], sc* combo, int p_maxcombo[DROP+1]);
int evaluate3(F_T field[ROW][COL], int flag, sc* combo, int p_maxcombo[DROP+1]);

ll zoblish_field[ROW][COL][DROP+1];

// === ★変更: 配列サイズを拡張 ===
// [色][相対位置1][相対位置2]
int data[15][REL_SIZE][REL_SIZE];

// --- 構造体定義 ---
struct node {
    T_T first_te;
    ll movei[(TRN/21)+1];
    int score;
    sc combo;
    sc nowC;
    sc nowR;
    sc prev;
    int prev_score;
    uc improving;
    ll hash;
    node() {
        this->score = 0;
        this->prev = -1;
    }
    bool operator < (const node& n)const {
        return score < n.score;
    }
}fff[NODE_SIZE];

struct Action {
    T_T first_te;
    int score;
    int maxcombo;
    ll moving[(TRN/21)+1];
    Action() {
        this->score = 0;
    }
};

Action BEAM_SEARCH(F_T f_field[ROW][COL]);
double part1 = 0, part2 = 0, part3 = 0, part4 = 0, MAXCOMBO = 0;

// --- ビームサーチ実装 ---
Action BEAM_SEARCH(F_T f_field[ROW][COL]) {
    int stop = 0;
    int p_maxcombo[DROP+1] = {0};
    int drop[DROP + 1] = { 0 };
    for (int row = 0; row < ROW; row++) {
        for (int col = 0; col < COL; col++) {
            if (1 <= f_field[row][col] && f_field[row][col] <= DROP) {
                drop[f_field[row][col]]++;
            }
        }
    }
    for (int i = 1; i <= DROP; i++) {
        stop += drop[i] / 3;
        p_maxcombo[i]=drop[i]/3;
    }
    MAXCOMBO += (double)stop;

    vector<node>dque;
    double start, st;
    dque.clear();
    for (int i = 0; i < ROW; i++) {
        for (int j = 0; j < COL; j++) {
            node cand;
            cand.nowR = i;
            cand.nowC = j;
            cand.prev = -1;
            cand.first_te = (T_T)YX(i, j);
            for (int trn = 0; trn <= TRN/21; trn++) {
                cand.movei[trn] = 0ll;
            }
            F_T ff_field[ROW][COL];
            memcpy(ff_field,f_field,sizeof(ff_field));
            sc cmb;
            ll ha;
            cand.prev_score=sum_e2(ff_field,&cmb,&ha,p_maxcombo);
            cand.improving=0;
            cand.hash=ha;
            dque.push_back(cand);
        }
    }
    int dx[DIR] = { -1, 0,0,1 }, dy[DIR] = { 0,-1,1,0 };
    Action bestAction;
    int maxValue = 0;
    bestAction.maxcombo = stop;
    unordered_map<ll, bool> checkNodeList[ROW*COL];

    for (int i = 0; i < MAX_TURN; i++) {
        int ks = (int)dque.size();
        start = omp_get_wtime();
#pragma omp parallel for private(st),reduction(+:part1,part4)
        for (int k = 0; k < ks; k++) {
            node temp = dque[k];
            F_T temp_field[ROW][COL];
            memcpy(temp_field, f_field, sizeof(temp_field));
            operation(temp_field, temp.first_te,temp.movei);
            for (int j = 0; j < DIR; j++) {
                node cand = temp;
                if (0 <= cand.nowC + dx[j] && cand.nowC + dx[j] < COL &&
                    0 <= cand.nowR + dy[j] && cand.nowR + dy[j] < ROW) {
                    if (cand.prev + j != 3) {
                        int ny=cand.nowR + dy[j];
                        int nx=cand.nowC + dx[j];
                        F_T field[ROW][COL];
                        memcpy(field,temp_field,sizeof(temp_field));
                        F_T tmp=field[cand.nowR][cand.nowC];
                        cand.hash^=(zoblish_field[cand.nowR][cand.nowC][tmp])^(zoblish_field[ny][nx][field[ny][nx]]);
                        cand.hash^=(zoblish_field[cand.nowR][cand.nowC][field[ny][nx]])^(zoblish_field[ny][nx][tmp]);
                        field[cand.nowR][cand.nowC]=field[ny][nx];
                        field[ny][nx]=tmp;
                        cand.nowC += dx[j];
                        cand.nowR += dy[j];
                        cand.movei[i/21] |= (((ll)(j+1))<<((3*i)%63));
                        st = omp_get_wtime();
                        sc cmb;
                        cand.score = evaluate3(field, EVAL_FALL | EVAL_COMBO, &cmb,p_maxcombo);
                        cand.combo = cmb;
                        part1 += omp_get_wtime() - st;
                        cand.prev = j;
                        st = omp_get_wtime();
                        fff[(4 * k) + j] = cand;
                        part4 += omp_get_wtime() - st;
                    } else {
                        cand.combo = -1;
                        fff[(4 * k) + j] = cand;
                    }
                } else {
                    cand.combo = -1;
                    fff[(4 * k) + j] = cand;
                }
            }
        }
        part2 += omp_get_wtime() - start;
        start = omp_get_wtime();
        dque.clear();
        vector<pair<int,int> >vec;
        int ks2 = 0;
        for (int j = 0; j < 4 * ks; j++) {
            if (fff[j].combo != -1) {
                if (fff[j].combo==stop) {
                    maxValue = stop;
                    bestAction.score = maxValue;
                    bestAction.first_te=fff[j].first_te;
                    memcpy(bestAction.moving, fff[j].movei, sizeof(fff[j].movei));
                    return bestAction;
                }
                if(fff[j].score>fff[j].prev_score){fff[j].improving=fff[j].improving+1;}
                fff[j].prev_score=fff[j].score;
                int sc=fff[j].score+(BONUS*fff[j].improving)+(fff[j].nowR*3);
                vec.push_back(make_pair(-sc,j));    
                ks2++;
            }
        }
        sort(vec.begin(),vec.end());
        int push_node=0;
        for (int j = 0; push_node < BEAM_WIDTH; j++) {
            if((int)vec.size()<=j){break;}
            int v=vec[j].second;
            node temp = fff[v];
            if (maxValue < temp.combo) {
                maxValue = temp.combo;
                bestAction.score = maxValue;
                bestAction.first_te=temp.first_te;
                memcpy(bestAction.moving, temp.movei, sizeof(temp.movei));
            }
            if (i < MAX_TURN - 1) {
                int pos=(temp.nowR*COL)+temp.nowC;
                if(!checkNodeList[pos][temp.hash]){
                    checkNodeList[pos][temp.hash]=true;
                    dque.push_back(temp);
                    push_node++;
                }
            }
        }
        part3 += omp_get_wtime() - start;
    }
    return bestAction;
}

void show_field(F_T field[ROW][COL]) {
    for (int i = 0; i < ROW; i++) {
        for (int j = 0; j < COL; j++) {
            printf("%d", field[i][j]);
        }
        printf("\n");
    }
}
void fall(int x,int h,F_T field[ROW][COL]) {
    int tgt;
    for (tgt = ROW - 1; tgt >= h && field[tgt][x] != 0; tgt--);
    for (int i = tgt - 1; i >= h; i--) {
        if (field[i][x] != 0) {
            F_T c = field[i][x];
            field[i][x] = 0;
            field[tgt][x] = c;
            tgt--;
        }
    }
}
void init(F_T field[ROW][COL]) { set(field, !0); }
void set(F_T field[ROW][COL], int force) {
    for (int i = 0; i < ROW; i++) {
        for (int j = 0; j < COL; j++) {
            if (field[i][j] == 0 || force) {
                field[i][j] = (F_T)rnd(force ? 0 : 1, DROP);
            }
        }
    }
}
int chain(int nrw, int ncl, F_T d, F_T field[ROW][COL],
    F_T chkflag[ROW][COL], F_T delflag[ROW][COL]) {
    int count = 0;
#define CHK_CF(Y,X) (field[Y][X] == d && chkflag[Y][X]==0 && delflag[Y][X] > 0)
    if (CHK_CF(nrw, ncl)) {
        ++count;
        chkflag[nrw][ncl]=1;
        if (0 < nrw && CHK_CF(nrw - 1, ncl)) count += chain(nrw - 1, ncl, d, field, chkflag, delflag);
        if (nrw < ROW - 1 && CHK_CF(nrw + 1, ncl)) count += chain(nrw + 1, ncl, d, field, chkflag, delflag);
        if (0 < ncl && CHK_CF(nrw, ncl - 1)) count += chain(nrw, ncl - 1, d, field, chkflag, delflag);
        if (ncl < COL - 1 && CHK_CF(nrw, ncl + 1)) count += chain(nrw, ncl + 1, d, field, chkflag, delflag);
    }
    return count;
}
int evaluate(F_T field[ROW][COL], int flag) {
    int combo = 0;
    while (1) {
        int cmb = 0;
        F_T chkflag[ROW][COL]={0};
        F_T delflag[ROW][COL]={0};
        F_T GetHeight[COL];
        for (int row = 0; row < ROW; row++) {
            for (int col = 0; col < COL; col++) {
                F_T num=field[row][col];
                if(row==0) GetHeight[col]=(F_T)ROW;
                if(num>0 && GetHeight[col]==(F_T)ROW) GetHeight[col]=(F_T)row;
                if (col <= COL - 3 && num == field[row][col + 1] && num == field[row][col + 2] && num > 0) {
                    delflag[row][col]=1; delflag[row][col+1]=1; delflag[row][col+2]=1;
                }
                if (row <= ROW - 3 && num == field[row + 1][col] && num == field[row + 2][col] && num > 0) {
                    delflag[row][col]=1; delflag[row+1][col]=1; delflag[row+2][col]=1;
                }
            }
        }
        for (int row = 0; row < ROW; row++) {
            for (int col = 0; col < COL; col++) {
                if (delflag[row][col] > 0) {
                    if (chain(row, col, field[row][col], field, chkflag, delflag) >= 3) cmb++;
                }
            }
        }
        combo += cmb;
        if (cmb == 0 || 0 == (flag & EVAL_COMBO)) break;
        for (int row = 0; row < ROW; row++) for (int col = 0; col < COL; col++) if (delflag[row][col]> 0) field[row][col] = 0;
        if (flag & EVAL_FALL) for(int x=0;x<COL;x++) fall(x,GetHeight[x],field);
        if (flag & EVAL_SET) set(field, 0);
    }
    return combo;
}
int evaluate2(F_T field[ROW][COL], int flag, sc* combo, ll* hash,int p_maxcombo[DROP+1]) {
    int ev = 0; *combo = 0; ll ha=0; int oti = 0; int d_maxcombo[DROP+1]={0};
    while (1) {
        int cmb = 0, cmb2 = 0;
        F_T chkflag[ROW][COL]={0}, delflag[ROW][COL]={0}, GetHeight[COL];
        int cnt_drop[DROP+1]={0}, right[DROP+1], left[DROP+1];
        for(int i=0;i<=DROP;i++){ right[i]=-1; left[i]=COL; }
        for (int row = 0; row < ROW; row++) {
            for (int col = 0; col < COL; col++) {
                F_T num = field[row][col];
                cnt_drop[(int)num]++;
                if(row==0) GetHeight[col]=(F_T)ROW;
                if(num>0 && GetHeight[col]==(F_T)ROW) GetHeight[col]=(F_T)row;
                if(oti==0) ha ^= zoblish_field[row][col][(int)num];
                if (col <= COL - 3 && num == field[row][col + 1] && num == field[row][col + 2] && num > 0) {
                    delflag[row][col]=1; delflag[row][col+1]=1; delflag[row][col+2]=1;
                }
                if (row <= ROW - 3 && num == field[row + 1][col] && num == field[row + 2][col] && num > 0) {
                    delflag[row][col]=1; delflag[row+1][col]=1; delflag[row+2][col]=1;
                }
            }
        }
        F_T erase_x[COL]={0};
        for (int row = 0; row < ROW; row++) {
            for (int col = 0; col < COL; col++) {
                if (delflag[row][col]>0) {
                    int c = chain(row, col, field[row][col], field, chkflag, delflag);
                    if (c >= 3) {
                        cmb++;
                        if (c == 3) cmb2 += 30; else cmb2 += 20;
                        d_maxcombo[(int)field[row][col]]++;
                    }
                    field[row][col]=0; erase_x[col]=1;
                } else {
                    right[(int)field[row][col]]=max(right[(int)field[row][col]],col);
                    left[(int)field[row][col]]=min(left[(int)field[row][col]],col);
                }
            }
        }
        for(int i=1;i<=DROP;i++){
            if(right[i]!=-1&&left[i]!=COL&&cnt_drop[i]>=3&&p_maxcombo[i]!=d_maxcombo[i]) cmb2-=right[i]-left[i];
        }
        *combo += cmb; ev += cmb2;
        if (cmb == 0 || 0 == (flag & EVAL_COMBO)) break;
        oti++;
        if (flag & EVAL_FALL) for(int x=0;x<COL;x++) if(erase_x[x]==1) fall(x,GetHeight[x],field);
        if (flag & EVAL_SET) set(field, 0);
    }
    ev += oti; *hash=ha; return ev;
}
int evaluate3(F_T field[ROW][COL], int flag, sc* combo, int p_maxcombo[DROP+1]) {
    int ev = 0; *combo = 0; int oti = 0; int d_maxcombo[DROP+1]={0};
    while (1) {
        int cmb = 0, cmb2 = 0;
        F_T chkflag[ROW][COL]={0}, delflag[ROW][COL]={0}, GetHeight[COL];
        int cnt_drop[DROP+1]={0}, right[DROP+1], left[DROP+1];
        for(int i=0;i<=DROP;i++){ right[i]=-1; left[i]=COL; }
        for (int row = 0; row < ROW; row++) {
            for (int col = 0; col < COL; col++) {
                F_T num = field[row][col];
                cnt_drop[(int)num]++;
                if(row==0) GetHeight[col]=(F_T)ROW;
                if(num>0 && GetHeight[col]==(F_T)ROW) GetHeight[col]=(F_T)row;
                if (col <= COL - 3 && num == field[row][col + 1] && num == field[row][col + 2] && num > 0) {
                    delflag[row][col]=1; delflag[row][col+1]=1; delflag[row][col+2]=1;
                }
                if (row <= ROW - 3 && num == field[row + 1][col] && num == field[row + 2][col] && num > 0) {
                    delflag[row][col]=1; delflag[row+1][col]=1; delflag[row+2][col]=1;
                }
            }
        }
        F_T erase_x[COL]={0};
        for (int row = 0; row < ROW; row++) {
            for (int col = 0; col < COL; col++) {
                if (delflag[row][col]>0) {
                    int c = chain(row, col, field[row][col], field, chkflag, delflag);
                    if (c >= 3) {
                        cmb++;
                        if (c == 3) cmb2 += 30; else cmb2 += 20;
                        d_maxcombo[(int)field[row][col]]++;
                    }
                    field[row][col]=0; erase_x[col]=1;
                } else {
                    right[(int)field[row][col]]=max(right[(int)field[row][col]],col);
                    left[(int)field[row][col]]=min(left[(int)field[row][col]],col);
                }
            }
        }
        for(int i=1;i<=DROP;i++){
            if(right[i]!=-1&&left[i]!=COL&&cnt_drop[i]>=3&&p_maxcombo[i]!=d_maxcombo[i]) cmb2-=right[i]-left[i];
        }
        *combo += cmb; ev += cmb2;
        if (cmb == 0 || 0 == (flag & EVAL_COMBO)) break;
        oti++;
        if (flag & EVAL_FALL) for(int x=0;x<COL;x++) if(erase_x[x]==1) fall(x,GetHeight[x],field);
        if (flag & EVAL_SET) set(field, 0);
    }
    ev += oti; return ev;
}
int sum_e3(F_T field[ROW][COL], sc* combo, int p_maxcombo[DROP+1]) { return evaluate3(field, EVAL_FALL | EVAL_COMBO, combo,p_maxcombo); }
int sum_e2(F_T field[ROW][COL], sc* combo, ll* hash,int p_maxcombo[DROP+1]) { return evaluate2(field, EVAL_FALL | EVAL_COMBO, combo,hash,p_maxcombo); }
int sum_e(F_T field[ROW][COL]) { return evaluate(field, EVAL_FALL | EVAL_COMBO); }
int sum_evaluate(F_T field[ROW][COL]) { return evaluate(field, EVAL_FS | EVAL_COMBO); }

void operation(F_T field[ROW][COL], T_T first_te,ll route[(TRN/21)+1]) {
    int prw = (int)YY(first_te), pcl = (int)XX(first_te), i,j;
    int dx[DIR] = { -1, 0,0,1 };
    int dy[DIR] = { 0,-1,1,0 };
    for (i = 0; i <= TRN/21; i++) {
        if (route[i] == 0ll) break;
        for(j=0;j<21;j++){
            int dir = (int)(7ll&(route[i]>>(3*j)));
            if(dir==0) break;
            int row=prw+dy[dir-1], col=pcl+dx[dir-1];
            F_T c = field[prw][pcl]; field[prw][pcl] = field[row][col]; field[row][col] = c;
            prw = row, pcl = col;
        }
    }
}
int rnd(int mini, int maxi) {
    static mt19937 mt((int)time(0));
    uniform_int_distribution<int> dice(mini, maxi);
    return dice(mt);
}
ll xor128() {
    static unsigned long long rx = 123456789, ry = 362436069, rz = 521288629, rw = 88675123;
    ll rt = (rx ^ (rx << 11)); rx = ry; ry = rz; rz = rw;
    return (rw = (rw ^ (rw >> 19)) ^ (rt ^ (rt >> 8)));
}

// === ★変更: memo関数 (相対座標ID化) ===
void memo(F_T field[ROW][COL]){
    vector<int>v[10];
    for(int i=0;i<ROW*COL;i++){
        int a = (int)(field[i/COL][i%COL]);
        v[a].push_back(i);
    }

    for(int i=0;i<10;i++){
        for(int j=0;j<(int)v[i].size();j++){
            if((int)v[i].size()<=j+2){break;}
            int p1 = v[i][j];
            int p2 = v[i][j+1];
            int p3 = v[i][j+2];
           
            // --- ここを変更 ---
            int r1 = p1 / COL; int c1 = p1 % COL;
            int r2 = p2 / COL; int c2 = p2 % COL;
            int r3 = p3 / COL; int c3 = p3 % COL;

            // 相対ベクトル (dy, dx)
            int dy1 = r2 - r1; int dx1 = c2 - c1;
            int dy2 = r3 - r1; int dx2 = c3 - c1;

            // ユニークID生成
            // dyは -ROW~+ROW, dxは -COL~+COL の範囲なのでオフセットを足して正の値にする
            int d1_idx = (dy1 + D_OFFSET_Y) * D_stride_X + (dx1 + D_OFFSET_X);
            int d2_idx = (dy2 + D_OFFSET_Y) * D_stride_X + (dx2 + D_OFFSET_X);

            // 安全チェック (範囲外アクセス防止)
            if (d1_idx >= 0 && d1_idx < REL_SIZE && d2_idx >= 0 && d2_idx < REL_SIZE) {
                data[i][d1_idx][d2_idx]++;
            }
        }
    }
}

void counting(F_T field[ROW][COL],string route){
    F_T f_field[ROW][COL];
    for(int i=0;i<ROW*COL;i++){f_field[i/COL][i%COL]=field[i/COL][i%COL];}
    int tgt=0;
    string top="";
    while(1){
        if(route[tgt]==','){tgt++;break;}
        top+=route[tgt];
        tgt++;
    }
    int pos;
    if((int)top.size()==2){int x=top[0]-'0';int y=(top[1]-'0')-5;pos=(y*COL)+x;}
    else{int x=top[0]-'0';int y=5;pos=(y*COL)+x;}
    //int tesuu=(int)route.size()-tgt;
    //int cnt=0;
    for(int j=tgt;j<(int)route.size();j++){
        memo(f_field);
        //cnt++;
        if(route[j]=='3'){swap(f_field[pos/COL][pos%COL],f_field[pos/COL][(pos%COL)-1]);pos--;}
        if(route[j]=='6'){swap(f_field[pos/COL][pos%COL],f_field[(pos/COL)-1][pos%COL]);pos-=COL;}
        if(route[j]=='1'){swap(f_field[pos/COL][pos%COL],f_field[(pos/COL)+1][pos%COL]);pos+=COL;}
        if(route[j]=='4'){swap(f_field[pos/COL][pos%COL],f_field[pos/COL][(pos%COL)+1]);pos++;}
    }
    memo(f_field);
}
double logN(double b, double x) {
    return log(x) / log(b);
}

int main() {
    int i, j, k;

    for(i=0;i<ROW;++i) for(j=0;j<COL;++j) for(k=0;k<=DROP;k++) zoblish_field[i][j][k]=xor128();

    int mistake=0;
    double avg = 0;
    double start;
    double t_sum = 0;
    double oti_avg = 0;
    int acc=0;
   
    bool start_test=false; // data.txtから読み込む場合
    if(start_test){
        ifstream myf ("data.txt");
        string ls;
        while(getline(myf,ls)){
            string parent="", child="";
            bool slash=false;
            for(i=0;i<(int)ls.size();i++){
                if(ls[i]=='/'){slash=true;continue;}
                if(slash){child+=ls[i];} else{parent+=ls[i];}
            }
            int counter=0;
            string xz[3]={"","",""};
            for(i=0;i<(int)parent.size();i++){
                if(parent[i]==','){counter++;continue;}
                xz[counter]+=parent[i];
            }
            // 配列サイズチェックを入れて安全に読み込む
            int d1 = stoi(xz[1]);
            int d2 = stoi(xz[2]);
            if(d1 >= 0 && d1 < REL_SIZE && d2 >= 0 && d2 < REL_SIZE) {
                data[stoi(xz[0])][d1][d2]=stoi(child);
            }
        }
        myf.close();
    }
   
    for (i = 0; i < PROBLEM; i++) {
        // --- ★変更: 保存ループを拡張されたサイズに対応 ---
        if(i%100==0){
            string mystr="";    
            for (int a1=0;a1<10;a1++){
                for(int a2=0;a2<REL_SIZE;a2++){ // ROW*COL (30) から REL_SIZE (256) に変更
                    for(int a3=0;a3<REL_SIZE;a3++){ // 同上
                        int value=data[a1][a2][a3];    
                        string ms=to_string(a1)+","+to_string(a2)+","+to_string(a3)+"/"+to_string(value);
                        if(value>0){
                            mystr+=ms+"\n";
                        }    
                    }
                }
            }
            ofstream fi("data.txt");
            fi<<mystr;
            fi.close();
        }
       
        F_T f_field[ROW][COL], field[ROW][COL], oti_field[ROW][COL];
        printf("input:No.%d/%d\n", i + 1, PROBLEM);
        init(f_field); set(f_field, 0);
        show_field(f_field);
       
        start = omp_get_wtime();
        Action tmp = BEAM_SEARCH(f_field);
        double diff = omp_get_wtime() - start;
        t_sum += diff;
       
        string layout="";
        for(int v=0;v<ROW;v++) for(int u=0;u<COL;u++) layout+=to_string(f_field[v][u]-1);
       
        string route="";
        int path_length=0;
        route+=to_string(XX(tmp.first_te))+to_string(YY(tmp.first_te)+5)+",";
        for (j = 0; j <= TRN/21; j++) {
            if (tmp.moving[j] == 0ll) break;
            for(k=0;k<21;k++){
                int dir = (int)(7ll&(tmp.moving[j]>>(3*k)));
                if (dir==0) break;
                if (dir==1) route+=to_string(3);
                if (dir==2) route+=to_string(6);
                if (dir==3) route+=to_string(1);
                if (dir==4) route+=to_string(4);
                path_length++;
            }
        }
        string url="http://serizawa.web5.jp/puzzdra_theory_maker/index.html?layout="+layout+"&route="+route+"&ctwMode=false";
        cout<<url<<endl;
        printf("\n");
        memcpy(field, f_field, sizeof(f_field));
        counting(field,route);        
        operation(field, tmp.first_te,tmp.moving);
        printf("output:No.%d/%d\n", i + 1, PROBLEM);
        show_field(field);
        memcpy(oti_field, field, sizeof(field));
        int combo = sum_e(field);
        int oti = sum_evaluate(oti_field);
        if(combo!=tmp.maxcombo){mistake++;} else{acc++;}
        printf("mistake=%d\n",mistake);
        printf("acc=%d\n",acc);
        printf("path_length=%d\n",path_length);
        printf("Normal:%d/%dCombo\n", combo, tmp.maxcombo);
        printf("Oti:%dCombo\n", oti);
        printf("Duration:%fSec\n", diff);
        printf("------------\n");
        avg += (double)combo;
        oti_avg += (double)oti;
    }
    printf("TotalDuration:%fSec\n", t_sum);
    printf("Avg.NormalCombo #:%f/%f\n", avg / (double)i, MAXCOMBO / (double)i);
    printf("Avg.OtiCombo #:%f\n", oti_avg / (double)i);
    printf("p1:%f,p2:%f,p3:%f,p4:%f\n", part1, part2, part3, part4);
    j = getchar();
    return 0;
}
