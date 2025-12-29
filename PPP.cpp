/*

Linux:

g++ -O2 -std=c++11 -fopenmp PPP.cpp loguru.cpp -o PPP -lpthread -ldl

Windows11:

g++ -O2 -std=c++11 -fopenmp -lpthread PPP.cpp loguru.cpp -o PPP

./PPP
*/
#pragma warning(disable:4710)
#pragma warning(disable:4711)
#pragma warning(disable:4820)
#include <tuple>
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
#include "hash_map.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;

// === パラメータ設定 ===
#define ROW 5
#define COL 6
#define DROP 8
#define TRN 150
#define MAX_TURN 150
#define BEAM_WIDTH 10000
#define PROBLEM 100
#define BONUS 10

#define NODE_SIZE MAX(500,4*BEAM_WIDTH)
#define DIR 4
#define DLT(ST,ED) ((double)((ED)-(ST))/CLOCKS_PER_SEC)
#define XX(PT)  ((PT)&15)
#define YY(PT)  XX((PT)>>4)
#define YX(Y,X) ((Y)<<4|(X))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

typedef char F_T;
typedef char T_T;
typedef signed char sc;
typedef unsigned char uc;
typedef unsigned long long ll;

enum { EVAL_NONE = 0, EVAL_FALL, EVAL_SET, EVAL_FS, EVAL_COMBO };

// --- 関数プロトタイプ宣言 ---
void init(F_T field[ROW][COL]);
void fall(int x,int h,F_T field[ROW][COL]);
void set(F_T field[ROW][COL], int force);
void show_field(F_T field[ROW][COL]);
int rnd(int mini, int maxi);
int chain(int nrw, int ncl, F_T d, F_T field[ROW][COL], F_T chkflag[ROW][COL], F_T delflag[ROW][COL]);
int sum_e(F_T field[ROW][COL]);
int sum_evaluate(F_T field[ROW][COL]);
int sum_e2(F_T field[ROW][COL], sc* combo, ll* hash,int p_maxcombo[DROP+1]);
void operation(F_T field[ROW][COL], T_T first_te,ll route[(TRN/21)+1]);
ll xor128();
int evaluate(F_T field[ROW][COL], int flag);

// --- グローバル変数 ---
ll zoblish_field[ROW][COL][DROP+1];

int data[15][ROW*COL][ROW*COL][ROW*COL]={0};

// 探索ノード
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
} fff[NODE_SIZE];

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

int NNUE_init_score(F_T board[ROW][COL]) {
    vector<int>v[10];
    for(int i=0;i<ROW*COL;i++){
        int a = (int)(board[i/COL][i%COL]);
        v[a].push_back(i);
    }

    int score=0;

    for(int i=0;i<10;i++){
        for(int j=0;j<(int)v[i].size();j+=3){
                if((int)v[i].size()<=j+2){break;}
                int p1 = v[i][j];
                int p2 = v[i][j+1];
                int p3 = v[i][j+2];
                score += data[i][p1][p2][p3];
        }
    }
   
    return score;
}

int NNUE_score(F_T board[ROW][COL],int c1) {
    int v[ROW*COL]={0};
    int cnty=0;
    for(int i=0;i<ROW*COL;i++){
        int a = (int)(board[i/COL][i%COL]);
        if(a==c1){
        v[cnty]=i;
        cnty++;
        }
    }

    int score=0;

    for(int j=0;j<cnty;j+=3){
    if(cnty<=j+2){break;}
    int p1 = v[j];
    int p2 = v[j+1];
    int p3 = v[j+2];
    score += data[c1][p1][p2][p3];
    }
   
    return score;
}

// --- ビームサーチ ---
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

    vector<node> dque;
    double st;
    dque.clear();
   
    // 1手目展開
    for (int i = 0; i < ROW; i++) {
        for (int j = 0; j < COL; j++) {
            node cand;
            cand.nowR = i;
            cand.nowC = j;
            cand.prev = -1;
            cand.first_te = (T_T)YX(i, j);
            for (int trn = 0; trn <= TRN/21; trn++) cand.movei[trn] = 0ll;
           
            F_T ff_field[ROW][COL];
            memcpy(ff_field,f_field,sizeof(ff_field));
            sc cmb;
            ll ha;
            //cand.prev_score=sum_e2(ff_field,&cmb,&ha,p_maxcombo);
            cand.score=NNUE_init_score(ff_field);
            cand.improving=0;
            cand.hash=ha;
            dque.push_back(cand);
        }
    }

    int dx[DIR] = { -1, 0,0,1 };
    int dy[DIR] = { 0,-1,1,0 };
    Action bestAction;
    int maxValue = 0;
    bestAction.maxcombo = stop;

    emilib::HashMap<ll, bool> checkNodeList[ROW*COL];

    // ビーム探索ループ
    for (int i = 0; i < MAX_TURN; i++) {
        int ks = (int)dque.size();
       
        #pragma omp parallel for private(st)
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
                        int c1=(int)field[cand.nowR][cand.nowC];
                        int c2=(int)field[ny][nx];
                        int OLD=NNUE_score(field,c1)+NNUE_score(field,c2);
                        cand.hash^=(zoblish_field[cand.nowR][cand.nowC][tmp])^(zoblish_field[ny][nx][field[ny][nx]]);
                        cand.hash^=(zoblish_field[cand.nowR][cand.nowC][field[ny][nx]])^(zoblish_field[ny][nx][tmp]);
                        field[cand.nowR][cand.nowC]=field[ny][nx];
                        field[ny][nx]=tmp;
                        cand.nowC += dx[j];
                        cand.nowR += dy[j];
                        cand.movei[i/21] |= (((ll)(j+1))<<((3*i)%63));
                        //cand.score += NNUE_score(field,c1,c2); // NNUE評価
                        int NEW=NNUE_score(field,c1)+NNUE_score(field,c2);
                        //cand.score = NNUE_init_score(field);
                        cand.score-=OLD;
                        cand.score+=NEW;
                        cand.combo = (sc)evaluate(field, EVAL_FALL | EVAL_COMBO);
                        cand.prev = j;
                        fff[(4 * k) + j] = cand;
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

        dque.clear();
        vector<tuple<int, int, int> >vec;
        int ks2 = 0;
        for (int j = 0; j < 4 * ks; j++) {
            if (fff[j].combo != -1) {
                if(fff[j].score>fff[j].prev_score){fff[j].improving=fff[j].improving+1;}
                fff[j].prev_score=fff[j].score;
                int sc=fff[j].score+(BONUS*fff[j].improving)+(fff[j].nowR*3);
				int cb=(int)(-fff[j].combo);
                vec.emplace_back(0, -sc, j);    
                ks2++;
            }
        }
        sort(vec.begin(),vec.end());
        int push_node=0;
        for (int j = 0; push_node < BEAM_WIDTH; j++) {
            if((int)vec.size()<=j){break;}
            int v=get<2>(vec[j]);
            node temp = fff[v];
           
            // 理論値コンボチェックは省略(NNUE学習中は評価値のみで走る)
           
            if (i < MAX_TURN - 1) {
                int pos=(temp.nowR*COL)+temp.nowC;
                if(!checkNodeList[pos][temp.hash]){
                    checkNodeList[pos][temp.hash]=true;
                    dque.push_back(temp);
                    push_node++;
                }
            }
            // 暫定ベスト更新
            if(maxValue<(int)temp.combo) {
                maxValue=(int)temp.combo;
                bestAction.first_te = temp.first_te;
                memcpy(bestAction.moving, temp.movei, sizeof(temp.movei));
                if(maxValue>=stop){return bestAction;}
            }
        }
        if(push_node==0) break;
    }
    return bestAction;
}

// --- その他の関数実装 ---
void show_field(F_T field[ROW][COL]) {
    for (int i = 0; i < ROW; i++) {
        for (int j = 0; j < COL; j++) printf("%d", field[i][j]);
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
            if (field[i][j] == 0 || force) field[i][j] = (F_T)rnd(force ? 0 : 1, DROP);
        }
    }
}
int chain(int nrw, int ncl, F_T d, F_T field[ROW][COL], F_T chkflag[ROW][COL], F_T delflag[ROW][COL]) {
    int count = 0;
    #define CHK_CF(Y,X) (field[Y][X] == d && chkflag[Y][X]==0 && delflag[Y][X] > 0)
    if (CHK_CF(nrw, ncl)) {
        ++count; chkflag[nrw][ncl]=1;
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
        F_T chkflag[ROW][COL]={0}, delflag[ROW][COL]={0}, GetHeight[COL];
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
int sum_e2(F_T field[ROW][COL], sc* combo, ll* hash,int p_maxcombo[DROP+1]) { return evaluate2(field, EVAL_FALL | EVAL_COMBO, combo,hash,p_maxcombo); }
int sum_e(F_T field[ROW][COL]) { return evaluate(field, EVAL_FALL | EVAL_COMBO); }
int sum_evaluate(F_T field[ROW][COL]) { return evaluate(field, EVAL_FS | EVAL_COMBO); }
void operation(F_T field[ROW][COL], T_T first_te,ll route[(TRN/21)+1]) {
    int prw = (int)YY(first_te), pcl = (int)XX(first_te), i,j;
    int dx[DIR] = { -1, 0,0,1 }, dy[DIR] = { 0,-1,1,0 };
    for (i = 0; i <= TRN/21; i++) {
        if (route[i] == 0ll) break;
        for(j=0;j<21;j++){
            int dir = (int)(7ll&(route[i]>>(3*j)));
            if(dir==0) break;
            int row=prw+dy[dir-1], col=pcl+dx[dir-1];
            F_T c = field[prw][pcl]; field[prw][pcl] = field[row][col]; field[row][col] = c;
            prw = row; pcl = col;
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
            data[i][p1][p2][p3]++;
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
int main() {
    for(int i1=0;i1<ROW;++i1){
    for(int i2=0;i2<COL;++i2){ 
    for(int i3=0;i3<=DROP;i3++){
    zoblish_field[i1][i2][i3]=xor128();
    }
    }
    }

    int i,j,k;

    bool start_test=true;
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
            string xz[4]={"","","",""};
            for(i=0;i<(int)parent.size();i++){
                if(parent[i]==','){counter++;continue;}
                xz[counter]+=parent[i];
            }
            data[stoi(xz[0])][stoi(xz[1])][stoi(xz[2])][stoi(xz[3])]=stoi(child);
        }
        myf.close();
    }

    int acc=0;
    int mistake=0;
    double avg = 0;//平均コンボ数
    double start;
    double t_sum = 0;
    double oti_avg = 0;//平均落ちコンボ数
    double MAXCOMBOT=0;
    for (i = 0; i < PROBLEM; i++) {//PROBLEM問解く

    F_T f_field[ROW][COL]; //スワイプ前の盤面
    F_T field[ROW][COL]; //盤面
    F_T oti_field[ROW][COL];//落ちコン用盤面
    printf("input:No.%d/%d\n", i + 1, PROBLEM);
    init(f_field); set(f_field, 0);//初期盤面生成
    /*
    string str="";
    cin>>str;
    for (j = 0; j < ROW; j++) {
    for (k = 0; k < COL; k++) {
    f_field[j][k] = (str[k+(COL*j)] - '0')+1;
    }
    }
    */
    show_field(f_field);//盤面表示
    start = omp_get_wtime();
    Action tmp = BEAM_SEARCH(f_field);//ビームサーチしてtmpに最善手を保存
    double diff = omp_get_wtime() - start;
    t_sum += diff;
    string layout="";

    for(int v=0;v<ROW;v++){
    for(int u=0;u<COL;u++){
    layout+=to_string(f_field[v][u]-1);
    }
    }
    string route="";
    //printf("(x,y)=(%d,%d)", XX(tmp.first_te), YY(tmp.first_te));
    int path_length=0;
    route+=to_string(XX(tmp.first_te))+to_string(YY(tmp.first_te)+5)+",";
    for (j = 0; j <= TRN/21; j++) {//y座標は下にいくほど大きくなる
    if (tmp.moving[j] == 0ll) { break; }
    for(k=0;k<21;k++){
    int dir = (int)(7ll&(tmp.moving[j]>>(3*k)));
    if (dir==0){break;}
    if (dir==1) { route+=to_string(3);}//printf("L"); } //"LEFT"); }
    if (dir==2) { route+=to_string(6);}//printf("U"); } //"UP"); }
    if (dir==3) { route+=to_string(1);}//printf("D"); } //"DOWN"); }
    if (dir==4) { route+=to_string(4);}//printf("R"); } //"RIGHT"); }
    path_length++;
    }
    }
    string url="http://serizawa.web5.jp/puzzdra_theory_maker/index.html?layout="+layout+"&route="+route+"&ctwMode=false";
    cout<<url<<endl;
    printf("\n");
    memcpy(field, f_field, sizeof(f_field));
    //counting(field,route);       
    operation(field, tmp.first_te,tmp.moving);
    printf("output:No.%d/%d\n", i + 1, PROBLEM);
    show_field(field);
    memcpy(oti_field, field, sizeof(field));
    int combo = sum_e(field);
    int oti = sum_evaluate(oti_field);
    if(combo!=tmp.maxcombo){mistake++;}
    else{acc++;}
    printf("mistake=%d\n",mistake);
    printf("acc=%d\n",acc);
    printf("path_length=%d\n",path_length);
    printf("Normal:%d/%dCombo\n", combo, tmp.maxcombo);
    printf("Oti:%dCombo\n", oti);
    printf("Duration:%fSec\n", diff);
    printf("------------\n");
    avg += (double)combo;
    oti_avg += (double)oti;
    MAXCOMBOT+=(double)tmp.maxcombo;
    }
    printf("TotalDuration:%fSec\n", t_sum);
    printf("Avg.NormalCombo #:%f/%f\n", avg / (double)i, MAXCOMBOT / (double)i);
    printf("Avg.OtiCombo #:%f\n", oti_avg / (double)i);
    printf("p1:%f,p2:%f,p3:%f,p4:%f\n", part1, part2, part3, part4);
    j = getchar();
	return 0;
}
