/*

Windows10,Linux ONLY

Linux導入手続き

//メモリ容量確認
free -h

//g++インストール

sudo apt install g++

//wgetインストール

sudo apt-get update

sudo apt-get install -y wget

//Laevateinn.cppをダウンロード

wget --no-check-certificate https://raw.githubusercontent.com/koduma/puzzdra_solver/master/Laevateinn.cpp

//hash_map.hpp,loguru.cpp,loguru.hppをダウンロード

wget --no-check-certificate https://raw.githubusercontent.com/koduma/puzzdra_solver/master/folder/hash_map.hpp

wget --no-check-certificate https://raw.githubusercontent.com/koduma/puzzdra_solver/master/folder/loguru.cpp

wget --no-check-certificate https://raw.githubusercontent.com/koduma/puzzdra_solver/master/folder/loguru.hpp

//ビーム幅調整

vi Laevateinn.cpp

//コンパイル

Windows10:g++ -O2 -std=c++11 -fopenmp -mbmi2 -lpthread Laevateinn.cpp loguru.cpp -o Laevateinn -mcmodel=large
Linux:g++ -O2 -std=c++11 -fopenmp -mbmi2 -lpthread -ldl Laevateinn.cpp loguru.cpp -o Laevateinn -mcmodel=large

//run

./Laevateinn

//input

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
#include "hash_map.hpp"
#include <immintrin.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;
#define DLT(ST,ED) ((double)((ED)-(ST))/CLOCKS_PER_SEC)//時間差分
#define XX(PT) ((PT)&15)
#define YY(PT) XX((PT)>>4)
#define YX(Y,X) ((Y)<<4|(X))
#define DIR 4//方向
#define ROW 5//縦//MAX6
#define COL 6//横//MAX7
#define DROP 8//ドロップの種類//MAX9
#define TRN 150//手数//MAX155
#define BEAM_WIDTH 2800000//ビーム幅//MAX200000
#define PROBLEM 1//問題数
#define BONUS 10//評価値改善係数
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define NODE_SIZE MAX(500,4*BEAM_WIDTH)
typedef char F_T;//盤面型
typedef char T_T;//手数型
typedef signed char sc;
typedef unsigned char uc;
typedef unsigned long long ll;
enum { EVAL_NONE = 0, EVAL_FALL, EVAL_SET, EVAL_FS, EVAL_COMBO };
void init(F_T field[ROW][COL]); //初期配置生成関数
void fall(int x,int h,F_T field[ROW][COL]); //ドロップの落下処理関数
void set(F_T field[ROW][COL], int force); //空マスを埋める関数
void show_field(F_T field[ROW][COL]); //盤面表示関数
int rnd(int mini, int maxi); //整数乱数
//上下左右に連結しているドロップを再帰的に探索していく関数
int chain(int nrw, int ncl, F_T d, F_T field[ROW][COL], F_T chkflag[ROW][COL], F_T delflag[ROW][COL]);
int evaluate(F_T field[ROW][COL], int flag); //コンボ数判定関数
int sum_e(F_T field[ROW][COL]);//落とし有り、落ちコン無しコンボ数判定関数
int sum_evaluate(F_T field[ROW][COL]);//落としも落ちコンも有りコンボ数判定関数
void operation(F_T field[ROW][COL], T_T first_te,ll route[(TRN/21)+1],ll dropBB[DROP+1]); //スワイプ処理関数

int evaluate2(F_T field[ROW][COL], int flag, sc* combo, ll* hash,int p_maxcombo[DROP+1]);//落とし減点評価関数
int sum_e2(F_T field[ROW][COL], sc* combo, ll* hash,int p_maxcombo[DROP+1]);//評価関数

ll xor128();//xorshift整数乱数
ll zoblish_field[ROW][COL][DROP+1];

ll sqBB[64];
int evaluate3(ll dropBB[DROP+1], int flag, sc* combo, ll* hash,int p_maxcombo[DROP+1]);//落とし減点評価関数
int sum_e3(ll dropBB[DROP+1], sc* combo, ll* hash,int p_maxcombo[DROP+1]);//評価関数
ll around(ll bitboard);
int table[64];
ll fill_64[64];
ll file_bb[COL];
ll calc_mask(ll bitboard);
ll fallBB(ll p,ll rest,ll mask);

int MSB64bit(ll v) {
   if(v == 0ll){return 0;}
   int out =63-__builtin_clzll(v);
   return out;
}

struct node {//どういう手かの構造体
	T_T first_te;
	ll movei[(TRN/21)+1];//スワイプ移動座標
	int score;//評価値
	sc combo;//コンボ数
	sc nowC;//今どのx座標にいるか
	sc nowR;//今どのy座標にいるか
	sc prev;//1手前は上下左右のどっちを選んだか
	int prev_score;//1手前の評価値
	uc improving;//評価値改善回数
	ll hash;//盤面のハッシュ値
	node() {//初期化
		this->score = 0;
		this->prev = -1;
		//memset(this->movei, STP, sizeof(this->movei));
	}
	bool operator < (const node& n)const {//スコアが高い方が優先される
		return score < n.score;
	}
}fff[NODE_SIZE];
struct Action {//最終的に探索された手
	T_T first_te;
	int score;//コンボ数
	int maxcombo;//理論コンボ数
	ll moving[(TRN/21)+1];//スワイプ移動座標
	Action() {//初期化
		this->score = 0;
		//memset(this->moving, STP, sizeof(this->moving));
	}
};
Action BEAM_SEARCH(F_T f_field[ROW][COL],int maxi,int MAX_TRN,int prev_dir,int now_pos); //ルート探索関数
double part1 = 0, part2 = 0, part3 = 0, MAXCOMBO = 0;
Action BEAM_SEARCH(F_T f_field[ROW][COL],int maxi,int MAX_TRN,int prev_dir,int now_pos) {

	int po=9+(8*(COL-1))+ROW-1;

	int stop = 0;//理論最大コンボ数

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
	//1手目を全通り探索する
	dque.clear();

	if(maxi==0){
	for (int i = 0; i < ROW; i++) {
		for (int j = 0; j < COL; j++) {
			node cand;
			cand.nowR = i;//y座標
			cand.nowC = j;//x座標
			cand.prev = -1;//1手前はない
			cand.first_te = (T_T)YX(i, j);//1手目のyx座標
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
	}// L, U,D,R //
	}
	else{
	node ca;
	ca.nowR = now_pos/COL;//y座標
	ca.nowC = now_pos%COL;//x座標
	if(maxi==1){
	ca.prev = -1;
	}
	else{
	ca.prev=prev_dir;
	}
	ca.first_te = (T_T)YX(ca.nowR,ca.nowC);
	for (int trn = 0; trn <= TRN/21; trn++) {
	ca.movei[trn] = 0ll;
	}
	F_T ff_field[ROW][COL];
	memcpy(ff_field,f_field,sizeof(ff_field));
	sc cmb;
	ll ha;
	ca.prev_score=sum_e2(ff_field,&cmb,&ha,p_maxcombo);
	ca.improving=0;
	ca.hash=ha;
	dque.push_back(ca);
	}

	int dx[DIR] = { -1, 0,0,1 },
		dy[DIR] = { 0,-1,1,0 };
	Action bestAction;//最善手
	int maxValue = 0;//最高スコア

	bestAction.maxcombo = stop;
	emilib::HashMap<ll, bool> checkNodeList[ROW*COL];

	ll rootBB[DROP+1]={0};

	for(int row=0;row<ROW;row++){
	for(int col=0;col<COL;col++){
	int pos=po-((8*col)+row);
	rootBB[f_field[row][col]]|=(1ll << (pos));
	}
	}

	//2手目以降をビームサーチで探索
	for (int i = 0; i < MAX_TRN; i++) {
		int ks = (int)dque.size();
		start = omp_get_wtime();
#pragma omp parallel for
		for (int k = 0; k < ks; k++) {
#ifdef _OPENMP
			if (i == 0 && k == 0) {
				printf("Threads[%d/%d]\n\n",
					omp_get_num_threads(),
					omp_get_max_threads());
			}
#endif
			node temp = dque[k];//que.front(); que.pop();
			F_T temp_field[ROW][COL];
			ll temp_dropBB[DROP+1]={0};
			memcpy(temp_field, f_field, sizeof(temp_field));
			memcpy(temp_dropBB,rootBB,sizeof(rootBB));
			operation(temp_field, temp.first_te,temp.movei,temp_dropBB);
			for (int j = 0; j < DIR; j++) {//上下左右の4方向が発生
				node cand = temp;
				if (0 <= cand.nowC + dx[j] && cand.nowC + dx[j] < COL &&
					0 <= cand.nowR + dy[j] && cand.nowR + dy[j] < ROW) {
					if (cand.prev + j != 3) {
						int ny=cand.nowR + dy[j];
						int nx=cand.nowC + dx[j];
						F_T field[ROW][COL];//盤面
						ll dropBB[DROP+1]={0};
						memcpy(field,temp_field,sizeof(temp_field));//盤面をもどす
						memcpy(dropBB,temp_dropBB,sizeof(temp_dropBB));
						F_T tmp=field[cand.nowR][cand.nowC];
						int pre_drop=(int)tmp;
						int pre_pos=po-((8*cand.nowC)+cand.nowR);
						int next_drop=(int)field[ny][nx];
						int next_pos=po-((8*nx)+ny);
						dropBB[pre_drop]^=(sqBB[pre_pos]|sqBB[next_pos]);
						dropBB[next_drop]^=(sqBB[pre_pos]|sqBB[next_pos]);
						field[cand.nowR][cand.nowC]=field[ny][nx];
						field[ny][nx]=tmp;
						cand.nowC += dx[j];
						cand.nowR += dy[j];
						cand.movei[i/21] |= (((ll)(j+1))<<((3*i)%63));
						//st = omp_get_wtime();
						sc cmb;
						ll ha;
						cand.score = sum_e3(dropBB, &cmb,&ha,p_maxcombo);
						cand.combo = cmb;
						cand.hash=ha;
						//part1 += omp_get_wtime() - st;
						cand.prev = j;
						//st = omp_get_wtime();
						//#pragma omp critical
											//{ pque.push(cand); }
						fff[(4 * k) + j] = cand;
						//part4 += omp_get_wtime() - st;
					}
					else {
						cand.combo = -1;
						fff[(4 * k) + j] = cand;
					}
				}
				else {
					cand.combo = -1;
					fff[(4 * k) + j] = cand;
				}
			}
		}
		//printf("depth=%d/%d\n",i+1,MAX_TRN);
		part1 += omp_get_wtime() - start;
		start = omp_get_wtime();
		dque.clear();
		deque<int>vec[3001];
		int ks2 = 0;
		for (int j = 0; j < 4 * ks; j++) {
			if (fff[j].combo != -1) {
			if (fff[j].combo == stop) {
				maxValue = stop;
				bestAction.score = maxValue;
				bestAction.first_te = fff[j].first_te;
				memcpy(bestAction.moving, fff[j].movei, sizeof(fff[j].movei));
				part2+=omp_get_wtime() - start;
				return bestAction;
			}
			if(fff[j].score>fff[j].prev_score){fff[j].improving=fff[j].improving+1;}
			fff[j].prev_score=fff[j].score;
			vec[fff[j].score+(BONUS*fff[j].improving)+(fff[j].nowR*3)+100].push_front(j);
			ks2++;
			}
		}
		part2+=omp_get_wtime() - start;
		if(i==MAX_TRN-1){return bestAction;}
		start = omp_get_wtime();
		int push_node=0;
		int possible_score=3000;
		for (int j = 0; push_node < BEAM_WIDTH ;j++) {
			if(possible_score<0){break;}
			if((int)vec[possible_score].size()==0){
			possible_score--;
			continue;
			}
			int v=vec[possible_score][0];
			node temp = fff[v];
			//swap(vec[possible_score][0], vec[possible_score].back());
			//vec[possible_score].pop_back();
			vec[possible_score].pop_front();
			if (maxValue < temp.combo) {//コンボ数が増えたらその手を記憶する
				maxValue = temp.combo;
				bestAction.score = maxValue;
				bestAction.first_te = temp.first_te;
				memcpy(bestAction.moving, temp.movei, sizeof(temp.movei));
			}
			if (i < MAX_TRN - 1) {
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
			if (field[i][j] == 0 || force) {//空マスだったらうめる
				field[i][j] = (F_T)rnd(force ? 0 : 1, DROP);//1-DROPの整数乱数
			}
		}
	}
}
int chain(int nrw, int ncl, F_T d, F_T field[ROW][COL],
	F_T chkflag[ROW][COL], F_T delflag[ROW][COL]) {
	int count = 0;
#define CHK_CF(Y,X) (field[Y][X] == d && chkflag[Y][X]==0 && delflag[Y][X] > 0)
	//連結している同じ色の消去ドロップが未探索だったら
	if (CHK_CF(nrw, ncl)) {
		++count; //連結ドロップ数の更新
		chkflag[nrw][ncl]=1;//探索済みにする
			//以下上下左右に連結しているドロップを再帰的に探索していく
		if (0 < nrw && CHK_CF(nrw - 1, ncl)) {
			count += chain(nrw - 1, ncl, d, field, chkflag, delflag);
		}
		if (nrw < ROW - 1 && CHK_CF(nrw + 1, ncl)) {
			count += chain(nrw + 1, ncl, d, field, chkflag, delflag);
		}
		if (0 < ncl && CHK_CF(nrw, ncl - 1)) {
			count += chain(nrw, ncl - 1, d, field, chkflag, delflag);
		}
		if (ncl < COL - 1 && CHK_CF(nrw, ncl + 1)) {
			count += chain(nrw, ncl + 1, d, field, chkflag, delflag);
		}
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
				if(row==0){
				GetHeight[col]=(F_T)ROW;
				}
				if(num>0 && GetHeight[col]==(F_T)ROW){
				GetHeight[col]=(F_T)row;
				}
				if (col <= COL - 3 && num == field[row][col + 1] && num == field[row][col + 2] && num > 0) {
					delflag[row][col]=1;
					delflag[row][col+1]=1;
					delflag[row][col+2]=1;
				}
				if (row <= ROW - 3 && num == field[row + 1][col] && num == field[row + 2][col] && num > 0) {
					delflag[row][col]=1;
					delflag[row+1][col]=1;
					delflag[row+2][col]=1;
				}
			}
		}
		for (int row = 0; row < ROW; row++) {
			for (int col = 0; col < COL; col++) {
				if (delflag[row][col] > 0) {
					if (chain(row, col, field[row][col], field, chkflag, delflag) >= 3) {
						cmb++;
					}
				}
			}
		}
		combo += cmb;
		//コンボが発生しなかったら終了
		if (cmb == 0 || 0 == (flag & EVAL_COMBO)) { break; }
		for (int row = 0; row < ROW; row++) {
			for (int col = 0; col < COL; col++) {
				//コンボになったドロップは空になる
				if (delflag[row][col]> 0) { field[row][col] = 0; }
			}
		}

		if (flag & EVAL_FALL){
		for(int x=0;x<COL;x++){
		fall(x,GetHeight[x],field);
		}
		}//落下処理発生
		if (flag & EVAL_SET){set(field, 0);}//落ちコン発生

	}
	return combo;
}
int evaluate2(F_T field[ROW][COL], int flag, sc* combo, ll* hash,int p_maxcombo[DROP+1]) {
	int ev = 0;
	*combo = 0;
	ll ha=0;
	int oti = 0;
	int d_maxcombo[DROP+1]={0};

	while (1) {
		int cmb = 0;
		int cmb2 = 0;
		F_T chkflag[ROW][COL]={0};
		F_T delflag[ROW][COL]={0};
		F_T GetHeight[COL];
		int cnt_drop[DROP+1]={0};
		int right[DROP+1];
		int left[DROP+1];
		for(int i=0;i<=DROP;i++){
		right[i]=-1;
		left[i]=COL;
		}	
		for (int row = 0; row < ROW; row++) {
			for (int col = 0; col < COL; col++) {
				F_T num = field[row][col];
				cnt_drop[(int)num]++;
				if(row==0){
				GetHeight[col]=(F_T)ROW;
				}
				if(num>0 && GetHeight[col]==(F_T)ROW){
				GetHeight[col]=(F_T)row;
				}
				if(oti==0){
					ha ^= zoblish_field[row][col][(int)num];
				}
				if (col <= COL - 3 && num == field[row][col + 1] && num == field[row][col + 2] && num > 0) {
					delflag[row][col]=1;
					delflag[row][col+1]=1;
					delflag[row][col+2]=1;
				}
				if (row <= ROW - 3 && num == field[row + 1][col] && num == field[row + 2][col] && num > 0) {
					delflag[row][col]=1;
					delflag[row+1][col]=1;
					delflag[row+2][col]=1;
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
						if (c == 3) { cmb2 += 30; }
						else { cmb2 += 20; }
						d_maxcombo[(int)field[row][col]]++;
					}
					field[row][col]=0;
					erase_x[col]=1;
				}
                                else{
					right[(int)field[row][col]]=max(right[(int)field[row][col]],col);
					left[(int)field[row][col]]=min(left[(int)field[row][col]],col);
                                }
			}
		}
		for(int i=1;i<=DROP;i++){
		if(right[i]!=-1&&left[i]!=COL&&cnt_drop[i]>=3&&p_maxcombo[i]!=d_maxcombo[i]){
		cmb2-=right[i]-left[i];
		}
		}
		*combo += cmb;
		ev += cmb2;
		//コンボが発生しなかったら終了
		if (cmb == 0 || 0 == (flag & EVAL_COMBO)) { break; }
		oti++;
		if (flag & EVAL_FALL){//落下処理発生
		for(int x=0;x<COL;x++){
		if(erase_x[x]==1){
		fall(x,GetHeight[x],field);
		}
		}
		}
		if (flag & EVAL_SET){set(field, 0);}//落ちコン発生

	}
	ev += oti;
	*hash=ha;
	return ev;
}
int evaluate3(ll dropBB[DROP+1], int flag, sc* combo, ll* hash,int p_maxcombo[DROP+1]) {
	int ev = 0;
	*combo = 0;
	ll ha=0;
	int oti = 0;
	ll occBB=0;

	for(int i=1;i<=DROP;i++){
	occBB|=dropBB[i];
	}

	int po=9+(8*(COL-1))+ROW-1;

	int d_maxcombo[DROP+1]={0};

	while (1) {
		int cmb = 0;
		int cmb2 = 0;

		ll linked[DROP+1]={0};

		for(int i=1;i<=DROP;i++){

		ll vert = (dropBB[i]) & (dropBB[i] << 1) & (dropBB[i] << 2);
		ll hori = (dropBB[i]) & (dropBB[i] << 8) & (dropBB[i] << 16);

		linked[i]=vert | (vert >> 1) | (vert >> 2) | hori | (hori >> 8) | (hori >> 16);

		}

		for(int i=1;i<=DROP;i++){
		long long tmp_linked=(long long)linked[i];
		while(1){
		long long chainBB=tmp_linked&(-tmp_linked);
		if(chainBB==0ll){break;}
		long long peek=chainBB;
		while (1) {
		long long tmp=tmp_linked&((long long)around(peek));
		if(peek==tmp){break;}
		peek=tmp;
		}
		int c=__builtin_popcountll(peek);
		tmp_linked^=peek;
		if(c>=3){
		cmb++;
		if(c==3){cmb2+=30;}
		else {cmb2+=20;}
		d_maxcombo[i]++;
		}
		}
		}


		for(int i=1;i<=DROP;i++){

		if(p_maxcombo[i]==d_maxcombo[i]){continue;}

		ll erased_dropBB=(dropBB[i])^(linked[i]);

		if(erased_dropBB==0ll){continue;}

		int c=__builtin_popcountll(erased_dropBB);

		if(c<3){continue;}

		long long tmp_drop=(long long)erased_dropBB;
		long long t=tmp_drop&(-tmp_drop);
		ll exist=(ll)t;
		if(exist==0ll){continue;}
		int h=( int ) ( ( exist * 0x03F566ED27179461ULL ) >> 58 );
		int pos=table[h];
		int LSB=(po-pos)/8;
		int MSB=MSB64bit(erased_dropBB);
		if(MSB==0){continue;}
		MSB=(po-MSB)/8;
		cmb2-=LSB-MSB;
		}

		if(oti==0){
		for(int i=1;i<=DROP;i++){

		long long tmp_drop=(long long)dropBB[i];

		int prev_x;
		int prev_y;

		for(int j=0;;j++){
		long long t=tmp_drop&(-tmp_drop);
		ll exist=(ll)t;
		if(exist==0ll){break;}
		int h=( int ) ( ( exist * 0x03F566ED27179461ULL ) >> 58 );
		int pos=table[h];
		int pos_x=(po-pos)/8;
		int pos_y=(po-pos)%8;
		ha ^= zoblish_field[pos_y][pos_x][i];
		prev_x=pos_x;
		prev_y=pos_y;
		tmp_drop=tmp_drop & ~(1ll<<(pos));
		}//j
		dropBB[i]^=linked[i];
		occBB^=linked[i];
		}//i
		}//if
		else{
		for(int i=1;i<=DROP;i++){
		dropBB[i]^=linked[i];
		occBB^=linked[i];
		}
		}

		*combo += cmb;
		ev += cmb2;
		//コンボが発生しなかったら終了
		if (cmb == 0 || 0 == (flag & EVAL_COMBO)) { break; }
		oti++;

		ll mask=calc_mask(occBB);

		for(int i=1;i<=DROP;i++){
		dropBB[i]=fallBB(dropBB[i],occBB,mask);
		}
		occBB=fallBB(occBB,occBB,mask);
	}
	ev += oti;
	*hash=ha;
	return ev;
}
int sum_e3(ll dropBB[DROP+1], sc* combo, ll* hash,int p_maxcombo[DROP+1]) {//落とし有り、落ちコン無し評価関数
	return evaluate3(dropBB, EVAL_FALL | EVAL_COMBO, combo,hash,p_maxcombo);
}
int sum_e2(F_T field[ROW][COL], sc* combo, ll* hash,int p_maxcombo[DROP+1]) {//落とし有り、落ちコン無し評価関数
	return evaluate2(field, EVAL_FALL | EVAL_COMBO, combo,hash,p_maxcombo);
}
int sum_e(F_T field[ROW][COL]) {//落とし有り、落ちコン無しコンボ数判定関数
	return evaluate(field, EVAL_FALL | EVAL_COMBO);
}
int sum_evaluate(F_T field[ROW][COL]) {//落としも落ちコンも有りコンボ数判定関数
	return evaluate(field, EVAL_FS | EVAL_COMBO);
}
//移動した後の盤面を生成する関数
void operation(F_T field[ROW][COL], T_T first_te,ll route[(TRN/21)+1],ll dropBB[DROP+1]) {
	int prw = (int)YY(first_te), pcl = (int)XX(first_te), i,j;
	int dx[DIR] = { -1, 0,0,1 };
	int dy[DIR] = { 0,-1,1,0 };
	int po=9+(8*(COL-1))+ROW-1;
	for (i = 0; i <= TRN/21; i++) {
		if (route[i] == 0ll) { break; }
		//移動したら、移動前ドロップと移動後ドロップを交換する
		for(j=0;j<21;j++){
		int dir = (int)(7ll&(route[i]>>(3*j)));
		if(dir==0){break;}
		int row=prw+dy[dir-1];
		int col=pcl+dx[dir-1];
		int pre_drop=(int)field[prw][pcl];
		int pre_pos=po-((8*pcl)+prw);
		int next_drop=(int)field[row][col];
		int next_pos=po-((8*col)+row);
		dropBB[pre_drop]^=(sqBB[pre_pos]|sqBB[next_pos]);
		dropBB[next_drop]^=(sqBB[pre_pos]|sqBB[next_pos]);
		F_T c = field[prw][pcl];
		field[prw][pcl] = field[row][col];
		field[row][col] = c;
		prw = row, pcl = col;
		}//j
	}//i
}
ll around(ll bitboard){
return bitboard | bitboard >> 1 | bitboard << 1 | bitboard >> 8 | bitboard << 8;
}
int rnd(int mini, int maxi) {
	static mt19937 mt((int)time(0));
	uniform_int_distribution<int> dice(mini, maxi);
	return dice(mt);
}
ll xor128() {//xorshift整数乱数
	static unsigned long long rx = 123456789, ry = 362436069, rz = 521288629, rw = 88675123;
	ll rt = (rx ^ (rx << 11));
	rx = ry; ry = rz; rz = rw;
	return (rw = (rw ^ (rw >> 19)) ^ (rt ^ (rt >> 8)));
}
ll calc_mask(ll bitboard){

ll ret=0;

for(int i=0;i<COL;i++){
ret|=fill_64[__builtin_popcountll((bitboard & file_bb[i]))] << (8*(COL-i));
}

return ret;

}
ll fallBB(ll p,ll rest, ll mask)
{
p = _pext_u64(p, rest);
p = _pdep_u64(p, mask);
return p;
}

int main() {
	
	int i, j, k;
	for(i=0;i<ROW;++i){
	for(j=0;j<COL;++j){
	for(k=0;k<=DROP;k++){
	zoblish_field[i][j][k]=xor128();
	}
	}
	}
	int po=9+(8*(COL-1))+ROW-1;
	for(i=0;i<ROW;i++){
	for(j=0;j<COL;j++){
	int pos=po-(8*j)-i;
	sqBB[pos]|=1ll<<pos;
	}
	}
	ll ha = 0x03F566ED27179461ULL;
	for (i = 0; i < 64; i++){
	table[ ha >> 58 ] = i;
	ha <<= 1;
	}

	ll res=2ll;
	fill_64[1]=res;
	for(i=2;i<64;i++){
	fill_64[i]=res+(1ll<<i);
	res=fill_64[i];
	}

	for(i=0;i<COL;i++){
	for(j=0;j<ROW;j++){
	file_bb[i] |= (1ll << (po-j));
	}
	po-=8;
	}

	string bestans="";
	string str="";
	string date="";


	double avg = 0;//平均コンボ数
	double start;
	double t_sum = 0;
	double oti_avg = 0;//平均落ちコンボ数
	for (i = 0; i < PROBLEM; i++) {//PROBLEM問解く
		F_T f_field[ROW][COL]; //スワイプ前の盤面
		F_T field[ROW][COL]; //盤面
		F_T oti_field[ROW][COL];//落ちコン用盤面
		printf("input:No.%d/%d\n", i + 1, PROBLEM);
		printf("date=");
		cin>>date;
		//init(f_field); set(f_field, 0);//初期盤面生成
		printf("layout=");
		cin>>str;
		for (j = 0; j < ROW; j++) {
			for (k = 0; k < COL; k++) {
				f_field[j][k] = (str[k+(COL*j)] - '0')+1;
			}
		}
		printf("\n");
		show_field(f_field);//盤面表示
		Action tmp,tmp2;
		double diff;
		string ans="";
		int now_pos,prev_dir;
		now_pos=0;
		prev_dir=0;
		int tesuu_min=TRN;
		for(int shots=0;shots<TRN;shots++){
		printf("\n-----search_start-----\n");
		printf("\nshots=%d/%d\n\n",shots+1,TRN);
		start = omp_get_wtime();
		tmp = BEAM_SEARCH(f_field,shots,TRN,prev_dir,now_pos);//ビームサーチしてtmpに最善手を保存
		diff = omp_get_wtime() - start;
		printf("\n-----search_end-----\n");
		t_sum += diff;
		int path_length;
		if(shots<=1){path_length=0;}
		else{path_length=shots-1;}
		for (j = 0; j <= TRN/21; j++) {//y座標は下にいくほど大きくなる
			if (tmp.moving[j] == 0ll) { break; }
			for(k=0;k<21;k++){
			int dir = (int)(7ll&(tmp.moving[j]>>(3*k)));
			if (dir==0){break;}
			path_length++;
			}
		}
		if(tesuu_min>path_length){tesuu_min=path_length;}
		tmp2=tmp;
		if(shots==0){
		ans+=to_string(XX(tmp.first_te))+to_string(YY(tmp.first_te)+5)+",";
		now_pos=((int)YY(tmp.first_te))*COL+(int)XX(tmp.first_te);
		}
		else{
		int dir = (int)(7ll&(tmp.moving[0]));
		if (dir==1) { ans+=to_string(3);swap(f_field[now_pos/COL][now_pos%COL],f_field[now_pos/COL][(now_pos%COL)-1]);now_pos--;prev_dir=dir-1;} //"LEFT"); }
		if (dir==2) { ans+=to_string(6);swap(f_field[now_pos/COL][now_pos%COL],f_field[(now_pos/COL)-1][now_pos%COL]);now_pos-=COL;prev_dir=dir-1;} //"UP"); }
		if (dir==3) { ans+=to_string(1);swap(f_field[now_pos/COL][now_pos%COL],f_field[(now_pos/COL)+1][now_pos%COL]);now_pos+=COL;prev_dir=dir-1;} //"DOWN"); }
		if (dir==4) { ans+=to_string(4);swap(f_field[now_pos/COL][now_pos%COL],f_field[now_pos/COL][(now_pos%COL)+1]);now_pos++;prev_dir=dir-1;} //"RIGHT"); }
		}
		memcpy(field,f_field,sizeof(f_field));
		int combo = sum_e(field);
		string url="";
		if(date=="null"){url="http://serizawa.web5.jp/puzzdra_theory_maker/index.html?layout="+str+"&route="+ans+"&ctwMode=false";}
		else{url="http://serizawa.web5.jp/puzzdra_theory_maker/index.html?layout="+str+"&route="+ans+"&date="+date+"&ctwMode=false";}
		if(combo==tmp.maxcombo){
		printf("\nResult=>{\n");
		printf("\ncombo=%d/%d\n",combo,tmp.maxcombo);
		printf("\npath_length=%d\n\n",shots);
		cout<<url<<endl;
		printf("\ncompleted\n");
		printf("\n}\n");
		break;
		}//if(combo
		else{
		int tesuu;
		if(shots<=1){tesuu=0;}
		else{tesuu=shots-1;}
		string route=ans;
		if(shots>0){
		route.pop_back();
		}
		for (j = 0; j <= TRN/21; j++) {//y座標は下にいくほど大きくなる
			if (tmp2.moving[j] == 0ll) { break; }
			for(k=0;k<21;k++){
			int dir = (int)(7ll&(tmp2.moving[j]>>(3*k)));
			if (dir==0){break;}
			if (dir==1) { route+=to_string(3); } //"LEFT"); }
			if (dir==2) { route+=to_string(6); } //"UP"); }
			if (dir==3) { route+=to_string(1); } //"DOWN"); }
			if (dir==4) { route+=to_string(4); } //"RIGHT"); }
			tesuu++;
			}
		}
		if(tesuu==tesuu_min){bestans=route;}
		if(date=="null"){url="http://serizawa.web5.jp/puzzdra_theory_maker/index.html?layout="+str+"&route="+route+"&ctwMode=false";}
		else{url="http://serizawa.web5.jp/puzzdra_theory_maker/index.html?layout="+str+"&route="+route+"&date="+date+"&ctwMode=false";}	
		printf("\nResult=>{\n");
		printf("\ncombo=%d/%d\n",tmp2.score,tmp2.maxcombo);
		printf("\npath_length=%d\n\n",tesuu);
		cout<<url<<endl;
		printf("\nplease wait...\n");
		printf("\n}\n");
		}
		}//shots
	}//i
	string url2="";
	if(date=="null"){url2="http://serizawa.web5.jp/puzzdra_theory_maker/index.html?layout="+str+"&route="+bestans+"&ctwMode=false";}
	else{url2="http://serizawa.web5.jp/puzzdra_theory_maker/index.html?layout="+str+"&route="+bestans+"&date="+date+"&ctwMode=false";}
	printf("\n");
	cout<<url2<<endl;
	printf("\n");
	printf("TotalDuration:%fSec\n", t_sum);
	printf("p1:%f,p2:%f,p3:%f\n", part1, part2, part3);
	cin>>i;
	cin>>j;
	cin>>k;
	return 0;
}
