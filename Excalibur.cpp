/*
Windows10,Windows11,Linux,MacOS

Linux導入手続き

//メモリ容量確認
free -h

//g++インストール
sudo apt install -y g++

//wgetインストール
sudo apt-get update
sudo apt-get install -y wget

//Excalibur.cppをダウンロード
wget --no-check-certificate https://raw.githubusercontent.com/koduma/puzzdra_solver/master/Excalibur.cpp

//hash_map.hpp,loguru.cpp,loguru.hppをダウンロード
wget --no-check-certificate https://raw.githubusercontent.com/koduma/puzzdra_solver/master/hash_map.hpp
wget --no-check-certificate https://raw.githubusercontent.com/koduma/puzzdra_solver/master/loguru.cpp
wget --no-check-certificate https://raw.githubusercontent.com/koduma/puzzdra_solver/master/loguru.hpp

//ビーム幅調整
vi Excalibur.cpp

//コンパイル
Linux:g++ -O2 -std=c++11 -fopenmp -mbmi2 -lpthread Excalibur.cpp loguru.cpp -o Excalibur -mcmodel=large -ldl
Windows10,Windows11:g++ -O2 -std=c++11 -fopenmp -mbmi2 -lpthread Excalibur.cpp loguru.cpp -o Excalibur -mcmodel=large
MacOS:g++ -std=c++11 -fopenmp -mbmi2 -lpthread Excalibur.cpp loguru.cpp -o Excalibur -ldl

//run
./Excalibur

//input

layout=367254402726710107527213362754
:path_length=52,10combo
	
layout=047631151072370164261053045210
:path_length=50,10combo
	
layout=242242100331023100110324132543
:path_length=26,9combo
	
layout=201053210251533425501353123221
:path_length=26,9combo
	
layout=015315151020442313510540210411
:path_length=27,9combo
	
layout=432015152244350331552132312515
:path_length=31,9combo
	
layout=323243441332042002331313014300
:path_length=19,8combo
	
layout=225530333313140355004550251403
:path_length=24,9combo
	
layout=224234425402054400304510125043
:path_length=30,8combo
	
layout=053241405407470557104053134522
:path_length=41,10combo

*/
#pragma warning(disable:4710)
#pragma warning(disable:4711)
#pragma warning(disable:4820)
#pragma GCC optimize("unroll-loops")
#include <vector>
#include <cfloat>
#include <cstdio>
#include <tuple>
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

#pragma GCC target ("sse4.2")

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
#define BEAM_WIDTH 2800000//MAX2800000
#define BEAM_WIDTH2 3//MAX30
#define PROBLEM 1//問題数
#define BONUS 10//評価値改善係数
#define DELTA_P 40//分解能
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define NODE_SIZE MAX(500,DIR*BEAM_WIDTH)
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
int evaluate3(ll dropBB[DROP+1], int flag, sc* combo, int p_maxcombo[DROP+1]);//落とし減点評価関数
int sum_e3(ll dropBB[DROP+1], sc* combo, int p_maxcombo[DROP+1]);//評価関数
ll around(ll bitboard);
int table[64];
ll fill_64[64];
ll file_bb[COL];
ll calc_mask(ll bitboard);
ll fallBB(ll p,ll rest,ll mask);
int shot=0;
multimap<ll, ll> visited;	
ll zoblish_field2[ROW*COL];

int MSB64bit(ll v) {
   if(v == 0ll){return 0;}
   int out =63-__builtin_clzll(v);
   return out;
}

int dfs(ll cur,int depth,emilib::HashMap<ll, bool>v){
if(v[cur]){return TRN;}
v[cur]=true;
auto p = visited.equal_range(cur);
int pl=TRN;
for (auto it = p.first; it != p.second; ++it) {
if((it->second)==(ll)1){return depth;}
}	
for (auto it = p.first; it != p.second; ++it) {
pl=min(pl,dfs(it->second,depth+1,v));
}
return pl;
}

ll c_hash(F_T board[ROW][COL]){
ll hash=0ll;
for (int row = 0; row < ROW; row++) {
for (int col = 0; col < COL; col++) {
F_T num = board[row][col];
hash ^= zoblish_field[row][col][(int)num];
}
}
return hash;
}

struct hash_chain{
	F_T field[ROW][COL];
	T_T first_te;
	ll movei[(TRN/21)+1];
	vector<ll>hashchain;
	ll check_hash(F_T board[ROW][COL]){
	return c_hash(board);
	}
	
	void calc_hashchain(){
	int pos=XX(first_te)+YY(first_te)*COL;
	hashchain.push_back(check_hash(field)^zoblish_field2[pos]);
	for (int i = 0; i <= TRN/21; i++) {//y座標は下にいくほど大きくなる
	if (movei[i] == 0ll) { break; }
	for(int k=0;k<21;k++){
	int dir = (int)(7ll&(movei[i]>>(3*k)));
	if (dir==0){break;}
	if (dir==1) { swap(field[pos/COL][pos%COL],field[(pos-1)/COL][(pos-1)%COL]);pos--; } //"LEFT"); }
	if (dir==2) { swap(field[pos/COL][pos%COL],field[(pos-COL)/COL][(pos-COL)%COL]);pos-=COL; } //"UP"); }
	if (dir==3) { swap(field[pos/COL][pos%COL],field[(pos+COL)/COL][(pos+COL)%COL]);pos+=COL; } //"DOWN"); }
	if (dir==4) { swap(field[pos/COL][pos%COL],field[(pos+1)/COL][(pos+1)%COL]);pos++; } //"RIGHT"); }
	hashchain.push_back(check_hash(field)^zoblish_field2[pos]);
	}
	}
	}
};
struct node {//どういう手かの構造体
	ll movei[(TRN/21)+1];//スワイプ移動座標
	ll hash;//盤面のハッシュ値
	int score;//評価値
	int prev_score;//1手前の評価値
	T_T first_te;
	uc improving;//評価値改善回数
	sc combo;//コンボ数
	sc pos;//今どの座標にいるか
	sc prev;//1手前は上下左右のどっちを選んだか
	node() {//初期化
		this->score = 0;
		this->prev = -1;
		this->prev_score = -1;
		//memset(this->movei, STP, sizeof(this->movei));
	}
	bool operator < (const node& n)const {//スコアが高い方が優先される
		return score < n.score;
	}
}fff[NODE_SIZE];

struct node2 {

	F_T field[ROW][COL];
	T_T first_te;
	ll movei[(TRN/21)+1];
	string path;
	int path_length;
	int pos;
	sc prev;
	ll hash;
	string true_path;
	int true_path_length;
	int second_score;

	void calc_path(){

	path_length=0;

	path=to_string(XX(first_te))+to_string(YY(first_te)+5)+",";

	for (int i = 0; i <= TRN/21; i++) {//y座標は下にいくほど大きくなる
	if (movei[i] == 0ll) { break; }
	for(int k=0;k<21;k++){
	int dir = (int)(7ll&(movei[i]>>(3*k)));
	if (dir==0){break;}
	if (dir==1) { path+=to_string(3); } //"LEFT"); }
	if (dir==2) { path+=to_string(6); } //"UP"); }
	if (dir==3) { path+=to_string(1); } //"DOWN"); }
	if (dir==4) { path+=to_string(4); } //"RIGHT"); }
	path_length++;
	}
	}
	}

	void calc_hash(){
	
	hash=c_hash(field);
	}
	
	int calc_pl(ll cur){
	emilib::HashMap<ll, bool>v;
	return dfs(cur,0,v);	
	}

}ff[DIR*BEAM_WIDTH2];

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

int adder(F_T field[ROW][COL]){
  
    int x_cnt[DROP+1][COL]={0};
    
    for(int r=0;r<ROW;r++){
    for(int c=0;c<COL;c++){
        x_cnt[(int)field[r][c]][c]++;
    }
    }
    
    int ret=0;
     
    for(int d=1;d<=DROP;d++){
    for(int c=0;c<COL;c++){       
    for(int pos=c+1;pos<COL;pos++){
    int dx=pos-c;
    int dd=x_cnt[d][pos];
    ret+=dx*dd*x_cnt[d][c];
    }
    }
    }

    return ret;
}
double logN(double b, double x) {
    return log(x) / log(b);
}
Action BEAM_SEARCH(F_T f_field[ROW][COL],int maxi,int MAX_TRN,int prev_dir,int now_pos,int stop); //ルート探索関数
double part1 = 0, part2 = 0, part3 = 0, MAXCOMBO = 0;
Action BEAM_SEARCH(F_T f_field[ROW][COL],int maxi,int MAX_TRN,int prev_dir,int now_pos,int stop) {
	shot++;if(shot%1000==0){cout<<"shot="<<(shot)<<endl;}
	int po=9+(8*(COL-1))+ROW-1;

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
			cand.pos = (i*COL)+j;
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
	ca.pos = now_pos;
	if(maxi==1){
	ca.prev = -1;
	}
	else{
	ca.prev=prev_dir;
	}
	ca.first_te = (T_T)YX((ca.pos)/COL,(ca.pos)%COL);
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
	if((int)cmb>=stop){
	Action accept;
	accept.first_te=ca.first_te;
	accept.score=stop;
	accept.maxcombo=stop;
	for (int trn = 0; trn <= TRN/21; trn++) {
	accept.moving[trn] = 0ll;
	}
	return accept;
	}
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

	int Q=(int)floor(logN(3.0,(double)BEAM_WIDTH))+1;

	//2手目以降をビームサーチで探索
	for (int i = 0; i < MAX_TRN; i++) {
		int ks = (int)dque.size();
		start = omp_get_wtime();
#pragma omp parallel for
		for (int k = 0; k < ks; k++) {
			node temp = dque[k];//que.front(); que.pop();
			F_T temp_field[ROW][COL];
			ll temp_dropBB[DROP+1]={0};
			memcpy(temp_field, f_field, sizeof(temp_field));
			memcpy(temp_dropBB,rootBB,sizeof(rootBB));
			operation(temp_field, temp.first_te,temp.movei,temp_dropBB);
			for (int j = 0; j < DIR; j++) {//上下左右の4方向が発生
				node cand = temp;
				int xx = (cand.pos)%COL;
				int yy = (cand.pos)/COL;
				if (0 <= xx + dx[j] && xx + dx[j] < COL &&
					0 <= yy + dy[j] && yy + dy[j] < ROW) {
					if (cand.prev + j != 3) {
						int ny=yy + dy[j];
						int nx=xx + dx[j];
						F_T field[ROW][COL];//盤面
						ll dropBB[DROP+1]={0};
						memcpy(field,temp_field,sizeof(temp_field));//盤面をもどす
						memcpy(dropBB,temp_dropBB,sizeof(temp_dropBB));
						F_T tmp=field[yy][xx];
						cand.hash^=(zoblish_field[yy][xx][tmp])^(zoblish_field[ny][nx][field[ny][nx]]);
						cand.hash^=(zoblish_field[yy][xx][field[ny][nx]])^(zoblish_field[ny][nx][tmp]);
						int pre_drop=(int)tmp;
						int pre_pos=po-((8*xx)+yy);
						int next_drop=(int)field[ny][nx];
						int next_pos=po-((8*nx)+ny);
						dropBB[pre_drop]^=(sqBB[pre_pos]|sqBB[next_pos]);
						dropBB[next_drop]^=(sqBB[pre_pos]|sqBB[next_pos]);
						field[yy][xx]=field[ny][nx];
						field[ny][nx]=tmp;
						cand.pos += (dy[j]*COL)+dx[j];
						cand.movei[i/21] |= (((ll)(j+1))<<((3*i)%63));
						//st = omp_get_wtime();
						//if(i<Q&&MAX_TRN>20&&CUT){
						//cand.score=0;
						//cand.combo=0;						
						//}
						//else{
						sc cmb;
						cand.score = evaluate3(dropBB, EVAL_FALL | EVAL_COMBO, &cmb,p_maxcombo);
						cand.score -= adder(field);
						cand.combo = cmb;
						//}
						//part1 += omp_get_wtime() - st;
						cand.prev = j;
						//st = omp_get_wtime();
						//#pragma omp critical
											//{ pque.push(cand); }
						fff[(DIR * k) + j] = cand;
						//part4 += omp_get_wtime() - st;
					}
					else {
						cand.combo = -1;
						fff[(DIR * k) + j] = cand;
					}
				}
				else {
					cand.combo = -1;
					fff[(DIR * k) + j] = cand;
				}
			}
		}
		//printf("depth=%d/%d\n",i+1,MAX_TRN);
		part1 += omp_get_wtime() - start;
		start = omp_get_wtime();
		dque.clear();
		deque<int>vec[5001];
		int ks2 = 0;
		bool congrats=false;
		for (int j = 0; j < DIR * ks; j++) {
			if (fff[j].combo != -1) {
			if (fff[j].combo >= stop) {
				maxValue = fff[j].combo;
				bestAction.score = maxValue;
				bestAction.first_te = fff[j].first_te;
				memcpy(bestAction.moving, fff[j].movei, sizeof(fff[j].movei));
				hash_chain hc;
				F_T abc[ROW][COL];
				memcpy(hc.field,f_field,sizeof(abc));
				hc.first_te=fff[j].first_te;
				memcpy(hc.movei, fff[j].movei, sizeof(fff[j].movei));
				hc.calc_hashchain();
				if((int)hc.hashchain.size()>0){
				for(int r=0;r<(int)hc.hashchain.size()-1;r++){
				ll cur=hc.hashchain[r];
				ll nexthash=hc.hashchain[r+1];
				bool find=false;
				auto p = visited.equal_range(cur);
				for (auto it = p.first; it != p.second; ++it) {
				if(it->second==nexthash){find=true;break;}
				}
				if(!find){
				visited.emplace(cur,nexthash);
				}
				if(r==(int)hc.hashchain.size()-2){
				find=false;
				p=visited.equal_range(nexthash);
				for (auto it = p.first; it != p.second; ++it) {
				if(it->second==(ll)1){find=true;break;}
				}
				if(!find){
				visited.emplace(nexthash,(ll)1);
				}
				}
				}
				}
				congrats=true;
			}
			if(fff[j].score>fff[j].prev_score){fff[j].improving=fff[j].improving+1;}
			fff[j].prev_score=fff[j].score;
			vec[fff[j].score+(BONUS*fff[j].improving)+(((fff[j].pos)/COL)*3)+2000].push_front(j);
			ks2++;
			}
		}
		part2+=omp_get_wtime() - start;
		if(congrats){return bestAction;}
		start = omp_get_wtime();
		int push_node=0;
		int possible_score=5000;
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
			int pos=temp.pos;
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
string BEAM_SEARCH2(F_T field[ROW][COL],int MAX_TRN); //ルート探索関数
string BEAM_SEARCH2(F_T field[ROW][COL],int MAX_TRN) {

/*
	T_T first_te;
	ll movei[(TRN/21)+1];//スワイプ移動座標
	int score;//評価値
	sc combo;//コンボ数
	sc pos;//今どの座標にいるか
	sc prev;//1手前は上下左右のどっちを選んだか
	int prev_score;//1手前の評価値
	uc improving;//評価値改善回数
	ll hash;//盤面のハッシュ値
	F_T field[ROW][COL];
	T_T first_te;
	ll movei[(TRN/21)+1];
	string path;
	int path_length;
	int pos;
	sc prev;
	ll hash;
	string true_path;
	int true_path_length;
*/

	int LIM=TRN;
	int stop=0;
	int drop[DROP + 1] = { 0 };
	int p_maxcombo[DROP + 1] = { 0 };
	for (int row = 0; row < ROW; row++) {
		for (int col = 0; col < COL; col++) {
			if (1 <= field[row][col] && field[row][col] <= DROP) {
				drop[field[row][col]]++;
			}
		}
	}
	for (int i = 1; i <= DROP; i++) {
		stop += drop[i] / 3;
		p_maxcombo[i] = drop[i] / 3;
	}

	double start = omp_get_wtime();

	printf("\n-----search_start_1/2-----\n");

	vector<node2>dque;

	deque<node2>pus[TRN+1];

	double avg=0;

	double path_length_array[ROW][COL];
	
	vector<int>pro_league;

	for (int i = 0; i < ROW; i++) {
	for (int j = 0; j < COL; j++) {
	node2 cand,cand2;
	int MLEN=cand2.calc_pl(c_hash(field)^zoblish_field2[(i*COL)+j]);
	if(i==0&&j==0){MLEN=TRN;}	
	int lim=LIM;
	if((int)pro_league.size()>=BEAM_WIDTH2){lim=pro_league[BEAM_WIDTH2-1];LIM=lim;}	
	Action tmp=BEAM_SEARCH(field,1,max(0,min(lim-1,MLEN-1)),-1,(i*COL)+j,stop);
	if(i==0&&j==0){stop=0;}
	stop=max(stop,tmp.score);
	F_T f_field[ROW][COL];
	memcpy(f_field,field,sizeof(f_field));
	if(sum_e(f_field)>=stop){return "05,";}	
	memcpy(cand.field,field,sizeof(f_field));
	cand.first_te = tmp.first_te;
	for (int trn = 0; trn <= TRN/21; trn++) {
	cand.movei[trn] = tmp.moving[trn];
	}
	cand.calc_path();
	cand.pos = (i*COL)+j;
	cand.prev = -1;
	cand.calc_hash();
	cand.true_path=to_string(j)+to_string(i+5)+",";
	cand.true_path_length=0;
	if(stop!=tmp.score){cand.path_length=TRN;}
	//printf("beam=%d,visited=%d\n",cand.path_length,MLEN);
	if(cand.path_length>MLEN){
		string layout="";
		for(int b=0;b<ROW*COL;b++){
		layout+=field[b/COL][b%COL]-1+'0';
		}
		//cout<<"layout="<<layout<<endl;
		//cout<<"prev="<<(int)cand.prev<<endl;
		//cout<<"pos="<<cand.pos+1<<endl;
	}
	cand.path_length=min(cand.path_length,MLEN);
	pus[cand.path_length].push_front(cand);
	cout<<"pos="<<cand.pos+1<<"/"<<ROW*COL<<endl;
	cout<<"path_length="<<cand.path_length<<endl;
	cout<<"combo="<<tmp.score<<"/"<<stop<<endl;
	if(cand.path_length<lim){
	avg+=(double)cand.path_length;
	pro_league.push_back(cand.path_length);
	sort(pro_league.begin(),pro_league.end());	
	}
	path_length_array[i][j]=(double)cand.path_length;	
	}
	}
	double delta_t = omp_get_wtime()-start;

	double variance=0;
	
	if((int)pro_league.size()>0){
	avg/=(double)(pro_league.size());
	}

	for (int i = 0; i < ROW; i++) {
	for (int j = 0; j < COL; j++) {
	if(path_length_array[i][j]<=149.0){
	variance+=pow(fabs(path_length_array[i][j]-avg),3.0);
	}
	}
	}

	if(variance<0.0001){printf("\ndifficulty=INF\n");}
	else{printf("\ndifficulty=%f\n",delta_t*(10000.0/variance)*(10000.0/variance));}

	printf("\n-----search_start_2/2-----\n");
	int cnt=0;

	for(int i=0;i<TRN;i++){
	if((int)pus[i].size()==0){continue;}
	while(1){
	if((int)pus[i].size()==0){break;}
	node2 cand=pus[i][0];
	pus[i].pop_front();
	if(cnt<BEAM_WIDTH2){
	dque.push_back(cand);
	cnt++;
	}
	else{break;}
	}
	}

	int dx[DIR] = { -1, 0,0,1 };
	int dy[DIR] = { 0,-1,1,0 };
	string bestAction;
	int maxValue = 0;

	emilib::HashMap<ll, bool> checkNodeList[ROW*COL];

	for (int i = 0; i < MAX_TRN; i++) {
	int ks = (int)dque.size();
	pro_league.clear();	
	for (int k = 0; k < ks; k++) {

	node2 temp = dque[k];
	for (int j = 0; j < DIR; j++) {
	node2 cand = temp;
	int x=cand.pos%COL;
	int y=cand.pos/COL;
	if (0 <= x + dx[j] && x + dx[j] < COL &&0 <= y + dy[j] && y + dy[j] < ROW) {
	if (cand.prev + j != 3) {
	F_T f_field[ROW][COL];
	memcpy(f_field,cand.field,sizeof(f_field));
	int nx=x + dx[j];
	int ny=y + dy[j];
	swap(f_field[y][x],f_field[ny][nx]);
	cand.pos = (ny*COL)+nx;
	if(j==0){cand.true_path+=to_string(3);}
	else if(j==1){cand.true_path+=to_string(6);}
	else if(j==2){cand.true_path+=to_string(1);}
	else{cand.true_path+=to_string(4);}
	cand.prev = j;
	int MLEN=cand.calc_pl(c_hash(f_field)^zoblish_field2[cand.pos]);
	int lim=LIM;
	if((int)pro_league.size()>=BEAM_WIDTH2){lim=pro_league[BEAM_WIDTH2-1];LIM=lim;}		
	Action tmp = BEAM_SEARCH(f_field,i+2,max(0,min(lim-1,MLEN-1)),cand.prev,cand.pos,stop);
	sc cmb;
	ll ha;
	F_T tmp_field[ROW][COL];
	memcpy(tmp_field,f_field,sizeof(tmp_field));
	cand.second_score = evaluate2(tmp_field, EVAL_FALL | EVAL_COMBO,&cmb ,&ha,p_maxcombo);
	cand.first_te = tmp.first_te;
	for (int trn = 0; trn <= TRN/21; trn++) {
	cand.movei[trn] = tmp.moving[trn];
	}
	memcpy(cand.field,f_field,sizeof(f_field));
	cand.calc_path();
	cand.calc_hash();
	if(stop!=tmp.score){cand.path_length=TRN;}
	//printf("beam=%d,visited=%d\n",cand.path_length,MLEN);
	if(cand.path_length>MLEN){
	string layout="";
	for(int b=0;b<ROW*COL;b++){
	layout+=f_field[b/COL][b%COL]-1+'0';
	}
	//cout<<"layout="<<layout<<endl;
	//cout<<"prev="<<(int)cand.prev<<endl;
	//cout<<"pos="<<cand.pos+1<<endl;
	}	
	cand.path_length=min(cand.path_length,MLEN);
	if(cand.path_length<lim){	
	pro_league.push_back(cand.path_length);
	sort(pro_league.begin(),pro_league.end());
	}	
	ff[(DIR * k) + j] = cand;
	}//if(cand.prev
	else {
	cand.path_length = -1;
	ff[(DIR * k) + j] = cand;
	}
	}//if(0<=x+dx[j]
	else {
	cand.path_length = -1;
	ff[(DIR * k) + j] = cand;
	}
	}//for(int j=0;
	}//for(int k=0;

	printf("depth=%d/%d\n",i+1,MAX_TRN);
	dque.clear();
	vector<tuple<int, int, int> > vec;
	for(int j=0;j<DIR*ks;j++){
	if(ff[j].path_length!=-1&&ff[j].path_length-i<10){
	F_T f_field[ROW][COL];
	memcpy(f_field,ff[j].field,sizeof(f_field));
	int combo = sum_e(f_field);
	if(combo>=stop){return ff[j].true_path;}
	}
	if(ff[j].path_length!=-1){	
	vec.push_back(make_tuple(ff[j].path_length,-ff[j].second_score,j));
	}
	}
	sort(vec.begin(), vec.end());	
	int push_node=0;
	int ps=get<0>(vec[0]);
	int opt1=TRN;	
	for (int j = 0; push_node < BEAM_WIDTH2 ;j++) {
	if(j>=(int)vec.size()){break;}
	int v=get<2>(vec[j]);
	if(push_node==0){
	printf("predict=%d\n",i+1+ps);	
	opt1=ps;
	}	
	node2 temp = ff[v];
	F_T f_field[ROW][COL];
	memcpy(f_field,temp.field,sizeof(f_field));
	int combo = sum_e(f_field);
	if (maxValue < combo) {//コンボ数が増えたらその手を記憶する
	maxValue = combo;
	bestAction=temp.true_path;
	}
	if (i < MAX_TRN - 1) {
	if(!checkNodeList[temp.pos][temp.hash]){
	checkNodeList[temp.pos][temp.hash]=true;
	if(ps<=opt1+DELTA_P){	
	dque.push_back(temp);
	push_node++;
	}	
	}
	}
	}
	}//for(int i=0;

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

		//cmb2*=4;
		cmb2+=cmb2;
		cmb2+=cmb2;

		for(int s=0;s<=COL-3;s++){
		int same_num[DROP+1]={0};
		for(int col=s;col<=s+2;col++){
		for(int row=0;row<ROW;row++){
		same_num[field[row][col]]++;
		}
		}
		for(int i=1;i<=DROP;i++){
		if(p_maxcombo[i]!=d_maxcombo[i]){
		cmb2+=same_num[i];
		}
		}
		}

		for(int col=0;col<COL;col++){
		int y_bonus[DROP+1]={0};
		for(int row=0;row<ROW;row++){
		y_bonus[field[row][col]]++;
		}
		for(int i=1;i<=DROP;i++){
		if(y_bonus[i]>=3){cmb2+=20;}
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
	
	int penalty=0;

	for(int i=1;i<=DROP;i++){
	penalty+=(p_maxcombo[i]-d_maxcombo[i])*10;
	}

	int alone=0;
	
	bool find=false;

	for(int x=0;x<COL;x++){
	if(field[ROW-1][x] == 0){
	if(!find){alone++;}
	find=true;	
	}
	else{find=false;}	
	}

	ev-=penalty*alone;
  
	*hash=ha;
	return ev;
}
int evaluate3(ll dropBB[DROP+1], int flag, sc* combo, int p_maxcombo[DROP+1]) {
	int ev = 0;
	*combo = 0;
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

		ll erased_dropBB=dropBB[i];

		if(erased_dropBB==0ll){continue;}

		int c=__builtin_popcountll(erased_dropBB);

		if(c<3){continue;}
			
		erased_dropBB^=linked[i];
		c=__builtin_popcountll(erased_dropBB);
		if(c<2){continue;}	

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

		for(int i=1;i<=DROP;i++){
		dropBB[i]^=linked[i];
		occBB^=linked[i];
		}

		//cmb2*=4;
		cmb2+=cmb2;
		cmb2+=cmb2;

		for(int s=0;s<=COL-3;s++){
		int same_num[DROP+1]={0};
		ll bp=0ll;
		for(int col=s;col<=s+2;col++){
		bp+=file_bb[col];
		}
		for(int i=1;i<=DROP;i++){
		if(p_maxcombo[i]!=d_maxcombo[i]){
		same_num[i]+=__builtin_popcountll(bp&dropBB[i]);
		cmb2+=same_num[i];
		}
		}
		}

		for(int col=0;col<COL;col++){
		ll bp=file_bb[col];
		for(int i=1;i<=DROP;i++){
		int yb=__builtin_popcountll(bp&dropBB[i]);
		if(yb>=3){cmb2+=20;}
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
	
	int penalty=0;

	ll board=occBB;

	for(int i=1;i<=DROP;i++){
	penalty+=(p_maxcombo[i]-d_maxcombo[i])*10;
	}

	int alone=0;
	
	bool find=false;

	for(int x=0;x<COL;x++){
	if(((board>>(po-((8*(x))+(ROW-1))))&1) == 0){
	if(!find){alone++;}
	find=true;	
	}
	else{find=false;}
	}

	ev-=penalty*alone;
  
	return ev;
}
int sum_e3(ll dropBB[DROP+1], sc* combo, int p_maxcombo[DROP+1]) {//落とし有り、落ちコン無し評価関数
	return evaluate3(dropBB, EVAL_FALL | EVAL_COMBO, combo,p_maxcombo);
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
	
	for(i=0;i<ROW*COL;i++){
	zoblish_field2[i]=xor128();
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
	string layout="";
	string date="";
	string url="";

	double avg = 0;//平均コンボ数
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
		cin>>layout;
		for (j = 0; j < ROW; j++) {
			for (k = 0; k < COL; k++) {
				f_field[j][k] = (layout[k+(COL*j)] - '0')+1;
			}
		}
		printf("\n");
		show_field(f_field);//盤面表示
		printf("\n");
		double start = omp_get_wtime();
		bestans=BEAM_SEARCH2(f_field,TRN);
		if(date=="null"){url="http://serizawa.web5.jp/puzzdra_theory_maker/index.html?layout="+layout+"&route="+bestans+"&ctwMode=false";}
		else{url="http://serizawa.web5.jp/puzzdra_theory_maker/index.html?layout="+layout+"&route="+bestans+"&date="+date+"&ctwMode=false";}
		double diff = omp_get_wtime() - start;
		t_sum += diff;
		int tgt=0;
		string top="";
		while(1){

		if(bestans[tgt]==','){tgt++;break;}
		top+=bestans[tgt];
		tgt++;

		}
		int pos;
		if((int)top.size()==2){int x=top[0]-'0';int y=(top[1]-'0')-5;pos=(y*COL)+x;}
		else{int x=top[0]-'0';int y=5;pos=(y*COL)+x;}

		for(j=tgt;j<(int)bestans.size();j++){
		if(bestans[j]=='3'){swap(f_field[pos/COL][pos%COL],f_field[pos/COL][(pos%COL)-1]);pos--;}
		if(bestans[j]=='6'){swap(f_field[pos/COL][pos%COL],f_field[(pos/COL)-1][pos%COL]);pos-=COL;}
		if(bestans[j]=='1'){swap(f_field[pos/COL][pos%COL],f_field[(pos/COL)+1][pos%COL]);pos+=COL;}
		if(bestans[j]=='4'){swap(f_field[pos/COL][pos%COL],f_field[pos/COL][(pos%COL)+1]);pos++;}
		}
		int maxcombo=0;
		int drop[DROP + 1] = { 0 };
		for (int row = 0; row < ROW; row++) {
		for (int col = 0; col < COL; col++) {
		if (1 <= f_field[row][col] && f_field[row][col] <= DROP) {
		drop[f_field[row][col]]++;
		}
		}
		}
		for (int loop = 1; loop <= DROP; loop++) {
		maxcombo += drop[loop] / 3;
		}
		int combo=sum_e(f_field);
		printf("\nResult=>{\n");
		printf("\ncombo=%d/%d\n",combo,maxcombo);
		int si=(int)bestans.size()-(int)top.size()-1;
		printf("\npath_length=%d\n\n",si);
		cout<<url<<endl;
		printf("\ncompleted\n");
		printf("\n}\n");

	}//i
	printf("TotalDuration:%fSec\n", t_sum);
	printf("p1:%f,p2:%f,p3:%f,p4:%f\n", part1, part2, part3,t_sum-(part1+part2+part3));
	cin>>i;
	cin>>j;
	cin>>k;
	return 0;
}
