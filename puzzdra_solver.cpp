/*
puzzdra_solver

パズドラソルバです。

printf("TotalDuration:%fSec\n", t_sum);
printf("Avg.Combo #:%lf\n", avg / (double)i);

これが改善されればpull request受け付けます。

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
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;
#define DLT(ST,ED) ((double)((ED)-(ST))/CLOCKS_PER_SEC)//時間差分
#define XX(PT)  ((PT)&15)
#define YY(PT)  XX((PT)>>4)
#define YX(Y,X) ((Y)<<4|(X))
#define DIR 4//方向
#define ROW 5//縦
#define COL 6//横
#define TRN  44//手数
#define STP YX(7,7)//無効手[無効座標]
#define MAX_TURN 40//最大ルート長
#define BEAM_WIDTH 50000//ビーム幅
typedef char F_T;//盤面型
typedef char T_T;//手数型
enum { EVAL_NONE = 0, EVAL_FALL, EVAL_SET, EVAL_FS, EVAL_COMBO };
void init(F_T field[ROW][COL]); //初期配置生成関数
void fall(F_T field[ROW][COL]); //ドロップの落下処理関数
void set(F_T field[ROW][COL], int force); //空マスを埋める関数
void show_field(F_T field[ROW][COL]); //盤面表示関数
unsigned int rnd(int mini, int maxi); //整数乱数
//上下左右に連結しているドロップを再帰的に探索していく関数
int chain(int nrw, int ncl, int d, F_T field[ROW][COL], F_T chkflag[ROW][COL], F_T delflag[ROW][COL]);
int evaluate(F_T field[ROW][COL], int flag); //コンボ数判定関数
int sum_e(F_T field[ROW][COL]);//落とし有り、落ちコン無しコンボ数判定関数
int sum_evaluate(F_T field[ROW][COL]);//落としも落ちコンも有りコンボ数判定関数
void operation(F_T field[ROW][COL], T_T route[TRN]); //スワイプ処理関数
struct member {//どういう手かの構造体
	T_T movei[TRN];//スワイプ移動座標
	int score;//コンボ数
	int nowC;//今どのx座標にいるか
	int nowR;//今どのy座標にいるか
	int prev;//1手前は上下左右のどっちを選んだか
	member() {//初期化
		this->score = 0;
		this->prev = -1;
		//memset(this->movei, STP, sizeof(this->movei));
	}
	bool operator < (const member& n)const {//スコアが高い方が優先される
		return score < n.score;
	}
}fff[BEAM_WIDTH * 4];
struct Action {//最終的に探索された手
	int score;//コンボ数
	T_T moving[TRN];//スワイプ移動座標
	Action() {//初期化
		this->score = 0;
		//memset(this->moving, STP, sizeof(this->moving));
	}
};
Action BEAM_SEARCH(F_T f_field[ROW][COL]); //ルート探索関数
double part1 = 0, part2 = 0, part3 = 0, part4 = 0;
Action BEAM_SEARCH(F_T f_field[ROW][COL]) {

	int stop = 0;//理論最大コンボ数

	int drop[7] = { 0 };
	for (int row = 0; row < ROW; row++) {
		for (int col = 0; col < COL; col++) {
			drop[f_field[row][col]]++;
		}
	}
	for (int i = 1; i <= 6; i++) {
		stop += drop[i] / 3;
	}

	deque<member>dque;
	clock_t start, st;
	//1手目を全通り探索する
	dque.clear();
	for (int i = 0; i < ROW; i++) {
		for (int j = 0; j < COL; j++) {
			member cand;
			cand.nowR = i;//y座標
			cand.nowC = j;//x座標
			cand.prev = -1;//1手前はない
			memset(cand.movei, STP, sizeof(cand.movei));
			cand.movei[0] = (T_T)YX(i, j);//1手目のyx座標
			dque.push_back(cand);
		}
	}         // L, U,D,R //
	int dx[DIR] = { -1, 0,0,1 },
		dy[DIR] = { 0,-1,1,0 };
	Action bestAction;//最善手
	int maxValue = 0;//最高スコア

	//2手目以降をビームサーチで探索
	for (int i = 1; i < MAX_TURN; i++) {
		//priority_queue<member, vector<member>, less<member> >pque;
		int ks = (int)dque.size();
		start = clock();
#pragma omp parallel for reduction(+:part1,part4)
		for (int k = 0; k < ks; k++) {
#ifdef _OPENMP
			if (i == 1 && k == 0) {
				printf("Threads[%d/%d]\n",
					omp_get_num_threads(),
					omp_get_max_threads());
			}
#endif
			member temp = dque[k];//que.front(); que.pop();
			for (int j = 0; j < DIR; j++) {//上下左右の4方向が発生
				member cand = temp;
				if (0 <= cand.nowC + dx[j] && cand.nowC + dx[j] < COL &&
					0 <= cand.nowR + dy[j] && cand.nowR + dy[j] < ROW) {
					if (cand.prev + j != 3) {
						cand.nowC += dx[j];
						cand.nowR += dy[j];
						cand.movei[i] = (T_T)YX(cand.nowR, cand.nowC);
						st = clock();
						F_T field[ROW][COL]; //盤面
						memcpy(field, f_field, sizeof(field));//盤面をもどす
						operation(field, cand.movei);
						cand.score = sum_e(field);
						memcpy(field, f_field, sizeof(field));//盤面をもどす
						operation(field, cand.movei);
						part1 += DLT(st, clock());
						cand.prev = j;
						st = clock();
						//#pragma omp critical
											//{ pque.push(cand); }
						fff[(4 * k) + j] = cand;
						part4 += DLT(st, clock());
					}
					else {
						cand.score = -114514;
						fff[(4 * k) + j] = cand;
					}
				}
				else {
					cand.score = -114514;
					fff[(4 * k) + j] = cand;
				}
			}
		}
		part2 += DLT(start, clock());
		start = clock();
		dque.clear();
		vector<pair<int, int> >vec;
		for (int j = 0; j < 4 * ks; j++) {
			vec.push_back(make_pair(fff[j].score, j));
		}
		sort(vec.begin(), vec.end());
		reverse(vec.begin(), vec.end());
		for (int j = 0; j < BEAM_WIDTH && j < 4 * ks; j++) {
			member temp = fff[vec[j].second];
			if (temp.score == -114514) { continue; }
			if (maxValue < temp.score) {//コンボ数が増えたらその手を記憶する
				maxValue = temp.score;
				bestAction.score = maxValue;
				memcpy(bestAction.moving, temp.movei, sizeof(temp.movei));
				//コンボ数が理論値になったらreturn
				if (temp.score == stop) { return bestAction; }
			}
			if (i < MAX_TURN - 1) {
				dque.push_back(temp);
			}
		}
		part3 += DLT(start, clock());
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
void fall(F_T field[ROW][COL]) {
	for (int j = 0; j < COL; j++) {
		int tgt;
		for (tgt = ROW - 1; tgt >= 0 && field[tgt][j] != 0; tgt--);
		for (int i = tgt - 1; i >= 0; i--) {
			if (field[i][j] != 0) {
				F_T c = field[i][j];
				field[i][j] = 0;
				field[tgt][j] = c;
				tgt--;
			}
		}
	}
}
void init(F_T field[ROW][COL]) { set(field, !0); }
void set(F_T field[ROW][COL], int force) {
	for (int i = 0; i < ROW; i++) {
		for (int j = 0; j < COL; j++) {
			if (field[i][j] == 0 || force) {//空マスだったらうめる
				field[i][j] = (F_T)rnd(force ? 0 : 1, 6);//1-6の整数乱数
			}
		}
	}
}
int chain(int nrw, int ncl, int d, F_T field[ROW][COL],
	F_T chkflag[ROW][COL], F_T delflag[ROW][COL]) {
	int count = 0;
#define CHK_CF(Y,X) (field[Y][X] == d && chkflag[Y][X] == 0)
	//連結している同じ色のドロップが未探索だったら
	if (CHK_CF(nrw, ncl)) {
		++count; //連結ドロップ数の更新
		chkflag[nrw][ncl] =   //探索済みにする
			delflag[nrw][ncl] = 1;//コンボがつながる可能性があるので、1に設定
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
		F_T chkflag[ROW][COL] = { 0 };
		F_T delflag[ROW][COL] = { 0 };
		F_T t_erase[ROW][COL] = { 0 };
		for (int row = 0; row < ROW; row++) {
			for (int col = 0; col < COL; col++) {
				if (col <= COL - 3 && field[row][col] == field[row][col + 1] && field[row][col] == field[row][col + 2] && field[row][col] > 0) {
					delflag[row][col] = field[row][col];
					delflag[row][col + 1] = field[row][col];
					delflag[row][col + 2] = field[row][col];
				}
				if (row <= ROW - 3 && field[row][col] == field[row + 1][col] && field[row][col] == field[row + 2][col] && field[row][col] > 0) {
					delflag[row][col] = field[row][col];
					delflag[row + 1][col] = field[row][col];
					delflag[row + 2][col] = field[row][col];
				}
			}
		}
		for (int row = 0; row < ROW; row++) {
			for (int col = 0; col < COL; col++) {
				if (delflag[row][col] > 0) {
					if (chain(row, col, field[row][col], field, chkflag, t_erase) >= 3) {
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
				//コンボになったドロップは空になる。
				if (t_erase[row][col] != 0 && delflag[row][col] > 0) { field[row][col] = 0; }
			}
		}

		if (flag & EVAL_FALL)fall(field);//落下処理発生
		if (flag & EVAL_SET)set(field, 0);//落ちコン発生

	}
	return combo;
}

int sum_e(F_T field[ROW][COL]) {//落とし有り、落ちコン無しコンボ数判定関数
	return evaluate(field, EVAL_FALL | EVAL_COMBO);
}
int sum_evaluate(F_T field[ROW][COL]) {//落としも落ちコンも有りコンボ数判定関数
	return evaluate(field, EVAL_FS | EVAL_COMBO);
}
//移動した後の盤面を生成する関数
void operation(F_T field[ROW][COL], T_T route[TRN]) {
	int prw = (int)YY(route[0]), pcl = (int)XX(route[0]), i;
	for (i = 1; i < MAX_TURN; i++) {
		if (route[i] == STP) { break; }
		//移動したら、移動前ドロップと移動後ドロップを交換する。
		int row = (int)YY(route[i]), col = (int)XX(route[i]);
		F_T c = field[prw][pcl];
		field[prw][pcl] = field[row][col];
		field[row][col] = c;
		prw = row, pcl = col;
	}
}
//double d_rnd(double mini, double maxi) {//xorshift実数乱数、おまじない
unsigned int rnd(int mini, int maxi) {//xorshift整数乱数、おまじない
	random_device rd;

	mt19937 mt(rd());

	uniform_int_distribution<int> dice(mini, maxi);
	return dice(mt);
	//return (((double)w / ((double)UINT_MAX + 1)) * (maxi - mini)) + mini;
}
int main() {

	int i, j;
	double avg = 0;//平均コンボ数
	clock_t start, end;
	double t_sum = 0;
	for (i = 0; i < 1000; i++) {//1000問解く
		F_T f_field[ROW][COL]; //スワイプ前の盤面
		F_T field[ROW][COL]; //盤面
		printf("problem:%d\n", i + 1);
		init(f_field); set(f_field, 0);//初期盤面生成
		show_field(f_field);//盤面表示
		start = clock();
		Action tmp = BEAM_SEARCH(f_field);//ビームサーチしてtmpに最善手を保存
		end = clock();
		printf("%fSec\n", DLT(start, end));
		t_sum += DLT(start, end);
		printf("(x,y)=(%d,%d)", XX(tmp.moving[0]), YY(tmp.moving[0]));
		for (j = 1; j < MAX_TURN; j++) {//y座標は下にいくほど大きくなる
			if (tmp.moving[j] == STP) { break; }
			if (XX(tmp.moving[j]) == XX(tmp.moving[j - 1]) + 1) { printf("R"); } //"RIGHT"); }
			if (XX(tmp.moving[j]) == XX(tmp.moving[j - 1]) - 1) { printf("L"); } //"LEFT"); }
			if (YY(tmp.moving[j]) == YY(tmp.moving[j - 1]) + 1) { printf("D"); } //"DOWN"); }
			if (YY(tmp.moving[j]) == YY(tmp.moving[j - 1]) - 1) { printf("U"); } //"UP"); }
		} printf("\n");
		memcpy(field, f_field, sizeof(f_field));
		operation(field, tmp.moving);
		int combo = sum_e(field);
		printf("%dCombo\n", combo);
		avg += (double)combo;
	}
	printf("TotalDuration:%fSec\n", t_sum);
	printf("Avg.Combo #:%lf\n", avg / (double)i);
	printf("p1:%f,p2:%f,p3:%f,p4:%f\n", part1, part2, part3, part4);
	//j = getchar();
	return 0;
}