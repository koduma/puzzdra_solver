/*
puzzdra_solver
パズドラのルート解析プログラムです。
なるべく少ない時間でなるべく大きいコンボを出したいです。
printf("TotalDuration:%fSec\n", t_sum);
printf("Avg.Combo #:%lf/%lf\n", avg / (double)i,MAXCOMBO/(double)i);
これらが改善されればpull request受け付けます。

914769
264812
379934
355886
951279

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
#define DROP 6//ドロップの種類//MAX9
#define TRN  44//手数//MAX305
#define STP YX(7,7)//無効手[無効座標]
#define MAX_TURN 40//最大ルート長//MAX300
#define BEAM_WIDTH 42000//ビーム幅//MAX500000
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

int evaluate2(F_T field[ROW][COL], int flag,int *combo);//落とし減点評価関数
int sum_e2(F_T field[ROW][COL],int *combo);//評価関数

struct member {//どういう手かの構造体
	T_T movei[TRN];//スワイプ移動座標
	int score;//評価値
	int combo;//コンボ数
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
double part1 = 0, part2 = 0, part3 = 0, part4 = 0, MAXCOMBO = 0;
Action BEAM_SEARCH(F_T f_field[ROW][COL]) {

	int stop = 0;//理論最大コンボ数

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
	}
	MAXCOMBO += (double)stop;

	deque<member>dque;
	double start, st;
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
		start = omp_get_wtime();
#pragma omp parallel for private(st),reduction(+:part1,part4)
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
						st = omp_get_wtime();
						F_T field[ROW][COL]; //盤面
						memcpy(field, f_field, sizeof(field));//盤面をもどす
						operation(field, cand.movei);
						int cmb;
						cand.score = sum_e2(field,&cmb);
						cand.combo = cmb;
						part1 += omp_get_wtime() - st;
						cand.prev = j;
						st = omp_get_wtime();
						//#pragma omp critical
											//{ pque.push(cand); }
						fff[(4 * k) + j] = cand;
						part4 += omp_get_wtime() - st;
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
		part2 += omp_get_wtime() - start;
		start = omp_get_wtime();
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
			if (maxValue < temp.combo) {//コンボ数が増えたらその手を記憶する
				maxValue = temp.combo;
				bestAction.score = maxValue;
				memcpy(bestAction.moving, temp.movei, sizeof(temp.movei));
				//コンボ数が理論値になったらreturn
				if (temp.combo == stop) { return bestAction; }
			}
			if (i < MAX_TURN - 1) {
				dque.push_back(temp);
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
				field[i][j] = (F_T)rnd(force ? 0 : 1, DROP);//1-DROPの整数乱数
			}
		}
	}
}
int chain(int nrw, int ncl, int d, F_T field[ROW][COL],
	F_T chkflag[ROW][COL], F_T delflag[ROW][COL]) {
	int count = 0;
#define CHK_CF(Y,X) (field[Y][X] == d && chkflag[Y][X] == 0 && delflag[Y][X] > 0)
	//連結している同じ色の消去ドロップが未探索だったら
	if (CHK_CF(nrw, ncl)) {
		++count; //連結ドロップ数の更新
		chkflag[nrw][ncl] = 1; //探索済みにする
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
				//コンボになったドロップは空になる。
				if (delflag[row][col] > 0) { field[row][col] = 0; }
			}
		}

		if (flag & EVAL_FALL)fall(field);//落下処理発生
		if (flag & EVAL_SET)set(field, 0);//落ちコン発生

	}
	return combo;
}

int evaluate2(F_T field[ROW][COL], int flag,int *combo) {
	int ev = 0;
	*combo = 0;
	while (1) {
		int cmb = 0;
		int cmb2 = 0;
		F_T chkflag[ROW][COL] = { 0 };
		F_T delflag[ROW][COL] = { 0 };
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
		char cnt[DROP + 1] = { 0 };
		char drop[DROP + 1][ROW * COL][2] = { 0 };

		for (int row = 0; row < ROW; row++) {
			for (int col = 0; col < COL; col++) {
				drop[field[row][col]][cnt[field[row][col]]][0] = (char)row;
				drop[field[row][col]][cnt[field[row][col]]][1] = (char)col;
				cnt[field[row][col]]++;
				if (delflag[row][col] > 0) {
					int c = chain(row, col, field[row][col], field, chkflag, delflag);
					if (c >= 3) {
						cmb++;
						if (c == 3) { cmb2 += 30; }
						else { cmb2 += 20; }
					}
				}
			}
		}
		for (int i = 1; i <= DROP; i++) {
			for (int j = 0; j < cnt[i] - 1; j++) {
				char add = max(drop[i][j + 1][0] - drop[i][j][0], drop[i][j][0] - drop[i][j + 1][0]) + max(drop[i][j + 1][1] - drop[i][j][1], drop[i][j][1] - drop[i][j + 1][1]);
				cmb2 -= (int)add;
			}
		}
		*combo += cmb;
		ev += cmb2;
		//コンボが発生しなかったら終了
		if (cmb == 0 || 0 == (flag & EVAL_COMBO)) { break; }
		for (int row = 0; row < ROW; row++) {
			for (int col = 0; col < COL; col++) {
				//コンボになったドロップは空になる。
				if (delflag[row][col] > 0) { field[row][col] = 0; }
			}
		}

		if (flag & EVAL_FALL)fall(field);//落下処理発生
		if (flag & EVAL_SET)set(field, 0);//落ちコン発生

	}
	return ev;
}
int sum_e2(F_T field[ROW][COL],int* combo) {//落とし有り、落ちコン無し評価関数
	return evaluate2(field, EVAL_FALL | EVAL_COMBO, combo);
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
}
int main() {

	int i, j, k;

	double avg = 0;//平均コンボ数
	double start;
	double t_sum = 0;
	for (i = 0; i < 1000; i++) {//1000問解く
		F_T f_field[ROW][COL]; //スワイプ前の盤面
		F_T field[ROW][COL]; //盤面
		printf("input:No.%d\n", i + 1);
		init(f_field); set(f_field, 0);//初期盤面生成
		show_field(f_field);//盤面表示
		/*
		for (j = 0; j < ROW; j++) {
			string s = "";
			cin >> s;
			for (k = 0; k < COL; k++) {
				f_field[j][k] = s[k] - '0';
			}
		}
		*/
		start = omp_get_wtime();
		Action tmp = BEAM_SEARCH(f_field);//ビームサーチしてtmpに最善手を保存
		double diff = omp_get_wtime() - start;
		printf("%fSec\n", diff);
		t_sum += diff;
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
		printf("output:No.%d\n",i+1);
		show_field(field);
		int combo = sum_e(field);
		printf("%dCombo\n", combo);
		avg += (double)combo;
	}
	printf("TotalDuration:%fSec\n", t_sum);
	printf("Avg.Combo #:%lf/%lf\n", avg / (double)i, MAXCOMBO / (double)i);
	printf("p1:%f,p2:%f,p3:%f,p4:%f\n", part1, part2, part3, part4);
	j = getchar();
	return 0;
}
