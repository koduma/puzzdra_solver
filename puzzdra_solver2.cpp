/*
puzzdra_solver
パズドラのルート解析プログラムです。
なるべく少ない時間でなるべく大きいコンボを出したいです。
printf("TotalDuration:%fSec\n", t_sum);
printf("Avg.Combo #:%lf/%lf\n", avg / (double)i,MAXCOMBO/(double)i);
これらが改善されればpull request受け付けます。
*/

#pragma warning(disable:4710)
#pragma warning(disable:4711)
#pragma warning(disable:4820)

#include <windows.h>
#include <fstream>
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
#include <functional>

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

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
#define DROP 6//ドロップの種類
#define TRN  44//手数
#define STP YX(7,7)//無効手[無効座標]
#define MAX_TURN 40//最大ルート長
#define BEAM_WIDTH 50000//ビーム幅
#define PICTURES 6//ドロップ画像の種類

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
double part1 = 0, part2 = 0, part3 = 0, part4 = 0, MAXCOMBO = 0;


class RowReso
{
private:
	cv::Mat* _org_img;
	cv::Mat* _reso_img;
	cv::Mat* _reso_org;

	int _reso;
	int _reso_width;
	int _reso_height;

public:
	RowReso()
	{
		_org_img = NULL;
		_reso_img = NULL;
		_reso_org = NULL;
	}
	~RowReso()
	{
		if (_reso_img != NULL) delete _reso_img;
		if (_reso_org != NULL) delete _reso_org;
	}

	// 初期化
	void Initialize(cv::Mat& img, int reso)
	{
		int width = img.cols / reso;
		int height = img.rows / reso;

		_org_img = &img;
		_reso_img = new cv::Mat(height, width, CV_MAKETYPE(img.depth(), img.channels()));
		_reso = reso;
		_reso_width = width;
		_reso_height = height;
	}
	// 低解像度を作成
	cv::Mat& Do()
	{
		for (int y = 0; y < _reso_height; ++y) {
			for (int x = 0; x < _reso_width; ++x) {
				int x1 = (_reso + 1) / 2 + _reso * x;
				int y1 = (_reso + 1) / 2 + _reso * y;
				cv::Vec3b& v = _org_img->at<cv::Vec3b>(y1, x1);
				// cout << x << "," << y << endl;
				_reso_img->at<cv::Vec3b>(y, x) = v;
			}
		}
		return *_reso_img;
	}
	// 確認用に元の画像の大きさに戻す
	cv::Mat& GetOriginalSize()
	{
		if (_reso_org == NULL) {
			_reso_org = new cv::Mat(
				_org_img->rows, _org_img->cols,
				CV_MAKETYPE(_org_img->depth(), _org_img->channels()));
		}
		for (int y = 0; y < _reso_height; ++y) {
			for (int x = 0; x < _reso_width; ++x) {
				cv::Vec3b& v = _reso_img->at<cv::Vec3b>(y, x);
				for (int y1 = 0; y1 < _reso; ++y1) {
					for (int x1 = 0; x1 < _reso; ++x1) {
						_reso_org->at<cv::Vec3b>(y * _reso + y1, x * _reso + x1) = v;
					}
				}
			}
		}
		return *_reso_org;
	}
};

int picture_analysis(int num, F_T field[ROW][COL])
{

	string str = "C:\\Users\\def\\puzzdra" + to_string(num) + ".png";

	// テンプレート画像
	cv::Mat img_koma[PICTURES];
	img_koma[0] = cv::imread("C:\\Users\\def\\Desktop\\gomi\\fire1.png", 1);
	img_koma[1] = cv::imread("C:\\Users\\def\\Desktop\\gomi\\water1.png", 1);
	img_koma[2] = cv::imread("C:\\Users\\def\\Desktop\\gomi\\wind1.png", 1);
	img_koma[3] = cv::imread("C:\\Users\\def\\Desktop\\gomi\\light1.png", 1);
	img_koma[4] = cv::imread("C:\\Users\\def\\Desktop\\gomi\\dark1.png", 1);
	img_koma[5] = cv::imread("C:\\Users\\def\\Desktop\\gomi\\cure1.png", 1);

	cv::Mat img = cv::imread(str, 1);//横(x)=img.cols,縦(y)=img.rows

	int reso = 3;
	RowReso Reso, ResoKoma[PICTURES];
	Reso.Initialize(img, reso);
	cv::Mat img_reso_komas[PICTURES];
	for (int i = 0; i < PICTURES; i++) {
		ResoKoma[i].Initialize(img_koma[i], reso);
		// 低解像度の教師画像
		img_reso_komas[i] = ResoKoma[i].Do();
	}

	// 枠線の色
	cv::Scalar cols[PICTURES];
	cols[0] = cv::Scalar(0, 0, 255);
	cols[1] = cv::Scalar(0, 255, 0);
	cols[2] = cv::Scalar(255, 0, 0);
	cols[3] = cv::Scalar(255, 0, 255);
	cols[4] = cv::Scalar(255, 255, 0);
	cols[5] = cv::Scalar(255, 255, 255);

	cv::Mat& img_reso = Reso.Do();
	cv::Mat& img_reso_org = Reso.GetOriginalSize();

	cv::Mat img_search, img_result;
	img_reso.copyTo(img_search);

	int cnt = 0;


	for (int j = 0; j < PICTURES; j++) {
		cv::Mat& img_reso_koma = img_reso_komas[j];
		// テンプレートマッチング
		cv::matchTemplate(img_search, img_reso_koma, img_result, CV_TM_CCOEFF_NORMED);
		// 50 個検出する
		for (int i = 0; i < 50; i++) {
			cv::Point max_pt;
			double maxVal;
			cv::minMaxLoc(img_result, NULL, &maxVal, NULL, &max_pt);
			// 一定スコア以下の場合は処理終了
			if (maxVal < 0.9) break;


			cv::Rect roi_rect(0, 0, img_reso_koma.cols, img_reso_koma.rows);
			roi_rect.x = max_pt.x;
			roi_rect.y = max_pt.y;
			cv::Rect roi_rect_org(roi_rect.x * reso, roi_rect.y * reso, img_reso_koma.cols * reso, img_reso_koma.rows * reso);

			cnt++;
			//std::cout << "(" << max_pt.x << ", " << max_pt.y << "), score=" << maxVal << std::endl;
			cv::rectangle(img_reso_org, roi_rect_org, cols[i], 3);

			for (int y = 0; y < ROW; y++) {
				for (int x = 0; x < COL; x++) {//79=(x*79/1440),372=(y*372/2560)
					if (79 * x <= max_pt.x && max_pt.x <= 79 * (x + 1) &&
						372 + (79 * y) <= max_pt.y && max_pt.y <= 372 + (79 * (y + 1))) {//(1274-(236/2))/3,236/3
						field[y][x] = j + 1;//同じのがh枚なら(j/h)+1
					}
				}

			}
			for (int y = 0; y < img_reso_koma.rows; y++) {
				for (int x = 0; x < img_reso_koma.cols; x++) {
					int xx = max_pt.x + x - img_reso_koma.cols / 2;
					int yy = max_pt.y + y - img_reso_koma.rows / 2;
					if (0 <= xx && xx < img_result.cols - 1) {
						if (0 <= yy && yy < img_result.rows - 1) {
							img_result.at<int>(yy, xx) = 0;
						}
					}
				}
			}
		}
	}

	std::cout << cnt << std::endl;

	//remove(str.c_str());

	return cnt;

}

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
void swipe_event(int num,int swipe_point[ROW][COL][2],Action ans) {

	int i;

	string add_x = "";
	string add_y = "";

	int stop = TRN;

	for (i = 0; i < TRN; i++) {
		if (ans.moving[i] == STP) {
			stop = i;
			break;
		}
	}

	for (i = 0; i < stop; i++) {
		int XXX = XX(ans.moving[i]);
		int YYY = YY(ans.moving[i]);
		add_x += to_string(swipe_point[YYY][XXX][0]);
		add_y += to_string(swipe_point[YYY][XXX][1]);
		if (i != stop - 1) {
			add_x += ",";
			add_y += ",";
		}
	}

	string str = "C:\\Users\\def\\Desktop\\gomi\\jaa" + to_string(num) + ".bat";

	string add2 = "adb shell am instrument -w -e x " + add_x + " -e y " + add_y + " -e s 6 -e class jp.android.uiautomator.autotest.CUI#swipe jp.android.uiautomator.autotest/.AutoTestRunner";

	string str2 = "@echo off\n" + add2 + "\n@echo on";

	ofstream ofs(str);

	ofs << str2 << endl;

	system(str.c_str());

}
void delete_android_screenshot(int num) {

	string str3 = "C:\\Users\\def\\Desktop\\gomi\\hoge" + to_string(num) + ".bat";
	string str4 = "@echo off\nadb shell rm -rf /sdcard/puzzdra" + to_string(num) + ".png\n@echo on";
	ofstream ofs(str3);

	ofs << str4 << endl;

	system(str3.c_str());
}
void delete_bat_file(int num) {
	string str = "C:\\Users\\def\\Desktop\\gomi\\ikuze" + to_string(num) + ".bat";
	string str3 = "C:\\Users\\def\\Desktop\\gomi\\hoge" + to_string(num) + ".bat";
	remove(str.c_str());
	remove(str3.c_str());
}
void get_screenshot(int num) {

	string str = "C:\\Users\\def\\Desktop\\gomi\\ikuze" + to_string(num) + ".bat";
	string str2 = "@echo off\nadb shell screencap -p /sdcard/puzzdra" + to_string(num) + ".png\nadb pull /sdcard/puzzdra" + to_string(num) + ".png\n@echo on";

	ofstream ofs(str);

	ofs << str2 << endl;

	system(str.c_str());

}
void delete_swipe_event(int num) {
	string str = "C:\\Users\\def\\Desktop\\gomi\\jaa" + to_string(num) + ".bat";
	remove(str.c_str());
}
void delete_picture(int num) {
	string str = "C:\\Users\\def\\puzzdra" + to_string(num) + ".png";
	remove(str.c_str());
}

Action solve(F_T field[ROW][COL]) {

	F_T f_field[ROW][COL];

	memcpy(f_field, field, sizeof(f_field));

	printf("input=\n");
	for (int row = 0; row < ROW; row++) {
		for (int col = 0; col < COL; col++) {
			printf("%d", field[row][col]);
		}
		printf("\n");
	}
	printf("\n");

	Action tmp = BEAM_SEARCH(field);

	operation(f_field,tmp.moving);

	printf("最低%dコンボ以上\n", sum_e(f_field));

	return tmp;
}

int main() {

	int swipe_point[ROW][COL][2];
	F_T field[ROW][COL] = {0};

	int j, k;

	for (int j = 0; j < ROW; j++) {
		for (int k = 0; k < COL; k++) {//126=(x*126/1440),236=(x*236/1440),1274=(y*1274/2560),xは横幅、yは縦幅
			swipe_point[j][k][0] = 126 + (236 * k);//x
			swipe_point[j][k][1] = 1274 + (236 * j);//y
		}
	}

	int i = 1;

	while (1) {
		Sleep(17000);
		get_screenshot(i);
		delete_android_screenshot(i);
		delete_bat_file(i);
		int cnt = picture_analysis(i, field);

		if (cnt == 30) {
			swipe_event(i, swipe_point, solve(field));
			delete_swipe_event(i);
			delete_picture(i);
		}
		else {
			break;
		}

		i++;
	}

	j = getchar();
	return 0;
}