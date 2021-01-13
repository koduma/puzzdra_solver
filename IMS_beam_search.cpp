/*
puzzdra_solver

�p�Y�h���̃��[�g��̓v���O�����ł�

�R���p�C����MinGW�𐄏����܂�

�Ȃ�ׂ����Ȃ����ԂłȂ�ׂ��傫���R���{���o�������ł�

printf("TotalDuration:%fSec\n", t_sum);
printf("Avg.NormalCombo #:%f/%f\n", avg / (double)i, MAXCOMBO / (double)i);

����炪���P������pull request�󂯕t���܂�

�`�F�b�N1�F�����10�R���{�ł��邩

962679
381515
489942
763852
917439

914769
264812
379934
355886
951279

�`�F�b�N2�F1000�Ֆʕ��ϗ����R���{����9.20�t�߂�

�`�F�b�N3�F1000�Ֆʕ��σR���{�������_�l�t�߂�

�S�`�F�b�N�B�������獇�i
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
#include <immintrin.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;
#define DLT(ST,ED) ((double)((ED)-(ST))/CLOCKS_PER_SEC)//���ԍ���
#define XX(PT) ((PT)&15)
#define YY(PT) XX((PT)>>4)
#define YX(Y,X) ((Y)<<4|(X))
#define DIR 4//����
#define ROW 5//�c//MAX6
#define COL 6//��//MAX7
#define DROP 8//�h���b�v�̎��//MAX9
#define TRN 55//�萔//MAX155
#define MAX_TURN 55//�ő僋�[�g��//MAX150
#define BEAM_WIDTH 2500000//�r�[����//MAX200000
#define PROBLEM 10000//��萔
#define BONUS 10//�]���l���P�W��
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define NODE_SIZE MAX(500,4*BEAM_WIDTH)
typedef char F_T;//�Ֆʌ^
typedef char T_T;//�萔�^
typedef unsigned long long ll;
enum { EVAL_NONE = 0, EVAL_FALL, EVAL_SET, EVAL_FS, EVAL_COMBO };
void init(F_T field[ROW][COL]); //�����z�u�����֐�
void fall(int x,int h,F_T field[ROW][COL]); //�h���b�v�̗��������֐�
void set(F_T field[ROW][COL], int force); //��}�X�𖄂߂�֐�
void show_field(F_T field[ROW][COL]); //�Ֆʕ\���֐�
int rnd(int mini, int maxi); //��������
//�㉺���E�ɘA�����Ă���h���b�v���ċA�I�ɒT�����Ă����֐�
int chain(int nrw, int ncl, F_T d, F_T field[ROW][COL], F_T chkflag[ROW][COL], F_T delflag[ROW][COL]);
int evaluate(F_T field[ROW][COL], int flag); //�R���{������֐�
int sum_e(F_T field[ROW][COL]);//���Ƃ��L��A�����R�������R���{������֐�
int sum_evaluate(F_T field[ROW][COL]);//���Ƃ��������R�����L��R���{������֐�
void operation(F_T field[ROW][COL], T_T first_te,ll route[(TRN/21)+1],ll dropBB[DROP+1]); //�X���C�v�����֐�

int evaluate2(F_T field[ROW][COL], int flag, int* combo, ll* hash);//���Ƃ����_�]���֐�
int sum_e2(F_T field[ROW][COL], int* combo, ll* hash);//�]���֐�

ll xor128();//xorshift��������
ll zoblish_field[ROW][COL][DROP+1];

ll sqBB[64];
int evaluate3(ll dropBB[DROP+1], int flag, int* combo, ll* hash);//���Ƃ����_�]���֐�
int sum_e3(ll dropBB[DROP+1], int* combo, ll* hash);//�]���֐�
ll around(ll bitboard);
int table[64];
ll fill_64[64];
ll file_bb[COL];
ll calc_mask(ll bitboard);
ll fallBB(ll p,ll rest,ll mask);

struct node {//�ǂ������肩�̍\����
	T_T first_te;
	ll movei[(TRN/21)+1];//�X���C�v�ړ����W
	int score;//�]���l
	int combo;//�R���{��
	int nowC;//���ǂ�x���W�ɂ��邩
	int nowR;//���ǂ�y���W�ɂ��邩
	int prev;//1��O�͏㉺���E�̂ǂ�����I�񂾂�
	int prev_score;//1��O�̕]���l
	int improving;//�]���l���P��
	ll hash;//�Ֆʂ̃n�b�V���l
	node() {//������
		this->score = 0;
		this->prev = -1;
		//memset(this->movei, STP, sizeof(this->movei));
	}
	bool operator < (const node& n)const {//�X�R�A�����������D�悳���
		return score < n.score;
	}
}fff[NODE_SIZE];
struct Action {//�ŏI�I�ɒT�����ꂽ��
	T_T first_te;
	int score;//�R���{��
	int maxcombo;//���_�R���{��
	ll moving[(TRN/21)+1];//�X���C�v�ړ����W
	Action() {//������
		this->score = 0;
		//memset(this->moving, STP, sizeof(this->moving));
	}
};
Action BEAM_SEARCH(F_T f_field[ROW][COL],int maxi); //���[�g�T���֐�
double part1 = 0, part2 = 0, part3 = 0, part4 = 0, MAXCOMBO = 0;
Action BEAM_SEARCH(F_T f_field[ROW][COL],int maxi) {

	int po=9+(8*(COL-1))+ROW-1;

	int stop = 0;//���_�ő�R���{��

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

	deque<node>dque;
	double start, st;
	//1��ڂ�S�ʂ�T������
	dque.clear();
	for (int i = 0; i < ROW; i++) {
		for (int j = 0; j < COL; j++) {
			node cand;
			cand.nowR = i;//y���W
			cand.nowC = j;//x���W
			cand.prev = -1;//1��O�͂Ȃ�
			cand.first_te = (T_T)YX(i, j);//1��ڂ�yx���W
			for (int trn = 0; trn <= TRN/21; trn++) {
				cand.movei[trn] = 0ll;
			}
			F_T ff_field[ROW][COL];
			memcpy(ff_field,f_field,sizeof(ff_field));
			int cmb;
			ll ha;
			cand.prev_score=sum_e2(ff_field,&cmb,&ha);
			cand.improving=0;
			cand.hash=ha;
			dque.push_back(cand);
		}
	}// L, U,D,R //
	int dx[DIR] = { -1, 0,0,1 },
		dy[DIR] = { 0,-1,1,0 };
	Action bestAction;//�őP��
	int maxValue = 0;//�ō��X�R�A

	bestAction.maxcombo = stop;
	unordered_map<ll, bool> checkNodeList[ROW*COL];

	ll rootBB[DROP+1]={0};

	for(int row=0;row<ROW;row++){
	for(int col=0;col<COL;col++){
	int pos=po-((8*col)+row);
	rootBB[f_field[row][col]]|=(1ll << (pos));
	}
	}

	//2��ڈȍ~���r�[���T�[�`�ŒT��
	for (int i = 0; i < MAX_TURN; i++) {
		printf("depth=%d/%d\n",i+1,MAX_TURN);
		int ks = (int)dque.size();
		start = omp_get_wtime();
#pragma omp parallel for private(st),reduction(+:part1,part4)
		for (int k = 0; k < ks; k++) {
#ifdef _OPENMP
			if (i == 0 && k == 0) {
				printf("Threads[%d/%d]\n",
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
			for (int j = 0; j < DIR; j++) {//�㉺���E��4����������
				node cand = temp;
				if (0 <= cand.nowC + dx[j] && cand.nowC + dx[j] < COL &&
					0 <= cand.nowR + dy[j] && cand.nowR + dy[j] < ROW) {
					if (cand.prev + j != 3) {
						int ny=cand.nowR + dy[j];
						int nx=cand.nowC + dx[j];
						F_T field[ROW][COL];//�Ֆ�
						ll dropBB[DROP+1]={0};
						memcpy(field,temp_field,sizeof(temp_field));//�Ֆʂ����ǂ�
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
						st = omp_get_wtime();
						int cmb;
						ll ha;
						cand.score = sum_e3(dropBB, &cmb,&ha);
						cand.combo = cmb;
						cand.hash=ha;
						part1 += omp_get_wtime() - st;
						cand.prev = j;
						st = omp_get_wtime();
						//#pragma omp critical
											//{ pque.push(cand); }
						fff[(4 * k) + j] = cand;
						part4 += omp_get_wtime() - st;
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
		part2 += omp_get_wtime() - start;
		start = omp_get_wtime();
		dque.clear();
		vector<pair<int, int> >vec;
		int ks2 = 0;
		for (int j = 0; j < 4 * ks; j++) {
			if (fff[j].combo != -1) {
			if (fff[j].combo == stop) {
				maxValue = stop;
				bestAction.score = maxValue;
				bestAction.first_te = fff[j].first_te;
				memcpy(bestAction.moving, fff[j].movei, sizeof(fff[j].movei));
				return bestAction;
			}
				if(fff[j].score>fff[j].prev_score){fff[j].improving=fff[j].improving+1;}
				fff[j].prev_score=fff[j].score;
				vec.push_back(make_pair(fff[j].score+(BONUS*fff[j].improving)+rnd(0,maxi), j));
				ks2++;
			}
		}
		sort(vec.begin(), vec.end());
		int push_node=0;
		for (int j = 0; push_node < BEAM_WIDTH && j < ks2; j++) {
			node temp = fff[vec[ks2-1-j].second];
			if (maxValue < temp.combo) {//�R���{�����������炻�̎���L������
				maxValue = temp.combo;
				bestAction.score = maxValue;
				bestAction.first_te = temp.first_te;
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
			if (field[i][j] == 0 || force) {//��}�X�������炤�߂�
				field[i][j] = (F_T)rnd(force ? 0 : 1, DROP);//1-DROP�̐�������
			}
		}
	}
}
int chain(int nrw, int ncl, F_T d, F_T field[ROW][COL],
	F_T chkflag[ROW][COL], F_T delflag[ROW][COL]) {
	int count = 0;
#define CHK_CF(Y,X) (field[Y][X] == d && chkflag[Y][X]==0 && delflag[Y][X] > 0)
	//�A�����Ă��铯���F�̏����h���b�v�����T����������
	if (CHK_CF(nrw, ncl)) {
		++count; //�A���h���b�v���̍X�V
		chkflag[nrw][ncl]=1;//�T���ς݂ɂ���
			//�ȉ��㉺���E�ɘA�����Ă���h���b�v���ċA�I�ɒT�����Ă���
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
		//�R���{���������Ȃ�������I��
		if (cmb == 0 || 0 == (flag & EVAL_COMBO)) { break; }
		for (int row = 0; row < ROW; row++) {
			for (int col = 0; col < COL; col++) {
				//�R���{�ɂȂ����h���b�v�͋�ɂȂ�
				if (delflag[row][col]> 0) { field[row][col] = 0; }
			}
		}

		if (flag & EVAL_FALL){
		for(int x=0;x<COL;x++){
		fall(x,GetHeight[x],field);
		}
		}//������������
		if (flag & EVAL_SET){set(field, 0);}//�����R������

	}
	return combo;
}
int evaluate2(F_T field[ROW][COL], int flag, int* combo, ll* hash) {
	int ev = 0;
	*combo = 0;
	ll ha=0;
	int oti = 0;
	while (1) {
		int cmb = 0;
		int cmb2 = 0;
		F_T chkflag[ROW][COL]={0};
		F_T delflag[ROW][COL]={0};
		F_T GetHeight[COL];
		for (int row = 0; row < ROW; row++) {
			for (int col = 0; col < COL; col++) {
				F_T num = field[row][col];
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

		F_T cnt[DROP + 1] = { 0 };
		F_T drop[DROP + 1][ROW * COL][2] = { 0 };

		for (int row = 0; row < ROW; row++) {
			for (int col = 0; col < COL; col++) {
				drop[field[row][col]][cnt[field[row][col]]][0] = (F_T)row;
				drop[field[row][col]][cnt[field[row][col]]][1] = (F_T)col;
				cnt[field[row][col]]++;
				if (delflag[row][col]>0) {
					int c = chain(row, col, field[row][col], field, chkflag, delflag);
					if (c >= 3) {
						cmb++;
						if (c == 3) { cmb2 += 30; }
						else { cmb2 += 20; }
					}
				}
			}
		}
		F_T erase_x[COL]={0};
		for (int i = 1; i <= DROP; i++) {
			for (int j = 0; j < cnt[i] - 1; j++) {
				int d1 = (int)drop[i][j][0];
				int d2 = (int)drop[i][j][1];
				int d3 = (int)drop[i][j + 1][0];
				int d4 = (int)drop[i][j + 1][1];
				int add = max(d1 - d3, d3 - d1) + max(d2 - d4, d4 - d2);
				add += add;
				add /= 3;
				cmb2 -= add;
				if (delflag[d1][d2]> 0) {
					field[d1][d2] = 0;
					erase_x[d2]=1;
				}
				if (delflag[d3][d4] > 0) {
					field[d3][d4] = 0;
					erase_x[d4]=1;
				}
			}
		}
		*combo += cmb;
		ev += cmb2;
		//�R���{���������Ȃ�������I��
		if (cmb == 0 || 0 == (flag & EVAL_COMBO)) { break; }
		oti++;
		if (flag & EVAL_FALL){//������������
		for(int x=0;x<COL;x++){
		if(erase_x[x]==1){
		fall(x,GetHeight[x],field);
		}
		}
		}
		if (flag & EVAL_SET){set(field, 0);}//�����R������

	}
	ev += oti;
	*hash=ha;
	return ev;
}
int evaluate3(ll dropBB[DROP+1], int flag, int* combo, ll* hash) {
	int ev = 0;
	*combo = 0;
	ll ha=0;
	int oti = 0;
	ll occBB=0;

	for(int i=1;i<=DROP;i++){
	occBB|=dropBB[i];
	}

	int po=9+(8*(COL-1))+ROW-1;

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
		else{cmb2+=20;}
		}
		}
		}


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
		if(oti==0){
		ha ^= zoblish_field[pos_y][pos_x][i];
		}
		if(j>0){
		int add=max(pos_x-prev_x,prev_x-pos_x)+max(pos_y-prev_y,prev_y-pos_y);
		add+=add;
		add/=3;
		cmb2-=add;
		}
		prev_x=pos_x;
		prev_y=pos_y;
		tmp_drop=tmp_drop & ~(1ll<<(pos));
		}
		dropBB[i]^=linked[i];
		occBB^=linked[i];
		}

		*combo += cmb;
		ev += cmb2;
		//�R���{���������Ȃ�������I��
		if (cmb == 0 || 0 == (flag & EVAL_COMBO)) { break; }
		oti++;

		ll mask=calc_mask(occBB);

		for(int i=1;i<=DROP;i++){
		dropBB[i]=fallBB(dropBB[i],occBB,mask);
		}
		occBB=fallBB(occBB,occBB,mask);

		/*
		F_T chkflag[ROW][COL]={0};
		F_T delflag[ROW][COL]={0};
		F_T GetHeight[COL];
		for (int row = 0; row < ROW; row++) {
			for (int col = 0; col < COL; col++) {
				F_T num = field[row][col];
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
		F_T cnt[DROP + 1] = { 0 };
		F_T drop[DROP + 1][ROW * COL][2] = { 0 };
		for (int row = 0; row < ROW; row++) {
			for (int col = 0; col < COL; col++) {
				drop[field[row][col]][cnt[field[row][col]]][0] = (F_T)row;
				drop[field[row][col]][cnt[field[row][col]]][1] = (F_T)col;
				cnt[field[row][col]]++;
				if (delflag[row][col]>0) {
					int c = chain(row, col, field[row][col], field, chkflag, delflag);
					if (c >= 3) {
						cmb++;
						if (c == 3) { cmb2 += 30; }
						else { cmb2 += 20; }
					}
				}
			}
		}
		F_T erase_x[COL]={0};
		for (int i = 1; i <= DROP; i++) {
			for (int j = 0; j < cnt[i] - 1; j++) {
				int d1 = (int)drop[i][j][0];
				int d2 = (int)drop[i][j][1];
				int d3 = (int)drop[i][j + 1][0];
				int d4 = (int)drop[i][j + 1][1];
				int add = max(d1 - d3, d3 - d1) + max(d2 - d4, d4 - d2);
				add += add;
				add /= 3;
				cmb2 -= add;
				if (delflag[d1][d2]> 0) {
					field[d1][d2] = 0;
					erase_x[d2]=1;
				}
				if (delflag[d3][d4] > 0) {
					field[d3][d4] = 0;
					erase_x[d4]=1;
				}
			}
		}
		*combo += cmb;
		ev += cmb2;
		//�R���{���������Ȃ�������I��
		if (cmb == 0 || 0 == (flag & EVAL_COMBO)) { break; }
		oti++;
		if (flag & EVAL_FALL){//������������
		for(int x=0;x<COL;x++){
		if(erase_x[x]==1){
		fall(x,GetHeight[x],field);
		}
		}
		}
		if (flag & EVAL_SET){set(field, 0);}//�����R������
*/

	}
	ev += oti;
	*hash=ha;
	return ev;
}
int sum_e3(ll dropBB[DROP+1], int* combo, ll* hash) {//���Ƃ��L��A�����R�������]���֐�
	return evaluate3(dropBB, EVAL_FALL | EVAL_COMBO, combo,hash);
}
int sum_e2(F_T field[ROW][COL], int* combo, ll* hash) {//���Ƃ��L��A�����R�������]���֐�
	return evaluate2(field, EVAL_FALL | EVAL_COMBO, combo,hash);
}
int sum_e(F_T field[ROW][COL]) {//���Ƃ��L��A�����R�������R���{������֐�
	return evaluate(field, EVAL_FALL | EVAL_COMBO);
}
int sum_evaluate(F_T field[ROW][COL]) {//���Ƃ��������R�����L��R���{������֐�
	return evaluate(field, EVAL_FS | EVAL_COMBO);
}
//�ړ�������̔Ֆʂ𐶐�����֐�
void operation(F_T field[ROW][COL], T_T first_te,ll route[(TRN/21)+1],ll dropBB[DROP+1]) {
	int prw = (int)YY(first_te), pcl = (int)XX(first_te), i,j;
	int dx[DIR] = { -1, 0,0,1 };
	int dy[DIR] = { 0,-1,1,0 };
	int po=9+(8*(COL-1))+ROW-1;
	for (i = 0; i <= TRN/21; i++) {
		if (route[i] == 0ll) { break; }
		//�ړ�������A�ړ��O�h���b�v�ƈړ���h���b�v����������
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
ll xor128() {//xorshift��������
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

//layout=053241405407470557104053134522
//No1:43��

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

	double avg = 0;//���σR���{��
	double start;
	double t_sum = 0;
	double oti_avg = 0;//���ϗ����R���{��
	for (i = 0; i < PROBLEM; i++) {//PROBLEM�����
		F_T f_field[ROW][COL]; //�X���C�v�O�̔Ֆ�
		F_T field[ROW][COL]; //�Ֆ�
		F_T oti_field[ROW][COL];//�����R���p�Ֆ�
		printf("input:No.%d/%d\n", i + 1, PROBLEM);
		string date="";
		printf("date=");
		cin>>date;
		//init(f_field); set(f_field, 0);//�����Ֆʐ���
		string str="";
		cin>>str;
		for (j = 0; j < ROW; j++) {
			for (k = 0; k < COL; k++) {
				f_field[j][k] = (str[k+(COL*j)] - '0')+1;
			}
		}
		show_field(f_field);//�Ֆʕ\��
		Action tmp,tmp2;
		double diff;
		int tesuu_min=10000;
		int tesuu;
		string ans;
		for(int loop=0;loop<10000;loop++){
		printf("loop=%d/%d\n",loop+1,10000);
		printf("tesuu_min=%d\n",tesuu_min);
		start = omp_get_wtime();
		tmp = BEAM_SEARCH(f_field,loop);//�r�[���T�[�`����tmp�ɍőP���ۑ�
		diff = omp_get_wtime() - start;
		t_sum += diff;
		//printf("(x,y)=(%d,%d)", XX(tmp.moving[0]), YY(tmp.moving[0]));
		tesuu=0;
		if(tmp.score==tmp.maxcombo){
		string route="";
		route+=to_string(XX(tmp.first_te))+to_string(YY(tmp.first_te)+5)+",";
		for (j = 0; j <= TRN/21; j++) {//y���W�͉��ɂ����قǑ傫���Ȃ�
			if (tmp.moving[j] == 0ll) { break; }
			for(k=0;k<21;k++){
			int dir = (int)(7ll&(tmp.moving[j]>>(3*k)));
			if (dir==0){break;}
			if (dir==1) { route+=to_string(3); } //"LEFT"); }
			if (dir==2) { route+=to_string(6); } //"UP"); }
			if (dir==3) { route+=to_string(1); } //"DOWN"); }
			if (dir==4) { route+=to_string(4); } //"RIGHT"); }
			tesuu++;
			}
		}
		if(tesuu_min>tesuu){
		tesuu_min=tesuu;
		tmp2=tmp;
		string url="http://serizawa.web5.jp/puzzdra_theory_maker/index.html?layout="+str+"&route="+route+"&date="+date+"&ctwMode=false";
		ans=url;
		}
		}//if(tmp
		cout<<ans<<endl;
		}//loop
	}//i
	printf("TotalDuration:%fSec\n", t_sum);
	printf("Avg.NormalCombo #:%f/%f\n", avg / (double)i, MAXCOMBO / (double)i);
	printf("Avg.OtiCombo #:%f\n", oti_avg / (double)i);
	printf("p1:%f,p2:%f,p3:%f,p4:%f\n", part1, part2, part3, part4);
	j = getchar();
	return 0;
}