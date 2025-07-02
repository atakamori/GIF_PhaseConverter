//Changelog

//2025/7/2
//gitにコミットするためのダミー。コードに変更は加えていない。

//2025/6/27
//dph計算のバグが取れた。正秒で値が狂っていたのを修正。実際に位相（fringe[])の値を決定してから差をとる方式に変えた。

//2025/6/26
// フィットの良否判定からdphを外す。
// 新谷氏のコメントを受けての変更だが、その部分のコードは残しておいても良いかもしれない。
// いずれにせよ、dph（隣接する位相の差）の計算にバグがあるので、デバグが必要。

//2025/6/3
//メモリリークの問題は未解決だが、2025/5/29 15:29の地震で一部dphに飛びが見られたので、フィットの良否判定にdphも加える試行を行う。


//2025/5/28, 29, 30
//長時間のデータ処理をするとメモリ使用量が多くなり、プロセスがkillされてしまう問題への対応を試みる。
//callocでのメモリ確保をcallocに変更-->多少改善?とはいえ本質的にはダメ
//短い配列（MM個）はcallocではなく、普通の配列として定義-->効果なし
//**************メモリリークの問題は未解決******************
//アイデア：細かくcalloc, freeするより、大きめのメモリを最初に確保しておいて、プログラム終了時に開放する方が良い?


//2025/5/27
//1分ファイル単位で、最初の区間では補間を行わないようにする（そもそもできない）
//2023/5/4 6:14開始で処理しようとすると、先頭でSegmentation faultが出ることから発覚
//将来的には外挿で処理するか?前の分の楕円パラメータも参照して補間?

//2025/5/21
//zabsの許容範囲を変更してテスト：
//th_zabsmin=0.3 -> 0.5, th_zabsmax=1.5 -> 1.3
// -> 条件を満たせない場合が見つかったので、元に戻す。
//th_zabsmin=0.5 -> 0.3, th_zabsmax=1.3 -> 1.5


//2025/5/9
// -i オプションの追加：指定すると、補間を行わない。

//2025/5/8
//バグ修正：
//関数fringe_detでzabsmaxが正しく返されていなかった（2025/5/7にコードを変えたときに、コピペミスをしていたため））
//1sファイルの書き出し機能追加：
//とりあえず、コメントアウトしていた箇所を有効化するなどしてファイルを書き出すようにした。問題無さそうだが、バグがないか検証は必要

//2025/5/7
//バグ修正：
//補間時のフラグ書き込み, 正秒でのll引き継ぎ（50 kHzの方。こちらはまだバグが残っているかも）
//フィッティングできなかったとき、強制的にzabs=0として補間を発生させるようにした（2023/5/5 12:05にこのようなケースがあった）

//2025/5/5
//20 Hz平均値計算のバグ修正（正秒でllが引き継がれず、2Piのずれが生じるケースがあった）

//2025/5/1
//フィット長を変更する前に、各段で必ず補間を行うことにする（テスト）
//位相変化が比較的ゆっくりで円弧が短いため楕円がフィットできないケースに対応できるか、テストするのが狙い。

//2025/4/30
//50 kHz出力をバイナリファイルで書き出すオプション（-H）追加

//2025/4/11 フィット楕円軸の傾き変化を記録
//補間に入る条件の候補とする。

//2025/4/7 補間処理完成
//補間に入る条件としてはzabsの範囲のみ参照。次期バージョンで楕円率（長短径の比）も条件として課してみる予定。



#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>
#include<time.h>
#include<unistd.h>

#define CC 299792458.e0				/* speed of light */
#define PI 3.1415926535897932384626
#define ADCFAC 5.525E-9

//default setting
#define SAMP 50000			// sampling frequency
#define fDS 20				// downsampling rate (Hz)
#define MM_ini 5000				// step of ellipse fitting
#define N 1				// data reduction number (50kHz/50kHz=1)
#define NN (60*SAMP/10) // number of data to process (1 min data=60*SAMP)

#define DTH 0.9			// threshold for discontinuity judgegment
#define ZTH 0.3			// threshold for zabs valuation (offset from 1.0)

#define test 0
#define debug1 0
#define debug2 0
#define singlestep 0
#define aut_range 0
#define moving_avg 0
#define interpolate 1

//long id0[MM],id1[MM],id2[MM],id3[MM];			/* fringe raw data for 50ms */
//int id0[MM_ini],id1[MM_ini],id2[MM_ini],id3[MM_ini];			/* ADC raw data for 50ms */
//int id0[60*SAMP],id1[60*SAMP],id2[60*SAMP],id3[60*SAMP];			/* ADC raw data for 1min */

//double fringe[MM_ini],zabs0[MM_ini]; // fringe and zabs for 50 ms
//double fringe[60*SAMP],zabs0[60*SAMP]; // fringe and zabs for 1min

//double px[NF];
struct binary_ph {
	float x; //位相絶対値から2piの整数倍を引いた残り
	long l; //位相絶対値を2piで割った時の整数
//} bph[M_ini];
};

struct ellip_param{
	int MM;
	double a00,b00,x0,y0,sp0,cp0;
}; // 楕円パラメータ構造体（位相決定に用いる）

struct fit_eval_param{
	double zabsmin, zabsmax, dfringemin, dfringemax, flatmin, flatmax, absdaa;
}; //eval_param[60*SAMP];

struct fringe_zabs{
	double fringe, zabs;
};// fr_zabs2[60*SAMP];

struct fringe_zabs_out{
	int numpnts;
	double fringe, zabs;
};

struct fit_param{
	int MM;
	double pin, amd, x0, y0, a00, b00, aa; //パラメータを追加すること（2024/4/15）
};

struct zabs_dph{
	double zabs, dph;
};



// 楕円フィット関数（その1）を定義
// 機能： ADCのデータ配列MM個を渡されると、フィット長MM個でフィットして、フィッティング結果を出力する。
// 入力: ADCデータ配列（iad0, iad1, iad2, iad3）, 配列長=フィット長(MM)（位相値はMM個の点について計算するが、出力は1点のみ）
// 出力： フィットパラメータの構造体（フィットに使った点数MM, 入射パワー値pin, 吸収線信号amd, 楕円原点座標x0,y0, 楕円半径a00,b00, 楕円軸傾斜角sp, cp） 

void elliptic_fit(int id0[], int id1[], int id2[], int id3[], int MM, struct ellip_param *elpara){

	int j,k,kk,ndph;
	long n,n1;
	//double pin,pinmax,pinmin,amd,amdmax,amdmin,fringemax,fringemin,fringe0;
	double x0max,x0min,y0max,y0min,a00max,a00min,b00max,b00min,aamax,aamin;
	//double xxmax,xxmin,xx1;
	double pp,ss,pp0,ss0,aa,xx,xx0,s0,s00;
	double x4,y4,x1y3,x2y2,x3y1,x3,y3,x1y2,x2y1,x2,y2,x1y1,x1,y1;
	double dsin,dsn,dcos,da,dab1,dab2,dab3,dab4,dab5,b1,b2,b3,b4,b5;
	double a00,b00,x0,y0,sp0,cp0;
	
	//double fringe[MM], zabs0[MM];
	
//MM個のデータを使っての楕円フィット開始（元プログラムを流用）
	/*	楕円フィット */
	x4=0.0;
	y4=0.0;
	x1y3=0.0;
	x2y2=0.0;
	x3y1=0.0;
	x3=0.0;
	y3=0.0;
	x1y2=0.0;
	x2y1=0.0;
	x2=0.0;
	y2=0.0;
	x1y1=0.0;
	x1=0.0;
	y1=0.0;
	pp0=0.0;
	ss0=0.0;
		
	for(j=0;j<MM;j+=N){ 
		pp0+=(double)id0[j]*ADCFAC/((double)id2[j]*ADCFAC);
		ss0+=(double)id1[j]*ADCFAC/((double)id2[j]*ADCFAC);
	}
	pp0/=(double)(MM/N); //(AD00/AD02) mean 桁落ち誤差を防ぐために平均値との差分で楕円フィットする
	ss0/=(double)(MM/N); //(AD01/AD02) mean
	for(j=0;j<MM;j+=N){
		pp=(double)id0[j]*ADCFAC/((double)id2[j]*ADCFAC)-pp0;
		ss=(double)id1[j]*ADCFAC/((double)id2[j]*ADCFAC)-ss0;
		x4+=pp*pp*pp*pp;
		y4+=ss*ss*ss*ss;
		x1y3+=pp*ss*ss*ss;
		x2y2+=pp*pp*ss*ss;
		x3y1+=pp*pp*pp*ss;
		x3+=pp*pp*pp;
		y3+=ss*ss*ss;
		x1y2+=pp*ss*ss;
		x2y1+=pp*pp*ss;
		x2+=pp*pp;
		y2+=ss*ss;
		x1y1+=pp*ss;
		x1+=pp;
		y1+=ss;
	}

	da= x1y1*x1y1*x1y3*x1y3*x4 - 2*x1y1*x1y1*x1y3*x2y2*x3y1 + x1y1*x1y1*x2y2*x2y2*x2y2 - 
	x1y1*x1y1*x2y2*x4*y4 + x1y1*x1y1*x3y1*x3y1*y4 - 2*x1y1*x1y2*x1y2*x1y3*x4 + 
	2*x1y1*x1y2*x1y2*x2y2*x3y1 + 2*x1y1*x1y2*x1y3*x2y1*x3y1 + 
	2*x1y1*x1y2*x1y3*x2y2*x3 - 4*x1y1*x1y2*x2y1*x2y2*x2y2 + 
	2*x1y1*x1y2*x2y1*x4*y4 + 2*x1y1*x1y2*x2y2*x4*y3 - 
	2*x1y1*x1y2*x3*x3y1*y4 - 2*x1y1*x1y2*x3y1*x3y1*y3 - 
	2*x1y1*x1y3*x1y3*x2y1*x3 + 2*x1y1*x1y3*x2y1*x2y1*x2y2 - 
	2*x1y1*x1y3*x2y1*x4*y3 + 2*x1y1*x1y3*x3*x3y1*y3 - 
	2*x1y1*x2y1*x2y1*x3y1*y4 + 2*x1y1*x2y1*x2y2*x3*y4 + 
	2*x1y1*x2y1*x2y2*x3y1*y3 - 2*x1y1*x2y2*x2y2*x3*y3 + x1y2*x1y2*x1y2*x1y2*x4 - 
	2*x1y2*x1y2*x1y2*x2y1*x3y1 - 2*x1y2*x1y2*x1y2*x2y2*x3 + 2*x1y2*x1y2*x1y3*x2y1*x3 + 
	x1y2*x1y2*x2*x2y2*x2y2 - x1y2*x1y2*x2*x4*y4 + 3*x1y2*x1y2*x2y1*x2y1*x2y2 - 
	2*x1y2*x1y2*x2y1*x4*y3 - x1y2*x1y2*x2y2*x4*y2 + x1y2*x1y2*x3*x3*y4 + 
	2*x1y2*x1y2*x3*x3y1*y3 + x1y2*x1y2*x3y1*x3y1*y2 - 2*x1y2*x1y3*x2*x2y1*x2y2 + 
	2*x1y2*x1y3*x2*x4*y3 - 2*x1y2*x1y3*x2y1*x2y1*x2y1 + 2*x1y2*x1y3*x2y1*x4*y2 - 
	2*x1y2*x1y3*x3*x3*y3 - 2*x1y2*x1y3*x3*x3y1*y2 + 2*x1y2*x2*x2y1*x3y1*y4 - 
	2*x1y2*x2*x2y2*x3y1*y3 - 2*x1y2*x2y1*x2y1*x3*y4 + 2*x1y2*x2y1*x2y1*x3y1*y3 - 
	2*x1y2*x2y1*x2y2*x3y1*y2 + 2*x1y2*x2y2*x2y2*x3*y2 + x1y3*x1y3*x2*x2y1*x2y1 - 
	x1y3*x1y3*x2*x4*y2 + x1y3*x1y3*x3*x3*y2 - 2*x1y3*x2*x2y1*x3y1*y3 + 
	2*x1y3*x2*x2y2*x3y1*y2 + 2*x1y3*x2y1*x2y1*x3*y3 - 2*x1y3*x2y1*x2y2*x3*y2 - 
	x2*x2y1*x2y1*x2y2*y4 + 2*x2*x2y1*x2y2*x2y2*y3 - x2*x2y2*x2y2*x2y2*y2 + 
	x2*x2y2*x4*y2*y4 - x2*x2y2*x4*y3*y3 - x2*x3y1*x3y1*y2*y4 + x2*x3y1*x3y1*y3*y3 + 
	x2y1*x2y1*x2y1*x2y1*y4 - 2*x2y1*x2y1*x2y1*x2y2*y3 + x2y1*x2y1*x2y2*x2y2*y2 - x2y1*x2y1*x4*y2*y4 + 
	x2y1*x2y1*x4*y3*y3 + 2*x2y1*x3*x3y1*y2*y4 - 2*x2y1*x3*x3y1*y3*y3 - 
	x2y2*x3*x3*y2*y4 + x2y2*x3*x3*y3*y3;

	dab1= x1y2*x1y2*(x1*(x1y3*x2y1 + x3*y4 + x3y1*y3) + x1y1*x1y1*x2y2 + 
	x1y1*(-2*x1y3*x2 + x3*y3 + 2*x3y1*y2) + x1y3*x3*y1 - x2*x2*y4 - 
	2*x2*x2y1*y3 + x2y1*x2y1*y2 + 2*x2y1*x2y2*y1) + 
	x1y2*(x1y1*(x1*x1y3*x2y2 - x1*x3y1*y4 + x1y3*x3y1*y1 + 3*x2*x2y1*y4 + 
	x2*x2y2*y3 + x2y1*x2y1*y3 - 4*x2y1*x2y2*y2 - x2y2*x2y2*y1) - 
	x1y3*(2*x1*x3*y3 + x1*x3y1*y2 - 2*x2*x2*y3 - x2*x2y1*y2 + x2*x2y2*y1 + 
	2*x2y1*x2y1*y1) - x1*x2y1*x2y1*y4 + x1*x2y2*x2y2*y2 + 
	x1y1*x1y1*(x1y3*x2y1 - x3*y4 - 2*x3y1*y3) + x2*x3y1*y1*y4 - 
	x2*x3y1*y2*y3 - x2y1*x3*y1*y4 + x2y1*x3*y2*y3 + x2y1*x3y1*y1*y3 - 
	x2y1*x3y1*y2*y2 - x2y2*x3*y1*y3 + x2y2*x3*y2*y2) - 
	x1y2*x1y2*x1y2*(x1*x2y2 + x1y1*x2y1 + x3*y2 + x3y1*y1) + 
	x1y1*(x1*((-x1y3*x1y3)*x2y1 + x1y3*x3y1*y3 + x2y2*(x2y1*y4 - x2y2*y3)) - 
	x1y3*x1y3*x3*y1 + x1y3*(-3*x2*x2y1*y3 + x2*x2y2*y2 + x2y1*x2y1*y2 + 
	x2y1*x2y2*y1) - x2*x3y1*y2*y4 + x2*x3y1*y3*y3 + x2y1*x3*y2*y4 - 
	x2y1*x3*y3*y3 - x2y1*x3y1*y1*y4 + x2y1*x3y1*y2*y3 + x2y2*x3*y1*y4 - 
	x2y2*x3*y2*y3) + x1*x1y3*x1y3*x3*y2 + x1*x1y3*x2y1*x2y1*y3 - 
	x1*x1y3*x2y1*x2y2*y2 + x1*x2y1*x3y1*y2*y4 - x1*x2y1*x3y1*y3*y3 - 
	x1*x2y2*x3*y2*y4 + x1*x2y2*x3*y3*y3 + x1y1*x1y1*x1y1*(x3y1*y4 - x1y3*x2y2) + 
	x1y1*x1y1*(x1y3*x1y3*x2 + x1y3*x3*y3 - x1y3*x3y1*y2 - x2*x2y2*y4 - 
	x2y1*x2y1*y4 + x2y1*x2y2*y3 + x2y2*x2y2*y2) + x1y2*x1y2*x1y2*x1y2*x2 - x1y3*x1y3*x2*x2*y2 + 
	x1y3*x1y3*x2*x2y1*y1 - x1y3*x2*x3y1*y1*y3 + x1y3*x2*x3y1*y2*y2 + 
	x1y3*x2y1*x3*y1*y3 - x1y3*x2y1*x3*y2*y2 + x2*x2*x2y2*y2*y4 - 
	x2*x2*x2y2*y3*y3 - x2*x2y1*x2y1*y2*y4 + x2*x2y1*x2y1*y3*y3 - 
	x2*x2y1*x2y2*y1*y4 + x2*x2y1*x2y2*y2*y3 + x2*x2y2*x2y2*y1*y3 - 
	x2*x2y2*x2y2*y2*y2 + x2y1*x2y1*x2y1*y1*y4 - x2y1*x2y1*x2y1*y2*y3 - x2y1*x2y1*x2y2*y1*y3 + 
	x2y1*x2y1*x2y2*y2*y2;

	dab2= x1y2*x1y2*(-2*x1*x2y1*x3y1 - x1*x2y2*x3 + x1y1*x1y1*(-x4) + 
	x2*(x1y1*x3y1 + x2y1*x2y1 - x4*y2) + x1y1*x2y1*x3 + x2*x2*x2y2 - 
	x2y1*x4*y1 + x3*x3*y2 + x3*x3y1*y1) + 
	x1y2*(x1y3*((-x1)*x1y1*x4 + x1*x2y1*x3 + x1y1*x2*x3 - x2*x2*x2y1 + 
	x2*x4*y1 - x3*x3*y1) - x1y1*((-x1)*x2y2*x3y1 + 4*x2*x2y1*x2y2 - 
	x2*x4*y3 + x2y1*x2y1*x2y1 - 3*x2y1*x4*y2 - x2y2*x4*y1 + x3*x3*y3 + 
	3*x3*x3y1*y2 + x3y1*x3y1*y1) + 2*x1*x2y1*x2y1*x2y2 - x1*x2y1*x4*y3 - 
	x1*x2y2*x4*y2 + x1*x3*x3y1*y3 + x1*x3y1*x3y1*y2 + 
	x1y1*x1y1*(x2y1*x3y1 + x2y2*x3) - x2*x2*x3y1*y3 + x2*x2y1*x3*y3 + 
	x2*x2y1*x3y1*y2 + x2*x2y2*x3*y2 - x2*x2y2*x3y1*y1 - 2*x2y1*x2y1*x3*y2 + 
	x2y1*x2y1*x3y1*y1) + 
	x1y1*(x1y3*(x1*x2y1*x3y1 + 2*x2*x2y1*x2y1 - x2*x4*y2 - x2y1*x4*y1 + 
	x3*x3*y2 + x3*x3y1*y1) - x1*(x2y1*x2y2*x2y2 - x2y2*x4*y3 + x3y1*x3y1*y3) - 
	x2*x2y2*x3*y3 + x2*x2y2*x3y1*y2 + x2y1*x2y1*x3*y3 - 2*x2y1*x2y1*x3y1*y2 + 
	x2y1*x2y2*x3*y2 + x2y1*x2y2*x3y1*y1 - x2y2*x2y2*x3*y1) + 
	x1y2*x1y2*x1y2*(x1*x4 - x2*x3) - x1*x1y3*x2y1*x2y1*x2y1 + x1*x1y3*x2y1*x4*y2 - 
	x1*x1y3*x3*x3y1*y2 + x1*x2y1*x2y1*x3y1*y3 - x1*x2y1*x2y2*x3*y3 - 
	x1*x2y1*x2y2*x3y1*y2 + x1*x2y2*x2y2*x3*y2 + x1y1*x1y1*x1y1*(x1y3*x4 - x2y2*x3y1) + 
	x1y1*x1y1*(x2*(x2y2*x2y2 - x1y3*x3y1) - x2y1*(2*x1y3*x3 + x4*y3) + 
	x2y1*x2y1*x2y2 - x2y2*x4*y2 + x3*x3y1*y3 + x3y1*x3y1*y2) + 
	x1y3*x2*x2*x3y1*y2 - x1y3*x2*x2y1*x3*y2 - x1y3*x2*x2y1*x3y1*y1 + 
	x1y3*x2y1*x2y1*x3*y1 + x2*x2*x2y1*x2y2*y3 - x2*x2*x2y2*x2y2*y2 - x2*x2y1*x2y1*x2y1*y3 + 
	x2*x2y1*x2y2*x2y2*y1 - x2*x2y2*x4*y1*y3 + x2*x2y2*x4*y2*y2 + 
	x2*x3y1*x3y1*y1*y3 - x2*x3y1*x3y1*y2*y2 + x2y1*x2y1*x2y1*x2y1*y2 - x2y1*x2y1*x2y1*x2y2*y1 + 
	x2y1*x2y1*x4*y1*y3 - x2y1*x2y1*x4*y2*y2 - 2*x2y1*x3*x3y1*y1*y3 + 
	2*x2y1*x3*x3y1*y2*y2 + x2y2*x3*x3*y1*y3 - x2y2*x3*x3*y2*y2;

	dab3= x1y2*(x1y1*((-x1)*x2y2*x2y2 + x1*x4*y4 + x1y3*x2*x2y1 - x1y3*x4*y1 - 
	x2*x3*y4 - 2*x2*x3y1*y3 - 2*x2y1*x3*y3 + x2y1*x3y1*y2 + 
	3*x2y2*x3*y2 + x2y2*x3y1*y1) + x1*((-x1y3)*x2y1*x2y1 + x1y3*x4*y2 - 
	x2y1*x3*y4 + x2y1*x3y1*y3 + x2y2*x3*y3 - x2y2*x3y1*y2) + 
	x1y1*x1y1*(2*x4*y3 - 2*x2y1*x2y2) - x1y3*x2*x3*y2 + x1y3*x2y1*x3*y1 + 
	x2*x2*x2y1*y4 - x2*x2*x2y2*y3 + x2*x2y1*x2y1*y3 - 2*x2*x2y1*x2y2*y2 + 
	x2*x2y2*x2y2*y1 - x2*x4*y1*y4 + x2*x4*y2*y3 + x2y1*x2y1*x2y1*(-y2) + 
	x2y1*x2y1*x2y2*y1 - x2y1*x4*y1*y3 + x2y1*x4*y2*y2 + x3*x3*y1*y4 - 
	x3*x3*y2*y3 + x3*x3y1*y1*y3 - x3*x3y1*y2*y2) + 
	x1y2*x1y2*(x1*x2y1*x2y2 - x1*x4*y3 + x1y1*(x2*x2y2 + x2y1*x2y1 - 2*x4*y2) + 
	x2*x3*y3 + x2*x3y1*y2 + x2y1*x3*y2 - x2y1*x3y1*y1 - 2*x2y2*x3*y1) + 
	x1y1*(x1*(x1y3*x2y1*x2y2 - x1y3*x4*y3 - x2y1*x3y1*y4 + x2y2*x3y1*y3) - 
	x2y1*(2*x1y3*x3*y2 - 3*x2*x2y2*y3 + x2y2*x2y2*y1 - x4*y1*y4 + 
	x4*y2*y3) + x1y3*x2*x3*y3 + x1y3*x2y2*x3*y1 + 
	x2y1*x2y1*(x2y2*y2 - 2*x2*y4) - x2*x2y2*x2y2*y2 + x2*x4*y2*y4 - 
	x2*x4*y3*y3 - x3*x3*y2*y4 + x3*x3*y3*y3 - x3*x3y1*y1*y4 + 
	x3*x3y1*y2*y3) + x1*x1y3*x2y1*x3*y3 - x1*x1y3*x2y2*x3*y2 + 
	x1*x2y1*x2y1*x2y1*y4 - 2*x1*x2y1*x2y1*x2y2*y3 + x1*x2y1*x2y2*x2y2*y2 - 
	x1*x2y1*x4*y2*y4 + x1*x2y1*x4*y3*y3 + x1*x3*x3y1*y2*y4 - 
	x1*x3*x3y1*y3*y3 + x1y1*x1y1*x1y1*(x2y2*x2y2 - x4*y4) + 
	x1y1*x1y1*((-x1y3)*x2*x2y2 + x1y3*x4*y2 + x2*x3y1*y4 + 2*x2y1*x3*y4 - 
	2*x2y2*x3*y3 - x2y2*x3y1*y2) + x1y2*x1y2*x1y2*(x4*y1 - x2*x2y1) - 
	x1y3*x2*x2*x2y1*y3 + x1y3*x2*x2*x2y2*y2 + x1y3*x2*x2y1*x2y1*y2 - 
	x1y3*x2*x2y1*x2y2*y1 + x1y3*x2*x4*y1*y3 - x1y3*x2*x4*y2*y2 - 
	x1y3*x3*x3*y1*y3 + x1y3*x3*x3*y2*y2 - x2*x2*x3y1*y2*y4 + x2*x2*x3y1*y3*y3 + 
	x2*x2y1*x3*y2*y4 - x2*x2y1*x3*y3*y3 + x2*x2y1*x3y1*y1*y4 - 
	x2*x2y1*x3y1*y2*y3 - x2*x2y2*x3y1*y1*y3 + x2*x2y2*x3y1*y2*y2 - 
	x2y1*x2y1*x3*y1*y4 + x2y1*x2y1*x3*y2*y3 + x2y1*x2y2*x3*y1*y3 - 
	x2y1*x2y2*x3*y2*y2;

	dab4= x1y2*(-2*x1*(x1y3*x2y1*x2y2 - x1y3*x4*y3 - x2y1*x3y1*y4 + x2y2*x3y1*y3) + 
	x1y1*x1y1*(x4*y4 - x2y2*x2y2) + x1y1*(x1y3*x2*x2y2 - x1y3*x2y1*x2y1 - 
	x2*x3y1*y4 - x2y1*x3*y4 + x2y1*x3y1*y3 + x2y2*x3*y3) - 
	2*x1y3*x2*x3*y3 - x1y3*x2*x3y1*y2 + x1y3*x2y1*x3*y2 + 
	x1y3*x2y1*x3y1*y1 + x1y3*x2y2*x3*y1 - x2*x2y1*x2y1*y4 + x2*x2y2*x2y2*y2 + 
	2*x2y1*x2y1*x2y2*y2 - 2*x2y1*x2y2*x2y2*y1 + x2y1*x4*y1*y4 - x2y1*x4*y2*y3 + 
	x2y2*x4*y1*y3 - x2y2*x4*y2*y2 - x3*x3y1*y1*y4 + x3*x3y1*y2*y3 - 
	x3y1*x3y1*y1*y3 + x3y1*x3y1*y2*y2) + 
	x1y2*x1y2*(x1*x2y2*x2y2 - x1*x4*y4 + x1y1*x2y1*x2y2 - x1y1*x4*y3 + 
	x1y3*x2*x2y1 - x1y3*x4*y1 + x2*x3*y4 + x2*x3y1*y3 - 2*x2y1*x3y1*y2 - 
	x2y2*x3*y2 + x2y2*x3y1*y1) + x1*x1y3*x1y3*x2y1*x2y1 - x1*x1y3*x1y3*x4*y2 - 
	2*x1*x1y3*x2y1*x3y1*y3 + 2*x1*x1y3*x2y2*x3y1*y2 - x1*x2y1*x2y1*x2y2*y4 + 
	2*x1*x2y1*x2y2*x2y2*y3 - x1*x2y2*x2y2*x2y2*y2 + x1*x2y2*x4*y2*y4 - 
	x1*x2y2*x4*y3*y3 - x1*x3y1*x3y1*y2*y4 + x1*x3y1*x3y1*y3*y3 + 
	x1y1*x1y1*(x1y3*x2y1*x2y2 - x1y3*x4*y3 - x2y1*x3y1*y4 + x2y2*x3y1*y3) + 
	x1y1*(x1y3*x1y3*(x4*y1 - x2*x2y1) + x1y3*(x2*x3y1*y3 + x2y1*x3*y3 + 
	x2y1*x3y1*y2 - x2y2*(x3*y2 + 2*x3y1*y1)) + 
	x2y2*(x2*x2y1*y4 - 2*x2y1*x2y1*y3 - x4*y1*y4 + x4*y2*y3) - 
	x2*x2y2*x2y2*y3 + x2y1*x2y1*x2y1*y4 - x2y1*x4*y2*y4 + x2y1*x4*y3*y3 + x2y2*x2y2*x2y2*y1 + 
	x3*x3y1*y2*y4 - x3*x3y1*y3*y3 + x3y1*x3y1*y1*y4 - x3y1*x3y1*y2*y3) + 
	x1y2*x1y2*x1y2*(x4*y2 - x2*x2y2) + x1y3*x1y3*x2*x3*y2 - x1y3*x1y3*x2y1*x3*y1 + 
	x1y3*x2*x2y1*x2y1*y3 - x1y3*x2*x2y1*x2y2*y2 - x1y3*x2y1*x2y1*x2y1*y2 + 
	x1y3*x2y1*x2y1*x2y2*y1 - x1y3*x2y1*x4*y1*y3 + x1y3*x2y1*x4*y2*y2 + 
	x1y3*x3*x3y1*y1*y3 - x1y3*x3*x3y1*y2*y2 + x2*x2y1*x3y1*y2*y4 - 
	x2*x2y1*x3y1*y3*y3 - x2*x2y2*x3*y2*y4 + x2*x2y2*x3*y3*y3 - 
	x2y1*x2y1*x3y1*y1*y4 + x2y1*x2y1*x3y1*y2*y3 + x2y1*x2y2*x3*y1*y4 - 
	x2y1*x2y2*x3*y2*y3 + x2y1*x2y2*x3y1*y1*y3 - x2y1*x2y2*x3y1*y2*y2 - 
	x2y2*x2y2*x3*y1*y3 + x2y2*x2y2*x3*y2*y2;

	dab5= x1y3*(x1y2*(x1*x2y1*x3y1 + x1*x2y2*x3 + x1y1*x1y1*(-x4) + 
	x2*(x1y1*x3y1 - 2*x2y1*x2y1 + x4*y2) + x1y1*x2y1*x3 - x2*x2*x2y2 + 
	2*x2y1*x4*y1 - x3*x3*y2 - 2*x3*x3y1*y1) + 
	x1*(-2*x1y1*x2y2*x3y1 + x2y1*x2y1*x2y2 - x2y1*x4*y3 + x3*x3y1*y3) + 
	x1y2*x1y2*(x2*x3 - x1*x4) + x1y1*x1y1*x2y2*x3 + x1y1*x2*x4*y3 - 
	x1y1*x2y1*x4*y2 - x1y1*x3*x3*y3 + x1y1*x3*x3y1*y2 - x2*x2*x3y1*y3 + 
	x2*x2y1*x3*y3 - x2*x2y1*x3y1*y2 + 2*x2*x2y2*x3y1*y1 + x2y1*x2y1*x3*y2 - 
	2*x2y1*x2y2*x3*y1) + 
	x1y2*(x1*(-2*x2y1*x2y2*x2y2 + x2y1*x4*y4 + x2y2*x4*y3 - x3*x3y1*y4 - 
	x3y1*x3y1*y3) + x1y1*x1y1*x2y2*x3y1 + 
	x1y1*((-x2)*x4*y4 + x2y1*x2y1*x2y2 - x2y1*x4*y3 + x2y2*x4*y2 + x3*x3*y4 + 
	x3*x3y1*y3 - x3y1*x3y1*y2) + x2*x2*x3y1*y4 - x2*x2y1*x3*y4 + 
	x2*x2y1*x3y1*y3 - x2*x2y2*x3*y3 - x2*x2y2*x3y1*y2 + x2y1*x2y1*x3y1*y2 - 
	2*x2y1*x2y2*x3y1*y1 + 2*x2y2*x2y2*x3*y1) + 
	x1y2*x1y2*(x1*x2y2*x3y1 - x1y1*(x2y1*x3y1 + 2*x2y2*x3) + 2*x2*x2y1*x2y2 - 
	x2y1*x4*y2 - x2y2*x4*y1 + x3*x3y1*y2 + x3y1*x3y1*y1) + 
	x1y3*x1y3*(x1*x1y1*x4 - x1*x2y1*x3 - x2*(x1y1*x3 + x4*y1) + x2*x2*x2y1 + 
	x3*x3*y1) + x1*x1y1*x2y2*x2y2*x2y2 - x1*x1y1*x2y2*x4*y4 + x1*x1y1*x3y1*x3y1*y4 - 
	x1*x2y1*x2y1*x3y1*y4 + x1*x2y1*x2y2*x3*y4 + x1*x2y1*x2y2*x3y1*y3 - 
	x1*x2y2*x2y2*x3*y3 - x1y1*x1y1*x2y1*x2y2*x2y2 + x1y1*x1y1*x2y1*x4*y4 - 
	x1y1*x1y1*x3*x3y1*y4 + x1y2*x1y2*x1y2*(x1y1*x4 - x2*x3y1) + x1y1*x2*x2y2*x3*y4 - 
	x1y1*x2*x2y2*x3y1*y3 - x1y1*x2y1*x2y1*x3*y4 + x1y1*x2y1*x2y2*x3*y3 + 
	x1y1*x2y1*x2y2*x3y1*y2 - x1y1*x2y2*x2y2*x3*y2 - x2*x2*x2y1*x2y2*y4 + 
	x2*x2*x2y2*x2y2*y3 + x2*x2y1*x2y1*x2y1*y4 - x2*x2y1*x2y1*x2y2*y3 + x2*x2y1*x2y2*x2y2*y2 - 
	x2*x2y2*x2y2*x2y2*y1 + x2*x2y2*x4*y1*y4 - x2*x2y2*x4*y2*y3 - x2*x3y1*x3y1*y1*y4 + 
	x2*x3y1*x3y1*y2*y3 - x2y1*x2y1*x2y1*x2y2*y2 + x2y1*x2y1*x2y2*x2y2*y1 - x2y1*x2y1*x4*y1*y4 + 
	x2y1*x2y1*x4*y2*y3 + 2*x2y1*x3*x3y1*y1*y4 - 2*x2y1*x3*x3y1*y2*y3 - 
	x2y2*x3*x3*y1*y4 + x2y2*x3*x3*y2*y3;

	b1=dab1/da;
	b2=dab2/da;
	b3=dab3/da;
	b4=dab4/da;
	b5=dab5/da;
	a00=sqrt(4*b2/(b3*b3 -4*b1*b2)*(-1+(b1*b5*b5+b2*b4*b4-b3*b4*b5)/(b3*b3-4*b1*b2))); //楕円軸半径?
	b00=sqrt(4*b1/(b3*b3 -4*b1*b2)*(-1+(b1*b5*b5+b2*b4*b4-b3*b4*b5)/(b3*b3-4*b1*b2))); //楕円軸半径?
	x0=(2*b2*b4 - b3*b5)/(b3*b3-4*b1*b2); //楕円中心座標x
	y0=(-b3*b4 + 2*b1*b5)/(b3*b3-4*b1*b2); //楕円中心座標y
	sp0=-b3/2/sqrt(b1*b2); 
	cp0=sqrt(1-sp0*sp0); 
	aa=atan2(sp0,cp0); //楕円の傾き
	x0+=pp0;
	y0+=ss0;

	if(x0max<x0) x0max=x0;
	if(x0min>x0) x0min=x0;
	if(y0max<y0) y0max=y0;
	if(y0min>y0) y0min=y0;
	if(a00max<a00) a00max=a00;
	if(a00min>a00) a00min=a00;
	if(b00max<b00) b00max=b00;
	if(b00min>b00) b00min=b00;
	if(aamax<aa) aamax=aa;
	if(aamin>aa) aamin=aa;

	// 戻り値
	(*elpara).MM=MM;
	(*elpara).a00=a00;
	(*elpara).b00=b00;
	(*elpara).x0=x0;
	(*elpara).y0=y0;
	(*elpara).sp0=sp0;
	(*elpara).cp0=cp0;
	

} // フィッティング関数（その1）の定義終わり


// 位相決定関数を定義
// 機能： ADCのデータ配列MM個と楕円パラメータを渡されると、位相値とzabsの入った構造体を計算する
// 入力: ADCデータ配列（iad0, iad1, iad2, iad3）, 配列長=フィット長(MM)、楕円パラメータ構造体（ellip_param）
// 出力： 結果の構造体（位相値&zabs）MM点

void fringe_det(int id0[], int id1[], int id2[], int id3[], int MM, struct ellip_param elpara, struct fringe_zabs *fz[], struct fit_eval_param *evalpara){

	int j,k,kk,ndph;
	int flg_nan;//nan（楕円が決まらなかった場合に生じる）検出フラグ
	long n,n1;
	double pin,pinmax,pinmin,amd,amdmax,amdmin,fringe0;
	double dfringemax,dfringemin,x0max,x0min,y0max,y0min;
	double a00max,a00min,b00max,b00min,aamax,aamin;
	//double xxmax,xxmin,xx1;
	double pp,ss,pp0,ss0,zabs,xx,xx0,s0,s00,zabsmax,zabsmin,flat,flatmin,flatmax,absdaa;
	double x4,y4,x1y3,x2y2,x3y1,x3,y3,x1y2,x2y1,x2,y2,x1y1,x1,y1;
	double dsin,dsn,dcos,da,dab1,dab2,dab3,dab4,dab5,b1,b2,b3,b4,b5;
	double a00,b00,x0,y0,sp0,cp0,aa;
	
	double fringe[MM], dfringe[MM], zabs0[MM], flat0[MM], aa0[MM];
	zabsmax=-1.e20;
	zabsmin=1.e20;
	dfringemax=-1.e20;
	dfringemin=1.e20;
	flatmax=-1.e20;
	flatmin=1.e20;
	absdaa=0.;
	
	a00=elpara.a00;
	b00=elpara.b00;
	x0=elpara.x0;
	y0=elpara.y0;
	sp0=elpara.sp0;
	cp0=elpara.cp0;
	aa=atan2(sp0,cp0);
	aa0[0]=0.;
	flg_nan=0;

	
	//楕円->単位円の変換式を当てはめて位相角決定（実際には理想的な単位円からずれる）
	//for(j=0;j<MM;j++){
	for(j=0;j<MM;j++){
	
		pp=(double)id0[j]*ADCFAC/((double)id2[j]*ADCFAC);
		ss=(double)id1[j]*ADCFAC/((double)id2[j]*ADCFAC);
		pin=(double)id2[j]*ADCFAC;
		amd=(double)id3[j]*ADCFAC;
		if(j%N==0){
			if(pinmax<pin) pinmax=pin;
			if(pinmin>pin) pinmin=pin;
			if(amdmax<amd) amdmax=amd;
			if(amdmin>amd) amdmin=amd;
		}
		dcos=(pp-x0)/a00;
		dsin=(ss-y0)/b00;
		dsn=(dsin-sin(aa)*dcos)/cos(aa);
		zabs=sqrt(dsn*dsn+dcos*dcos);
		
		if(isnan(zabs)){
			flg_nan=1;
			//printf("nan detected during fitting\n");
			//exit(0);
		}
		
		//printf("j: %d, zabs: %f\n",j,zabs);
		xx=atan2(dsn,dcos);
//				fprintf(feout,"%25.16e %25.16e %25.16e %25.16e %25.16e %25.16e %25.16e %25.16e \n",pp,ss,a00*cos((double)k)+x0,b00*sin((double)k+aa)+y0,dcos,dsn,cos((double)k),sin((double)k));
		flat=b00/a00;
		fringe[j]=xx;
		absdaa=fabs(aa-aa0[j]);
		
		zabs0[j]=zabs;
		flat0[j]=flat;
		aa0[j]=aa;
		if(j==0){
			dfringe[0]==0;
			fringe0=xx;
		}else{
			dfringe[j]=xx-fringe0;
			fringe0=xx;
		}
				
		if(zabsmax<zabs0[j]) zabsmax=zabs0[j];
		if(zabsmin>zabs0[j]) zabsmin=zabs0[j];
		if(dfringemax<dfringe[j]) dfringemax=dfringe[j];
		if(dfringemin>dfringe[j]) dfringemin=dfringe[j];
		if(flatmax<flat0[j]) flatmax=flat0[j];
		if(flatmin>flat0[j]) flatmin=flat0[j];
		
		// 戻り値代入
		(*fz)[j].fringe=xx;
		(*fz)[j].zabs=zabs;
		//printf("xx:%f\n",xx);
		//printf("zabs:%f\n",zabs);
		//printf("zabsmax:%f\n",zabsmax);	
	}
		
		// 戻り値代入（max, min）
		
		if(flg_nan==1){//フィッティングできずnanが発生した場合は、強制的に区間分割、補間を行うため、敢えてzabsに異常値を代入する。
			(*evalpara).zabsmin=0.;			
			(*evalpara).zabsmax=0.;
		}else{
			(*evalpara).zabsmax=zabsmax;
			(*evalpara).zabsmin=zabsmin;
		}
		
		//(*evalpara).zabsmin=zabsmin;
		//(*evalpara).zabsmin=zabsmin;
		(*evalpara).dfringemin=dfringemin;			
		(*evalpara).dfringemax=dfringemax;			
		(*evalpara).flatmin=flatmin;			
		(*evalpara).flatmax=flatmax;
		(*evalpara).absdaa=absdaa;

} // 位相決定関数の定義終わり



/////////////////////////////////////////////////////////////////////////////////
// 段階的フィッティングのmain
/////////////////////////////////////////////////////////////////////////////////

int main(int argc, char*argv[]){
	int opt; // to enable "-x" style options
	int opt50ka; // option to flag 50 kHz file output (ascii)
	int opt50kb; // option to flag 50 kHz files output (binary)
	int optbelp; // option to flag belp file output
	int optint; // option to enable interpolation
	
	FILE *fin0,*fin1,*fin2,*fin3,*fout,*fbph,*fzdph, *fbelp, *fdiag, *fraw, *fflen, *ffringe, *fzabs;
	char fpath[256],fnin[256],fnin1[256],fnout0[256],fnout[512],fn[256],fn1[512];

	int UNIT_LEN=25; // unit (minimum) data point for elliptic fitting (50 = 1 msec)
//	int MMtab[Nstep]={50,10,5,2,1};	// 各段階でフィッティングに使う点の数。単位はUNIT_LEN で定義（デフォルト50点=1 msec）
//	int MMtab[]={100,20,10,4,2,1};	// 各段階でフィッティングに使う点の数。単位はUNIT_LEN で定義（デフォルト50点=1 msec）
//	int MMtab[]={100,50,20,10,5,4,2,1};	// 各段階でフィッティングに使う点の数。単位はUNIT_LEN で定義（デフォルト25点=0.5msec）
//	int MMtab[]={100,50};	// 各段階でフィッティングに使う点の数。単位はUNIT_LEN で定義（デフォルト25点=0.5msec）
	int MMtab[]={100,50,25,20,10,5,4,2,1};	// 各段階でフィッティングに使う点の数。単位はUNIT_LEN で定義（デフォルト25点=0.5msec）
//	int MMtab[]={20,10,5,4,2,1};	// 各段階でフィッティングに使う点の数。単位はUNIT_LEN で定義（デフォルト25点=0.5msec）
//	int MMtab[]={1};	// 各段階でフィッティングに使う点の数。単位はUNIT_LEN で定義（デフォルト25点=0.5msec）




	int Nstep=sizeof(MMtab)/sizeof(MMtab[0]); // number of gradual fitting

	int M=20; //average output length in Hz
	int MM; // elliptic fit length ### MUST BE EVEN ###
	int j,k, ndph, inp, pnt;
	
	int *flg,*flg0,*flg1; // 各ブロックのフィッティング結果の合否を入れるフラグ（インデックスはUNIT_LENを単位としたブロックの番号（0スタート）。
	// フラグの意味 0:フィット不良（または未定。デフォルト）, 1:フィット良好
	// flg0: 前段のフラグ, flg1: 処理中の段のフラグ
	// ブロック総数が多すぎるため、すべてのブロックについてフラグを作ることはできない。2段階分だけメモリを割り当て、不要になったら解放する。
	
	int flgblk; // フィット処理するブロックに含まれるフラグの和を代入する変数
	
	int istep=0; // フィッティングの段階を示すインデックス fitting stage (0 to 3)
	
	int Nblock_step[Nstep]; // 各段階に含まれるブロックの数
	for(j=0;j<Nstep;j++){
		Nblock_step[j]=60*SAMP/MMtab[j]/UNIT_LEN;
		//printf("Nblock_step[%d]:%d\n",j,Nblock_step[j]);
	};
	
	int iMM=0; // index for MM
	
	int iflg_s, iflg_e; //チェックするフラグのインデックス（先頭と最後）
		
	int *id0,*id1,*id2,*id3; // .AD0*ファイルから読み込んだデータを格納する配列のポインタ（デフォルト1分長）

	int *idfit0,*idfit1,*idfit2,*idfit3; // フィッティングに使うデータを格納する配列のポインタ
	
	struct fringe_zabs *fz_res;
	struct fringe_zabs_out *fz_res_out;
//	struct fringe_zabs fz_res;
	struct fit_param fpara_res;
	struct binary_ph *bph;
	struct ellip_param *belp;
	struct ellip_param elpara; // フィットするブロックのフィッティング結果（楕円パラメータを入れる変数）
	struct ellip_param *elpara_blk; //楕円パラメータを格納する配列用のポインタ（補間処理に用いる）

	struct zabs_dph zdph;
	struct fit_eval_param evalpara, evalpara_int; // フィットするブロックの良否判定に用いるパラメータを入れる変数（暫定的にzabsmin, zabsmax, fringemin, fringemaxとする）
	
	double th_zabsmax, th_zabsmin; // フィッティング結果の合否に用いるzabsの範囲を決める閾値
	th_zabsmax=1.3;
	th_zabsmin=0.5;
	double th_dphmax, th_dphmin; // フィッティング良否判定に用いるdphの範囲を決める閾値
	th_dphmax=PI*3.;
	th_dphmin=-PI*3.;
	
	
	//int NN; //total samples to process
	long n, n1;
	unsigned int isy,ism,isd,ish,ismin,iey,iem,ied,ieh,iemin;
	unsigned int iy,im,id,ih,imin,isec;
	int ims,ids,ihs,imins,ime,ide,ihe,imine;
	//double pin,pinmax,pinmin,amd,amdmax,amdmin,fringe,fringemax,fringemin,fringe0;
	//double pin,pinmax,pinmin,amd,amdmax,amdmin,fringemax,fringemin,fringe0;
	double pinmax,pinmin,amdmax,amdmin,fringemax,fringemin,fringe0;
	double dfringemax,dfringemin,dfringe,dphmax,dphmin,x0max,x0min,y0max,y0min;
	double a00max,a00min,b00max,b00min,aamax,aamin;
	double xxmax,xxmin,xx1;
	//double pp,ss,pp0,ss0,aa,xx,xx0,s0,s00;
	double pp,ss,pp0,ss0,xx,xx0,s0,s00;
	//double zabs,zabsmax,zabsmin;
	double zabsmax,zabsmin;
	double x4,y4,x1y3,x2y2,x3y1,x3,y3,x1y2,x2y1,x2,y2,x1y1,x1,y1;
	//double p00=1.e20,p01=1.e20;
	double dsin,dsn,dcos,da,dab1,dab2,dab3,dab4,dab5,b1,b2,b3,b4,b5;
	//double a00,b00,x0,y0,sp0,cp0,x00,y00,dxy;
	double sp0,cp0,x00,y00,dxy;
	
	double *fringe, *zabs, *pin, *amd, *a00, *b00, *x0, *y0, *aa; // pointers for array variables 
	int *MMp; // MM for record
	
	double *fringe_out, *zabs_out, *pin_out, *amd_out, *a00_out, *b00_out, *x0_out, *y0_out, *aa_out; //最終結果を入れる変数
	int *MMp_out;
	
	int ll,ll0,llmax,llmin,ll1,dl;
	//int nmin1;
	//double sx,sxx,sy,sxy,g1,g2,g3,gav,gsd,a,b,q1,qmin1,nn;
//	double PI=3.1415926535897932384626;

	//補間処理に関わる変数
	int flg_int; //補間中か否かを区別するフラグ（0:補間していない, 1:補間処理中）
	int int_s, int_e; //補間開始&終了ブロック番号
	struct ellip_param ellp_s, ellp_e; //補間開始&終了ブロックの楕円パラメータ
	struct ellip_param ellp_int; //補間点の楕円パラメータ

	char num[10][2]={"0","1","2","3","4","5","6","7","8","9"};
	char out_rate[12], fit_len[12];

	// print start time
	time_t t_s=time(NULL);
	printf("start time: %s",ctime(&t_s));


	snprintf(out_rate, 12, "%d", M*MM/1000);
	snprintf(fit_len, 12, "%d", 1000*MM/SAMP);

	//printf("sizeof bph: %d\n",sizeof(struct binary_ph));
	//printf("sizeof bzabs: %d\n\n",sizeof(double));

	s0=0.0;	// 位相平均値算出用変数（20Hz, 50ms）
	s00=0.0;	// 位相平均値算出用変数（1 Hz, 1s） total phase for 1s
	n=0;	// data number
	n1=0;	// data number (high)
	zabsmax=-1.e20;
	zabsmin=1.e20;
	fringemax=-1.e20;
	fringemin=1.e20;
	dfringemax=-1.e20;
	dfringemin=1.e20;
	dphmax=-1.e20;
	dphmin=1.e20;
	pinmax=-1.e20;
	pinmin=1.e20;
	amdmax=-1.e20;
	amdmin=1.e20;
	x0max=-1.e20;
	x0min=1.e20;
	y0max=-1.e20;
	y0min=1.e20;
	a00max=-1.e20;
	a00min=1.e20;
	b00max=-1.e20;
	b00min=1.e20;
	aamax=-1.e20;
	aamin=1.e20;
	
	x00=0.0;
	y00=0.0;
	dxy=0.0;
	
	while((opt = getopt(argc, argv, "Hhei")) != -1){
		switch(opt){
			case 'h':
				opt50ka =1;
				break;
			case 'H':
				opt50kb =1;
				break;
			case 'e':
				optbelp =1;
				break;
			case 'i':
				optint =1;
				break;
			default:
				printf("error! \'%c\' \'%c\'\n", opt, optopt);
				return 1;
		}
		
	}

#if test // test mode
	strcpy(fpath,"/home/takamori/GIF/phase/data");
	//strcpy(fpath,"C:\\Users\\takamori\\Desktop\\phase_devel\\data");
	////strcpy(fpath,"V:\\w_backup\\GIF\\phase\\data");
	//strcpy(fpath,"C:\\Users\\takamori\\Desktop\\GIF\\phase\\data");
	isy=23;
	ism=5;
	isd=5;
	ish=14;
	ismin=42;
	iey=isy;
	iem=ism;
	ied=isd;
	ieh=ish;
	iemin=ismin;


#else // normal mode
	//printf("input data path (e.g. g:\\GIF\\)? \n");
	//scanf("%s",fpath);
	//strcpy(fpath,"h:\\GIF\\");
	//strcpy(fpath,"180111\\");
	//strcpy(fpath,"V:\\w_backup\\GIF\\phase\\data");
	//strcpy(fpath,"C:\\Users\\takamori\\Desktop\\GIF\\phase\\data");
	strcpy(fpath,"/home/takamori/GIF/phase/data");
	//strcpy(fpath,"C:\\Users\\takamori\\Desktop\\GIF\\phase\\data");
	//strcat(fpath,"2024\\01\\");
	//strcat(fpath,"20");
	printf("start year (2 digits)? \n");
	scanf("%d",&isy);
	//isy=24;
	printf("start month ? \n");
	scanf("%d",&ism);
	//ism=1;
	printf("start day ? \n");
	scanf("%d",&isd);
	printf("start hour ? \n");
	scanf("%d",&ish);
	printf("start minute ? \n");
	scanf("%d",&ismin);

	//printf("stop year ? \n");
	//scanf("%d",&iey);
	//iey=24;
	iey=isy;
	//printf("stop month ? \n");
	//scanf("%d",&iem);
	//iem=1;
	iem=ism;
	printf("stop day ? \n");
	scanf("%d",&ied);
	printf("stop hour ? \n");
	scanf("%d",&ieh);
	printf("stop minute ? \n");
	scanf("%d",&iemin);

//	printf("output file name ? [*_1s.txt] [*.bph] [*.zdph]\n");
//	scanf("%s",fnout);

#endif


	strcpy(fnout0,"dev");
	//sprintf(fnout,"%s%02d%02d%02d%02d%02d",fnout,(unsigned int)isy,(unsigned int)ism,(unsigned int)isd,(unsigned int)ish,(unsigned int)ismin);
	sprintf(fnout,"%s_20%02d%02d%02d%02d%02d_%02d%02d%02d",fnout0,isy,ism,isd,ish,ismin,ied,ieh,iemin);

	//strcpy(fnout,"fit50msblk_dev");
	//sprintf(fn,"%s_20%s%s%s%s",fnout,isy,isd,ish,ism);

	strcpy(fn,fnout);
	//sprintf(fn, "%s_%skHz_fit%sms", fn, out_rate, fit_len);
	strcat(fn,"_1s.txt");
	if((fout=fopen(fn,"w"))==NULL){
		printf("file open error (%s) !! \n",fn);
		return(-1);
	}
	strcpy(fn,fnout);
	//sprintf(fn, "%s_%skHz_fit%sms", fn, out_rate, fit_len);
	strcat(fn,".bph");
	if((fbph=fopen(fn,"wb"))==NULL){
		puts("output file open error (.bph) !!");
		return(-1);
	}

	strcpy(fn,fnout);
	//sprintf(fn, "%s_%skHz_fit%sms", fn, out_rate, fit_len);
	strcat(fn,".zdph");
	if((fzdph=fopen(fn,"wb"))==NULL){
		puts("output file open error (.zdph) !!");
		return(-1);
	}
	
	if(optbelp==1){
		strcpy(fn,fnout);
		//sprintf(fn, "%s_%skHz_fit%sms", fn, out_rate, fit_len);
		strcat(fn,".belp");
		if((fbelp=fopen(fn,"wb"))==NULL){
			puts("output file open error (.belp) !!");
			return(-1);
		}
	}

	strcpy(fn,fnout);
	//sprintf(fn, "%s_%skHz_fit%sms", fn, out_rate, fit_len);
	strcat(fn,".diag");
	if((fdiag=fopen(fn,"wb"))==NULL){
		puts("output file open error (.diag) !!");
		return(-1);
	}
	
	if(opt50ka==1){
		strcpy(fn,fnout);
		sprintf(fn1, "%s_50k", fn);
		strcat(fn1,".dat");
		if((fraw=fopen(fn1,"w"))==NULL){
			puts("output file open error (.dat) !!");
			return(-1);
		}
	}

	if(opt50kb==1){
		strcpy(fn,fnout);
		sprintf(fn1, "%s_50k", fn);
		strcat(fn1,".bfl");
		if((fflen=fopen(fn1,"wb"))==NULL){
			puts("output file open error (.bfl) !!");
			return(-1);
		}
			
		strcpy(fn,fnout);
		sprintf(fn1, "%s_50k", fn);
		strcat(fn1,".bfr");
		if((ffringe=fopen(fn1,"wb"))==NULL){
			puts("output file open error (.bfr) !!");
			return(-1);
		}

		strcpy(fn,fnout);
		sprintf(fn1, "%s_50k", fn);
		strcat(fn1,".bzabs");
		if((fzabs=fopen(fn1,"wb"))==NULL){
			puts("output file open error (.bzabs) !!");
			return(-1);
		}
	}
	

	//fprintf(fout,"# zabsmax zabsmin phase_1s.x phase_1s.l phasemax.x phasemax.l phasemin.x phasemin.l dphasemax dphasemin pinmax pinmin amdmax amdmin x0max x0min y0max y0min a00max a00min b00max b00min aamax aamin \n");
	//fprintf(fout,"# zabsmax zabsmin phase_1s.x phase_1s.l phasemax.x phasemax.l phasemin.x phasemin.l dphasemax dphasemin pinmax pinmin amdmax amdmin x0max x0min y0max y0min a00max a00min b00max b00min aamax aamin \n");
	fprintf(fout,"year month day hour minute second zabsmax zabsmin phase_1s.x phase_1s.l phasemax.x phasemax.l phasemin.x phasemin.l dphasemax dphasemin pinmax pinmin amdmax amdmin\n");

	for(iy=isy;iy<=iey;iy++){ // year count 
		if(iy==isy)	ims=ism;
		else	ims=1;
		if(iy==iey)	ime=iem;
		else	ime=12;

		for(im=ims;im<=ime;im++){ // month count 
			if((iy==isy)&&(im==ism))	ids=isd;
			else	ids=1;
			if((iy==iey)&&(im==iem))	ide=ied;
			else{
				if((im==4)||(im==6)||(im==9)||(im==11)) ide=30;
				else if((im==2)&&((iy%4)==0)) ide=29;
				else if(im==2) ide=28;
				else ide=31;
			}

			for(id=ids;id<=ide;id++){ // day counter 
				if((iy==isy)&&(im==ism)&&(id==isd))	ihs=ish;
				else	ihs=0;
				if((iy==iey)&&(im==iem)&&(id==ied))	ihe=ieh;
				else	ihe=23;

				for(ih=ihs;ih<=ihe;ih++){ // hour counter 
					if((iy==isy)&&(im==ism)&&(id==isd)&&(ih==ish))	imins=ismin;
					else	imins=0;
					if((iy==iey)&&(im==iem)&&(id==ied)&&(ih==ieh))	imine=iemin;
					else	imine=59;

					for(imin=imins;imin<=imine;imin++){ // minute counter 
					printf("\nprocessing data of: 20%02d%02d%02d%02d%02d\n", iy,im,id,ih,imin);
					//printf("n:%ld, n1:%ld\n", n, n1);

					// make file name 
						strcpy(fnin,fpath);
					//	strcat(fnin,"50000Hz\\20");
						strcat(fnin,"/20");
						strcat(fnin,num[iy/10]);
						strcat(fnin,num[iy-(iy/10)*10]);
						strcat(fnin,"/");
						strcat(fnin,num[im/10]);
						strcat(fnin,num[im-(im/10)*10]);
						strcat(fnin,"/");
						strcat(fnin,num[id/10]);
						strcat(fnin,num[id-(id/10)*10]);
						strcat(fnin,"/");
						strcat(fnin,num[ih/10]);
						strcat(fnin,num[ih-(ih/10)*10]);
						strcat(fnin,"/");
						strcat(fnin,num[iy/10]);
						strcat(fnin,num[iy-(iy/10)*10]);
						strcat(fnin,num[im/10]);
						strcat(fnin,num[im-(im/10)*10]);
						strcat(fnin,num[id/10]);
						strcat(fnin,num[id-(id/10)*10]);
						strcat(fnin,num[ih/10]);
						strcat(fnin,num[ih-(ih/10)*10]);
						strcat(fnin,num[imin/10]);
						strcat(fnin,num[imin-(imin/10)*10]);

						if(imin==0)	printf("%s \n",fnin);

						strcpy(fnin1,fnin);
						strcat(fnin1,".AD00");
						if((fin0=fopen(fnin1,"rb"))==NULL){
							printf("file open error (%s) !! \n",fnin1);
							return(-1);
						}
						strcpy(fnin1,fnin);
						strcat(fnin1,".AD01");
						if((fin1=fopen(fnin1,"rb"))==NULL){
							printf("file open error (%s) !! \n",fnin1);
							return(-1);
						}
						strcpy(fnin1,fnin);
						strcat(fnin1,".AD02");
						if((fin2=fopen(fnin1,"rb"))==NULL){
							printf("file open error (%s) !! \n",fnin1);
							return(-1);
						}
						strcpy(fnin1,fnin);
						strcat(fnin1,".AD03");
						if((fin3=fopen(fnin1,"rb"))==NULL){
							printf("file open error (%s) !! \n",fnin1);
							return(-1);
						}
					//	printf("output data start time [s] ? \n");
					//	scanf("%d",&ntest);

					//printf("n= ? \n");
					//scanf("%d",&jk);

					/*	if(nf==1){
							strcpy(fnin,fn);
							strcat(fnin,".txt");
							if((fpout=fopen(fnin,"w"))==NULL){
								puts("output file open error (.txt) !!");
								return(-1);
							}
							if((fraw=fopen("test2.txt","w"))==NULL){
								puts("output file open error (.txt) !!");
								return(-1);
							}
						}
					*/
					
					// データ読み込み用配列を確保
						id0=calloc(60*SAMP,sizeof(int)); if (id0==NULL)exit(0);
						id1=calloc(60*SAMP,sizeof(int)); if (id1==NULL)exit(0);
						id2=calloc(60*SAMP,sizeof(int)); if (id2==NULL)exit(0);
						id3=calloc(60*SAMP,sizeof(int)); if (id3==NULL)exit(0);

					// データファイル読み込み（1分単位）
						//printf("reading: %s\n", fnin1);
						if(fread(id0,sizeof(int),60*SAMP,fin0)!=60*SAMP){
							puts("data reading error ! ");
							return(-1);
						}
						if(fread(id1,sizeof(int),60*SAMP,fin1)!=60*SAMP){
							puts("data reading error ! ");
							return(-1);
						}
						if(fread(id2,sizeof(int),60*SAMP,fin2)!=60*SAMP){
							puts("data reading error ! ");
							return(-1);
						}
						if(fread(id3,sizeof(int),60*SAMP,fin3)!=60*SAMP){
							puts("data reading error ! ");
							return(-1);
						}
						
						/*
						for(j=0;j<MM;j++){
							printf("data read\n");
							printf("j=%d: id0 =%d \n",j,id0[j]);
						}
						*/
					

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// フィッティング処理
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

						fz_res_out=calloc(60*SAMP,sizeof(struct fringe_zabs_out)); // 良好なフィッティング結果を入れる変数（デフォルト1分ぶん）
						fringe=calloc(60*SAMP,sizeof(double)); // 位相の2PI指数を決める際に一時的に使う変数（デフォルト1分ぶん）

						//zabs_out = calloc(sizeof(double)*60*SAMP); // 生データ（50 kHz）のzabs出力用変数
						
						if(optbelp==1){
							belp = calloc(60*SAMP,sizeof(struct ellip_param)); //生データ（50 kHz）の楕円パラメータ出力用構造体
						}
						
						bph = calloc((int)(60*SAMP)/(SAMP/fDS),sizeof(struct binary_ph)); // bph（平均値）出力用構造体
						//printf("size of bph:%ld\n",sizeof(bph));

						flg=calloc(60*SAMP/UNIT_LEN,sizeof(int)); // フィッティング良否フラグ。単位データ長（UNIT_LEN）毎に設定
						for(int i=0;i<60*SAMP/UNIT_LEN;i++){
							flg[i]=0;
						};
						
						elpara_blk=calloc(60*SAMP/UNIT_LEN,sizeof(struct ellip_param)); // 楕円パラメータ。単位データ長（50点）毎に設定
						
						for(istep=0;istep<Nstep;istep++){
							//printf("istep: %d\n",istep);
							MM=MMtab[istep]*UNIT_LEN; // フィッティングに用いるデータ点数=各ブロックに含まれるデータ点数
							//printf("MM: %d\n",MM);

							//printf("number of points for fitting:%d\n",MM);
							//printf("number of blocks at this step: %d\n",Nblock_step[istep]);
							
							//flg1=calloc(sizeof(int)*Nblock_step[istep]); // 処理中の段のフィッティング良否を入れるフラグ
							//printf("flg1 for istep=%d is ready. number of elements:%d\n", istep, Nblock_step[istep]);
							
							k=0; // フィット長切り替えの段階中で何番目のブロックを処理しているかを示すインデックス
							
							while(k < Nblock_step[istep]){
								/*
								if(istep>0){
									printf("current block:%d, parent block:%d flag[%d]:%d\n",k, k/(MMtab[istep-1]/MMtab[istep]),k/(MMtab[istep-1]/MMtab[istep]),flg0[k/(MMtab[istep-1]/MMtab[istep])]);
								}else{
									printf("current block:%d, parent block:na flag[na]:na\n",k);
								};
								*/
								
								flgblk=0;
								
								for(int i=k*MMtab[istep];i<(k+1)*MMtab[istep];i++){
									//printf("flg[%d]:%d\n",i,flg[i]);
									flgblk+=flg[i]; //現在処理中のブロックに含まれるフラグの和（すべて正常にフィットできていれば、ブロックに含まれる単位データ長の個数に一致する）
									//printf("flg for current block:%d\n",flgblk);
								};
								
								//printf("flag check: done\n");
								
								//flg1=calloc(sizeof(int)*Nblock_step[istep]); // 処理中の段のフィッティング良否を入れるフラグ
								if(k==0){
									//printf("istep:%d, MM:%d, flgblk:%d, MMtab[istep]:%d,\n",istep,MM,flgblk,MMtab[istep]);
									printf("step:%d, fit length:%3.1f msec\n",istep,(float)MMtab[istep]*UNIT_LEN/SAMP*1000.);
									
								}
																
								if((istep==0)||(flgblk!=MMtab[istep])){ // 処理中の段が最初の段にあたる or ひとつ上の段の対応するブロックが不良の場合のみフィッティングを行う。
									//printf("execute fitting\n");
									// 単純にMM点で区切った区間を使ってフィッティング（端付近などの違いは考慮しない）
									// フィッティング関数に渡すデータ点を決定
									
									idfit0=calloc(MM,sizeof(int));
									idfit1=calloc(MM,sizeof(int));
									idfit2=calloc(MM,sizeof(int));
									idfit3=calloc(MM,sizeof(int));
					
									for(int j=0;j<MM;j++){
										idfit0[j]=id0[k*MM+j];
										idfit1[j]=id1[k*MM+j];
										idfit2[j]=id2[k*MM+j];
										idfit3[j]=id3[k*MM+j];
									}
									
									fz_res=calloc(MM,sizeof(struct fringe_zabs)); // ブロックのフィッティング結果（fringe, zabs）を入れる変数

									//pnt=k-(int)floor(k/MM)*MM;
									
									elliptic_fit(idfit0, idfit1, idfit2, idfit3, MM, &elpara); //フィッティング実行
									

									//printf("elpara.b00: %f\n",elpara.b00);
									//printf("fitting done\n");
									fringe_det(idfit0, idfit1, idfit2, idfit3, MM, elpara, &fz_res, &evalpara); //位相決定
									// for debugging
									//printf("fringe: %f\n",fz_res[MM-1].fringe);
									//printf("zabsmax: %f\n",evalpara.zabsmax);
									//printf("zabsmin: %f\n",evalpara.zabsmin);
									
									
										// 処理中ブロックのフィッティング結果の良否判定
										//printf("zabsmax:%f\n",evalpara.zabsmax);
										//printf("zabsmin:%f\n",evalpara.zabsmin);
										if(isnan(zabsmax)){printf("nan !!\n"); exit(0);};
										
										/*
										if(evalpara.absdaa>PI/4.){
											printf("absdaa:%lf at block no. %d\n",evalpara.absdaa,k);
										}
										*/
										
										if((evalpara.zabsmax<th_zabsmax)&&(evalpara.zabsmin>=th_zabsmin)){//&&(evalpara.dfringemax<th_dphmax)&&(evalpara.dfringemin>=th_dphmin)){//&&(evalpara.zabsmax>0)&&(evalpara.absdaa<PI/4.)){
											//printf("good fit!\n");
											for(int i=k*MMtab[istep];i<(k+1)*MMtab[istep];i++){
												elpara_blk[i]=elpara;//楕円パラメータの書き込み（フィッティングの良否にかかわらず実行）
												flg[i]=1; //現在処理中のブロックに含まれるフラグに1を書き込む（良好）
											};
											
											for(int j=0;j<MM;j++){
												fz_res_out[k*MM+j].numpnts=MM;
												fz_res_out[k*MM+j].fringe=fz_res[j].fringe;
												fz_res_out[k*MM+j].zabs=fz_res[j].zabs;
												
												//fringe_out[k*MM+j]=fz_res[j].fringe;
												//zabs_out[k*MM+j]=fz_res[j].zabs;
												
												// binary elliptic parameter data set for saving
												if(optbelp==1){
													belp[k*MM+j].MM=elpara.MM;
													belp[k*MM+j].a00=elpara.a00;
													belp[k*MM+j].b00=elpara.b00;
													belp[k*MM+j].x0=elpara.x0;
													belp[k*MM+j].y0=elpara.y0;
													belp[k*MM+j].sp0=elpara.sp0;
													belp[k*MM+j].cp0=elpara.cp0;
												}
									//			printf("index: %d\n", k%M);
												
											}
										}else{
											for(int i=k*MMtab[istep];i<(k+1)*MMtab[istep];i++){
												elpara_blk[i]=elpara;//楕円パラメータの書き込み（フィッティングの良否にかかわらず実行）
												flg[i]=0; //現在処理中のブロックに含まれるフラグに0を書き込む（不良）
												//printf("flg=0\n");
											};
											
											// 最後の段階でも良好なフィット結果が得られなかった場合の書き出し（段階番号を-1として書き出す）。
											if(istep==Nstep-1){
												for(j=0;j<MM;j++){
												//printf("fringe_out[k*MM]=%f\n",fz_res[j].fringe);
												fz_res_out[k*MM+j].numpnts=-1;
												fz_res_out[k*MM+j].fringe=fz_res[j].fringe;
												fz_res_out[k*MM+j].zabs=fz_res[j].zabs;
												
												// binary elliptic parameter data set for saving
												if(optbelp==1){
													belp[k*MM+j].MM=elpara.MM;
													belp[k*MM+j].a00=elpara.a00;
													belp[k*MM+j].b00=elpara.b00;
													belp[k*MM+j].x0=elpara.x0;
													belp[k*MM+j].y0=elpara.y0;
													belp[k*MM+j].sp0=elpara.sp0;
													belp[k*MM+j].cp0=elpara.cp0;
												}
												
												}
											}
										}
									//printf("k:%d, flag for current block:%d\n",k, flg1[k]);
						
									free(idfit0); free(idfit1); free(idfit2); free(idfit3);
									free(fz_res);
																		
									
								}else{
									//flg1[k]=2; // 前段までで良好なフィット結果が得られている場合、フィッティングは行わない（フラグ=2）
									//printf("no need to fit\n");
								}
												
							k++;
							//if(istep>0&&k>20)exit(0); // for debug
							}
							/*
							//フィット良否判定フラグの引き継ぎ
							free(flg0);
							flg0=calloc(sizeof(int)*Nblock_step[istep]);
							for(j=0;j<Nblock_step[istep];j++){
								flg0[j]=flg1[j];
								//printf("substitute flg0[%d]=flg1[%d]:%d\n",j,j,flg1[j]);
							};
							free(flg1);
							//printf("flg1 freed\n");
							*/
						//}
#if interpolate			
							
							//////////////////////////////////////////////////////////////////
							//////////////////////////////////////////////////////////////////
							//////////////////////////////////////////////////////////////////
							//補間処理
							//すべてのステップで補間を行う
							//先頭からフラグをチェックして、ゼロ（フィット不良）が現れたら補間処理に入る
							// 使用変数
							// int flg_int; //補間中か否かを区別するフラグ（0:補間していない, 1:補間処理中）
							// int int_s, int_e; //補間開始&終了ブロック番号
							// struct ellip_param ellp_s, ellp_e; //補間開始&終了ブロックの楕円パラメータ
							// struct ellip_param ellp_int; //補間点の楕円パラメータ
							//////////////////////////////////////////////////////////////////
							//////////////////////////////////////////////////////////////////
							//////////////////////////////////////////////////////////////////
							//printf("MM:%d\n",MM);
							//printf("Nblock_step:%d\n",Nblock_step[Nstep-1]);
							if((optint==1)&&((n!=0)&&(n1!=0))){
								printf("attempting interpolation...\n");
								printf("istep:%d\n",istep);
								printf("Nblock_step[istep]:%d\n",Nblock_step[istep]);
								printf("MM:%d\n",MM);
								
								flg_int=0;
								int_s=0;
								int_e=Nblock_step[istep];
								
								fz_res=calloc(MM,sizeof(struct fringe_zabs)); // ブロックのフィッティング結果（fringe, zabs）を入れる変数（再定義）
								
								for(k=0;k < Nblock_step[istep];k++){
									/*
									if(flg[k]==0){
										printf("flg[%d]:%d\n",k*MM,flg[k*MM]);
									}
									*/
									if(flg_int==0&&flg[k*MM/UNIT_LEN]!=1){
										printf("found a bad fit\n");
										printf("block no. (k):%d\n",k);
										flg_int=1;
										int_s=k;
										ellp_s=elpara_blk[k*MM/UNIT_LEN-1];
										printf("start of a00:%lf\n",ellp_s.a00);
									}
									
									if(flg_int==1&&flg[k*MM/UNIT_LEN]!=0){
										int_e=k;
										printf("found a good fit\n");
										//printf("k:%d, k*MM:%d, flg[k*MM/UNIT_LEN]:%d\n",k, k*MM, flg[k*MM/UNIT_LEN]);
										printf("block no. (k):%d\n",k);
										printf("no. of bad blocks:%d\n",int_e-int_s);
										ellp_e=elpara_blk[k*MM/UNIT_LEN];
										printf("end of a00:%lf\n",ellp_e.a00);
										
										if(int_s*MM/UNIT_LEN>1){//1分単位で先頭の区間については補間を行わない（行えない）
											printf("executing interpolation\n");
											for(int i=int_s;i<int_e;i++){
												//楕円補間
												ellp_int.a00=(ellp_e.a00-ellp_s.a00)/(int_e-int_s+1)*(i-int_s+1)+ellp_s.a00;
												ellp_int.b00=(ellp_e.b00-ellp_s.b00)/(int_e-int_s+1)*(i-int_s+1)+ellp_s.b00;
												ellp_int.x0=(ellp_e.x0-ellp_s.x0)/(int_e-int_s+1)*(i-int_s+1)+ellp_s.x0;
												ellp_int.y0=(ellp_e.y0-ellp_s.y0)/(int_e-int_s+1)*(i-int_s+1)+ellp_s.y0;
												ellp_int.sp0=(ellp_e.sp0-ellp_s.sp0)/(int_e-int_s+1)*(i-int_s+1)+ellp_s.sp0;
												ellp_int.cp0=(ellp_e.cp0-ellp_s.cp0)/(int_e-int_s+1)*(i-int_s+1)+ellp_s.cp0;
												
												printf("i:%d, ellp_int.a00:%lf\n",i,ellp_int.a00);
												
												//補間楕円を用いた位相決定
												idfit0=calloc(MM,sizeof(int));
												idfit1=calloc(MM,sizeof(int));
												idfit2=calloc(MM,sizeof(int));
												idfit3=calloc(MM,sizeof(int));
								
												for(int j=0;j<MM;j++){
													idfit0[j]=id0[i*MM+j];
													idfit1[j]=id1[i*MM+j];
													idfit2[j]=id2[i*MM+j];
													idfit3[j]=id3[i*MM+j];
												}
												//printf("first point of idfit0:%d\n",i*MM);
												
												fringe_det(idfit0, idfit1, idfit2, idfit3, MM, ellp_int, &fz_res, &evalpara_int); //位相決定

												///////////////////良否判定をフラグに、楕円パラメータを構造体に書き込む/////////////////////////////
												if((evalpara_int.zabsmax<th_zabsmax)&&(evalpara_int.zabsmin>=th_zabsmin)){
													for(int j=i*MM/UNIT_LEN;j<(i+1)*MM/UNIT_LEN;j++){
														elpara_blk[j]=ellp_int;//楕円パラメータの書き込み（フィッティングの良否にかかわらず実行）
														flg[j]=1; //現在処理中のブロックに含まれるフラグに1を書き込む（良）
													}
												}else{
													for(int j=i*MM/UNIT_LEN;j<(i+1)*MM/UNIT_LEN;j++){
														elpara_blk[j]=ellp_int;//楕円パラメータの書き込み（フィッティングの良否にかかわらず実行）
														flg[j]=0; //現在処理中のブロックに含まれるフラグに0を書き込む（不良）
														printf("MM:%d\n",MM);
														printf("i*MM/UNIT_LEN:%d\n", i*MM/UNIT_LEN);
														printf("(i+1)*MM/UNIT_LEN:%d\n", (i+1)*MM/UNIT_LEN);
														//printf("exit!\n");exit(0);
														//printf("bad interpolation\n");
														//printf("istep:%d, block:%d\n",istep,j);
													}
												}	
											
											
											
												///////////////////50 kHzデータの代入/////////////////////////////
												for(int j=0;j<MM;j++){
													fz_res_out[i*MM+j].numpnts=-MM;//補間によって決まった場合は負の値にする
													fz_res_out[i*MM+j].fringe=fz_res[j].fringe;
													fz_res_out[i*MM+j].zabs=fz_res[j].zabs;
													//printf("fz_res[%d].fringe:%f\n",j,fz_res[j].fringe);
													//printf("fz_res_out[%d].fringe:%f\n",i*MM+j,fz_res_out[i*MM+j].fringe);
													//printf("fz_res_out[%d].zabs:%f\n",i*MM+j,fz_res[j].zabs);
													
													//fringe_out[k*MM+j]=fz_res[j].fringe;
													//zabs_out[k*MM+j]=fz_res[j].zabs;
													
													// binary elliptic parameter data set for saving
													if(optbelp==1){
														belp[i*MM+j].MM=-MM;//補間によって決まった場合は負の値にする
														belp[i*MM+j].a00=ellp_int.a00;
														belp[i*MM+j].b00=ellp_int.b00;
														belp[i*MM+j].x0=ellp_int.x0;
														belp[i*MM+j].y0=ellp_int.y0;
														belp[i*MM+j].sp0=ellp_int.sp0;
														belp[i*MM+j].cp0=ellp_int.cp0;
													}
												}
												free(idfit0);free(idfit1);free(idfit2);free(idfit3);
											}
										}
										flg_int=0;
										
									}
								}
							} // 補間処理終わり
#endif


						} // 1分ぶんの位相決定終わり
						
						
						// 位相の2PI指数決め
						for(j=0;j<60*SAMP;j++){
							xx=fz_res_out[j].fringe;
							if((n==0)&&(n1==0)){	//最初のデータ－＞ll0=0, ll=0で変換
								ll0=0;
								ll=0;
								xx0=fz_res_out[0].fringe;
								fringe[j]=fz_res_out[0].fringe;
								fringe0=fz_res_out[0].fringe;
								dfringe=0.0;
								n++;
							}
							else{
								/*
								if((j+1)%SAMP==0){	//正秒
									dl=ll-ll0;	//前の1s間のΔl
									ll0=ll;	//1サンプル前のl
									ll=0;
								}
								*/
								
								if(xx-xx0>PI){
									ll--;
								}else{
									if(xx-xx0<-PI){
									ll++;
									}
								}
							
								xx0=xx;
								fringe[j]=fz_res_out[j].fringe+2.*PI*(double)ll;	//lは現正秒基準
								if((j+1)%SAMP==0) fringe0-=2.*PI*(double)dl;	//正秒の場合のfringe0は1s前の正秒基準のためdlで補正
								
								if(j==0){
									dfringe=fringe[0]-fringe0;
								}else{
									dfringe=fringe[j]-fringe[j-1];
								}
								
								//fringe0=fringe[j];
								if(dphmax<dfringe) dphmax=dfringe;
								if(dphmin>=dfringe) dphmin=dfringe;
								//if(j%1000==0)printf("j:%d, dphmax:%f, dphmin:%f \n",j,dphmax,dphmin);
								
								n++;
								if(n==0) n1++;
								//printf("j:%d, fringe[%d]:%lf, ll:%d, dfringemax:%f, dfringemin:%f \n",j,j,fringe[j],ll,dfringemax,dfringemin);
							}
						//}
								

								//位相のダウンサンプリング（デフォルト20 Hz）
								s0+=fringe[j];
								
								if(j%(SAMP/fDS)==(SAMP/fDS)-1){
									s0/=(double)(SAMP/fDS);
									//bph[(int)j/(SAMP/fDS)].l=(long)floor((s0+PI)/2./PI);
									bph[(int)j/(SAMP/fDS)].l=(long)fmod((s0+PI),(2.*PI));
									bph[(int)j/(SAMP/fDS)].x=(double)(s0-2.*PI*(double)bph[(int)j/(SAMP/fDS)].l);
									bph[(int)j/(SAMP/fDS)].l+=ll0;
									s0=0.0;
									//printf("bph[%d].x=%f, bph[%d].l=%d\n",j/(SAMP/fDS),bph[j/(SAMP/fDS)].x, j/(SAMP/fDS),bph[j/(SAMP/fDS)].l);
									//fprintf(fraw,"%f, %d\n",bph[j/(SAMP/fDS)].x,bph[j/(SAMP/fDS)].l);
								}
								
								//
								// 1秒毎の処理
								//
								
								s00+=fringe[j];//位相の1秒平均値準備
								/*
								if(j==0){
									fringe0=fringe[0];
									dfringe=0.;
								}else{
									dfringe=fringe[j]-fringe0;
									fringe0=fringe[j];
								}
								*/
								
								if(zabsmax<fz_res_out[j].zabs) zabsmax=fz_res_out[j].zabs;
								if(zabsmin>fz_res_out[j].zabs) zabsmin=fz_res_out[j].zabs;
								if(fringemax<fringe[j]) fringemax=fringe[j];
								if(fringemin>fringe[j]) fringemin=fringe[j];
								//if(dfringemax<dfringe) dfringemax=dfringe;
								//if(dfringemin>=dfringe) dfringemin=dfringe;
								if(pinmax<(double)id2[j]*ADCFAC) pinmax=(double)id2[j]*ADCFAC;
								if(pinmin>(double)id2[j]*ADCFAC) pinmin=(double)id2[j]*ADCFAC;
								if(amdmax<(double)id3[j]*ADCFAC) amdmax=(double)id3[j]*ADCFAC;
								if(amdmin>(double)id3[j]*ADCFAC) amdmin=(double)id3[j]*ADCFAC;
								
								
								/*
								//modified on 2025/6/4
								if(dphmax>=DTH*2.*PI) dphmax-=2.*PI;
								if(dphmin<=-DTH*2.*PI) dphmin+=2.*PI;
								*/
									
								if((j+1)%SAMP==0){
									isec=(int)j/SAMP; // 処理中のデータの秒数（hh:mm:ssのss部分）
									s00/=(double)SAMP;
									ll1=(int)floor((s00+PI)/2./PI);
									xx1=(float)(s00-2.*PI*(double)ll1);
									llmax=(int)floor((fringemax+PI)/2./PI);
									xxmax=(float)(fringemax-2.*PI*(double)llmax);
									llmin=(int)floor((fringemin+PI)/2./PI);
									xxmin=(float)(fringemin-2.*PI*(double)llmin);
									
									printf("dphmax:%f, dphmin:%f \n",dphmax,dphmin);
									
									
									//fprintf(fout,"%d %d %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",zabsmax, zabsmin, )
									fprintf(fout,"%d %d %d %d %d %d %f %f %f %d %f %d %f %d %f %f %f %f %f %f\n",\
									2000+iy,im,id,ih,imin,isec,(float)zabsmax,(float)zabsmin,(float)xx1,ll1+ll0,(float)xxmax,llmax+ll0,(float)xxmin,llmin+ll0,\
									(float)dphmax,(float)dphmin,(float)pinmax,(float)pinmin,(float)amdmax,(float)amdmin);
									/*
									if((dfringemax>=DTH*2.*PI)||(dfringemin<=-DTH*2.*PI)) ndph=1;
									if(dphmax<dfringemax) dphmax=dfringemax;
									if(dphmax<-dfringemin) dphmax=-dfringemin;
									*/
									
									
									ll+=ll0;
									zabsmax=-1.e20;
									zabsmin=1.e20;
									s00=0.0;
									fringemax=-1.e20;
									fringemin=1.e20;
									dphmax=-1.e20;
									dphmin=1.e20;
									pinmax=-1.e20;
									pinmin=1.e20;
									amdmax=-1.e20;
									amdmin=1.e20;
								}
								
								
								//printf("istep = %d\n", fz_res_out[j].numpnts);
								//printf("fringe[%d] = %f\n", j, fringe[j]);
								//printf("zabs[%d] = %f\n", j, fz_res_out[j].zabs);
								
								//50 kHz (ascii)データ書き出し: フィットに用いたデータ点数, 位相 [rad], zabs
								if(opt50ka==1) fprintf(fraw,"%d, %lf, %lf\n",fz_res_out[j].numpnts,fringe[j],fz_res_out[j].zabs);
								//if(isnan(fringe_out[j]))exit(0);
								
								//50 kHz (binary)データ書き出し: フィットに用いたデータ点数(.bfl), 位相 [rad](.bfr), zabs (.bzabs)
								if(opt50kb==1){
									fwrite(&fz_res_out[j].numpnts,sizeof(int),1,fflen);
									fwrite(&fringe[j],sizeof(double),1,ffringe);
									fwrite(&fz_res_out[j].zabs,sizeof(double),1,fzabs);
								}
								
								
							}
							
							// 楕円パラメータのバイナリファイル書き出し（50 kHz）
							if(optbelp==1){
								if(fwrite(belp,sizeof(struct ellip_param),60*SAMP,fbelp)!=60*SAMP){
									puts("data write error (.belp) !!");
									return(-1);
								}
							}

							
														
							// 生データ書き出し（50 kHz）： フィット使用点数, fringe[rad], zabs
							//if(opt50ka==1)fprintf(fraw,"%d, %f, %f\n",fz_res_out[j].numpnts,fringe[j],fz_res_out[j].zabs);
						//}
							// 1分間分のbph ダウンサンプルデータ書き出し（デフォルト20 Hz）
							if(fwrite(bph,sizeof(struct binary_ph),(int)(60*SAMP)/(SAMP/fDS),fbph)!=(int)(60*SAMP)/(SAMP/fDS)){
							puts("data write error (.bph) !!");
							return(-1);
							}
							
						
						
						
						//1sファイル書き出し
						//フォーマット："# year month day hour minute second zabsmax zabsmin phase_1s.x phase_1s.l phasemax.x phasemax.l phasemin.x phasemin.l\
						dphasemax dphasemin pinmax pinmin amdmax amdmin\n"
							/*
							for(j=0;j<60*SAMP;j++){
								/*
								isec=(int)j/SAMP; // 処理中のデータの秒数（hh:mm:ssのss部分）
								s0+=fringe[j];
								//s00+=fringe[j];
								if(zabsmax<fz_res_out[j].zabs) zabsmax=fz_res_out[j].zabs;
								if(zabsmin>fz_res_out[j].zabs) zabsmin=fz_res_out[j].zabs;
								if(fringemax<fringe[j]) fringemax=fringe[j];
								if(fringemin>fringe[j]) fringemin=fringe[j];
								if(pinmax<(double)id2[j]*ADCFAC) pinmax=(double)id2[j]*ADCFAC;
								if(pinmin>(double)id2[j]*ADCFAC) pinmin=(double)id2[j]*ADCFAC;
								if(amdmax<(double)id3[j]*ADCFAC) amdmax=(double)id3[j]*ADCFAC;
								if(amdmin>(double)id3[j]*ADCFAC) amdmin=(double)id3[j]*ADCFAC;
								
								if(j%SAMP==0){
									//fprintf(fout,"%d %d %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",zabsmax, zabsmin, )
									fprintf(fout,"%d %d %d %d %d %d %f %f %f %d %f %d %f %d %f %f %f %f %f %f\n",\
									2000+iy,im,id,ih,imin,isec,\
									(double)zabsmax,(float)zabsmin,(float)xx1,ll1+ll0,(float)xxmax,llmax+ll0,(float)xxmin,llmin+ll0,\
									(float)dfringemax,(float)dfringemin,(float)pinmax,(float)pinmin,(float)amdmax,(float)amdmin);
								}
								//printf("istep = %d\n", fz_res_out[j].step);
								//printf("fringe[%d] = %f\n", j, fringe[j]);
								//printf("zabs[%d] = %f\n", j, fz_res_out[j].zabs);
								fprintf(fraw,"%d, %f, %f\n",fz_res_out[j].numpnts,fringe[j],fz_res_out[j].zabs);
								//if(isnan(fringe_out[j]))exit(0);
							}
							
							/*
								s0/=(double)(MM/N); //50 ms平均
								bph[j%(MM/N)].l=(long)floor((s0+PI)/2./PI);
								bph[j%(MM/N)].x=(float)(s0-2.*PI*(double)bph[j].l);
								bph[j%(MM/N)].l+=ll0;
								printf("j%(MM/N)=%d\n",j%(MM/N));
								//printf("bph[%d].x=%f, bph[%d].l=%f\n",(int)j%(MM/n), bph[j%(MM/N)].x, (int)j%(MM/n), bph[j%(MM/N)].l);
								s0=0.0;
							*/
					
						
				// 終了処理		
						//free(fringe); free(zabs); free(pin); free(amd); free(a00); free(b00); free(x0); free(y0); free(aa); free(MMp);
						//free(fringe_out); free(zabs_out); free(pin_out); free(amd_out); free(a00_out); free(b00_out); free(x0_out); free(y0_out); free(aa_out); free(MMp_out);
						//free(fz_res);
						if(fclose(fin0)==EOF){
							puts("file close error !! (.AD00)");
							return(-1);
						}
						if(fclose(fin1)==EOF){
							puts("file close error !! (.AD01)");
							return(-1);
						}
						if(fclose(fin2)==EOF){
							puts("file close error !! (.AD02)");
							return(-1);
						}
						if(fclose(fin3)==EOF){
							puts("file close error !! (.AD03)");
							return(-1);
						}

						free(id0); free(id1); free(id2); free(id3);
						free(bph);
						free(fz_res_out);
						if(optbelp==1) free(belp);
						//free(zabs_out);
						free(fringe); free(flg);

					} // minute count
					


				} // hour count
			} // day count

			if(fclose(fbph)==EOF){
				puts("file close error !! (.bph)");
				return(-1);
			}

			if(optbelp==1){
				if(fclose(fbelp)==EOF){
					puts("file close error !! (.belp)");
					return(-1);
				}
			}

			if(fclose(fzdph)==EOF){
				puts("file close error !! (.zdph)");
				return(-1);
			}
			
			if(fclose(fdiag)==EOF){
				puts("file close error !! (.diag)");
				return(-1);
			}

			if(opt50ka==1){
				if(fclose(fraw)==EOF){
					puts("file close error !! (.dat)");
					return(-1);
				}
			}
			
			if(opt50kb==1){
				if(fclose(fflen)==EOF){
					puts("file close error !! (.bfl)");
					return(-1);
				}
				
				if(fclose(ffringe)==EOF){
					puts("file close error !! (.bfr)");
					return(-1);
				}
				
				if(fclose(fzabs)==EOF){
					puts("file close error !! (.bzabs)");
					return(-1);
				}

			}
										
						

		} // month count
	} // year count

	//print end time
	time_t t_e=time(NULL);
	printf("\nend time: %s",ctime(&t_e));
	printf("process time: %ld second(s)\n",t_e-t_s);
		
	exit(0);
}



/*
//データの書き出し
							//if((k+1)%SAMP==0){	//1s間の最終処理。元のプログラムを改変（元はSAMP->M）
							if((k+1)%SAMP==0){	//1s間の最終処理。元のプログラムを改変（元はSAMP->M）
//								if(fwrite(bph,sizeof(struct binary_ph),SAMP,fbph)!=SAMP){ //平均しないので、元のプログラムを改変（元はSAMP->MM）
//									puts("data write error (.bph) !!");
//									return(-1);
//								}
								s00/=(double)(MM*M/N);
								ll1=(int)floor((s00+PI)/2./PI);
								xx1=(float)(s00-2.*PI*(double)ll1);
								llmax=(int)floor((fringemax+PI)/2./PI);
								xxmax=(float)(fringemax-2.*PI*(double)llmax);
								llmin=(int)floor((fringemin+PI)/2./PI);
								xxmin=(float)(fringemin-2.*PI*(double)llmin);
								// # zabsmax zabsmin phase_1s.x phase_1s.l phasemax.x phasemax.l phasemin.x phasemin.l dphasemax dphasemin pinmax pinmin amdmax amdmin x0max x0min y0max y0min a00max a00min b00max b00min aamax aamin 
					//			if((n==0)&&(n1==0)) fprintf(fout,"# zabsmax zabxmin phase_1s.x phase_1s.l phasemax.x phasemax.l phasemin.x phasemin.l dphasemax dphasemin pinmax pinmin amdmax amdmin x0max x0min y0max y0min a00max a00min b00max b00min aamax aamin \n");
					//			fprintf(fout,"%25.16e %25.16e %25.16e %d %25.16e %d %25.16e %d %25.16e %25.16e %25.16e %25.16e %25.16e %25.16e %25.16e %25.16e %25.16e %25.16e %25.16e %25.16e %25.16e %25.16e %25.16e %25.16e %25.16e %25.16e \n",(float)zabsmax,(float)zabsmin,(float)xx,ll+ll0,(float)xxmax,llmax+ll0,(float)xxmin,llmin+ll0,(float)dfringemax,(float)dfringemin,(float)pinmax,(float)pinmin,(float)amdmax,(float)amdmin,(float)x0max,(float)x0min,(float)y0max,(float)y0min,(float)a00max,(float)a00min,(float)b00max,(float)b00min,(float)aamax,(float)aamin);
								fprintf(fout,"%d %d %d %d %d %d %e %e %e %d %e %d %e %d %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n",2000+iy,im,id,ih,imin,k/M,(float)zabsmax,(float)zabsmin,(float)xx1,ll1+ll0,(float)xxmax,llmax+ll0,(float)xxmin,llmin+ll0,(float)dfringemax,(float)dfringemin,(float)pinmax,(float)pinmin,(float)amdmax,(float)amdmin,(float)x0max,(float)x0min,(float)y0max,(float)y0min,(float)a00max,(float)a00min,(float)b00max,(float)b00min,(float)aamax,(float)aamin);
								if((dfringemax>=DTH*2.*PI)||(dfringemin<=-DTH*2.*PI)) ndph=1;
								if(dphmax<dfringemax) dphmax=dfringemax;
								if(dphmax<-dfringemin) dphmax=-dfringemin;
								ll+=ll0;
								zabsmax=-1.e20;
								zabsmin=1.e20;
								s00=0.0;
								fringemax=-1.e20;
								fringemin=1.e20;
								dfringemax=-1.e20;
								dfringemin=1.e20;
								pinmax=-1.e20;
								pinmin=1.e20;
								amdmax=-1.e20;
								amdmin=1.e20;
								x0max=-1.e20;
								x0min=1.e20;
								y0max=-1.e20;
								y0min=1.e20;
								a00max=-1.e20;
								a00min=1.e20;
								b00max=-1.e20;
								b00min=1.e20;
								aamax=-1.e20;
								aamin=1.e20;
							}
						
*/