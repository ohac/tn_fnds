#include "world.h"

#include <stdio.h> // for debug
#include <stdlib.h>

// spectrum, cepstrumは毎回malloc, freeするのが面倒だから．
/*
int getOneFrameSegment(double *f0, int tLen, double **specgram, double **aperiodicity, int fftl, double framePeriod, double currentTime, int fs, double defaultF0,
						fftw_complex *spectrum, fftw_complex *cepstrum, 
						double *response, int xLen);
*/
// 特定時刻の応答を取得する．
void getOneFrameSegment(double *f0, int tLen, double **specgram, int fftl, double **aperiodicity, int cLen, double targetF0, double framePeriod, double currentTime, int fs, double defaultF0,
						fftw_complex *spectrum, fftw_complex *cepstrum, 
						double *response, int xLen)
{
	int i, offset, noiseSize;
	double real, imag, tmp;
	fftw_plan	inverseFFT_RP;				// FFTセット
	fftw_plan	inverseFFT_RA;				// FFTセット
	fftw_plan	forwardFFT_R;				// FFTセット

	int currentFrame, currentPosition;
	double *periodicSpec, *aperiodicSpec, *aperiodicRatio;
	periodicSpec   = (double *)malloc(sizeof(double)* fftl);
	aperiodicSpec  = (double *)malloc(sizeof(double)* fftl);
	aperiodicRatio = (double *)malloc(sizeof(double)* fftl);

	fftw_complex *noiseSpec;
	noiseSpec   = (fftw_complex *)malloc(sizeof(fftw_complex)* fftl);

	double *aperiodicResponse, *periodicResponse;
	aperiodicResponse = (double *)malloc(sizeof(double)* fftl);
	periodicResponse  = (double *)malloc(sizeof(double)* fftl);

	forwardFFT_R = fftw_plan_dft_r2c_1d(fftl, aperiodicResponse, noiseSpec, FFTW_ESTIMATE);
	inverseFFT_RP = fftw_plan_dft_c2r_1d(fftl, spectrum, periodicResponse ,  FFTW_ESTIMATE);
	inverseFFT_RA = fftw_plan_dft_c2r_1d(fftl, spectrum, aperiodicResponse,  FFTW_ESTIMATE);

	currentFrame = (int)(currentTime/(framePeriod/1000.0) + 0.5);	
	currentPosition = (int)(currentTime*(double)fs);

	tmp = currentTime + 1.0/(f0[currentFrame] == 0.0 ? defaultF0 : f0[currentFrame]);
	noiseSize = (int)(tmp*(double)fs) - currentPosition;

	if(f0[currentFrame] == 0.0)
	{ // 無声音として扱う(波形の再利用)
		offset = currentPosition - (int)((double)currentFrame*framePeriod/1000.0*(double)fs);
		if(offset < 0)
		{
			currentFrame-=1;
			offset = currentPosition - (int)((double)currentFrame*framePeriod/1000.0*(double)fs);
		}
		for(i = 0;i < (int)((double)fs/defaultF0);i++)
		{
			if(i+currentPosition >= xLen || i+offset > fftl/2) break;
			response[i] = specgram[currentFrame][i+offset]; 
		}
	}
	else
	{ // 有声音として扱う

		// 非周期性指標の計算
		calculateAperiodicity(aperiodicity[currentFrame], cLen, fftl, f0[currentFrame], fs, targetF0, aperiodicRatio);
//		calculateAperiodicity(aperiodicity[currentFrame], fftl, fs, aperiodicRatio);

		// 値を取り出す
		for(i = 0;i <= fftl/2;i++)
		{
			periodicSpec[i] = specgram[currentFrame][i] * aperiodicRatio[i];
		}
		getMinimumPhaseSpectrum(periodicSpec, spectrum, cepstrum, fftl);

		fftw_execute(inverseFFT_RP);

		// ここから非周期音の合成
		for(i = 0;i <= fftl/2;i++)
		{
			aperiodicSpec[i] = specgram[currentFrame][i] * ((1-aperiodicRatio[i])+0.000000000000001);
		}
		getMinimumPhaseSpectrum(aperiodicSpec, spectrum, cepstrum, fftl);
		for(i = 0;i < noiseSize;i++)
			aperiodicResponse[i] = randn();
		for(;i < fftl;i++)
			aperiodicResponse[i] = 0.0;
		fftw_execute(forwardFFT_R);

		for(i = 0;i <= fftl/2;i++)
		{
			real = spectrum[i][0]*noiseSpec[i][0] - spectrum[i][1]*noiseSpec[i][1];
			imag = spectrum[i][0]*noiseSpec[i][1] + spectrum[i][1]*noiseSpec[i][0];
			spectrum[i][0] = real;
			spectrum[i][1] = imag;
		}
		fftw_execute(inverseFFT_RA);

		// 1.633はエネルギー調整用のマジックナンバー
		// 基本周期の３倍のハニング窓で波形を切り出した際のエネルギー損失を補償している．
		// *8.0は完全にマジックナンバー
		// ここは，将来リリースでの調整が必須です．
		for(i = 0;i < fftl;i++)
		{
			response[i] = (periodicResponse[i]*sqrt((double)noiseSize) + aperiodicResponse[i])*1.633/(double)fftl*12.0/sqrt((double)noiseSize);
//			response[i] = (periodicResponse[i]/sqrt((double)noiseSize))*1.633/(double)fftl*12.0;
//				aperiodicResponse[i]/(double)noiseSize;
//			response[i] = aperiodicResponse[i]/(double)fftl/sqrt((double)noiseSize);
		}

	}

	fftw_destroy_plan(inverseFFT_RP);
	fftw_destroy_plan(inverseFFT_RA);
	fftw_destroy_plan(forwardFFT_R);
	free(periodicSpec);
	free(aperiodicSpec);
	free(periodicResponse);
	free(aperiodicResponse);
	free(aperiodicRatio);
	free(noiseSpec);
}



void synthesis_ap(double *f0, int tLen, double **specgram, int fftl, double **aperiodicity, int cLen, double targetF0, double framePeriod, int fs, 
			   double *synthesisOut, int xLen)
{
	int i,j;
//	double defaultF0 = 300.0; // うーん．
	double defaultF0 = 160.0; // うーん．
	double *impulseResponse;
	impulseResponse = (double *)malloc(sizeof(double) * fftl);
	fftw_complex		*cepstrum, *spectrum;	// ケプストラムとスペクトル
	cepstrum = (fftw_complex *)malloc(sizeof(fftw_complex) * fftl);
	spectrum = (fftw_complex *)malloc(sizeof(fftw_complex) * fftl);

	double currentTime = 0.0;
	int currentPosition = 0;//currentTime / framePeriod;
	int currentFrame = 0;
	for(i = 0;;i++)
	{
		for(j = 0;j < fftl;j++) impulseResponse[j] = 0.0; // 配列は毎回初期化

		getOneFrameSegment(f0, tLen, specgram, fftl, aperiodicity, cLen, targetF0, framePeriod, currentTime, fs, defaultF0,
						spectrum, cepstrum, impulseResponse, xLen);

		currentPosition = (int)(currentTime*(double)fs);

		for(j = 0;j < fftl/2;j++)
		{
			if(j+currentPosition >= xLen) break;
			synthesisOut[j+currentPosition] += impulseResponse[j];
		}

		// 更新
		currentTime += 1.0/(f0[currentFrame] == 0.0 ? defaultF0 : f0[currentFrame]);
		currentFrame = (int)(currentTime/(framePeriod/1000.0) + 0.5);
		currentPosition = (int)(currentTime*(double)fs);
		if(j+currentPosition >= xLen || currentFrame >= tLen) break;
	}

	free(cepstrum); free(spectrum);
	free(impulseResponse);
	return;
}
