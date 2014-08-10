//tn_fnds v0.0.4   2012/3/10
//追加されているコメントには誤りがあるかもしれません。
#include <windows.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "wavread.h"

#pragma warning(disable:4996)

/* wavread関数の移植 */
double * wavread(char* filename, int *fs, int *Nbit, int *waveLength, int *offset, int *endbr)
{
	FILE *fp;
	char dataCheck[5]; // 少し多めに
	unsigned char forIntNumber[4];
	double tmp, signBias, zeroLine;
	short int channel;
	int quantizationByte;
	double *waveForm;

	dataCheck[4] = '\0'; // 文字列照合のため，最後に終了文字を入れる．
//	fp = fopen(filename, "rb");
	fp = fopen(filename, "rb");
	if(NULL == fp) 
	{
		printf("ファイルのロードに失敗\n");
		return NULL;
	}
	//ヘッダのチェック
size_t result =
	fread(dataCheck, sizeof(char), 4, fp); // "RIFF"
assert(result == 4);
	if(0 != strcmp(dataCheck,"RIFF"))
	{
		fclose(fp);
		printf("ヘッダRIFFが不正\n");
		return NULL;
	}
	fseek(fp, 4, SEEK_CUR); // 4バイト飛ばす
result =
	fread(dataCheck, sizeof(char), 4, fp); // "WAVE"
assert(result == 4);
	if(0 != strcmp(dataCheck,"WAVE"))
	{
		fclose(fp);
		printf("ヘッダWAVEが不正\n");
		return NULL;
	}
result =
	fread(dataCheck, sizeof(char), 4, fp); // "fmt "
assert(result == 4);
	if(0 != strcmp(dataCheck,"fmt "))
	{
		fclose(fp);
		printf("ヘッダfmt が不正\n");
		return NULL;
	}
result =
	fread(dataCheck, sizeof(char), 4, fp); //1 0 0 0
assert(result == 4);
	if(!(16 == dataCheck[0] && 0 == dataCheck[1] && 0 == dataCheck[2] && 0 == dataCheck[3]))
	{
		fclose(fp);
		printf("ヘッダfmt (2)が不正\n");
		return NULL;
	}
result =
	fread(dataCheck, sizeof(char), 2, fp); //1 0
assert(result == 2);
	if(!(1 == dataCheck[0] && 0 == dataCheck[1]))
	{
		fclose(fp);
		printf("フォーマットIDが不正\n");
		return NULL;
	}
	/*
result =
	fread(dataCheck, sizeof(char), 2, fp); //1 0
assert(result == 2);
	if(!(1 == dataCheck[0] && 0 == dataCheck[1]))
	{
		fclose(fp);
		printf("ステレオには対応していません\n");
		return NULL;
	}
	*/
	//チャンネル
result =
	fread(&channel, sizeof(short int), 1, fp); 
assert(result == 1);

	// サンプリング周波数
result =
	fread(forIntNumber, sizeof(char), 4, fp);
assert(result == 4);
	*fs = 0;
	for(int i = 3;i >= 0;i--)
	{
		*fs = *fs*256 + forIntNumber[i];
	}
	// 量子化ビット数
	fseek(fp, 6, SEEK_CUR); // 6バイト飛ばす
result =
	fread(forIntNumber, sizeof(char), 2, fp);
assert(result == 2);
	*Nbit = forIntNumber[0];
	// ヘッダ
	int dummy;
result =
	fread(dataCheck, sizeof(char), 4, fp); // "data"
assert(result == 4);
	while(0 != strcmp(dataCheck,"data"))
	{
result =
		fread(&dummy, sizeof(char), 4, fp);
assert(result == 4);
		fseek(fp, dummy, SEEK_CUR); // 無関係なチャンクを読み飛ばす
result =
		fread(dataCheck, sizeof(char), 4, fp); // "data"
assert(result == 4);
//		fclose(fp);
//		printf("ヘッダdataが不正\n");
//		return NULL;
	}
	// サンプル点の数
result =
	fread(forIntNumber, sizeof(char), 4, fp); // "data"
assert(result == 4);
	*waveLength = 0;
	for(int i = 3;i >= 0;i--)
	{
		*waveLength = *waveLength*256 + forIntNumber[i];
	}
	*waveLength /= (*Nbit/8 * channel);

	if(*endbr < 0) // 負の場合はoffsetからの距離
	{
		*endbr = (*waveLength * 1000 / *fs) - (*offset-*endbr);
	}

	int st, ed;
	st = max(0, min(*waveLength-1, (int)((*offset-100) * *fs / 1000))); 
	ed = max(0, min(*waveLength-1, *waveLength - (int)(max(0, *endbr - 100) * *fs / 1000)));
	*endbr = (ed*1000 / *fs) - ((*waveLength * 1000 / *fs) - *endbr);
	*offset = *offset - (st*1000 / *fs);
	*waveLength = (ed - st + 1);

	// 波形を取り出す
	waveForm = (double *)malloc(sizeof(double) * *waveLength);
	if(waveForm == NULL) return NULL;

	quantizationByte = *Nbit/8;
	zeroLine = pow(2.0,*Nbit-1);
//	for(int i = 0;i < *waveLength;i++)

	fseek(fp, st * quantizationByte * channel, SEEK_CUR);  //スタート位置まで読み飛ばす

	unsigned char *wavbuff;
	wavbuff = (unsigned char *) malloc(sizeof(char) * *waveLength * quantizationByte * channel);
result =
	fread(wavbuff, sizeof(char), *waveLength * quantizationByte * channel, fp); // 全部メモリに読み込む
assert(result == *waveLength * quantizationByte * channel);
	int seekindex;

	for(int i = 0;i < *waveLength;i++)
	{
		seekindex = i * quantizationByte * channel;
		signBias = 0.0;
		tmp = 0.0;
		// 符号の確認
		if(wavbuff[seekindex + quantizationByte-1] >= 128)
		{
			signBias = pow(2.0,*Nbit-1);
			wavbuff[seekindex + quantizationByte-1] = wavbuff[seekindex + quantizationByte-1] & 0x7F;
		}
		// データの読み込み
		for(int j = quantizationByte-1;j >= 0;j--)
		{
			tmp = tmp*256.0 + (double)(wavbuff[seekindex + j]);
		}
		waveForm[i] = (double)((tmp - signBias) / zeroLine);

	}
	// 成功
	free(wavbuff);
	fclose(fp);
	return waveForm;
}

