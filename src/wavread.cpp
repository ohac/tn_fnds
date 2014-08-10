//tn_fnds v0.0.4   2012/3/10
//�ǉ�����Ă���R�����g�ɂ͌�肪���邩������܂���B
#include <windows.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "wavread.h"

#pragma warning(disable:4996)

/* wavread�֐��̈ڐA */
double * wavread(char* filename, int *fs, int *Nbit, int *waveLength, int *offset, int *endbr)
{
	FILE *fp;
	char dataCheck[5]; // �������߂�
	unsigned char forIntNumber[4];
	double tmp, signBias, zeroLine;
	short int channel;
	int quantizationByte;
	double *waveForm;

	dataCheck[4] = '\0'; // ������ƍ��̂��߁C�Ō�ɏI������������D
//	fp = fopen(filename, "rb");
	fp = fopen(filename, "rb");
	if(NULL == fp) 
	{
		printf("�t�@�C���̃��[�h�Ɏ��s\n");
		return NULL;
	}
	//�w�b�_�̃`�F�b�N
	fread(dataCheck, sizeof(char), 4, fp); // "RIFF"
	if(0 != strcmp(dataCheck,"RIFF"))
	{
		fclose(fp);
		printf("�w�b�_RIFF���s��\n");
		return NULL;
	}
	fseek(fp, 4, SEEK_CUR); // 4�o�C�g��΂�
	fread(dataCheck, sizeof(char), 4, fp); // "WAVE"
	if(0 != strcmp(dataCheck,"WAVE"))
	{
		fclose(fp);
		printf("�w�b�_WAVE���s��\n");
		return NULL;
	}
	fread(dataCheck, sizeof(char), 4, fp); // "fmt "
	if(0 != strcmp(dataCheck,"fmt "))
	{
		fclose(fp);
		printf("�w�b�_fmt ���s��\n");
		return NULL;
	}
	fread(dataCheck, sizeof(char), 4, fp); //1 0 0 0
	if(!(16 == dataCheck[0] && 0 == dataCheck[1] && 0 == dataCheck[2] && 0 == dataCheck[3]))
	{
		fclose(fp);
		printf("�w�b�_fmt (2)���s��\n");
		return NULL;
	}
	fread(dataCheck, sizeof(char), 2, fp); //1 0
	if(!(1 == dataCheck[0] && 0 == dataCheck[1]))
	{
		fclose(fp);
		printf("�t�H�[�}�b�gID���s��\n");
		return NULL;
	}
	/*
	fread(dataCheck, sizeof(char), 2, fp); //1 0
	if(!(1 == dataCheck[0] && 0 == dataCheck[1]))
	{
		fclose(fp);
		printf("�X�e���I�ɂ͑Ή����Ă��܂���\n");
		return NULL;
	}
	*/
	//�`�����l��
	fread(&channel, sizeof(short int), 1, fp); 

	// �T���v�����O���g��
	fread(forIntNumber, sizeof(char), 4, fp);
	*fs = 0;
	for(int i = 3;i >= 0;i--)
	{
		*fs = *fs*256 + forIntNumber[i];
	}
	// �ʎq���r�b�g��
	fseek(fp, 6, SEEK_CUR); // 6�o�C�g��΂�
	fread(forIntNumber, sizeof(char), 2, fp);
	*Nbit = forIntNumber[0];
	// �w�b�_
	int dummy;
	fread(dataCheck, sizeof(char), 4, fp); // "data"
	while(0 != strcmp(dataCheck,"data"))
	{
		fread(&dummy, sizeof(char), 4, fp);
		fseek(fp, dummy, SEEK_CUR); // ���֌W�ȃ`�����N��ǂݔ�΂�
		fread(dataCheck, sizeof(char), 4, fp); // "data"
//		fclose(fp);
//		printf("�w�b�_data���s��\n");
//		return NULL;
	}
	// �T���v���_�̐�
	fread(forIntNumber, sizeof(char), 4, fp); // "data"
	*waveLength = 0;
	for(int i = 3;i >= 0;i--)
	{
		*waveLength = *waveLength*256 + forIntNumber[i];
	}
	*waveLength /= (*Nbit/8 * channel);

	if(*endbr < 0) // ���̏ꍇ��offset����̋���
	{
		*endbr = (*waveLength * 1000 / *fs) - (*offset-*endbr);
	}

	int st, ed;
	st = max(0, min(*waveLength-1, (int)((*offset-100) * *fs / 1000))); 
	ed = max(0, min(*waveLength-1, *waveLength - (int)(max(0, *endbr - 100) * *fs / 1000)));
	*endbr = (ed*1000 / *fs) - ((*waveLength * 1000 / *fs) - *endbr);
	*offset = *offset - (st*1000 / *fs);
	*waveLength = (ed - st + 1);

	// �g�`�����o��
	waveForm = (double *)malloc(sizeof(double) * *waveLength);
	if(waveForm == NULL) return NULL;

	quantizationByte = *Nbit/8;
	zeroLine = pow(2.0,*Nbit-1);
//	for(int i = 0;i < *waveLength;i++)

	fseek(fp, st * quantizationByte * channel, SEEK_CUR);  //�X�^�[�g�ʒu�܂œǂݔ�΂�

	unsigned char *wavbuff;
	wavbuff = (unsigned char *) malloc(sizeof(char) * *waveLength * quantizationByte * channel);
	fread(wavbuff, sizeof(char), *waveLength * quantizationByte * channel, fp); // �S���������ɓǂݍ���
	int seekindex;

	for(int i = 0;i < *waveLength;i++)
	{
		seekindex = i * quantizationByte * channel;
		signBias = 0.0;
		tmp = 0.0;
		// �����̊m�F
		if(wavbuff[seekindex + quantizationByte-1] >= 128)
		{
			signBias = pow(2.0,*Nbit-1);
			wavbuff[seekindex + quantizationByte-1] = wavbuff[seekindex + quantizationByte-1] & 0x7F;
		}
		// �f�[�^�̓ǂݍ���
		for(int j = quantizationByte-1;j >= 0;j--)
		{
			tmp = tmp*256.0 + (double)(wavbuff[seekindex + j]);
		}
		waveForm[i] = (double)((tmp - signBias) / zeroLine);

	}
	// ����
	free(wavbuff);
	fclose(fp);
	return waveForm;
}

