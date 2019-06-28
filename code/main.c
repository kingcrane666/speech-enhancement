#include "downsample.h"
#include "aec_srp_mvdr.h"
#define BYTE char
#define HEAD 44

extern frame;
void main(void)
{
	Word16 in_sp[NUM][FRM_LEN2], mic_sp[NUM][FRM_LEN], spk[FRM_LEN2], spk_sp[FRM_LEN], out_sp[6][FRM_LEN];
	Word32 i, flag, locate_mid_num, valid_mid_num, locate_moment_num, valid_moment_num;
	Float32 error_mid, error_moment;
	AEC_SRP_ST my_str;
    DOWNSAMPLE_STR down_str[NUM+1];
	BYTE wav_head[HEAD];

	clock_t start, finish, start_file, finish_file;
	double total_time;

	FILE *fq[NUM], *fr, *fs[6], *fp, *fn;

    for(i=0; i<NUM+1; i++)
    {
        down_sample_init(&down_str[i]);   
    }
	
    aec_srp_gsc_init(&my_str, 0.0f);

    fp = fopen("azimuth_mid2.txt", "w");
	fn = fopen("azimuth_moment.txt", "w");

	//fq[0]=fopen("0.pcm", "rb");
	//fq[1]=fopen("1.pcm", "rb");
	//fq[2]=fopen("2.pcm", "rb");
	//fq[3]=fopen("3.pcm", "rb");
	//fq[4]=fopen("4.pcm", "rb");
	//fq[5]=fopen("5.pcm", "rb");
	fq[0]=fopen("s001.ch01.wav", "rb");
	fq[1]=fopen("s001.ch02.wav", "rb");
	fq[2]=fopen("s001.ch03.wav", "rb");
	fq[3]=fopen("s001.ch04.wav", "rb");
	fq[4]=fopen("s001.ch05.wav", "rb");
	fq[5]=fopen("s001.ch06.wav", "rb");

	for(i=0; i<NUM; i++)
	{		
    	if(fq[i]==NULL) {printf("Cann't open the No.%ld microphone file!!!\n", i);  exit(0);}
	}

	fr=fopen("s001.ch07.wav", "rb");
	//fr=fopen("spk.pcm", "rb");
   	if(fr==NULL) {printf("Cann't open spk_sp.pcm file!!!\n");  exit(0);}

	fs[0] = fopen("out_0.pcm", "wb");
	fs[1] = fopen("out_1.pcm", "wb");
	fs[2] = fopen("out_2.pcm", "wb");
	fs[3] = fopen("out_3.pcm", "wb");
	fs[4] = fopen("out_4.pcm", "wb");
	fs[5] = fopen("out_5.pcm", "wb");
	for (i = 0; i < NUM; i++)
	{
		if (fs[i] == NULL) { printf("Cann't open the No.%ld out file!!!\n", i);  exit(0); }
	}
	
	for(i=0; i<NUM; i++)
	{
		fread(wav_head, sizeof(BYTE), HEAD, fq[i]);
	}
	fread(wav_head, sizeof(BYTE), HEAD, fr);

	start_file = clock();
    for(total_time=0, error_mid=error_moment=0,locate_mid_num=valid_mid_num=locate_moment_num=valid_moment_num=0,frame=0; ;frame++)
	{
        for(flag=i=0; i<NUM; i++)
		{
			//if( fread(in_sp[i], sizeof(short), FRM_LEN2, fq[i]) != FRM_LEN2 )  //48k
			if( fread(mic_sp[i], sizeof(short), FRM_LEN, fq[i]) != FRM_LEN )  //16k
			{flag=1; break;}
			
            //down_sample_run(&down_str[i], in_sp[i], mic_sp[i]);  //48k->16k
		}
		if(flag==1) break;

		if( fread(spk_sp, sizeof(short), FRM_LEN, fr) != FRM_LEN ) break;//16k
		//if( fread(spk, sizeof(short), FRM_LEN2, fr) != FRM_LEN2 ) break;//48k
        //down_sample_run(&down_str[NUM], spk, spk_sp);//48k->16k

	    if((frame&2047L)==0) printf("frame=%ld\n", frame);
		if(frame*160L>=350752)
		{
			i=i;
			//set_position(&my_str, -1, -1, HUMAN_FACE);
			//get_position(&my_str);
		}

		start = clock();
		aec_srp_gsc(&my_str, mic_sp, spk_sp, out_sp);
		finish = clock();
		total_time += (double)(finish-start)/CLOCKS_PER_SEC;
		
		if(get_position(&my_str)>0)
		{
			locate_mid_num++;
			if((fabs(get_position(&my_str)-(12*PI/6))<=PI/18)||(fabs(get_position(&my_str))<=PI/18))
			//if(fabs(get_position(&my_str)-(1*PI/6+PI/18))<=PI/36)
			{
				valid_mid_num++;
				//error_mid += fabs(get_position(&my_str)-(1*PI/6+PI/18));
				if(fabs(get_position(&my_str))<=PI/36)
				{
					error_mid += (Float32)fabs(get_position(&my_str));
				}
				else
				{
					error_mid += (Float32)fabs(get_position(&my_str)-12*PI/6+PI/36);
				}
			}
		}

		if(my_str.audio_az[MED_NUM-1]>0)
		{
			locate_moment_num++;
			if((fabs(my_str.audio_az[MED_NUM-1]-(12*PI/6))<=PI/18)||(fabs(my_str.audio_az[MED_NUM-1])<=PI/18))
			//if(fabs(my_str.audio_az[MED_NUM-1]-(1*PI/6+PI/18))<=PI/36)
			{
				valid_moment_num++;
				error_moment += (Float32)fabs(my_str.audio_az[MED_NUM-1]-1*PI/6);
			}
		}

		fprintf(fp, "%f\n", get_position(&my_str)/PI*180);
		fprintf(fn, "%f\n", my_str.audio_az[MED_NUM-1]/PI*180);

		for (i = 0; i < 6; i++)
		{
			fwrite(out_sp[i], sizeof(short), FRM_LEN, fs[i]);
		}
	}

	finish_file = clock();

	printf("中值滤波的定位成功率为%f\n", ((Float32)valid_mid_num)/((Float32)locate_mid_num));
	printf("中值滤波的平均绝对误差为%f\n", error_mid/valid_mid_num);
	printf("中值滤波的平均相对误差为%f\n", error_mid/valid_mid_num/(12*PI/6));

	printf("瞬时值的定位成功率为%f\n", ((Float32)valid_moment_num)/((Float32)locate_moment_num));
	printf("瞬时值的平均绝对误差为%f\n", error_moment/valid_moment_num);
	printf("瞬时值的平均相对误差为%f\n", error_moment/valid_moment_num/(1*PI/6));

	printf("噪声帧数为%d\n", my_str.frame_noise);
	printf("语音帧数为%d\n", my_str.frame_voice);
	printf("大于20dB的语音帧数为%d\n", my_str.frame_success);
	printf("定位成功的帧数为%d\n", my_str.frame_locate);
	printf("真实定位成功率为%f\n", (double)(my_str.frame_locate)/my_str.frame_success);

	for (i = 0; i < 6; i++) fclose(fs[i]);
	fclose(fr);

	for(i=0; i<NUM; i++) fclose(fq[i]);
	fclose(fp);
	fclose(fn);

	printf("程序运行时间为%f秒\n", (double)(finish_file-start_file)/CLOCKS_PER_SEC);
	printf("程序运行时间为%f秒\n", total_time);
	printf("AEC运行时间为%f秒\n", my_str.aec_time);
	printf("SRP运行时间为%f秒\n", my_str.srp_time);
	printf("MVDR运行时间为%f秒\n", my_str.mvdr_time);

	system("pause");
}