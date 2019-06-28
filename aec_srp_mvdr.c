#include "fft320.h"
#include "aec_srp_mvdr.h"
#include "inverse.h"

const Float32 ptrGCoh[2]={0.8f, 0.2f};
const Float32 kTargetSupp[3] = {-6.9f, -11.5f, -18.4f};
const Float32 min_overdrive[3] = {1.0f, 2.0f, 5.0f};
const Float32 weightCurve[160] = {
    0.0000f, 0.1000f, 0.1239f, 0.1338f, 0.1413f, 0.1477f,
	0.1534f, 0.1585f, 0.1631f, 0.1675f, 0.1716f, 0.1755f,
	0.1792f, 0.1827f, 0.1861f, 0.1893f, 0.1924f, 0.1955f,
	0.1984f, 0.2013f, 0.2040f, 0.2067f, 0.2094f, 0.2119f,
	0.2145f, 0.2169f, 0.2193f, 0.2217f, 0.2240f, 0.2263f,
	0.2285f, 0.2307f, 0.2329f, 0.2350f, 0.2371f, 0.2392f,
	0.2412f, 0.2432f, 0.2452f, 0.2471f, 0.2490f, 0.2509f,
	0.2528f, 0.2547f, 0.2565f, 0.2583f, 0.2601f, 0.2619f,
	0.2636f, 0.2654f, 0.2671f, 0.2688f, 0.2704f, 0.2721f,
	0.2738f, 0.2754f, 0.2770f, 0.2786f, 0.2802f, 0.2818f,
	0.2833f, 0.2849f, 0.2864f, 0.2879f, 0.2894f, 0.2909f,
	0.2924f, 0.2939f, 0.2954f, 0.2968f, 0.2983f, 0.2997f,
	0.3011f, 0.3025f, 0.3039f, 0.3053f, 0.3067f, 0.3081f,
	0.3094f, 0.3108f, 0.3121f, 0.3135f, 0.3148f, 0.3161f,
	0.3174f, 0.3187f, 0.3200f, 0.3213f, 0.3226f, 0.3239f,
	0.3252f, 0.3264f, 0.3277f, 0.3289f, 0.3302f, 0.3314f,
	0.3326f, 0.3338f, 0.3351f, 0.3363f, 0.3375f, 0.3387f,
	0.3399f, 0.3410f, 0.3422f, 0.3434f, 0.3446f, 0.3457f,
	0.3469f, 0.3480f, 0.3492f, 0.3503f, 0.3515f, 0.3526f,
	0.3537f, 0.3548f, 0.3559f, 0.3571f, 0.3582f, 0.3593f,
	0.3604f, 0.3614f, 0.3625f, 0.3636f, 0.3647f, 0.3658f,
	0.3668f, 0.3679f, 0.3690f, 0.3700f, 0.3711f, 0.3721f,
	0.3732f, 0.3742f, 0.3752f, 0.3763f, 0.3773f, 0.3783f,
	0.3794f, 0.3804f, 0.3814f, 0.3824f, 0.3834f, 0.3844f,
	0.3854f, 0.3864f, 0.3874f, 0.3884f, 0.3894f, 0.3904f,
	0.3913f, 0.3923f, 0.3933f, 0.3942f, 0.3952f, 0.3962f,
	0.3971f, 0.3981f, 0.3990f, 0.4000f
};
const Word16 ch_tbl[NUM_CHAN][2] = {{ 2,  3}, { 4,  5}, { 6,  8}, { 9,  11}, {12, 14}, {15, 17}, {18, 20}, {21, 24}, {25, 28}, {29, 33}, {34, 38}, {39, 44}, {45, 51}, {52, 60}, {61, 69}, {70, 79}};
const Word16 vm_tbl[90] = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 7, 7, 7, 8, 8, 9, 9, 10, 10, 11, 12, 12, 13, 13, 14, 15, 15, 16, 17, 17, 18, 19, 20, 20, 21, 22, 23, 24, 24, 25, 26, 27, 28, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50};
enum SOUND_SOURCE_CATEGORY
{
	HUMAN_FACE,
	NOISE_SOURCE
};

Word32 frame;
Float32 prototype_filter[ORD2];                 //coefficient of prototype filter
Float32 cos_tab0[D][2*D], sin_tab0[D][2*D];
Float32 wr[AZ_NUM][NUM][D], wi[AZ_NUM][NUM][D]; //AZ_NUM*NUM*FFT_LEN/2
int inMicLevel = 0, outMicLevel = 0;
uint8_t saturationWarning;

static int CmpFloat(const void* a, const void* b) {
  const float* da = (const float*)a;
  const float* db = (const float*)b;

  return (*da > *db) - (*da < *db);
}

void AGC_init(Agc_t *st)
{
	int minLevel, maxLevel, agcMode, fs;
	AGC_config_t agcConfig;

	minLevel = 0;
	maxLevel = 255;
	agcMode = kAgcModeFixedDigital;
	fs = 16000;
	st->lastError = 0;

	AGC_Init(st, minLevel, maxLevel, agcMode, fs);

	agcConfig.compressionGaindB = 6;
	agcConfig.limiterEnable = 1;
	agcConfig.targetLevelDbfs = 3;
	AGC_set_config(st, agcConfig);
}

void set_position(AEC_SRP_ST *st,          
	              Float32 video_pit,       
	              Float32 video_az,        
				  Word32 source_category  
			     )
{
	st->video_pit = video_pit;
	st->video_az = video_az;
	st->source_category = source_category;
}

Float32 get_position(AEC_SRP_ST *st        
		            )
{
	Word32 i, j;
	Float32 tmp[MED_NUM], acc2;
	
	for(j=0; j<MED_NUM; j++) tmp[j]=st->audio_az[j];

	for(i=0; i<MED_NUM-1; i++)
	{
		for(j=MED_NUM-1; j>i; j--)
		{
			if(tmp[j]<tmp[j-1])
			{
				acc2 = tmp[j];
				tmp[j]=tmp[j-1];
				tmp[j-1]=acc2;
			}
		}
	}
		
	return(tmp[MED_NUM/2]);
}

void aec_srp_gsc_init(AEC_SRP_ST *st,            
				      Float32 appointed_pitch   
		             )
{
	Word32 i, j, k, m, n;
	Float32 theta, tao;

	for(j=0; j<ORD2; j++) prototype_filter[j] = (Float32)(0.54-0.46*cos((2.0*j+1)/ORD2)) * (Float32)sin(PI*(2*j-ORD2+1)/(4*D)) / (Float32)(PI*(2*j-ORD2+1)/2.0);

	for(k=0; k<D; k++)
	{
		//cos_tab2[k]=(Float32)cos(2*PI*k/(2*D));
		//sin_tab2[k]=(Float32)sin(2*PI*k/(2*D));
		for(j=0; j<2*D; j++) cos_tab0[k][j]=(Float32)cos(2*PI*k*j/(2*D));
		for(j=0; j<2*D; j++) sin_tab0[k][j]=(Float32)sin(2*PI*k*j/(2*D));
	}

	fft_320_init();

	//for(k=0; k<5; k++)
	//{
	//	for(j=0; j<5; j++) cos_tab1[k][j]=(Float32)cos(2*PI*k*j/5);
	//	for(j=0; j<5; j++) sin_tab1[k][j]=(Float32)sin(2*PI*k*j/5);
	//}

	for(n=0; n<AZ_NUM; n++)
	{
		theta=(Float32)(n*2*PI/AZ_NUM);
		for(i=0; i<NUM; i++)
		{
    		tao=(Float32)(PI/D)*(Float32)(FS*DIAMETER/2/SPEED)*(Float32)cos(appointed_pitch)*(Float32)cos(theta-2*PI*i/NUM);
			for(k=0; k<D; k++)
			{
			    wr[n][i][k] = (Float32)cos(k*tao);
			    wi[n][i][k] = (Float32)sin(k*tao);
			}
		}
	}

	for(j=0; j<ORD2; j++) st->buf_spk[j] = 0;
	for(m=0; m<NUM; m++) {for(j=0; j<ORD2; j++) st->buf_mic[m][j] = 0;}
	for(k=0; k<D; k++) {for(j=0; j<ORD3; j++) st->spk_ana_re[k][j] = st->spk_ana_im[k][j] = 0;}

	for(m=0; m<NUM; m++) 
	{
		for(k=0; k<D; k++) {for(j=0; j<ORD3; j++) st->h_r[m][k][j] = st->h_i[m][k][j] = 0;}
	}
	for(k=0; k<D; k++) {for(j=0; j<ORD4; j++) st->echo_r[k][j] = st->echo_i[k][j] = 0;}
	for(k=0; k<D; k++) {for(j=0; j<ORD3; j++) st->spk_spec[k][j] = 0;}
	for(j=0; j<ORD3; j++) st->spk_ener[j] = 0;

	for(k=0; k<D; k++) st->sub_az[k]=0;
	st->audio_pit=0;
	for(j=0; j<MED_NUM; j++) st->audio_az[j]=0;
	st->video_pit=-1;
	st->video_az=-1;
	st->source_category=HUMAN_FACE;

	for (i = 0; i < 6; i++)
	{
		for (j = 0; j < ORD2 * 2; j++) st->syn[i][j] = 0;
	}

	st->init_farme_cnt = 0;
	st->voice_frame_cnt = 0;
	st->alpha = 1;
	
	NS_Init(&(st->my_str_anc), 16000, 2);
	AGC_init(&(st->my_str_agc));

    for(k=0; k<D; k++)
    {
        st->mic0_power[k]=0;
        st->spk_power[k]=0;
        st->e0_power[k]=0;
        st->mic0_e0_power_re[k]=0;
        st->mic0_e0_power_im[k]=0;
        st->mic0_spk_power_re[k]=0;
        st->mic0_spk_power_im[k]=0;
    }
    st->divergeState=0;
    st->micState=0;
    st->echoState=0;
    st->hNlMic0SpkAvgMin=1;
    st->nlp_mode = 1;
    st->hNlFbMin = 1;
    st->hNlFbLocalMin = 1;
    st->hNlNewMin = 0;
    st->hNlMinCtr = 0;
    st->overDrive = 2;
    st->overDriveSm = 2;
    st->mult = 2;

	st->pn = 0.0f;

	st->frame_success = 0;
	st->frame_noise = 0;
	st->frame_voice = 0;
	st->frame_locate = 0;
	///VAD
    st->vad_first = TRUE;
    st->hyster_cnt=0;
    st->last_update_cnt=0;
    st->update_cnt=0; 
    st->frame_cnt=0;
	for (i=0; i<NUM_CHAN; i++) st->ch_enrg_long_db[i]=st->ch_enrg[NUM_CHAN]=st->ch_noise[NUM_CHAN]=0;
	
	for(i = 0; i < D; i++)///
	{
        for(m = 0; m < NUM; m++)
		{
            for(n = 0; n < NUM; n++) st->qxx_re[i][m][n] = st->qxx_im[i][m][n] = st->qnn_re[i][m][n] = st->qnn_im[i][m][n] = 0.0;
            st->qnn_re[i][m][m] = 1.0;
        }
	}

	st->aec_time = 0;
	st->srp_time = 0;
	st->mvdr_time = 0;
}

void aec_srp_gsc(AEC_SRP_ST *st,        
            Word16 mic_sp[NUM][FRM_LEN], 
			Word16 spk_sp[FRM_LEN],      
			Word16 out_sp[6][FRM_LEN]      
            )
{
	Float32 mic_ana_re[NUM][D], mic_ana_im[NUM][D];    
	Float32 tmp[2*D], tmp2[2*D], e_r[NUM][D], e_i[NUM][D], ener[D], mic0_spec[D], spectrum[NUM][D], mic0_spec1[D];
	Float32 re[D], im[D], re2[NUM][D], im2[NUM][D], re5[2*D], im5[2*D], re6[NUM][D], im6[NUM][D];
	Float32 acc2, acc3, pearson, aec_step, mic0_ener, cxx[NUM], peak, att[D], power_sub[D], peak_sub[D];
    Float32 mic0_power_sum, e0_power_sum;
    Float32 coh_mic0_e0[D], coh_mic0_spk[D];
    Float32 hNl_mic0_e0_Avg, hNl_mic0_spk_Avg;
    Float32 hNlFb=0, hNlFblow=0;
    Float32 hNl[D]; 
    Float32 hNlPref[60];
    Float32 prefBandQuant = 0.75f, prefBandQuantLow = 0.5f;
	Float32 pn,ps;
	Float32 snr;
	Float32 ch_enrg_db[NUM_CHAN], tne, tce, vv, ch_enrg_dev,alpha;
	Float32 tr, ti, matrix_re[NUM], matrix_im[NUM];
	Float64 inv1_re[NUM][NUM], inv1_im[NUM][NUM], inv_re[NUM][NUM], inv_im[NUM][NUM];

	Word16 minPrefBand = 2, prefBandSize = 60;
	Word16 update_flag, vm_sum, ch_snr[NUM_CHAN];
	Word32  i, j, k, m, n, tot_az, sub_azimuth, num_sub[D];

	clock_t start1, start2, start3, finish1, finish2, finish3;

	start1 = clock();
	//update speaker buffer
	for(j=0; j<ORD2-FRM_LEN; j++) 
	{
		st->buf_spk[j] = st->buf_spk[j+FRM_LEN];
	}
	for(j=0; j<FRM_LEN; j++) 
	{
		st->buf_spk[ORD2-FRM_LEN+j] = spk_sp[j];
	}
    
	//update mic buffer
	for(m=0; m<NUM; m++)
	{
    	for(j=0; j<ORD2-FRM_LEN; j++) 
		{
			st->buf_mic[m][j] = st->buf_mic[m][j+FRM_LEN];
		}
	    for(j=0; j<FRM_LEN; j++)
		{
			st->buf_mic[m][ORD2-FRM_LEN+j] = mic_sp[m][j];
		}
	}

	//update spk subband buffer
	for(k=1; k<D; k++) 
	{
		for(j=0; j<ORD3-1; j++) 
		{
			st->spk_ana_re[k][j] = st->spk_ana_re[k][j+1];
			st->spk_ana_im[k][j] = st->spk_ana_im[k][j+1];
			st->spk_spec[k][j] = st->spk_spec[k][j+1];
		}
	}
	for(j=0; j<ORD3-1; j++)
	{
		st->spk_ener[j] = st->spk_ener[j+1];
	}

    //Decomposing spk-signal into subbands
    for(i=0; i<2*D; i++)
	{
	    for(acc2=0.0f, j=0; j<ORD2/(2*D); j++)	
		{
			acc2 += st->buf_spk[j*2*D+i]*prototype_filter[j*2*D+i];
		}
	    tmp[2*D-1-i]=acc2;
        tmp2[i] = 0;
	}
    fft_320(tmp, tmp2);

    //Decompose mic-signal into subbands
    for(k=0; k<D; k++)
    {
        st->spk_ana_re[k][ORD3-1] = tmp[k];
        st->spk_ana_im[k][ORD3-1] = tmp2[k];
    }
    for(m=0; m<NUM; m++)
	{
    	for(i=0; i<2*D; i++)
		{
		    for(acc2=0.0f, j=0; j<ORD2/(2*D); j++)
			{
				acc2 += st->buf_mic[m][j*2*D+i]*prototype_filter[j*2*D+i];
			}
		    tmp[2*D-1-i]=acc2;
            tmp2[i]=0;
		}
        fft_320(tmp, tmp2);
	    for(k=1; k<D; k++)
		{
			mic_ana_re[m][k] = tmp[k];
		    mic_ana_im[m][k] = tmp2[k];			
		}
	}

	//update magnitude spectra and energy
	for(acc3=0.0f, k=1; k<D; k++) 
	{
		acc2 = st->spk_ana_re[k][ORD3-1]*st->spk_ana_re[k][ORD3-1] + st->spk_ana_im[k][ORD3-1]*st->spk_ana_im[k][ORD3-1];
		acc3 += acc2;
		st->spk_spec[k][ORD3-1] = (Float32)sqrt(acc2);
	}
	st->spk_ener[ORD3-1]=acc3;
	for(mic0_ener=0.0f, k=1; k<D; k++) 
	{
		acc2 = mic_ana_re[0][k]*mic_ana_re[0][k] + mic_ana_im[0][k]*mic_ana_im[0][k];
		mic0_ener +=acc2;
		mic0_spec[k] = (Float32)sqrt(acc2);
	}

	//Calculating Pearson correlation coefficient
	for(pearson=0.0f, j=0; j<ORD3; j++)
	{
    	for(acc2=0.0f, k=1; k<D; k++) 
		{
			acc2 += st->spk_spec[k][j]*mic0_spec[k];
		}
		acc3 = (Float32)sqrt(mic0_ener*st->spk_ener[j]);
		acc2 = (acc3>1024.0f) ? acc2/acc3 : 0;
		if(pearson<acc2) 
		{
			pearson=acc2;
		}
	}
	pearson = (Float32)pow(pearson, 2);

	//choosing step size
	if(pearson>0.9)      aec_step = 1.0f;
	else if(pearson>0.7) aec_step = 0.4f;
	else if(pearson>0.5) aec_step = 0.16f;
	else if(pearson>0.3) aec_step = 0.006f;
    //else if(pearson>0.3) aec_step = 0.8f;
	else                 aec_step = 0.0001f;

	//update the echo
	for(k=1; k<D; k++)
	{
		for(j=0; j<ORD4-1; j++) 
		{
			st->echo_r[k][j] = st->echo_r[k][j+1];
			st->echo_i[k][j] = st->echo_i[k][j+1];
		}
	}

	//filtering spk-signal
	for(k=1; k<D; k++)
	{
		//calculating subband energy
		for(acc2=0.0f, j=0; j<ORD3; j++) 
		{
			acc2 += st->spk_ana_re[k][j]*st->spk_ana_re[k][j] + st->spk_ana_im[k][j]*st->spk_ana_im[k][j];
		}
	    ener[k]=acc2;

		//filtering
		for(acc3=acc2=0.0f, j=0; j<ORD3; j++) 
		{
			acc2 += st->spk_ana_re[k][j]*st->h_r[0][k][j] - st->spk_ana_im[k][j]*st->h_i[0][k][j];
			acc3 += st->spk_ana_im[k][j]*st->h_r[0][k][j] + st->spk_ana_re[k][j]*st->h_i[0][k][j];
		}
		st->echo_r[k][ORD4-1]=acc2;		
		st->echo_i[k][ORD4-1]=acc3;

		//calculating error signal of the first channel
		e_r[0][k] = mic_ana_re[0][k] - acc2;
		e_i[0][k] = mic_ana_im[0][k] - acc3;
	}

	//update the cofficient of adaptive filter by NLMS algorithm
   	for(k=1; k<D; k++)
	{
		if(ener[k]>8.0f)
		{
			acc2 = aec_step/ener[k]; 
     	    for(j=0; j<ORD3; j++) 
			{
				st->h_r[0][k][j] += acc2*(st->spk_ana_re[k][j]*e_r[0][k] + st->spk_ana_im[k][j]*e_i[0][k]);
				st->h_i[0][k][j] += acc2*(st->spk_ana_re[k][j]*e_i[0][k] - st->spk_ana_im[k][j]*e_r[0][k]);
			}
		}
	}

	//calculating energy of echo
   	for(k=1; k<D; k++)
	{
	    
   		for(acc2=0.0f, j=0; j<ORD4; j++)
		{
			acc2 += st->echo_r[k][j]*st->echo_r[k][j] + st->echo_i[k][j]*st->echo_i[k][j];
		}
        ener[k]=acc2;
	}

	//calculating error signals of 2-6 channels
	for(m=1; m<NUM; m++)
	{
    	for(k=1; k<D; k++)
		{
	    	for(acc3=acc2=0.0f, j=0; j<ORD4; j++) 
			{
				acc2 += st->echo_r[k][j]*st->h_r[m][k][j] - st->echo_i[k][j]*st->h_i[m][k][j];
				acc3 += st->echo_i[k][j]*st->h_r[m][k][j] + st->echo_r[k][j]*st->h_i[m][k][j];
			}

    		e_r[m][k] = mic_ana_re[m][k] - acc2;
	    	e_i[m][k] = mic_ana_im[m][k] - acc3;
		}
     	for(k=1; k<D; k++)
		{
		    if(ener[k]>8.0f)
			{
    			acc2 = aec_step/ener[k];
         	    for(j=0; j<ORD4; j++)
				{
					st->h_r[m][k][j] += acc2*(st->echo_r[k][j]*e_r[m][k] + st->echo_i[k][j]*e_i[m][k]);
					st->h_i[m][k][j] += acc2*(st->echo_r[k][j]*e_i[m][k] - st->echo_i[k][j]*e_r[m][k]);
				}
			}
		}
	}

    for(mic0_power_sum=e0_power_sum=0, k=1; k<D; k++)
    {
		//calculating the PSD and CPSD
        st->mic0_power[k] = ptrGCoh[0]*st->mic0_power[k]+ptrGCoh[1]*(mic_ana_re[0][k]*mic_ana_re[0][k]+mic_ana_im[0][k]*mic_ana_im[0][k]);
        st->spk_power[k] = ptrGCoh[0]*st->spk_power[k]+ptrGCoh[1]*(st->spk_ana_re[k][ORD3-1]*st->spk_ana_re[k][ORD3-1] + st->spk_ana_im[k][ORD3-1]*st->spk_ana_im[k][ORD3-1]);
        st->e0_power[k] = ptrGCoh[0]*st->e0_power[k]+ptrGCoh[1]*(e_r[0][k]*e_r[0][k]+e_i[0][k]*e_i[0][k]);
        st->mic0_e0_power_re[k] = ptrGCoh[0]*st->mic0_e0_power_re[k]+ptrGCoh[1]*(mic_ana_re[0][k]*e_r[0][k]+mic_ana_im[0][k]*e_i[0][k]);
        st->mic0_e0_power_im[k] = ptrGCoh[0]*st->mic0_e0_power_im[k]+ptrGCoh[1]*(mic_ana_re[0][k]*e_i[0][k]-mic_ana_im[0][k]*e_r[0][k]);
        st->mic0_spk_power_re[k] = ptrGCoh[0]*st->mic0_spk_power_re[k]+ptrGCoh[1]*(mic_ana_re[0][k]*st->spk_ana_re[k][ORD3-1]+mic_ana_im[0][k]*st->spk_ana_im[k][ORD3-1]);
        st->mic0_spk_power_im[k] = ptrGCoh[0]*st->mic0_spk_power_im[k]+ptrGCoh[1]*(mic_ana_re[0][k]*st->spk_ana_im[k][ORD3-1]-mic_ana_im[0][k]*st->spk_ana_re[k][ORD3-1]);

        mic0_power_sum += st->mic0_power[k];
        e0_power_sum += st->e0_power[k];
    }

	//Judging divergence
    if(st->divergeState == 0)
    {
        if(e0_power_sum > mic0_power_sum) 
		{ 
			st->divergeState = 1;
		}
    }
    else
    {
        if(e0_power_sum * 1.05f < mic0_power_sum)
		{
			st->divergeState = 0;
		}
    }

    if(st->divergeState == 1)
    {
        memcpy(e_r[0], mic_ana_re[0], D*sizeof(float));
        memcpy(e_i[0], mic_ana_im[0], D*sizeof(float));
    }
    
	//calculating correlation coefficients
    for(k=1; k<D; k++)
    {
        coh_mic0_e0[k]=(st->mic0_e0_power_re[k]*st->mic0_e0_power_re[k]+st->mic0_e0_power_im[k]*st->mic0_e0_power_im[k])/(st->mic0_power[k]*st->e0_power[k]+1e-10f);
        coh_mic0_spk[k]=(st->mic0_spk_power_re[k]*st->mic0_spk_power_re[k]+st->mic0_spk_power_im[k]*st->mic0_spk_power_im[k])/(st->mic0_power[k]*st->spk_power[k]+1e-10f);
    }

	//calculating the average correlation coefficient between mic0 and spk in preferred bands
    hNl_mic0_spk_Avg=0;
    for(k=minPrefBand; k<minPrefBand+prefBandSize; k++)
    {
        hNl_mic0_spk_Avg += coh_mic0_spk[k];
    }
    hNl_mic0_spk_Avg /= prefBandSize;
    hNl_mic0_spk_Avg = 1 - hNl_mic0_spk_Avg;

	//calculating the average correlation coefficient between mic0 and e0 in preferred bands
    hNl_mic0_e0_Avg=0;
    for(k=minPrefBand; k<minPrefBand+prefBandSize; k++)
    {
        hNl_mic0_e0_Avg += coh_mic0_e0[k];
    }
    hNl_mic0_e0_Avg /= prefBandSize;

	//Judging the proximal signal and the distal signal
    if(hNl_mic0_spk_Avg < 0.75f && hNl_mic0_spk_Avg < st->hNlMic0SpkAvgMin)
    {
        st->hNlMic0SpkAvgMin = hNl_mic0_spk_Avg;
    }

    if(hNl_mic0_e0_Avg > 0.98f && hNl_mic0_spk_Avg > 0.9f)
    {
        st->micState = 1;
    }
    else if(hNl_mic0_e0_Avg < 0.95f || hNl_mic0_spk_Avg < 0.8f)
    {
        st->micState = 0;
    }

	//calculating weighs
    if(st->hNlMic0SpkAvgMin == 1)
    {
        st->echoState = 0;
		st->overDrive = min_overdrive[st->nlp_mode];
        if(st->micState == 1)
        {
            memcpy(hNl, coh_mic0_e0, D*sizeof(float));
            hNlFb = hNl_mic0_e0_Avg;
            hNlFblow = hNl_mic0_e0_Avg;
        }
        else
        {
            for(k=1; k<D; k++) 
			{ 
				hNl[k] = 1 - coh_mic0_spk[k];
			}
            hNlFb = hNl_mic0_spk_Avg;
            hNlFblow = hNl_mic0_spk_Avg;
        }
    }
    else
    {
        if(st->micState == 1)
        {
            st->echoState = 0;
            memcpy(hNl, coh_mic0_e0, D*sizeof(float));
            hNlFb = hNl_mic0_e0_Avg;
            hNlFblow = hNl_mic0_e0_Avg;
        }
        else
        {
            st->echoState = 1;
            for(k=1; k<D; k++) 
			{
				hNl[k] = min(coh_mic0_e0[k], 1-coh_mic0_spk[k]);
			}
            memcpy(hNlPref, &hNl[minPrefBand], sizeof(float)*prefBandSize);
            qsort(hNlPref, prefBandSize, sizeof(float), CmpFloat);
            hNlFb = hNlPref[(int)floor(prefBandQuant*(prefBandSize-1))];
            hNlFblow = hNlPref[(int)floor(prefBandQuantLow*(prefBandSize-1))];
        }
    }

    if(hNlFblow < 0.6f && hNlFblow < st->hNlFbLocalMin)
    {
        st->hNlFbLocalMin = hNlFblow;
        st->hNlFbMin = hNlFblow;
        st->hNlNewMin = 1;
        st->hNlMinCtr = 0;
    }

    st->hNlFbLocalMin = min(st->hNlFbLocalMin + 0.0008f / st->mult, 1);
    st->hNlMic0SpkAvgMin = min(st->hNlMic0SpkAvgMin + 0.0006f / st->mult, 1);

    if(st->hNlNewMin == 1)
    {
        st->hNlMinCtr++;
    }

    if(st->hNlMinCtr == 2)
    {
        st->hNlNewMin = 0;
        st->hNlMinCtr = 0;
        st->overDrive = max(kTargetSupp[st->nlp_mode]/((float)log(st->hNlFbMin + 1e-10f) + 1e-10f), min_overdrive[st->nlp_mode]);
    }

    //Smooth the overdrive
    if(st->overDrive < st->overDriveSm)
    {
        st->overDriveSm = 0.99f * st->overDriveSm + 0.01f * st->overDrive;
    }
    else
    {
        st->overDriveSm = 0.9f * st->overDriveSm + 0.1f * st->overDrive;
    }

    for(k=1; k<D; k++)
    {
		//smooth the weigh
        if(hNl[k] > hNlFb)
        {
            hNl[k] = weightCurve[k] * hNlFb + (1 - weightCurve[k]) * hNl[k];
        }
		hNl[k] *= 1.0f;

		//weighing the error signals
		for(j=0; j<NUM; j++)
		{
	    	e_r[j][k] *= (hNl[k]);
	    	e_i[j][k] *= (hNl[k]);
		}
    }
    
    finish1 = clock();
	st->aec_time += (double)(finish1-start1)/CLOCKS_PER_SEC;

	start2 = clock();
	//calculating the magnitude spectra of error signals
	for(m=0; m<NUM; m++)
	{
		for(cxx[m]=0, k=1; k<D; k++)
		{
			acc2 = e_r[m][k]*e_r[m][k]+e_i[m][k]*e_i[m][k];
			cxx[m] += acc2;
			spectrum[m][k]=(Float32)sqrt(acc2);		    
		}
	}
	
	//ALA:Automatic Level Alignment
	for(acc2=0, m=0; m<NUM; m++) 
	{
		acc2 += cxx[m];
	}
	for(m=0; m<NUM; m++) 
	{
		if(cxx[m] > 0.0f)
		{
			acc3 = (Float32)sqrt(acc2/NUM/cxx[m]);
		    for(k=1; k<D; k++) 
			{
				re6[m][k] = e_r[m][k] * acc3;
				im6[m][k] = e_i[m][k] * acc3;
			}
		}
		else
		{
			for(k=1; k<D; k++) 
			{
				re6[m][k] = 0;
				im6[m][k] = 0;
			}
		}
	}

	//PHAT
	for(m=0; m<NUM; m++)
	{
		for(k=1; k<D; k++) 
		{
			acc2 = (Float32)sqrt(re6[m][k]*re6[m][k] + im6[m][k]*im6[m][k]);
			if (acc2 > 0.0f)
			{
				re2[m][k] = re6[m][k]  / acc2;
				im2[m][k] = im6[m][k]  / acc2;
			}
		}
	}

	//calculating the power-weigh
	for (peak = 0.0f, k = 1; k < D; k++)
	{
		acc3 = re6[0][k] * re6[0][k] + im6[0][k] * im6[0][k];
		mic0_spec1[k] = (Float32)sqrt(acc3);
		if (peak < mic0_spec1[k]) peak = mic0_spec1[k];
	}

	if (peak > 0.1f)
	{
		//for (k = 1; k < D; k++) mic0_spec[k] = max(2.0f + 1.0f*(Float32)log10(mic0_spec1[k] / peak), 0);
		//for (peak = 0.0f, k = 1; k < D; k++) { if (peak < mic0_spec[k]) peak = mic0_spec[k]; }
		//for (k = 1; k < D; k++) mic0_spec[k] = mic0_spec[k] / peak;
		for(k=1; k<D; k++)
		{
			acc3 = mic0_spec1[k]/peak;
			if(acc3>=0.01f) mic0_spec[k] = 1;
			else mic0_spec[k] = 0;
		}
	}

	for (i = ORD2 * 2 - 1; i >= 2 * D; i--) st->syn[0][i] = st->syn[0][i - 2 * D];
	for (i = 0; i < 2 * D; i++)
	{
		for(acc2=0.0f, k=1; k<D; k++) acc2 += e_r[0][k]*cos_tab0[k][i] - e_i[0][k]*sin_tab0[k][i];
		st->syn[0][i]=acc2;	    
	}
	for (i = 0; i < FRM_LEN; i++)
	{
		for(acc2=0.0f, j=0; j<ORD2/D; j++) acc2 += prototype_filter[j*D+i]*st->syn[0][j*2*D+(j&1)*D+i];
		out_sp[0][i]=(short)(acc2*D*32);
	}

	//VAD
	alpha = (st->vad_first==TRUE) ? 1.0f : 0.55f;
    for (i=LO_CHAN; i<=HI_CHAN; i++)
	{
        for(vv=0.0f, j=ch_tbl[i][0]; j<=ch_tbl[i][1]; j++) vv += re6[0][j]*re6[0][j] + im6[0][j]*im6[0][j];
        st->ch_enrg[i] = max((1.0f-alpha)*st->ch_enrg[i]+alpha*vv/(ch_tbl[i][1]-ch_tbl[i][0]+1), MIN_CHAN_ENRG);
	}
	
	st->frame_cnt++;    
    /* Initialize channel noise estimate to channel energy of vad_first four frames */
    if (st->frame_cnt<6) {for(i=LO_CHAN; i<=HI_CHAN; i++) st->ch_noise[i] = max(st->ch_enrg[i], INE);}    
    
    /* Compute the channel SNR indices */
    for (i=LO_CHAN; i<=HI_CHAN; i++) ch_snr[i] = (int)((max(10.0f*log10(st->ch_enrg[i]/st->ch_noise[i]), 0.0f) + 0.1875f)/0.375f);
    
    /* Compute the sum of voice metrics */
    for (vm_sum=0, i=LO_CHAN; i<=HI_CHAN; i++) vm_sum += vm_tbl[min(ch_snr[i], 89)];
    /* Compute the total noise estimate (tne) and total channel energy estimate (tce) */
    for (tne=0.0f, i=LO_CHAN; i<=HI_CHAN; i++) tne += st->ch_noise[i];
    for (tce=0.0f, i=LO_CHAN; i<=HI_CHAN; i++) tce += st->ch_enrg[i];
    /* Calculate log spectral deviation */
    for (i=LO_CHAN; i<=HI_CHAN; i++) ch_enrg_db[i] = 10.0f*(Float32)log10(st->ch_enrg[i]);
   
    if (st->vad_first==TRUE) {for (i=LO_CHAN; i<=HI_CHAN; i++) st->ch_enrg_long_db[i] = ch_enrg_db[i];}

    for (ch_enrg_dev=0.0f, i=LO_CHAN; i<=HI_CHAN; i++) ch_enrg_dev += (Float32)fabs(st->ch_enrg_long_db[i] - ch_enrg_db[i]);
    /* Calculate long term integration constant as a function of total channel energy (tce) *//* (i.e., high tce (-40 dB) -> slow integration (alpha = 0.99), low tce (-60 dB) -> fast integration (alpha = 0.50) */
    alpha = max(min(0.99f - (0.49f/TCE_RANGE)*(HIGH_TCE_DB - 10.0f*(Float32)log10(tce)), 0.9f), 0.5f);
    /* Calc long term log spectral energy */
    for (i=LO_CHAN; i<=HI_CHAN; i++) st->ch_enrg_long_db[i] = alpha*st->ch_enrg_long_db[i] + (1.0f-alpha)*ch_enrg_db[i];
    
    /* Set or reset the update flag */
    update_flag = FALSE;
    if (vm_sum<=UPDATE_THLD)
    {
		st->update_cnt = 0;
        update_flag = TRUE;
    }
    else if ((tce>NOISE_FLOOR)&&(ch_enrg_dev<DEV_THLD))
    {
        st->update_cnt++;
        if (st->update_cnt>=UPDATE_CNT_THLD) update_flag = TRUE;
    }
    
    if (st->update_cnt==st->last_update_cnt) st->hyster_cnt++;
    else                                     st->hyster_cnt = 0;
    
	st->last_update_cnt = st->update_cnt;
    if (st->hyster_cnt>HYSTER_CNT_THLD) st->update_cnt = 0;

    /* Update the channel noise estimates */
    if (update_flag==TRUE) {for(i=LO_CHAN; i<=HI_CHAN; i++) st->ch_noise[i] = max(0.9f*st->ch_noise[i]+0.1f*st->ch_enrg[i], MIN_CHAN_ENRG);}    
    st->vad_first = FALSE;
    //alpha = (update_flag == TRUE) ? AA : 1.0f;


	if(update_flag == 1) for(i=0; i<FRM_LEN; i++) out_sp[1][i] = 0;
	else         for(i=0; i<FRM_LEN; i++) out_sp[1][i] = 4000;


	//if(vad_cnt == 6) for(i=0; i<FRM_LEN; i++) out_sp[2][i] = 4000;
	//else         for(i=0; i<FRM_LEN; i++) out_sp[2][i] = 0;

	if(st->frame_cnt >= 6)
	{
		if(update_flag == 1)
		{
			//Estimating noise power
			snr=0.0f;
			st->frame_noise++;
			for(pn=0.0f,k=1; k<D; k++)
			{
				pn += re6[0][k]*re6[0][k] + im6[0][k]*im6[0][k];
			}

			if(st->frame_noise>30)
			{
				if((pn/st->pn<=1.2f)&&(st->pn/pn<=1.2f))
				{
					st->pn = (Float32)(0.9*st->pn + 0.1*pn);
				}
			}
			else
			{
				st->pn = (Float32)(0.9*st->pn + 0.1*pn);
			}

			//st->pn = 0.9*st->pn + 0.1*pn;
			//st->pn = pn;
		}
		else
		{
			st->frame_voice++;
			for(ps=0.0f,k=1; k<D; k++)
			{
				ps += re6[0][k]*re6[0][k] + im6[0][k]*im6[0][k];
			}
			//calculating the SNR
			if(st->pn>0.0f)
			{
				snr = 10.0f * (Float32)log10(ps/st->pn);
			}
			else
			{
				snr = 0.0f;
			}
			

			if(snr > 15.0f)
			{
				if(frame>=178)
				{
					i=i;
				}
				st->frame_success++;

				//approximate search of azimuth
				for(k=1; k<D; k++) { peak_sub[k] = 0; }
		        for(peak=0.0f, j=n=0; n<AZ_NUM; n+=AZ_STEP)
				{
					//calculating subband power and total power
			        for(acc2=0.0f, k=1; k<D; k++)			
					{ 			
				        for(re[k]=0, m=0; m<NUM; m++) 
						{
							re[k] += re2[m][k]*wr[n][m][k] + im2[m][k]*wi[n][m][k];
						}
				        for(im[k]=0, m=0; m<NUM; m++) 
						{
							im[k] += im2[m][k]*wr[n][m][k] - re2[m][k]*wi[n][m][k];
						}
				        power_sub[k] = re[k]*re[k]+im[k]*im[k];
                        if(peak_sub[k]<power_sub[k])
						{
                            peak_sub[k] = power_sub[k];
                            num_sub[k] = n;
						}
				        acc2 += power_sub[k]*mic0_spec[k];
					}
			        if(peak<acc2)
					{
				        peak=acc2;
				        j=n; 
					}
				}

				//elaborate search of subband azimuths
                for(k=1; k<D; k++)
				{
                    for(sub_azimuth=num_sub[k], i=num_sub[k]-AZ_STEP/2; i<num_sub[k]+AZ_STEP/2; i++)			
					{
	    		        n=i;
		    	        if(n<0) n += AZ_NUM;
			            if(n>AZ_NUM) n -= AZ_NUM;      			
				        for(re[k]=0, m=0; m<NUM; m++)
						{
							re[k] += re2[m][k]*wr[n][m][k] + im2[m][k]*wi[n][m][k];
						}
				        for(im[k]=0, m=0; m<NUM; m++)
						{
							im[k] += im2[m][k]*wr[n][m][k] - re2[m][k]*wi[n][m][k];
						}

				        acc2 = re[k]*re[k]+im[k]*im[k];
				        if(peak_sub[k]<acc2) 
						{
							peak_sub[k]=acc2;
							sub_azimuth=n;
						}
					}
                    st->sub_az[k] = sub_azimuth*2.0f*(Float32)PI/AZ_NUM;
				}

				//elaborate search of fullband azimuths
		        for(tot_az=j, i=j-AZ_STEP/2; i<=j+AZ_STEP/2; i++)
				{
			        n=i;
			        if(n<0) n += AZ_NUM;
			        if(n>AZ_NUM) n -= AZ_NUM;
			
			        for(acc2=0.0f, k=1; k<D; k++)
					{  			
				        for(re[k]=0, m=0; m<NUM; m++) 
						{
							re[k] += re2[m][k]*wr[n][m][k] + im2[m][k]*wi[n][m][k];
						}
				        for(im[k]=0, m=0; m<NUM; m++) 
						{
							im[k] += im2[m][k]*wr[n][m][k] - re2[m][k]*wi[n][m][k];
						}

				        acc2 += (re[k]*re[k]+im[k]*im[k])*mic0_spec[k];
					}
			        if(peak<acc2)
					{
				        peak=acc2;
				        tot_az=n;
					}
				}

		        for(j=0; j<MED_NUM-1; j++) 
				{
					st->audio_az[j]=st->audio_az[j+1];
				}
		        st->audio_az[MED_NUM-1] = tot_az*2.0f*(Float32)PI/AZ_NUM;
				//if(fabs(get_position(st)-PI/2)<=PI/18)
				if(fabs(st->audio_az[MED_NUM-1]-PI/2)<=PI/18)
				//if((fabs(get_position(st)-2*PI)<=PI/18)||(fabs(get_position(st))<=PI/18))
				//if((fabs(st->audio_az[MED_NUM-1]-2*PI)<=PI/36)||(fabs(st->audio_az[MED_NUM-1])<=PI/36))
				{
					st->frame_locate++;
				}
			}
		}
	}

	finish2 = clock();
	st->srp_time += (double)(finish2-start2)/CLOCKS_PER_SEC;
	
	if( (0<=st->video_az)&&(st->video_az<2*PI)&&(-PI/2.0<st->video_pit)&&(st->video_pit<PI/2) )
	{
		peak = st->video_az;
	}
	else
	{
		peak = get_position(st);
	}

   	n = (Word32)floor(peak*AZ_NUM/2.0f/PI+0.5);
    if(n<0) n += AZ_NUM;
    if(n>AZ_NUM) n -= AZ_NUM;
    
    start3 = clock();
    //MVDR
    alpha = (update_flag == TRUE) ? AA : 1.0f;
    for(k=1; k<D; k++)
    {
		//Estimating the signal and noise covariance matrix
        for(m=0; m<NUM; m++)
		{
            for(n=0; n<=m; n++)
            {
                tr = re6[m][k]*re6[n][k] + im6[m][k]*im6[n][k];
                ti = im6[m][k]*re6[n][k] - re6[m][k]*im6[n][k];
                st->qxx_re[k][m][n] = 0.2f*tr + 0.8f*st->qxx_re[k][m][n];
                st->qxx_im[k][m][n] = 0.2f*ti + 0.8f*st->qxx_im[k][m][n];      
				st->qnn_re[k][m][n] = (1-alpha)*tr + alpha*st->qnn_re[k][m][n];
                st->qnn_im[k][m][n] = (1-alpha)*ti + alpha*st->qnn_im[k][m][n];
			    inv1_re[m][n] = st->qxx_re[k][m][n] + (BB - 1)*st->qnn_re[k][m][n];
                inv1_im[m][n] = st->qxx_im[k][m][n] + (BB - 1)*st->qnn_im[k][m][n];
				//inv1_re[m][n] = st->qxx_re[k][m][n] ;
				//inv1_im[m][n] = st->qxx_im[k][m][n] ;
			}              
        }
        for(m=0; m<NUM; m++)
		{
            for(n=m+1; n<NUM; n++)
            {
                st->qxx_re[k][m][n] = st->qxx_re[k][n][m];
				st->qxx_im[k][m][n] = -st->qxx_im[k][n][m];
                st->qnn_re[k][m][n] = st->qnn_re[k][n][m];
                st->qnn_im[k][m][n] = -st->qnn_im[k][n][m];
				inv1_re[m][n] = inv1_re[n][m];
                inv1_im[m][n] = -inv1_im[n][m];
            }
		}
        inverse(inv1_re, inv1_im, inv_re, inv_im, NUM);

        for(m=0; m<NUM; m++)
        {
            matrix_re[m] = st->qxx_re[k][m][0] - st->qnn_re[k][m][0];
            matrix_im[m] = st->qxx_im[k][m][0] - st->qnn_im[k][m][0];
        }

		//for(m=0; m<NUM; m++)
		//{
		//	acc3 = sqrt(matrix_re[m]*matrix_re[m] + matrix_im[m]*matrix_im[m]);
		//	if(acc3>0.0f)
		//	{
		//		matrix_re[m] /= acc3;
		//		matrix_im[m] /= acc3;
		//	}
		//	else
		//	{
		//		matrix_re[m] = 0.0f;
		//		matrix_im[m] = 0.0f;
		//	}
		//}
        e_r[0][k]=e_i[0][k]=0.0f;
        for(m=0; m<NUM; m++)
        {
            tr=ti=0.0f;
			for (n=0; n<NUM; n++)
			{
    			tr += (Float32)(inv_re[m][n] * matrix_re[n] - inv_im[m][n] * matrix_im[n]);
    			ti += (Float32)(inv_re[m][n] * matrix_im[n] + inv_im[m][n] * matrix_re[n]);
				//tr += (Float32)(inv_re[m][n] * re2[n][k] - inv_im[m][n] * im2[n][k]);
				//ti += (Float32)(inv_re[m][n] * im2[n][k] + inv_im[m][n] * re2[n][k]);
			}
			e_r[0][k] += tr*re6[m][k] + ti*im6[m][k];
            e_i[0][k] += tr*im6[m][k] - ti*re6[m][k];
        }
    }

	finish3 = clock();
	st->mvdr_time += (double)(finish3-start3)/CLOCKS_PER_SEC;
    
    for (i = ORD2 * 2 - 1; i >= 2 * D; i--) st->syn[2][i] = st->syn[2][i - 2 * D];
	for (i = 0; i < 2 * D; i++)
	{
		for (acc2 = 0.0f, k = 1; k < D; k++) acc2 += e_r[0][k] * cos_tab0[k][i] - e_i[0][k] * sin_tab0[k][i];
		st->syn[2][i] = acc2;
	}
	for (i = 0; i < FRM_LEN; i++)
	{
		for (acc2 = 0.0f, j = 0; j < ORD2 / D; j++) acc2 += prototype_filter[j*D + i] * st->syn[2][j * 2 * D + (j & 1)*D + i];
		out_sp[2][i] = (short)(acc2*D * 32);
	}

	//Geometric blind separation
	for(k=1; k<D; k++)
	{
		att[k] = (Float32)pow(0.5 + 0.5*cos(st->sub_az[k] - peak), 6);
		e_r[0][k] *= att[k];
		e_i[0][k] *= att[k]; 
	}

	for(i=ORD2*2-1; i>=2*D; i--) st->syn[3][i] = st->syn[3][i-2*D];

    for(i=1; i<D; i++)
	{
		re5[2*D-i] = re5[i] = e_r[0][i];
		im5[2*D-i] = -e_i[0][i];
		im5[i] = e_i[0][i];
	}
	re5[D] = re5[0] = im5[D] = im5[0] = 0;

	fft_320(re5, im5);

    for(i=0; i<2*D; i++)
	{
		st->syn[3][i] = re5[i];
	}

	for(i=0; i<FRM_LEN; i++)
	{
	    for(acc2=0.0f, j=0; j<ORD2/D; j++) acc2 += prototype_filter[j*D+i]*st->syn[3][j*2*D+(j&1)*D+i];
		//out_sp[i]=(short)(acc2*D*32);
		out_sp[3][i]=(short)(acc2*D*16);
		//temp1[i]=(short)(acc2*D*32);
	}
	
	//ANC
	NS_run(&(st->my_str_anc), out_sp[3], out_sp[4]);
	//AGC
	AGC_Process(&(st->my_str_agc), out_sp[4], NULL, FRM_LEN, out_sp[5], NULL, inMicLevel, &outMicLevel, 0, (uint8_t*)(&saturationWarning));
    
	st->init_farme_cnt = (st->init_farme_cnt<INIT_LEN) ? (st->init_farme_cnt+1) : INIT_LEN;
}
