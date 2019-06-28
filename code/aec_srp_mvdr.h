//#include "stdio.h"
//#include "stdlib.h"
#include "math.h"
#include "time.h"
//#include "fft320.h"
#include "typedefs.h"
#include "typedefs1.h"

#define Word16 short
#define Word32 long
#define Float32 float
#define Float64 double

#define PI           3.1415926535897932
#define FS           16000           //sampling rate
#define FRM_LEN      160             //length of a frame(16k)
#define NUM          6               //number of microphone
#define D            FRM_LEN         //half of the number of subband
#define ORD2         (6*2*D)         //order of prototype low pass filter
#define T60          200             //reverberation time
#define ORD3         (T60*FS/1000/D) //order of subband filter in first channel
#define ORD4         6               //order of subband filter in other channels
#define ORD6         10

#define DIAMETER     0.08f           //diameter of microphone array(cm)
#define SPEED        343             //velocity of sound
#define AZ_STEP      18              //step size in search
#define AZ_NUM       (20*AZ_STEP)    //number of azimuth
#define INIT_LEN     64              //Initialization number of noise spectrum
#define MED_NUM      7               //order of median filter 
#define AA           0.95f 
#define BB           11    

///VAD
#define TRUE			1
#define FALSE			0
#define	NUM_CHAN		16
#define	LO_CHAN			0
#define	MID_CHAN		5
#define	HI_CHAN			15
#define	UPDATE_THLD		35
#define	METRIC_THLD		45
#define	INDEX_THLD		12
#define	SETBACK_THLD	12
#define	SNR_THLD		6
#define	INDEX_CNT_THLD	5
#define	UPDATE_CNT_THLD	50
#define	NORM_ENRG		1.0f	     // use (32768.0 * 32768.0) for fractional 
#define	NOISE_FLOOR		1.0f
#define	MIN_CHAN_ENRG	0.0625f
#define	INE			    16.0f
#define	MIN_GAIN		-13.0f
#define	HYSTER_CNT_THLD	6		     // forced update constants... 
#define	HIGH_TCE_DB		50.0f	     // 50 - 10log10(NORM_ENRG) 
#define	LOW_TCE_DB		30.0f	     // 30 - 10log10(NORM_ENRG) 
#define	TCE_RANGE		(HIGH_TCE_DB-LOW_TCE_DB)
#define	DEV_THLD		28.0f

typedef struct
{
	Word16  buf_mic[NUM][ORD2];       //buffer of mic-signal
	Word16  buf_spk[ORD2];            //buffer of speaker-signal
	Word16 divergeState;  
    Word16 micState;      
    Word16 echoState;
    Word16 hNlNewMin, hNlMinCtr;
    Word16 nlp_mode;
    Word16 mult;
    Word16  init_farme_cnt;
	Word16  voice_frame_cnt;
	Word16 frame_success;
	Word16 frame_noise;
	Word16 frame_voice;
	Word16 frame_locate;
	Word16 vad_first;
    Word16 hyster_cnt;
    Word16 last_update_cnt;
    Word16 update_cnt; 
    Word32 frame_cnt;
	Word32  source_category;

	Float32 spk_ana_re[D][ORD3];      //buffer of real part of speaker subband 
	Float32 spk_ana_im[D][ORD3];      //buffer of imaginary part of speaker subband 
	Float32 spk_spec[D][ORD3];        //buffer of magnitude spectra of speaker  
	Float32 spk_ener[ORD3];           //buffer of energy of speaker  
    Float32 mic0_power[D];            //PSD of the first mic-signal
    Float32 spk_power[D];             //PSD of speaker signal
    Float32 e0_power[D];              //PSD of error signal
    Float32 mic0_e0_power_re[D];      //real part of CPSD of mic0 and e0
    Float32 mic0_e0_power_im[D];      //imaginary part of CPSD of mic0 and e0 
    Float32 mic0_spk_power_re[D];     //real part of CPSD of mic0 and spk
    Float32 mic0_spk_power_im[D];     //imaginary part of CPSD of mic0 and spk 
    Float32 h_r[NUM][D][ORD3];        //real part of adaptive filter 
	Float32 h_i[NUM][D][ORD3];        //imaginary part of adaptive filter
	Float32 echo_r[D][ORD4];          //real part fo simulative echo
	Float32 echo_i[D][ORD4];	      //imaginary part of simulative echo
	Float32 noise_spectrum[NUM][D];   //spectrum of noise
	Float32 syn[6][ORD2*2];
	Float32 sub_az[D];                //azimuth of subband
    Float32 audio_az[MED_NUM];        //azimuth of fullband
	Float32 ch_enrg[NUM_CHAN];
    Float32 ch_noise[NUM_CHAN];
    Float32 ch_enrg_long_db[NUM_CHAN];
    Float32 qxx_re[D][NUM][NUM];      //covariance matrix of input signals 
	Float32 qxx_im[D][NUM][NUM];      
    Float32 qnn_re[D][NUM][NUM];      //covariance matrix of noise 
	Float32 qnn_im[D][NUM][NUM];

    Float32 pn;
    Float32 hNlMic0SpkAvgMin;  
    Float32 hNlFbMin, hNlFbLocalMin;
    Float32 overDrive, overDriveSm;
	Float32 audio_pit;
	Float32 video_pit;
	Float32 video_az;
	Float32 pearson;
	Float32 alpha;

	Float64 aec_time;
	Float64 srp_time;
	Float64 mvdr_time;
	
	NS_STRUCT my_str_anc;
	Agc_t my_str_agc;
} AEC_SRP_ST;

static int CmpFloat(const void* a, const void* b);
void aec_srp_gsc_init(AEC_SRP_ST *st,Float32 appointed_pitch);
void set_position(AEC_SRP_ST *st, Float32 video_pit,Float32 video_az,Word32 source_category);
Float32 get_position(AEC_SRP_ST *st);
void aec_srp_gsc(AEC_SRP_ST *st,Word16 mic_sp[NUM][FRM_LEN], Word16 spk_sp[FRM_LEN],Word16 out_sp[6][FRM_LEN]);