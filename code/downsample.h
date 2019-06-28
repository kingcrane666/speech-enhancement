typedef short Word16;
typedef long  Word32;

#define FRM_LEN2     480             //length of a frame(48k)
#define FRM_LEN3     160
#define ORD5         72              //order of anti-alias filter
#define DELAY_LEN    (ORD5-1)       

typedef struct
{
    Word16 x[DELAY_LEN];
} DOWNSAMPLE_STR;

void down_sample_init(DOWNSAMPLE_STR *st);
void down_sample_run(DOWNSAMPLE_STR *st, Word16 in_sp[FRM_LEN2], Word16 out_sp[FRM_LEN3]);