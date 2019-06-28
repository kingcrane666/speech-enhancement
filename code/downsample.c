#include "downsample.h"

const Word16 lpf3[ORD5]={-63, -112, -119, -22, 179, 399, 497, 373, 63, -248, -336, -109, 290, 548, 410, -93, -612, -709, -215, 572, 1039, 714, -300, -1287, -1387, -294, 1352, 2263, 1415, -1031, -3492, -3746, -339, 6262, 13598, 18408, 18408, 13598, 6262, -339, -3746, -3492, -1031, 1415, 2263, 1352, -294, -1387, -1287, -300, 714, 1039, 572, -215, -709, -612, -93, 410, 548, 290, -109, -336, -248, 63, 373, 497, 399, 179, -22, -119, -112, -63};

void down_sample_init(DOWNSAMPLE_STR *st)      
{
	Word32 i;

	for(i=0; i<DELAY_LEN; i++) st->x[i]=0;	
}

void down_sample_run(
			DOWNSAMPLE_STR *st,      
			Word16 in_sp[FRM_LEN2], 
            Word16 out_sp[FRM_LEN3]  
            )
{
	Word16 delay[FRM_LEN2+DELAY_LEN];
	Word32 acc2, i, j;

	for(j=0; j<DELAY_LEN; j++) delay[j]=st->x[j];
	for(j=0; j<FRM_LEN2; j++) delay[DELAY_LEN+j]=in_sp[j];
	for(j=0; j<DELAY_LEN; j++) st->x[j]=in_sp[FRM_LEN2-DELAY_LEN+j];

	for(i=0; i<FRM_LEN3; i++)
	{
		for(acc2=j=0; j<ORD5; j++) acc2 += delay[3*i+j]*lpf3[j];
		out_sp[i] = (Word16)(acc2>>16);
	}
}

