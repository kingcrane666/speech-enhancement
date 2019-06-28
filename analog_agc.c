/* analog_agc.c  Using a feedback system, determines an appropriate analog volume level given an input signal and current volume level. Targets a conservative signal level and is intended for use with a digital AGC to apply additional gain. */
#include "typedefs1.h"

#define kMuteGuardTimeMs 8000
/* Default settings if config is not used */
#define AGC_DEFAULT_TARGET_LEVEL 3
#define AGC_DEFAULT_COMP_GAIN 9
/* This is the target level for the analog part in ENV scale. To convert to RMS scale you have to add OFFSET_ENV_TO_RMS. */
#define ANALOG_TARGET_LEVEL 11
#define ANALOG_TARGET_LEVEL_2 5 // ANALOG_TARGET_LEVEL / 2
/* Offset between RMS scale (analog part) and ENV scale (digital part). This value actually varies with the FIXED_ANALOG_TARGET_LEVEL, hence we should in the future replace it with a table. */
#define OFFSET_ENV_TO_RMS 9
/* The reference input level at which the digital part gives an output of targetLevelDbfs (desired level) if we have no compression gain. This level should be set high enough not to compress the peaks due to the dynamics. */
#define DIGITAL_REF_AT_0_COMP_GAIN 4
/* Speed of reference level decrease. */
#define DIFF_REF_TO_ANALOG 5

/* Table for target energy levels. Values in Q(-7)   targetLevelTable = fprintf('%d,\t%d,\t%d,\t%d,\n', round((32767*10.^(-(0:63)'/20)).^2*16/2^7) */
static const int32_t kTargetLevelTable[64] = {134209536, 106606424, 84680493, 67264106, 53429779, 42440782, 33711911, 26778323, 21270778, 16895980, 13420954, 10660642, 8468049, 6726411, 5342978, 4244078, 3371191, 2677832, 2127078, 1689598, 1342095, 1066064, 846805, 672641, 534298, 424408, 337119, 267783, 212708, 168960, 134210, 106606, 84680, 67264, 53430, 42441, 33712, 26778, 21271, 16896, 13421, 10661, 8468, 6726, 5343, 4244, 3371, 2678, 2127, 1690, 1342, 1066, 847, 673, 534, 424, 337, 268, 213, 169, 134, 107, 85, 67};
static const int16_t kSlope1[8] = {21793, 12517, 7189, 4129, 2372, 1362, 472, 78};  // in Q13
static const int16_t kOffset1[8] = {25395, 23911, 22206, 20737, 19612, 18805, 17951, 17367};  // in Q14
static const int16_t kSlope2[8] = {2063, 1731, 1452, 1218, 1021, 857, 597, 337};  // in Q13
static const int16_t kOffset2[8] = {18432, 18379, 18290, 18177, 18052, 17920, 17670, 17286};  // in Q14

void AGC_ZeroCtrl(Agc_t *st, int32_t *inMicLevel, int32_t *env)
{
    int16_t i;
    int32_t midVal, tmp32;

    /* Is the input signal zero? */
    for (tmp32=i=0; i<10; i++) tmp32 += env[i];

    /* Each block is allowed to have a few non-zero samples. */
    if (tmp32 < 500) st->msZero += 10;
    else             st->msZero = 0;

    if (st->muteGuardMs > 0) st->muteGuardMs -= 10;
    if (st->msZero > 500)
    {
        st->msZero = 0;
        /* Increase microphone level only if it's less than 50% */
        midVal = (st->maxAnalog + st->minLevel + 1)>>1;
        if (*inMicLevel < midVal) st->micVol = *inMicLevel = min((*inMicLevel)*1126L>>10, st->zeroCtrlMax); /* *inMicLevel *= 1.1; *//* Reduces risk of a muted mic repeatedly triggering excessive levels due to zero signal detection. */
        st->activeSpeech = 0;
        st->Rxx16_LPw32Max = 0;
        /* The AGC has a tendency (due to problems with the VAD parameters), to vastly increase the volume after a muting event. This timer prevents upwards adaptation for a short period. */
        st->muteGuardMs = kMuteGuardTimeMs;
    }
}

void AGC_ExpCurve(int16_t volume, int16_t *index) // volume in Q14  // index in [0-7]    8 different curves
{
    if (volume > 5243)
    {
        if (volume > 7864)
        {
            if (volume > 12124) *index = 7;
            else                *index = 6;
        }
		else
        {
            if (volume > 6554)  *index = 5;
            else                *index = 4;
        }
    }
	else
    {
        if (volume > 2621)
        {
            if (volume > 3932) *index = 3;
			else               *index = 2;
        }
		else
        {
            if (volume > 1311) *index = 1;
			else               *index = 0;
        }
    }
}


static const uint16_t kGenFuncTable[kGenFuncTableSize] = {256, 485, 786, 1126, 1484, 1849, 2217, 2586, 2955, 3324, 3693, 4063, 4432, 4801, 5171, 5540, 5909, 6279, 6648,  7017,  7387,  7756,  8125,  8495,  8864,  9233,  9603,  9972, 10341, 10711, 11080, 11449, 11819, 12188, 12557, 12927, 13296, 13665, 14035, 14404, 14773, 15143, 15512, 15881, 16251, 16620, 16989, 17359, 17728, 18097, 18466, 18836, 19205, 19574, 19944, 20313, 20682, 21052, 21421, 21790, 22160, 22529, 22898, 23268, 23637, 24006, 24376, 24745, 25114, 25484, 25853, 26222, 26592, 26961, 27330, 27700, 28069, 28438, 28808, 29177, 29546, 29916, 30285, 30654, 31024, 31393, 31762, 32132, 32501, 32870, 33240, 33609, 33978, 34348, 34717, 35086, 35456, 35825, 36194, 36564, 36933, 37302, 37672, 38041, 38410, 38780, 39149, 39518, 39888, 40257, 40626, 40996, 41365, 41734, 42104, 42473, 42842, 43212, 43581, 43950, 44320, 44689, 45058, 45428, 45797, 46166, 46536, 46905};

// decimator
void AGC_DownsampleBy2(const int16_t* in, int16_t len, int16_t* out, int32_t* filtState) 
{
  int32_t tmp1, tmp2, diff, in32;
  int16_t i;

  for(i=0; i<(len>>1); i++)
  {
    // lower allpass filter
    in32 = (int32_t)in[2*i] << 10;
    diff = in32 - filtState[1];
    tmp1 = MUL_ACCUM_1(12199, diff, filtState[0]);
    filtState[0] = in32;
    diff = tmp1 - filtState[2];
    tmp2 = MUL_ACCUM_1(37471, diff, filtState[1]);
    filtState[1] = tmp1;
    diff = tmp2 - filtState[3];
    filtState[3] = MUL_ACCUM_1(60255, diff, filtState[2]);
    filtState[2] = tmp2;

    // upper allpass filter
    in32 = (int32_t)in[2*i+1] << 10;
    diff = in32 - filtState[5];
    tmp1 = MUL_ACCUM_1(3284, diff, filtState[4]);
    filtState[4] = in32;
    diff = tmp1 - filtState[6];
    tmp2 = MUL_ACCUM_1(24441, diff, filtState[5]);
    filtState[5] = tmp1;
    diff = tmp2 - filtState[7];
    filtState[7] = MUL_ACCUM_1(49528, diff, filtState[6]);
    filtState[6] = tmp2;

    // add two allpass outputs, divide by two and round // limit amplitude to prevent wrap-around, and write to output array  //*out++ = AGC_SatW32ToW16(out32);
	out[i] = (int16_t) min(max((filtState[3] + filtState[7] + 1024) >> 11, -32768), 32767);
  }
}

/* Algorithm:
Six term Taylor Series is used here to compute the square root of a number
y^0.5 = (1+x)^0.5 where x = y-1
= 1+(x/2)-0.5*((x/2)^2+0.5*((x/2)^3-0.625*((x/2)^4+0.875*((x/2)^5)
0.5 <= x < 1

Example of how the algorithm works, with ut=sqrt(in), and
with in=73632 and ut=271 (even shift value case):

in=73632
y= in/131072
x=y-1
t = 1 + (x/2) - 0.5*((x/2)^2) + 0.5*((x/2)^3) - 0.625*((x/2)^4) + 0.875*((x/2)^5)
ut=t*(1/sqrt(2))*512

or:

in=73632
in2=73632*2^14
y= in2/2^31
x=y-1
t = 1 + (x/2) - 0.5*((x/2)^2) + 0.5*((x/2)^3) - 0.625*((x/2)^4) + 0.875*((x/2)^5)
ut=t*(1/sqrt(2))
ut2=ut*2^9

which gives:

in  = 73632
in2 = 1206386688
y   = 0.56176757812500
x   = -0.43823242187500
t   = 0.74973506527313
ut  = 0.53014274874797
ut2 = 2.714330873589594e+002

or:

in=73632
in2=73632*2^14
y=in2/2
x=y-2^30
x_half=x/2^31
t = 1 + (x_half) - 0.5*((x_half)^2) + 0.5*((x_half)^3) - 0.625*((x_half)^4)
    + 0.875*((x_half)^5)
ut=t*(1/sqrt(2))
ut2=ut*2^9

which gives:

in  = 73632
in2 = 1206386688
y   = 603193344
x   = -470548480
x_half =  -0.21911621093750
t   = 0.74973506527313
ut  = 0.53014274874797
ut2 = 2.714330873589594e+002
*/

int32_t AGC_Sqrt(int32_t value)
{
    int16_t x_norm, nshift, sh, x_half;
    int32_t A, AA, B, x2;

    if (value == 0)
	{
		return (int32_t)0;
	}
	else
	{
        sh = AGC_NormW32(value);
        A = value<<sh;
        if (A < ((int32_t)0x7fffffff - 32767)) A = A + 32768L; // Round off bit
        else                                   A = (int32_t)0x7fffffff;

        x_norm = (int16_t)(A>>16); // x_norm = AH
        nshift = -(sh>>1); // nshift = sh>>1 // Negate the power for later de-normalization

        A = labs((int32_t)x_norm<<16); // A = abs(x_norm<<16)

        //A = AGC_SqrtLocal(A); // A = sqrt(A)
        /* The following block performs:
         y=A/2
         x=y-2^30
         x_half=x/2^31
         t = 1 + (x_half) - 0.5*((x_half)^2) + 0.5*((x_half)^3) - 0.625*((x_half)^4) + 0.875*((x_half)^5)  */

        B = (A>>1) - ((int32_t)0x40000000); // B = A/2 - 1/2
        x_half = (int16_t)(B>>16);// x_half = x/2 = (A-1)/2
        // B = 1 + x/2
        B = B + 0x40000000 + 0x40000000; // Add 0.5 twice (since 1.0 does not exist A Q31)

        x2 = (int32_t)x_half*x_half; // AA = (x/2)^2
        B = B - x2; // B = 1 + x/2 - 0.5*(x/2)^2

        AA = (-x2>>15)*(-x2>>15); // AA = (x/2)^4
        // B = B - 0.625*AA  // After this, B = 1 + x/2 - 0.5*(x/2)^2 - 0.625*(x/2)^4
        // AA = (x/2)^5
        // B = B + 0.875*AA    // After this, B = 1 + x/2 - 0.5*(x/2)^2 - 0.625*(x/2)^4 + 0.875*(x/2)^5
        // B = B + 0.5*AA  // After this, B = 1 + x/2 - 0.5*(x/2)^2 + 0.5*(x/2)^3 - 0.625*(x/2)^4 + 0.875*(x/2)^5
        A = B + -20480L*(AA>>15)*2 + 28672L*(x_half*(AA>>15)>>15)*2 + x_half*(x2>>15) + 32768; // Round off bit

    	// Even shift value case
        if (-2*nshift == sh) A = ((11585L*(A>>16)+8192) & 0x1fffc000)>>13; // Round off  // A = 1/sqrt(2)*t16  // 1/sqrt2 (==5a82)
	    else                 A = A>>16;
        A = (A&65535L)>>(-nshift); // De-normalize the result
        return A;
	}
}


// This function generates the compressor gain table used in the fixed digital part.
int32_t AGC_CalculateGainTable(int32_t *gainTable, int16_t digCompGaindB, int16_t targetLevelDbfs, uint8_t limiterEnable, int16_t analogTarget)
{
    uint32_t tmpU32no1, tmpU32no2, absInLevel, logApprox;
    int32_t inLevel, limiterLvl, tmp32, tmp32no1, tmp32no2, numFIX, den, y32;
    uint16_t constMaxGain, tmpU16, intPart, fracPart;
    int16_t limiterOffset = 0; // Limiter offset
    int16_t limiterIdx, constLinApprox, zeroGainLvl, maxGain, diffGain, i, tmp16no1;
    int zeros, zerosScale;

    // Calculate maximum digital gain and zero gain level
    tmp32no1 = (int32_t)(digCompGaindB - analogTarget)*2;
    tmp16no1 = analogTarget - targetLevelDbfs;
    tmp16no1 += (tmp32no1 + 1)/3;
    maxGain = max(tmp16no1, (analogTarget - targetLevelDbfs));
    tmp32no1 = maxGain*3;
    zeroGainLvl = digCompGaindB;
    zeroGainLvl -= (tmp32no1+1)>>1;
    if ((digCompGaindB <= analogTarget) && (limiterEnable))
    {
        zeroGainLvl += (analogTarget - digCompGaindB + kSoftLimiterLeft);
        limiterOffset = 0;
    }

    // Calculate the difference between maximum gain and gain at 0dB0v:  diffGain = maxGain + (compRatio-1)*zeroGainLvl/compRatio = (compRatio-1)*digCompGaindB/compRatio
    tmp32no1 = digCompGaindB*2;
    diffGain = (tmp32no1 + 1)/3;
    if (diffGain < 0 || diffGain >= kGenFuncTableSize) {assert(0);  return -1;}

    // Calculate the limiter level and index:
    //  limiterLvlX = analogTarget - limiterOffset
    //  limiterLvl  = targetLevelDbfs + limiterOffset/compRatio
    limiterIdx = 2 + (short)((analogTarget - limiterOffset)*2048L/6165);  //49321 = 10*log10(2)  in Q14
    tmp16no1 = (limiterOffset + 1)/ 3;
    limiterLvl = targetLevelDbfs + tmp16no1;

    // Calculate (through table lookup): constMaxGain = log2(1+2^(log2(e)*diffGain)); (in Q8)
    constMaxGain = kGenFuncTable[diffGain]; // in Q8

    // Calculate a parameter used to approximate the fractional part of 2^x with a piecewise linear function in Q14: constLinApprox = round(3/2*(4*(3-2*sqrt(2))/(log(2)^2)-0.5)*2^14);
    constLinApprox = 22817; // in Q14

    // Calculate a denominator used in the exponential part to convert from dB to linear scale: den = 20*constMaxGain (in Q8)
    den = 20L*constMaxGain; // in Q8

    for (i = 0; i < 32; i++)
    {
        // Calculate scaled input level (compressor): inLevel = fix((-constLog10_2*(compRatio-1)*(1-i)+fix(compRatio/2))/compRatio)
        inLevel = (2*(i-1)*49321L + 1)/3; // Q14 //49321 = 10*log10(2)  in Q14

        // Calculate diffGain-inLevel, to map using the genFuncTable
        inLevel = ((int32_t)diffGain<<14) - inLevel; // Q14
        // Make calculations on abs(inLevel) and compensate for the sign afterwards.
        absInLevel = labs(inLevel); // Q14

        // LUT with interpolation
        intPart = (uint16_t)(absInLevel>>14);
        fracPart = (uint16_t)(absInLevel & 0x00003FFF); // extract the fractional part
        tmpU16 = kGenFuncTable[intPart+1] - kGenFuncTable[intPart]; // Q8
        tmpU32no1 = (uint32_t)tmpU16*fracPart + ((uint32_t)kGenFuncTable[intPart]<<14); // Q22
        logApprox = tmpU32no1>>8; // Q14

        // Compensate for negative exponent using the relation:  //  log2(1 + 2^-x) = log2(1 + 2^x) - x
        if (inLevel < 0)
        {
            zeros = AGC_NormU32(absInLevel);
            zerosScale = 0;
            if (zeros < 15)
            {
                // Not enough space for multiplication
                tmpU32no2 = (absInLevel>>(15-zeros))*23637L; // Q(zeros+13)   //23637 = log2(e)      in Q14
                if (zeros < 9)
                {
                    tmpU32no1 = tmpU32no1>>(9-zeros); // Q(zeros+13)
                    zerosScale = 9-zeros;
                }
				else
                {
                    tmpU32no2 = tmpU32no2>>(zeros-9); // Q22
                }
            }
			else
            {
                tmpU32no2 = (uint32_t)absInLevel*23637L>>6; // Q22    //23637 = log2(e)      in Q14
            }
            logApprox = 0;
            if (tmpU32no2 < tmpU32no1) logApprox = (tmpU32no1 - tmpU32no2)>>(8 - zerosScale); //Q14
        }
        numFIX = ((int32_t)maxGain*constMaxGain<<6) - (int32_t)logApprox*diffGain; // Q14

        // Calculate ratio
        // Shift |numFIX| as much as possible. // Ensure we avoid wrap-around in |den| as well.
        if (numFIX > (den >> 8)) zeros = AGC_NormW32(numFIX); // |den| is Q8.
        else                     zeros = AGC_NormW32(den) + 8;
        numFIX = (numFIX<<zeros); // Q(14+zeros)

        // Shift den so we end up in Qy1
		tmp32no1 = (zeros>=8) ? (den<<(zeros-8)) : (den>>(8-zeros));  // Q(zeros)
        if (numFIX < 0) numFIX -= tmp32no1>>1;
        else            numFIX += tmp32no1>>1;

        y32 = numFIX/tmp32no1; // in Q14
        if (limiterEnable && (i < limiterIdx)) y32 = ((i-1)*49321L - limiterLvl*16384L + 10)/20;  //49321 = 10*log10(2)  in Q14
        if (y32 > 39000) tmp32 = ((y32>>1)*27213 + 2048)>>12; // in Q14   //54426 = log2(10)     in Q14
		else             tmp32 = (y32*27213 + 4096)>>13; // in Q14        //54426 = log2(10)     in Q14
        tmp32 += 262144L; //16 in Q14 (Make sure final output is in Q16)

        // Calculate power
        if (tmp32 > 0)
        {
            intPart = (int16_t)(tmp32>>14);
            fracPart = (tmp32 & 16383); // in Q14
            if (fracPart>>13) tmp32no2 = 16384 - ((16384-fracPart)*(32768L-constLinApprox)>>13);
			else              tmp32no2 = (fracPart*(constLinApprox-16384)>>13);
            fracPart = (uint16_t)tmp32no2;
            gainTable[i] = (1L<<intPart) + ((intPart>=14) ? (fracPart<<(intPart-14)) : (fracPart>>(14-intPart)));
        }
		else
        {
            gainTable[i] = 0;
        }
    }

    return 0;
}

int16_t AGC_ProcessVad(AgcVad_t *state, // (i) VAD state
                                   const int16_t *in, // (i) Speech signal
                                   int16_t nrSamples) // (i) number of samples
{
    int16_t k, subfr, buf1[8], buf2[4], zeros, dB;
	int32_t out, nrg, tmp32;    

    // process in 10 sub frames of 1 ms (to save on memory)
    for (nrg = subfr = 0; subfr < 10; subfr++)
    {        
        if (nrSamples == 160) // downsample to 4 kHz
        {
            for (k = 0; k < 8; k++) buf1[k] = (int16_t)(((int32_t)in[16*subfr+2*k] + in[16*subfr+2*k+1])>>1);
            AGC_DownsampleBy2(buf1, 8, buf2, state->downState);
        }
		else
        {
            AGC_DownsampleBy2(in+8*subfr, 8, buf2, state->downState);
            in += 8;
        }

        // high pass filter and compute energy
        for (k = 0; k < 4; k++)
        {
            out = buf2[k] + state->HPstate;
            state->HPstate = (int16_t)((75L*out>>7) - buf2[k]);
            nrg += (int32_t)out*out>>6;
        }
    }

    // find number of leading zeros
    if (!(0xFFFF0000 & nrg)) zeros = 16;
	else                     zeros = 0;

    if (!(0xFF000000 & (nrg << zeros))) zeros += 8;
    if (!(0xF0000000 & (nrg << zeros))) zeros += 4;
    if (!(0xC0000000 & (nrg << zeros))) zeros += 2;
    if (!(0x80000000 & (nrg << zeros))) zeros += 1;

    // energy level (range {-32..30}) (Q10)
    dB = (15-zeros)<<11;
    // Update statistics // decay time = AvgDecTime * 10 ms
    if (state->counter < kAvgDecayTime) state->counter++;
    // update short-term estimate of mean energy level (Q10)
    state->meanShortTerm = (int16_t)((state->meanShortTerm*15L + dB)>>4);
    // update short-term estimate of variance in energy level (Q8)
    state->varianceShortTerm = (((int32_t)dB*dB>>12) + state->varianceShortTerm*15L)>>4;
    // update short-term estimate of standard deviation in energy level (Q10)
    state->stdShortTerm = (int16_t)AGC_Sqrt((state->varianceShortTerm<<12) - (int32_t)state->meanShortTerm*state->meanShortTerm);

    // update long-term estimate of mean energy level (Q10)
    tmp32 = (int32_t)state->meanLongTerm*state->counter + dB;
    state->meanLongTerm = tmp32/(short)min(state->counter+1L, 32767);

    // update long-term estimate of variance in energy level (Q8)
    tmp32 = ((int32_t)dB*dB>>12) + state->varianceLongTerm*state->counter;
    state->varianceLongTerm = tmp32/(short)min(state->counter+1L, 32767);

    // update long-term estimate of standard deviation in energy level (Q10)
    state->stdLongTerm = (int16_t)AGC_Sqrt((state->varianceLongTerm<<12) - (int32_t)state->meanLongTerm*state->meanLongTerm);

    // update voice activity measure (Q10)
    if (state->stdLongTerm != 0) tmp32 = 12288L*(dB - state->meanLongTerm)/state->stdLongTerm;
    else                         tmp32 = 0x7FFFFFFF;
    tmp32 += state->logRatio*52L;
    // limit
    state->logRatio = (int16_t)min(max(tmp32>>6, -2048), 2048);

    return state->logRatio; // Q10
}

int32_t AGC_ProcessDigital(DigitalAgc_t *stt, const int16_t *in_near, const int16_t *in_near_H, int16_t *out, int16_t *out_H, uint32_t FS, int16_t lowlevelSignal)
{
    // array for gains (one value per ms, incl start & end)
    int32_t gains[11], env[10];
    int32_t tmp32, nrg, max_nrg, cur_level, gain32, delta;
    int16_t logratio, lower_thr, upper_thr;
    int16_t zeros, zeros_fast, frac, decay;
    int16_t gate, gain_adj;
    int16_t k, n;
    int16_t L, L2; // samples/subframe

    // determine number of samples per ms
    if (FS == 8000)       {L = 8;        L2 = 3;}
	else if (FS == 16000) {L = 16;       L2 = 4;}
	else if (FS == 32000) {L = 16;       L2 = 4;}
	else                  {return -1;}

    // TODO(andrew): again, we don't need input and output pointers...
    if (in_near != out) memcpy(out, in_near, 10*L*sizeof(int16_t)); // Only needed if they don't already point to the same place.
    if ( (FS == 32000)&&(in_near_H != out_H) ) memcpy(out_H, in_near_H, 10 * L * sizeof(int16_t));

    // VAD for near end
    logratio = AGC_ProcessVad(&stt->vadNearend, out, (int16_t)(L*10));

    // Account for far end VAD
    if (stt->vadFarend.counter > 10) logratio = (int16_t)((3L*logratio - stt->vadFarend.logRatio)>>2);

    // Determine decay factor depending on VAD  //  upper_thr = 1.0f;    //  lower_thr = 0.25f;
    upper_thr = 1024; // Q10
    lower_thr = 0; // Q10
    if (logratio > upper_thr)      decay = -65; // decay = -2^17 / DecayTime;  ->  -65
	else if (logratio < lower_thr) decay = 0;
	else                           decay = (int16_t)((lower_thr - logratio)*65L>>10); // decay = (int16_t)(((lower_thr - logratio) * (2^27/(DecayTime*(upper_thr-lower_thr)))) >> 10);     // SUBSTITUTED: 2^27/(DecayTime*(upper_thr-lower_thr))  ->  65

    // adjust decay factor for long silence (detected as low standard deviation)
    // This is only done in the adaptive modes
    if (stt->agcMode != kAgcModeFixedDigital)
    {
        if (stt->vadNearend.stdLongTerm < 4000)      decay = 0;
		else if (stt->vadNearend.stdLongTerm < 8096) decay = (int16_t)((int32_t)(stt->vadNearend.stdLongTerm - 4000)*decay>>12);
        if (lowlevelSignal != 0) decay = 0;
    }
    // Find max amplitude per sub frame   // iterate over sub frames
    for (k = 0; k < 10; k++)
    {
        for (max_nrg =n = 0; n < L; n++) // iterate over samples
        {
            nrg = (int32_t)out[k*L+n]*out[k*L+n];
            if (nrg > max_nrg) max_nrg = nrg;
        }
        env[k] = max_nrg;
    }

    // Calculate gain per sub frame
    gains[0] = stt->gain;
    for (k = 0; k < 10; k++)
    {
        // Fast envelope follower
        //  decay time = -131000 / -1000 = 131 (ms)
        stt->capacitorFast = max(stt->capacitorFast + (stt->capacitorFast>>16)*-1000 + ((0x0000FFFF&stt->capacitorFast)*-1000>>16), env[k]);
        // Slow envelope follower
        if (env[k] > stt->capacitorSlow) stt->capacitorSlow = AGC_SCALEDIFF32(500, (env[k] - stt->capacitorSlow), stt->capacitorSlow);  // increase capacitorSlow
		else                             stt->capacitorSlow = stt->capacitorSlow + (stt->capacitorSlow>>16)*decay + ((0x0000FFFF&stt->capacitorSlow)*decay>>16); // decrease capacitorSlow
        // use maximum of both capacitors as current level
        cur_level = max(stt->capacitorFast, stt->capacitorSlow);
        // Translate signal level into gain, using a piecewise linear approximation
        // find number of leading zeros
        zeros = AGC_NormU32((uint32_t)cur_level);
        if (cur_level == 0) zeros = 31;
        tmp32 = (cur_level<<zeros) & 0x7FFFFFFF;
        frac = (int16_t)(tmp32>>19); // Q12
        tmp32 = (int32_t)(stt->gainTable[zeros-1] - stt->gainTable[zeros])*frac;
        gains[k + 1] = stt->gainTable[zeros] + (tmp32>>12);
    }

    // Gate processing (lower gain during absence of speech)
    zeros = (zeros<<9) - (frac>>3);
    // find number of leading zeros
    zeros_fast = AGC_NormU32((uint32_t)stt->capacitorFast);
    if (stt->capacitorFast == 0) zeros_fast = 31;
    tmp32 = (stt->capacitorFast<<zeros_fast) & 0x7FFFFFFF;
    zeros_fast = zeros_fast<<9;
    zeros_fast -= (int16_t)(tmp32>>22);

    gate = 1000 + zeros_fast - zeros - stt->vadNearend.stdShortTerm;

    if (gate < 0) stt->gatePrevious = 0;
	else          stt->gatePrevious = gate = (int16_t)((gate + stt->gatePrevious*7L)>>3);

    //gate<0  -> no gate   //gate>2500  -> max gate
    if (gate > 0)
    {
        if (gate < 2500) gain_adj = (2500 - gate)>>5;
		else             gain_adj = 0;

        for (k = 0; k < 10; k++)
        {
            if ((gains[k+1] - stt->gainTable[0]) > 8388608) tmp32 = ((gains[k+1] - stt->gainTable[0])>>8)*(178 + gain_adj); // To prevent wraparound
			else                                            tmp32 = (gains[k+1] - stt->gainTable[0])*(178 + gain_adj)>>8;
            gains[k+1] = stt->gainTable[0] + tmp32;
        }
    }

    // Limit gain to avoid overload distortion
    for (k = 0; k < 10; k++)
    {
        // To prevent wrap around
        zeros = 10;
        if (gains[k+1] > 47453132L) zeros = 16 - AGC_NormW32(gains[k+1]);
        gain32 = (gains[k+1]>>zeros) + 1;
        gain32 = (int32_t)gain32*gain32;
        // check for overflow
        while ( (gain32>>13)*((env[k]>>12)+1) + (((8191&gain32)*((env[k]>>12)+1))>>13) > ((11>=zeros) ? (32767L<<2*(11-zeros)) : (32767L>>2*(zeros-11))) )
        {
            // multiply by 253/256 ==> -0.1 dB // Prevent wrap around
            if (gains[k+1] > 8388607) gains[k+1] = (gains[k+1]>>8)*253L; 
			else                      gains[k+1] = (gains[k+1]*253L>>8);
            gain32 = (gains[k+1]>>zeros) + 1;
            gain32 = gain32*gain32;
        }
    }
    // gain reductions should be done 1 ms earlier than gain increases
    for (k = 1; k < 10; k++)
    {
        if (gains[k] > gains[k+1]) gains[k] = gains[k+1];
    }
    // save start gain for next frame
    stt->gain = gains[10];

    // Apply gain
    // handle first sub frame separately
    delta = (gains[1] - gains[0])<<(4 - L2);
    gain32 = gains[0]<<4;
    // iterate over samples
    for (n=0; n<L; n++)
    {
        // For lower band
        out[n] = (int16_t)max(min((out[n]*(gain32>>4)>>16), 32767), -32768);
        // For higher band
        if (FS == 32000) out_H[n] = (int16_t)max(min((out_H[n]*(gain32>>4)>>16), 32767), -32768);
        gain32 += delta;
    }
    // iterate over subframes
    for (k = 1; k < 10; k++)
    {
        delta = (gains[k+1]-gains[k])<<(4 - L2);
        gain32 = (gains[k]<<4);
        // iterate over samples
        for (n = 0; n < L; n++)
        {
            out[k*L+n] = (int16_t)((int32_t)out[k*L+n]*(gain32>>4)>>16); // For lower band
            if (FS == 32000) out_H[k*L+n] = (int16_t)((int32_t)out_H[k*L+n]*(gain32>>4)>>16); // For higher band
            gain32 += delta;
        }
    }

    return 0;
}

void AGC_InitVad(AgcVad_t *state)
{
    int16_t k;

    state->HPstate = 0; // state of high pass filter
    state->logRatio = 0; // log( P(active) / P(inactive) )
    // average input level (Q10)
    state->meanLongTerm = 15<<10;
    // variance of input level (Q8)
    state->varianceLongTerm = 500<<8;
    state->stdLongTerm = 0; // standard deviation of input level in dB
    // short-term average input level (Q10)
    state->meanShortTerm = 15<<10;
    // short-term variance of input level (Q8)
    state->varianceShortTerm = 500<<8;
    state->stdShortTerm = 0; // short-term standard deviation of input level in dB
    state->counter = 3; // counts updates
    for (k = 0; k < 8; k++) state->downState[k] = 0; // downsampling filter
}

int32_t AGC_ProcessAnalog(void *state, int32_t inMicLevel, int32_t *outMicLevel, int16_t vadLogRatio, int16_t echo, uint8_t *saturationWarning)
{
    int32_t Rxx16w32, inMicLevelTmp, lastMicVol;
    int16_t i, tmpW16, vadThresh;
    uint8_t saturated = 0;
    Agc_t *st;

    st = (Agc_t *)state;
    inMicLevelTmp = (inMicLevel<<st->scale);

    if (inMicLevelTmp > st->maxAnalog)     return -1;
	else if (inMicLevelTmp < st->minLevel) return -1;

    if (st->firstCall == 0)
    {
        int32_t tmpVol;
        st->firstCall = 1;
        tmpVol = st->minLevel + ((st->maxLevel - st->minLevel)*51L>>9);
        /* If the mic level is very low at start, increase it! */
        if ((inMicLevelTmp < tmpVol) && (st->agcMode == kAgcModeAdaptiveAnalog)) inMicLevelTmp = tmpVol;
        st->micVol = inMicLevelTmp;
    }

    /* Set the mic level to the previous output value if there is digital input gain */
    if ((inMicLevelTmp == st->maxAnalog) && (st->micVol > st->maxAnalog)) inMicLevelTmp = st->micVol;
    /* If the mic level was manually changed to a very low value raise it! */
    if ((inMicLevelTmp != st->micVol) && (inMicLevelTmp < st->minOutput)) st->micVol = inMicLevelTmp = st->minLevel + ((st->maxLevel - st->minLevel)*51L>>9);
    if (inMicLevelTmp != st->micVol) st->micVol = inMicLevelTmp; // Incoming level mismatch; update our level. This could be the case if the volume is changed manually, or if the sound device has a low volume resolution.
    if (inMicLevelTmp > st->maxLevel) st->maxLevel = inMicLevelTmp; // Always allow the user to raise the volume above the maxLevel.

    // Store last value here, after we've taken care of manual updates etc.
    lastMicVol = st->micVol;

    /* Checks if the signal is saturated. Also a check if individual samples are larger than 12000 is done. If they are the counter for increasing the volume level is set to -100ms */
    //AGC_SaturationCtrl(st, &saturated, st->env[0]);
    for (i = 0; i < 10; i++)
    {
        tmpW16 = (int16_t)(st->env[0][i]>>20);
        if (tmpW16 > 875) st->envSum += tmpW16;
    }
    if (st->envSum > 25000)
    {
        saturated = 1;
        st->envSum = 0;
    }    
    st->envSum = (int16_t)(st->envSum*4055L>>12); /* st->envSum *= 0.99; */

    /* The AGC is always allowed to lower the level if the signal is saturated */
    if (saturated == 1)
    {
        /* Lower the recording level
		 * Rxx160_LP is adjusted down because it is so slow it could cause the AGC to make wrong decisions. */
        st->Rxx160_LPw32 = (st->Rxx160_LPw32>>3)*7L; // st->Rxx160_LPw32 *= 0.875;
        st->zeroCtrlMax = st->micVol;
        /* st->micVol *= 0.903; */
        inMicLevelTmp = st->micVol = min((29591L*(inMicLevelTmp - st->minLevel)>>15) + st->minLevel, lastMicVol - 2);
        if (st->micVol < st->minOutput) *saturationWarning = 1;

        /* Reset counter for decrease of volume level to avoid decreasing too much. The saturation control can still lower the level if needed. */
        st->msTooHigh = -100;
        /* Enable the control mechanism to ensure that our measure, Rxx160_LP, is in the correct range. This must be done since the measure is very slow. */
        st->activeSpeech = 0;
        st->Rxx16_LPw32Max = 0;

        /* Reset to initial values */
        st->msecSpeechInnerChange = kMsecSpeechInner;
        st->msecSpeechOuterChange = kMsecSpeechOuter;
        st->changeToSlowMode = 0;
        st->muteGuardMs = 0;
        st->upperLimit = st->startUpperLimit;
        st->lowerLimit = st->startLowerLimit;
    }

    /* Check if the input speech is zero. If so the mic volume is increased. On some computers the input is zero up as high level as 17% */
    AGC_ZeroCtrl(st, &inMicLevelTmp, st->env[0]);

    /* Check if the near end speaker is inactive. If that is the case the VAD threshold is increased since the VAD speech model gets more sensitive to any sound after a long silence. */
    //AGC_SpeakerInactiveCtrl(st);
    if (st->vadMic.stdLongTerm < 2500)
    {
        st->vadThreshold = 1500;
    }
	else
    {
        vadThresh = kNormalVadThreshold;
        if (st->vadMic.stdLongTerm < 4500) vadThresh += (4500 - st->vadMic.stdLongTerm)>>1; /* Scale between min and max threshold */
        st->vadThreshold = (int16_t)((vadThresh + 31L*st->vadThreshold)>>5); /* st->vadThreshold = (31 * st->vadThreshold + vadThresh) / 32; */
    }

    for (i = 0; i < 5; i++)
    {
        /* Computed on blocks of 16 samples */
        Rxx16w32 = st->Rxx16w32_array[0][i];

        /* Rxx160w32 in Q(-7) */
        st->Rxx160w32 = st->Rxx160w32 + ((Rxx16w32 - st->Rxx16_vectorw32[st->Rxx16pos])>>3);
        st->Rxx16_vectorw32[st->Rxx16pos] = Rxx16w32;

        /* Circular buffer */
        st->Rxx16pos++;
        if (st->Rxx16pos == RXX_BUFFER_LEN) st->Rxx16pos = 0;

        /* Rxx16_LPw32 in Q(-4) */
        st->Rxx16_LPw32 = st->Rxx16_LPw32 + ((Rxx16w32 - st->Rxx16_LPw32)>>kAlphaShortTerm);

        if (vadLogRatio > st->vadThreshold)
        {
            /* Speech detected! */
            /* Check if Rxx160_LP is in the correct range. If it is too high/low then we set it to the maximum of Rxx16_LPw32 during the first 200ms of speech. */
            if (st->activeSpeech < 250)
            {
                st->activeSpeech += 2;
                st->Rxx16_LPw32Max = max(st->Rxx16_LPw32Max, st->Rxx16_LPw32);
            }
			else if (st->activeSpeech == 250)
            {
                st->activeSpeech += 2;
                st->Rxx160_LPw32 = (int32_t)(st->Rxx16_LPw32Max>>3)*RXX_BUFFER_LEN;
            }

            st->Rxx160_LPw32 = st->Rxx160_LPw32 + ((st->Rxx160w32 - st->Rxx160_LPw32)>>kAlphaLongTerm);

            if (st->Rxx160_LPw32 > st->upperSecondaryLimit)
            {
                st->msTooHigh += 2;
                st->msTooLow = 0;
                st->changeToSlowMode = 0;

                if (st->msTooHigh > st->msecSpeechOuterChange)
                {
                    st->msTooHigh = 0;

                    /* Lower the recording level */
                    /* Multiply by 0.828125 which corresponds to decreasing ~0.8dB */
                    st->Rxx160_LPw32 = (st->Rxx160_LPw32>>6)*53L;

                    /* Reduce the max gain to avoid excessive oscillation (but never drop below the maximum analog level).
                     * st->maxLevel = (15 * st->maxLevel + st->micVol) / 16; */
                    st->maxLevel = max((15*st->maxLevel + st->micVol)>>4, st->maxAnalog);
                    st->zeroCtrlMax = st->micVol;
                    /* 0.95 in Q15 */
                    inMicLevelTmp = st->micVol = min((31130L*(inMicLevelTmp - st->minLevel)>>15) + st->minLevel, lastMicVol - 1);
                    /* Enable the control mechanism to ensure that our measure, Rxx160_LP, is in the correct range.  */
                    st->activeSpeech = 0;
                    st->Rxx16_LPw32Max = 0;
                }
            }
			else if (st->Rxx160_LPw32 > st->upperLimit)
            {
                st->msTooHigh += 2;
                st->msTooLow = 0;
                st->changeToSlowMode = 0;

                if (st->msTooHigh > st->msecSpeechInnerChange)
                {
                    /* Lower the recording level */
                    st->msTooHigh = 0;                    
                    st->Rxx160_LPw32 = (st->Rxx160_LPw32>>6)*53L; /* Multiply by 0.828125 which corresponds to decreasing ~0.8dB */
                    /* Reduce the max gain to avoid excessive oscillation (but never drop below the maximum analog level). st->maxLevel = (15 * st->maxLevel + st->micVol) / 16; */
                    st->maxLevel = max((15*st->maxLevel + st->micVol)>>4, st->maxAnalog);
                    st->zeroCtrlMax = st->micVol;
                    /* 0.965 in Q15 */
                    inMicLevelTmp = st->micVol = min((31621L*(inMicLevelTmp - st->minLevel)>>15) + st->minLevel, lastMicVol - 1);
                }
            } 
			else if (st->Rxx160_LPw32 < st->lowerSecondaryLimit)
            {
                st->msTooHigh = 0;
                st->changeToSlowMode = 0;
                st->msTooLow += 2;

                if (st->msTooLow > st->msecSpeechOuterChange)
                {
                    /* Raise the recording level */
                    int16_t index, weightFIX;
                    int16_t volNormFIX = 16384; // =1 in Q14.
                    st->msTooLow = 0;
                    /* Normalize the volume level */
                    if (st->maxInit != st->minLevel) volNormFIX = (int16_t)(((inMicLevelTmp - st->minLevel)<<14)/(st->maxInit - st->minLevel));
                    /* Find correct curve */
                    AGC_ExpCurve(volNormFIX, &index);

                    /* Compute weighting factor for the volume increase, 32^(-2*X)/2+1.05 */
                    weightFIX = kOffset1[index] - (int16_t)((int32_t)kSlope1[index]*volNormFIX>>13);
                    /* st->Rxx160_LPw32 *= 1.047 [~0.2 dB]; */
                    st->Rxx160_LPw32 = (uint32_t)(st->Rxx160_LPw32>>6)*67L;
                    inMicLevelTmp = st->micVol = (int32_t)min(((uint32_t)weightFIX*(inMicLevelTmp - st->minLevel)>>14) 
						+ st->minLevel, lastMicVol + 2);
                }
            } 
			else if (st->Rxx160_LPw32 < st->lowerLimit)
            {
                st->msTooHigh = 0;
                st->changeToSlowMode = 0;
                st->msTooLow += 2;

                if (st->msTooLow > st->msecSpeechInnerChange)
                {
                    /* Raise the recording level */
                    int16_t index, weightFIX, volNormFIX;

                    st->msTooLow = 0;
					volNormFIX = 16384; // =1 in Q14.
                    /* Normalize the volume level */
                    if (st->maxInit != st->minLevel) volNormFIX = (int16_t)(((inMicLevelTmp - st->minLevel)<<14)/(st->maxInit - st->minLevel));
                    /* Find correct curve */
                    AGC_ExpCurve(volNormFIX, &index);
                    weightFIX = kOffset2[index] - (int16_t)((int32_t)kSlope2[index]*volNormFIX>>13); /* Compute weighting factor for the volume increase, (3.^(-2.*X))/8+1 */
                    /* st->Rxx160_LPw32 *= 1.047 [~0.2 dB]; */
                    st->Rxx160_LPw32 = (st->Rxx160_LPw32>>6)*67L;
                    inMicLevelTmp = st->micVol = max(((uint32_t)weightFIX*(inMicLevelTmp - st->minLevel)>>14) + st->minLevel, lastMicVol + 1);
                }
            } 
			else
            {
                if (st->changeToSlowMode > 4000) /* The signal is inside the desired range which is: lowerLimit < Rxx160_LP/640 < upperLimit */
                {
                    st->msecSpeechInnerChange = 1000;
                    st->msecSpeechOuterChange = 500;
                    st->upperLimit = st->upperPrimaryLimit;
                    st->lowerLimit = st->lowerPrimaryLimit;
                } 
				else
                {
                    st->changeToSlowMode += 2; // in milliseconds
                }
                st->msTooLow = 0;
                st->msTooHigh = 0;
                st->micVol = inMicLevelTmp;
            }
        }
    }

    /* Ensure gain is not increased in presence of echo or after a mute event
     * (but allow the zeroCtrl() increase on the frame of a mute detection).
     */
    if (echo == 1 || (st->muteGuardMs > 0 && st->muteGuardMs < kMuteGuardTimeMs))
    {
        if (st->micVol > lastMicVol) st->micVol = lastMicVol;
    }

    /* limit the gain */
    if (st->micVol > st->maxLevel)       st->micVol = st->maxLevel;
    else if (st->micVol < st->minOutput) st->micVol = st->minOutput;

    *outMicLevel = st->micVol>>st->scale;
    if (*outMicLevel > (st->maxAnalog>>st->scale)) *outMicLevel = st->maxAnalog>>st->scale;

    return 0;
}

/*
 This function processes a 10/20ms frame and adjusts (normalizes) the gain both analog and digitally. 
 The gain adjustments are done only during active periods of speech. 
 The input speech length can be either 10ms or 20ms and the output is of the same length. 
 The length of the speech vectors must be given in samples (80/160 when FS=8000, and 160/320 when FS=16000 or FS=32000). 
 The echo parameter can be used to ensure the AGC will not adjust upward in the presence of echo.
 
 This function should be called after processing the near-end microphone signal, in any case after any echo cancellation.
 
 * Input:
 *      - agcInst           : AGC instance
 *      - inNear            : Near-end input speech vector (10 or 20 ms) for
 *                            L band
 *      - inNear_H          : Near-end input speech vector (10 or 20 ms) for
 *                            H band
 *      - samples           : Number of samples in input/output vector
 *      - inMicLevel        : Current microphone volume level
 *      - echo              : Set to 0 if the signal passed to add_mic is almost certainly free of echo; otherwise set to 1. If you have no information regarding echo set to 0.
 *
 * Output:
 *      - outMicLevel       : Adjusted microphone volume level
 *      - out               : Gain-adjusted near-end speech vector (L band) May be the same vector as the input.
 *      - out_H             : Gain-adjusted near-end speech vector (H band)
 *      - saturationWarning : A returned value of 1 indicates a saturation event has occurred and the volume cannot be further reduced. Otherwise will be set to 0.
 *
 * Return value:
 *                          :  0 - Normal operation.
 *                          : -1 - Error
 */
int AGC_Process(Agc_t *st, const int16_t *in_near, const int16_t *in_near_H, int16_t samples, int16_t *out, int16_t *out_H, int32_t inMicLevel, int32_t *outMicLevel, int16_t echo, uint8_t *saturationWarning)
{
    int32_t inMicLevelTmp;
    int16_t subFrames, i;
    uint8_t satWarningTmp = 0;

    if (st->fs == 8000)       subFrames = 80;
	else if (st->fs == 16000) subFrames = 160;
	else if (st->fs == 32000) subFrames = 160;
	else                      return -1;

    /* Check for valid pointers based on sampling rate */
    if (st->fs == 32000 && in_near_H == NULL) return -1;

    *saturationWarning = 0;
    //TODO: PUT IN RANGE CHECKING FOR INPUT LEVELS
    *outMicLevel = inMicLevelTmp = inMicLevel;

    // TODO(andrew): clearly we don't need input and output pointers...
    //   Change the interface to take a shared input/output.
    if (in_near != out) memcpy(out, in_near, samples * sizeof(int16_t)); // Only needed if they don't already point to the same place.
    if ( (st->fs == 32000)&&(in_near_H != out_H) ) memcpy(out_H, in_near_H, samples * sizeof(int16_t));
    for (i = 0; i < samples; i += subFrames)
    {
        if (AGC_ProcessDigital(&st->digitalAgc, &in_near[i], &in_near_H[i], &out[i], &out_H[i], st->fs, st->lowLevelSignal) == -1) return -1;
        if ((st->agcMode < kAgcModeFixedDigital) && ((st->lowLevelSignal == 0) || (st->agcMode != kAgcModeAdaptiveDigital)))
        {
            if (AGC_ProcessAnalog(st, inMicLevelTmp, outMicLevel, st->vadMic.logRatio, echo, saturationWarning) == -1) return -1;
        }
        /* update queue */
        if (st->inQueue > 1)
        {
            memcpy(st->env[0], st->env[1], 10*sizeof(int32_t));
            memcpy(st->Rxx16w32_array[0], st->Rxx16w32_array[1], 5*sizeof(int32_t));
        }

        if (st->inQueue > 0) st->inQueue--;
        /* If 20ms frames are used the input mic level must be updated so that the analog AGC does not think that there has been a manual volume change. */
        inMicLevelTmp = *outMicLevel;
        /* Store a positive saturation warning. */
        if (*saturationWarning == 1) satWarningTmp = 1;
    }
    /* Trigger the saturation warning if displayed by any of the frames. */
    *saturationWarning = satWarningTmp;

    return 0;
}

int AGC_set_config(Agc_t *st, AGC_config_t agcConfig)
{
    if ((agcConfig.limiterEnable != 0) && (agcConfig.limiterEnable != 1)) {st->lastError = AGC_BAD_PARAMETER_ERROR;   return -1;}
    st->limiterEnable = agcConfig.limiterEnable;
    st->compressionGaindB = agcConfig.compressionGaindB;
    if ((agcConfig.targetLevelDbfs < 0) || (agcConfig.targetLevelDbfs > 31)) {st->lastError = AGC_BAD_PARAMETER_ERROR;  return -1;}
    st->targetLevelDbfs = agcConfig.targetLevelDbfs;

    if (st->agcMode == kAgcModeFixedDigital) st->compressionGaindB += agcConfig.targetLevelDbfs; /* Adjust for different parameter interpretation in FixedDigital mode */

    /* Update threshold levels for analog adaptation */
    //AGC_UpdateAgcThresholds(st);
    /* Set analog target level in envelope dBOv scale */
    st->analogTarget = DIGITAL_REF_AT_0_COMP_GAIN + (DIFF_REF_TO_ANALOG * st->compressionGaindB + ANALOG_TARGET_LEVEL_2)/ANALOG_TARGET_LEVEL;
    if (st->analogTarget < DIGITAL_REF_AT_0_COMP_GAIN) st->analogTarget = DIGITAL_REF_AT_0_COMP_GAIN;
    if (st->agcMode == kAgcModeFixedDigital) st->analogTarget = st->compressionGaindB; /* Adjust for different parameter interpretation in FixedDigital mode */
    
	/* Since the offset between RMS and ENV is not constant, we should make this into a table, but for now, we'll stick with a constant, tuned for the chosen analog target level. */
    st->targetIdx = ANALOG_TARGET_LEVEL + OFFSET_ENV_TO_RMS;
    /* Analog adaptation limits *//* analogTargetLevel = round((32767*10^(-targetIdx/20))^2*16/2^7) */
    st->analogTargetLevel = RXX_BUFFER_LEN * kTargetLevelTable[st->targetIdx]; /* ex. -20 dBov */
    st->upperLimit = st->startUpperLimit = RXX_BUFFER_LEN * kTargetLevelTable[st->targetIdx - 1];/* -19 dBov */
    st->lowerLimit = st->startLowerLimit = RXX_BUFFER_LEN * kTargetLevelTable[st->targetIdx + 1];/* -21 dBov */
    st->upperPrimaryLimit = RXX_BUFFER_LEN * kTargetLevelTable[st->targetIdx - 2];/* -18 dBov */
    st->lowerPrimaryLimit = RXX_BUFFER_LEN * kTargetLevelTable[st->targetIdx + 2];/* -22 dBov */
    st->upperSecondaryLimit = RXX_BUFFER_LEN * kTargetLevelTable[st->targetIdx - 5];/* -15 dBov */
    st->lowerSecondaryLimit = RXX_BUFFER_LEN * kTargetLevelTable[st->targetIdx + 5];/* -25 dBov */

    /* Recalculate gain table */
    if (AGC_CalculateGainTable(&(st->digitalAgc.gainTable[0]), st->compressionGaindB, st->targetLevelDbfs, st->limiterEnable, st->analogTarget) == -1)
    {
        return -1;
    }
    /* Store the config in a AGC_config_t */
    st->usedConfig.compressionGaindB = agcConfig.compressionGaindB;
    st->usedConfig.limiterEnable = agcConfig.limiterEnable;
    st->usedConfig.targetLevelDbfs = agcConfig.targetLevelDbfs;

    return 0;
}

int AGC_Init(Agc_t *st, int32_t minLevel, int32_t maxLevel, int16_t agcMode, int32_t fs)
{
    int32_t i, max_add;

	//AGC_InitDigital(&st->digitalAgc, agcMode);
    if (agcMode == kAgcModeFixedDigital) st->digitalAgc.capacitorSlow = 0; // start at minimum to find correct gain faster
    else                                 st->digitalAgc.capacitorSlow = 134217728; // start out with 0 dB gain  // (int32_t)(0.125f * 32768.0f * 32768.0f);
    st->digitalAgc.capacitorFast = 0;
    st->digitalAgc.gain = 65536;
    st->digitalAgc.gatePrevious = 0;
    st->digitalAgc.agcMode = agcMode;    
	// initialize VADs
    AGC_InitVad(&st->digitalAgc.vadNearend);
    AGC_InitVad(&st->digitalAgc.vadFarend);	

    /* Analog AGC variables */
    st->envSum = 0;

    /* agcMode  = 0 - Only saturation protection
     *            1 - Analog Automatic Gain Control [-targetLevelDbfs (default -3 dBOv)]
     *            2 - Digital Automatic Gain Control [-targetLevelDbfs (default -3 dBOv)]
     *            3 - Fixed Digital Gain [compressionGaindB (default 8 dB)] */
    st->agcMode = agcMode;
    st->fs = fs;

    /* initialize input VAD */
    AGC_InitVad(&st->vadMic);

    st->scale = 0;
    /* Make minLevel and maxLevel static in AdaptiveDigital */
    if (st->agcMode == kAgcModeAdaptiveDigital)
    {
        minLevel = 0;
        maxLevel = 255;
    }
    // The maximum supplemental volume range is based on a vague idea of how much lower the gain will be than the real analog gain.
    max_add = (maxLevel - minLevel)>>2;

    /* Minimum/maximum volume level that can be set */
    st->minLevel = minLevel;
    st->maxAnalog = maxLevel;
    st->maxLevel = maxLevel + max_add;
    st->maxInit = st->maxLevel;
    st->zeroCtrlMax = st->maxAnalog;

    /* Initialize micVol parameter */
    st->micVol = st->maxAnalog;
    if (st->agcMode == kAgcModeAdaptiveDigital) st->micVol = 127; /* Mid-point of mic level */
    st->micRef = st->micVol;
    st->micGainIdx = 127;
    /* Minimum output volume is 4% higher than the available lowest volume level */
    st->minOutput = st->minLevel + ((st->maxLevel - st->minLevel)*5L>>7);

    st->msTooLow = 0;
    st->msTooHigh = 0;
    st->changeToSlowMode = 0;
    st->firstCall = 0;
    st->msZero = 0;
    st->muteGuardMs = 0;
    st->gainTableIdx = 0;
    st->msecSpeechInnerChange = kMsecSpeechInner;
    st->msecSpeechOuterChange = kMsecSpeechOuter;
    st->activeSpeech = 0;
    st->Rxx16_LPw32Max = 0;
    st->vadThreshold = kNormalVadThreshold;
    st->inActive = 0;

    for (i = 0; i < RXX_BUFFER_LEN; i++) st->Rxx16_vectorw32[i] = 1000; /* -54dBm0 */
    st->Rxx160w32 = 125*RXX_BUFFER_LEN; /* (st->Rxx16_vectorw32[0]>>3) = 125 */

    st->Rxx16pos = 0;
    st->Rxx16_LPw32 = 16284; /* Q(-4) */

    for (i=0; i<5; i++) st->Rxx16w32_array[0][i] = 0;
    for (i=0; i<10; i++) st->env[1][i] = st->env[0][i] = 0;
    st->inQueue = 0;
    for (i=0; i<8; i++) st->filterState[i] = 0;

    // Default config settings.
    st->defaultConfig.limiterEnable = 1;
    st->defaultConfig.targetLevelDbfs = AGC_DEFAULT_TARGET_LEVEL;
    st->defaultConfig.compressionGaindB = AGC_DEFAULT_COMP_GAIN;

    if (AGC_set_config(st, st->defaultConfig) == -1)
    {
        st->lastError = AGC_UNSPECIFIED_ERROR;
        return -1;
    }
    st->Rxx160_LPw32 = st->analogTargetLevel; // Initialize rms value
    st->lowLevelSignal = 0;

    /* Only positive values are allowed that are not too large */
    if ((minLevel >= maxLevel) || (maxLevel & 0xFC000000)) return -1;
	else                                                   return 0;
}
