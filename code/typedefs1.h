#ifndef AGC_TYPEDEFS_H_
#define AGC_TYPEDEFS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

typedef signed short        int16_t;
typedef signed int          int32_t;
typedef unsigned char       uint8_t;
typedef unsigned short      uint16_t;
typedef unsigned int        uint32_t;

#define MUL_ACCUM_1(a, b, c) (c + (b >> 16) * a + (((uint32_t)(0x0000FFFF & b) * a) >> 16))
#define AGC_SCALEDIFF32(A, B, C)    ((C) + ((B)>>16)*(A) + ( ((0x0000FFFF & (B))*(A)) >> 16 ))

// Errors
#define AGC_UNSPECIFIED_ERROR           18000
#define AGC_UNSUPPORTED_FUNCTION_ERROR  18001
#define AGC_UNINITIALIZED_ERROR         18002
#define AGC_NULL_POINTER_ERROR          18003
#define AGC_BAD_PARAMETER_ERROR         18004

#define kAgcModeUnchanged        0 
#define kAgcModeAdaptiveAnalog   1
#define kAgcModeAdaptiveDigital  2
#define kAgcModeFixedDigital     3

/* Analog Automatic Gain Control variables:
 * Constant declarations (inner limits inside which no changes are done)
 * In the beginning the range is narrower to widen as soon as the measure 'Rxx160_LP' is inside it. 
 Currently the starting limits are -22.2+/-1dBm0 and the final limits -22.2+/-2.5dBm0. 
 These levels makes the speech signal go towards -25.4dBm0 (-31.4dBov). 
 Tuned with wbfile-31.4dBov.pcm
 The limits are created by running the AGC with a file having the desired signal level and thereafter plotting Rxx160_LP in the dBm0-domain defined by out=10*log10(in/260537279.7); 
 Set the target level to the average level of our measure Rxx160_LP.
 Remember that the levels are in blocks of 16 in Q(-7). (Example matlab code: round(db2pow(-21.2)*16/2^7) )
 */
#define RXX_BUFFER_LEN  10
#define kMsecSpeechInner  520
#define kMsecSpeechOuter  340
#define kNormalVadThreshold  400
#define kAlphaShortTerm  6       // 1 >> 6 = 0.0156
#define kAlphaLongTerm  10       // 1 >> 10 = 0.000977

// To generate the gaintable, copy&paste the following lines to a Matlab window:
// MaxGain = 6; MinGain = 0; CompRatio = 3; Knee = 1;
// zeros = 0:31; lvl = 2.^(1-zeros);
// A = -10*log10(lvl) * (CompRatio - 1) / CompRatio;
// B = MaxGain - MinGain;
// gains = round(2^16*10.^(0.05 * (MinGain + B * ( log(exp(-Knee*A)+exp(-Knee*B)) - log(1+exp(-Knee*B)) ) / log(1/(1+exp(Knee*B))))));
// fprintf(1, '\t%i, %i, %i, %i,\n', gains);
// % Matlab code for plotting the gain and input/output level characteristic (copy/paste the following 3 lines):
// in = 10*log10(lvl); out = 20*log10(gains/65536);
// subplot(121); plot(in, out); axis([-30, 0, -5, 20]); grid on; xlabel('Input (dB)'); ylabel('Gain (dB)');
// subplot(122); plot(in, in+out); axis([-30, 0, -30, 5]); grid on; xlabel('Input (dB)'); ylabel('Output (dB)');
// zoom on;
// Generator table for y=log2(1+e^x) in Q8.
#define kGenFuncTableSize 128
#define kAvgDecayTime 250    // frames; < 3000
#define kSoftLimiterLeft 1

typedef struct
{
    int32_t downState[8];
    int16_t HPstate;
    int16_t counter;
    int16_t logRatio; // log( P(active) / P(inactive) ) (Q10)
    int16_t meanLongTerm; // Q10
    int32_t varianceLongTerm; // Q8
    int16_t stdLongTerm; // Q10
    int16_t meanShortTerm; // Q10
    int32_t varianceShortTerm; // Q8
    int16_t stdShortTerm; // Q10
} AgcVad_t; // total = 54 bytes

typedef struct
{
    int32_t  capacitorSlow;
    int32_t  capacitorFast;
    int32_t  gain;
    int32_t  gainTable[32];
    int16_t  gatePrevious;
    int16_t  agcMode;
    AgcVad_t vadNearend;
    AgcVad_t vadFarend;
} DigitalAgc_t;

typedef struct
{
    int16_t targetLevelDbfs;   // default 3 (-3 dBOv)
    int16_t compressionGaindB; // default 9 dB
    uint8_t limiterEnable;     // default kAgcTrue (on)
} AGC_config_t;


typedef struct
{
    // Configurable parameters/variables
    uint32_t            fs;                 // Sampling frequency
    int16_t             compressionGaindB;  // Fixed gain level in dB
    int16_t             targetLevelDbfs;    // Target level in -dBfs of envelope (default -3)
    int16_t             agcMode;            // Hard coded mode (adaptAna/adaptDig/fixedDig)
    uint8_t             limiterEnable;      // Enabling limiter (on/off (default off))
    AGC_config_t  defaultConfig;
    AGC_config_t  usedConfig;

    // General variables
    int16_t             lastError;

    // Target level parameters
    // Based on the above: analogTargetLevel = round((32767*10^(-22/20))^2*16/2^7)
    int32_t             analogTargetLevel;  // = RXX_BUFFER_LEN * 846805;       -22 dBfs
    int32_t             startUpperLimit;    // = RXX_BUFFER_LEN * 1066064;      -21 dBfs
    int32_t             startLowerLimit;    // = RXX_BUFFER_LEN * 672641;       -23 dBfs
    int32_t             upperPrimaryLimit;  // = RXX_BUFFER_LEN * 1342095;      -20 dBfs
    int32_t             lowerPrimaryLimit;  // = RXX_BUFFER_LEN * 534298;       -24 dBfs
    int32_t             upperSecondaryLimit;// = RXX_BUFFER_LEN * 2677832;      -17 dBfs
    int32_t             lowerSecondaryLimit;// = RXX_BUFFER_LEN * 267783;       -27 dBfs
    uint16_t            targetIdx;          // Table index for corresponding target level
    int16_t             analogTarget;       // Digital reference level in ENV scale

    // Analog AGC specific variables
    int32_t             filterState[8];     // For downsampling wb to nb
    int32_t             upperLimit;         // Upper limit for mic energy
    int32_t             lowerLimit;         // Lower limit for mic energy
    int32_t             Rxx160w32;          // Average energy for one frame
    int32_t             Rxx16_LPw32;        // Low pass filtered subframe energies
    int32_t             Rxx160_LPw32;       // Low pass filtered frame energies
    int32_t             Rxx16_LPw32Max;     // Keeps track of largest energy subframe
    int32_t             Rxx16_vectorw32[RXX_BUFFER_LEN];// Array with subframe energies
    int32_t             Rxx16w32_array[2][5];// Energy values of microphone signal
    int32_t             env[2][10];         // Envelope values of subframes

    int16_t             Rxx16pos;           // Current position in the Rxx16_vectorw32
    int16_t             envSum;             // Filtered scaled envelope in subframes
    int16_t             vadThreshold;       // Threshold for VAD decision
    int16_t             inActive;           // Inactive time in milliseconds
    int16_t             msTooLow;           // Milliseconds of speech at a too low level
    int16_t             msTooHigh;          // Milliseconds of speech at a too high level
    int16_t             changeToSlowMode;   // Change to slow mode after some time at target
    int16_t             firstCall;          // First call to the process-function
    int16_t             msZero;             // Milliseconds of zero input
    int16_t             msecSpeechOuterChange;// Min ms of speech between volume changes
    int16_t             msecSpeechInnerChange;// Min ms of speech between volume changes
    int16_t             activeSpeech;       // Milliseconds of active speech
    int16_t             muteGuardMs;        // Counter to prevent mute action
    int16_t             inQueue;            // 10 ms batch indicator

    // Microphone level variables
    int32_t             micRef;             // Remember ref. mic level for virtual mic
    uint16_t            gainTableIdx;       // Current position in virtual gain table
    int32_t             micGainIdx;         // Gain index of mic level to increase slowly
    int32_t             micVol;             // Remember volume between frames
    int32_t             maxLevel;           // Max possible vol level, incl dig gain
    int32_t             maxAnalog;          // Maximum possible analog volume level
    int32_t             maxInit;            // Initial value of "max"
    int32_t             minLevel;           // Minimum possible volume level
    int32_t             minOutput;          // Minimum output volume level
    int32_t             zeroCtrlMax;        // Remember max gain => don't amp low input

    int16_t             scale;              // Scale factor for internal volume levels
    // Structs for VAD and digital_agc
    AgcVad_t            vadMic;
    DigitalAgc_t        digitalAgc;
    int16_t             lowLevelSignal;
} Agc_t;

static __inline int AGC_NormW32(int32_t a)
{
  int zeros;

  if (a == 0) return 0;
  else if (a < 0) a = ~a;
  if (!(0xFFFF8000 & a)) zeros = 16;
  else                   zeros = 0;
  if (!(0xFF800000 & (a << zeros))) zeros += 8;
  if (!(0xF8000000 & (a << zeros))) zeros += 4;
  if (!(0xE0000000 & (a << zeros))) zeros += 2;
  if (!(0xC0000000 & (a << zeros))) zeros += 1;

  return zeros;
}

static __inline int AGC_NormU32(uint32_t a) {
  int zeros;

  if (a == 0) return 0;

  if (!(0xFFFF0000 & a)) zeros = 16;
  else                   zeros = 0;
  if (!(0xFF000000 & (a << zeros))) zeros += 8;
  if (!(0xF0000000 & (a << zeros))) zeros += 4;
  if (!(0xC0000000 & (a << zeros))) zeros += 2;
  if (!(0x80000000 & (a << zeros))) zeros += 1;

  return zeros;
}

int AGC_Init(void *agcInst, int32_t minLevel, int32_t maxLevel, int16_t agcMode, int32_t fs);
int AGC_set_config(Agc_t *st, AGC_config_t config);
int AGC_Process(void* agcInst, const int16_t* inNear, const int16_t* inNear_H, int16_t samples, int16_t* out, int16_t* out_H, int32_t inMicLevel, int32_t* outMicLevel, int16_t echo, uint8_t* saturationWarning);

#endif  // AGC_TYPEDEFS_H_
