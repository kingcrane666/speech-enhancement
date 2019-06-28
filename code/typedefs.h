#ifndef TYPEDEFS_H_
#define TYPEDEFS_H_

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

typedef short Word16;
typedef long  Word32;
typedef unsigned short UWord16;
typedef unsigned long  UWord32;

#define FRM_LEN                 160      // ‰»Î–≈∫≈÷°≥§10∫¡√Î
#define STAGES     8
#define ANA_LEN   (1<<STAGES)
#define ANA_LEN2  (1<<(STAGES-1))
#define HALF_ANAL_BLOCKL        ((1<<(STAGES-1))+1) /* Half max analysis block length + 1 */
#define SIMULT                  3
/* Probability threshold for noise state in speech/noise likelihood. */
#define HIST_PAR_EST            1000 /* Histogram size for estimation of parameters */

typedef struct
{
  UWord32      fs;
  Word16       analysisBuffer[ANA_LEN];
  Word16       synthesisBuffer[ANA_LEN];
  UWord16      noiseSupFilter[HALF_ANAL_BLOCKL];
  UWord16      overdrive; /* Q8 */
  UWord16      denoiseBound; /* Q14 */
  const Word16 *factor2Table;
  Word16       noiseEstLogQuantile[SIMULT* HALF_ANAL_BLOCKL];
  Word16       noiseEstDensity[SIMULT* HALF_ANAL_BLOCKL];
  Word16       noiseEstCounter[SIMULT];
  Word16       noiseEstQuantile[HALF_ANAL_BLOCKL];

  Word32       aggrMode;
  Word32       gainMap;

  Word32       maxLrt;
  Word32       minLrt;
  // Log LRT factor with time-smoothing in Q8.
  Word32       logLrtTimeAvgW32[HALF_ANAL_BLOCKL];
  Word32       featureLogLrt;
  Word32       thresholdLogLrt;
  Word16       weightLogLrt;

  UWord32      featureSpecDiff;
  UWord32      thresholdSpecDiff;
  Word16       weightSpecDiff;

  UWord32      featureSpecFlat;
  UWord32      thresholdSpecFlat;
  Word16       weightSpecFlat;

  // Conservative estimate of noise spectrum.
  Word32       avgMagnPause[HALF_ANAL_BLOCKL];
  UWord32      magnEnergy;
  UWord32      sumMagn;
  UWord32      curAvgMagnEnergy;
  UWord32      timeAvgMagnEnergy;
  UWord32      timeAvgMagnEnergyTmp;

  UWord32      whiteNoiseLevel;  // Initial noise estimate.
  // Initial magnitude spectrum estimate.
  UWord32      initMagnEst[HALF_ANAL_BLOCKL];
  // Pink noise parameters:
  Word32       pinkNoiseNumerator;  // Numerator.
  Word32       pinkNoiseExp;  // Power of freq.
  Word32       minNorm;  // Smallest normalization factor.
  Word32       zeroInputSignal;  // Zero input signal flag.

  // Noise spectrum from previous frame.
  UWord32      prevNoiseU32[HALF_ANAL_BLOCKL];
  // Magnitude spectrum from previous frame.
  UWord16      prevMagnU16[HALF_ANAL_BLOCKL];
  // Prior speech/noise probability in Q14.
  Word16       priorNonSpeechProb;

  Word32       blockIndex;  // Frame index counter.
  // Parameter for updating or estimating thresholds/weights for prior model.
  Word32       modelUpdate;
  Word32       cntThresUpdate;

  // Histograms for parameter estimation.
  Word16       histLrt[HIST_PAR_EST];
  Word16       histSpecFlat[HIST_PAR_EST];
  Word16       histSpecDiff[HIST_PAR_EST];

  // Quantities for high band estimate.
  Word16       dataBufHBFX[ANA_LEN];  // Q0

  Word16       qNoise;
  Word32       prevQNoise;
  Word32       prevQMagn;

  Word16       real[ANA_LEN];
  Word16       imag[ANA_LEN];
  Word32       energyIn;
  Word32       scaleEnergyIn;
  Word16       normData;
} NS_STRUCT;

void NS_Init(NS_STRUCT* st, Word32 fs, Word32 mode);
void NS_run(NS_STRUCT* st, short* speechFrame, short* outFrame);

#endif  // TYPEDEFS_H_
