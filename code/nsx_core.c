#include "typedefs.h"

/* Probability threshold for noise state in speech/noise likelihood. */
#define ONE_MINUS_PROB_RANGE_Q8 205 /* 205 ~= Q8(0.8) */

#define END_STARTUP_LONG        200
#define END_STARTUP_SHORT       50
#define FACTOR_Q16              2621440 /* 40 in Q16 */
#define FACTOR_Q7               5120 /* 40 in Q7 */
#define FACTOR_Q7_STARTUP       1024 /* 8 in Q7 */
#define WIDTH_Q8                3 /* 0.01 in Q8 (or 25 ) */

/* PARAMETERS FOR NEW METHOD */
#define DD_PR_SNR_Q11           2007 /* ~= Q11(0.98) DD update of prior SNR */
#define ONE_MINUS_DD_PR_SNR_Q11 41 /* DD update of prior SNR */
#define SPECT_FLAT_TAVG_Q14     4915 /* (0.30) tavg parameter for spectral flatness measure */
#define SPECT_DIFF_TAVG_Q8      77 /* (0.30) tavg parameter for spectral flatness measure */
#define PRIOR_UPDATE_Q14        1638 /* Q14(0.1) Update parameter of prior model */
#define NOISE_UPDATE_Q8         26 /* 26 ~= Q8(0.1) Update parameter for noise */

/* FEATURE EXTRACTION CONFIG  */
/* Bin size of histogram */
#define BIN_SIZE_LRT            10
/* Scale parameters: multiply dominant peaks of the histograms by scale factor to obtain. */
/* Thresholds for prior model */
#define FACTOR_1_LRT_DIFF       6 /* For LRT and spectral difference (5 times bigger) */
/* For spectral_flatness: used when noise is flatter than speech (10 times bigger). */
#define FACTOR_2_FLAT_Q10       922
/* Peak limit for spectral flatness (varies between 0 and 1) */
#define THRES_PEAK_FLAT         24 /* * 2 * BIN_SIZE_FLAT_FX */
/* Limit on spacing of two highest peaks in histogram: spacing determined by bin size. */
#define LIM_PEAK_SPACE_FLAT_DIFF    4 /* * 2 * BIN_SIZE_DIFF_FX */
/* Limit on relevance of second peak */
#define LIM_PEAK_WEIGHT_FLAT_DIFF   2
#define THRES_FLUCT_LRT         10240 /* = 20 * st->modelUpdate; fluctuation limit of LRT feat. */
/* Limit on the max and min values for the feature thresholds */
#define MAX_FLAT_Q10            38912 /* 2 * BIN_SIZE_FLAT_FX */
#define MIN_FLAT_Q10            4096  /* 2 * BIN_SIZE_FLAT_FX */
#define MAX_DIFF                100   /* 2 * BIN_SIZE_DIFF_FX */
#define MIN_DIFF                16    /* 2 * BIN_SIZE_DIFF_FX */
/* Criteria of weight of histogram peak  to accept/reject feature */
#define THRES_WEIGHT_FLAT_DIFF  154 /*(int)(0.3*(st->modelUpdate)) for flatness and difference */

#define STAT_UPDATES            9 /* Update every 512 = 1 << 9 block */
#define ONE_MINUS_GAMMA_PAUSE_Q8    13    /* ~= Q8(0.05) Update for conservative noise estimate */
#define GAMMA_NOISE_TRANS_AND_SPEECH_Q8 3 /* ~= Q8(0.01) Update for transition and noise region */
#define KSTARTBAND    5    // Skip first frequency bins during estimation. (0 <= value < 64)

static const Word16 NS_kLogTable[9] = {0, 177, 355, 532, 710, 887, 1065, 1242, 1420};
static const Word16 NS_kCounterDiv[201] = {32767, 16384, 10923, 8192, 6554, 5461, 4681, 4096, 3641, 3277, 2979, 2731,  2521, 2341, 2185, 2048, 1928, 1820, 1725, 1638, 1560, 1489, 1425, 1365, 1311,  1260, 1214, 1170, 1130, 1092, 1057, 1024, 993, 964, 936, 910, 886, 862, 840,  819, 799, 780, 762, 745, 728, 712, 697, 683, 669, 655, 643, 630, 618, 607,  596, 585, 575, 565, 555, 546, 537, 529, 520, 512, 504, 496, 489, 482, 475,  468, 462, 455, 449, 443, 437, 431, 426, 420, 415, 410, 405, 400, 395, 390,  386, 381, 377, 372, 368, 364, 360, 356, 352, 349, 345, 341, 338, 334, 331,  328, 324, 321, 318, 315, 312, 309, 306, 303, 301, 298, 295, 293, 290, 287,  285, 282, 280, 278, 275, 273, 271, 269, 266, 264, 262, 260, 258, 256, 254,  252, 250, 248, 246, 245, 243, 241, 239, 237, 236, 234, 232, 231, 229, 228,  226, 224, 223, 221, 220, 218, 217, 216, 214, 213, 211, 210, 209, 207, 206,  205, 204, 202, 201, 200, 199, 197, 196, 195, 194, 193, 192, 191, 189, 188,  187, 186, 185, 184, 183, 182, 181, 180, 179, 178, 177, 176, 175, 174, 173,  172, 172, 171, 170, 169, 168, 167, 166, 165, 165, 164, 163};
static const Word16 NS_kLogTableFrac[256] = {0,   1,   3,   4,   6,   7,   9,  10,  11,  13,  14,  16,  17,  18,  20,  21,  22,  24,  25,  26,  28,  29,  30,  32,  33,  34,  36,  37,  38,  40,  41,  42,  44,  45,  46,  47,  49,  50,  51,  52,  54,  55,  56,  57,  59,  60,  61,  62,  63,  65,  66,  67,  68,  69,  71,  72,  73,  74,  75,  77,  78,  79,  80,  81,  82,  84,  85,  86,  87,  88,  89,  90,  92,  93,  94,  95,  96,  97,  98,  99,  100, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 116,  117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131,  132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146,  147, 148, 149, 150, 151, 152, 153, 154, 155, 155, 156, 157, 158, 159, 160,  161, 162, 163, 164, 165, 166, 167, 168, 169, 169, 170, 171, 172, 173, 174,  175, 176, 177, 178, 178, 179, 180, 181, 182, 183, 184, 185, 185, 186, 187,  188, 189, 190, 191, 192, 192, 193, 194, 195, 196, 197, 198, 198, 199, 200,  201, 202, 203, 203, 204, 205, 206, 207, 208, 208, 209, 210, 211, 212, 212,  213, 214, 215, 216, 216, 217, 218, 219, 220, 220, 221, 222, 223, 224, 224,  225, 226, 227, 228, 228, 229, 230, 231, 231, 232, 233, 234, 234, 235, 236,  237, 238, 238, 239, 240, 241, 241, 242, 243, 244, 244, 245, 246, 247, 247,  248, 249, 249, 250, 251, 252, 252, 253, 254, 255, 255};

static const Word16 kIndicatorTable[17] = {0, 2017, 3809, 5227, 6258, 6963, 7424, 7718,  7901, 8014, 8084, 8126, 8152, 8168, 8177, 8183, 8187};
// hybrib Hanning & flat window
static const Word16 kBlocks160w256x[256] = {0,   268,   536,   804,  1072,  1339,  1606,  1872,  2139,  2404,  2669,  2933,  3196,  3459,  3720,  3981,  4240,  4499,  4756,  5012,  5266,  5520,  5771,  6021,  6270,  6517,  6762,  7005,  7246,  7486,  7723,  7959,  8192,  8423,  8652,  8878,  9102,  9324,  9543,  9760,  9974, 10185, 10394, 10600, 10803, 11003, 11200, 11394,  11585, 11773, 11958, 12140, 12318, 12493, 12665, 12833,  12998, 13160, 13318, 13472, 13623, 13770, 13913, 14053,  14189, 14321, 14449, 14574, 14694, 14811, 14924, 15032,  15137, 15237, 15334, 15426, 15515, 15599, 15679, 15754,  15826, 15893, 15956, 16015, 16069, 16119, 16165, 16207,  16244, 16277, 16305, 16329, 16349, 16364, 16375, 16382,  16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384,  16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384,  16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384,  16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384,  16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384,  16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384,  16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384,  16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384,  16384, 16382, 16375, 16364, 16349, 16329, 16305, 16277,  16244, 16207, 16165, 16119, 16069, 16015, 15956, 15893,  15826, 15754, 15679, 15599, 15515, 15426, 15334, 15237,  15137, 15032, 14924, 14811, 14694, 14574, 14449, 14321,  14189, 14053, 13913, 13770, 13623, 13472, 13318, 13160,  12998, 12833, 12665, 12493, 12318, 12140, 11958, 11773,  11585, 11394, 11200, 11003, 10803, 10600, 10394, 10185,  9974,  9760,  9543,  9324,  9102,  8878,  8652,  8423,  8192,  7959,  7723,  7486,  7246,  7005,  6762,  6517,  6270,  6021,  5771,  5520,  5266,  5012,  4756,  4499,  4240,  3981,  3720,  3459,  3196,  2933,  2669,  2404,  2139,  1872,  1606,  1339,  1072,   804,   536,   268};
static const Word16 kFactor1Table[257] = {8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8233, 8274, 8315, 8355, 8396, 8436, 8475, 8515, 8554, 8592, 8631, 8669,  8707, 8745, 8783, 8820, 8857, 8894, 8931, 8967, 9003, 9039, 9075, 9111, 9146, 9181,  9216, 9251, 9286, 9320, 9354, 9388, 9422, 9456, 9489, 9523, 9556, 9589, 9622, 9655,  9687, 9719, 9752, 9784, 9816, 9848, 9879, 9911, 9942, 9973, 10004, 10035, 10066,  10097, 10128, 10158, 10188, 10218, 10249, 10279, 10308, 10338, 10368, 10397, 10426,  10456, 10485, 10514, 10543, 10572, 10600, 10629, 10657, 10686, 10714, 10742, 10770,  10798, 10826, 10854, 10882, 10847, 10810, 10774, 10737, 10701, 10666, 10631, 10596,  10562, 10527, 10494, 10460, 10427, 10394, 10362, 10329, 10297, 10266, 10235, 10203,  10173, 10142, 10112, 10082, 10052, 10023, 9994, 9965, 9936, 9908, 9879, 9851, 9824,  9796, 9769, 9742, 9715, 9689, 9662, 9636, 9610, 9584, 9559, 9534, 9508, 9484, 9459,  9434, 9410, 9386, 9362, 9338, 9314, 9291, 9268, 9245, 9222, 9199, 9176, 9154, 9132,  9110, 9088, 9066, 9044, 9023, 9002, 8980, 8959, 8939, 8918, 8897, 8877, 8857, 8836,  8816, 8796, 8777, 8757, 8738, 8718, 8699, 8680, 8661, 8642, 8623, 8605, 8586, 8568,  8550, 8532, 8514, 8496, 8478, 8460, 8443, 8425, 8408, 8391, 8373, 8356, 8339, 8323,  8306, 8289, 8273, 8256, 8240, 8224, 8208, 8192};
static const Word16 kFactor2Aggressiveness1[257] = {7577, 7577, 7577, 7577, 7577, 7577,  7577, 7577, 7577, 7577, 7577, 7577, 7577, 7577, 7577, 7577, 7577, 7596, 7614, 7632,  7650, 7667, 7683, 7699, 7715, 7731, 7746, 7761, 7775, 7790, 7804, 7818, 7832, 7845,  7858, 7871, 7884, 7897, 7910, 7922, 7934, 7946, 7958, 7970, 7982, 7993, 8004, 8016,  8027, 8038, 8049, 8060, 8070, 8081, 8091, 8102, 8112, 8122, 8132, 8143, 8152, 8162,  8172, 8182, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192};
static const Word16 kFactor2Aggressiveness2[257] = {7270, 7270, 7270, 7270, 7270, 7306,  7339, 7369, 7397, 7424, 7448, 7472, 7495, 7517, 7537, 7558, 7577, 7596, 7614, 7632,  7650, 7667, 7683, 7699, 7715, 7731, 7746, 7761, 7775, 7790, 7804, 7818, 7832, 7845,  7858, 7871, 7884, 7897, 7910, 7922, 7934, 7946, 7958, 7970, 7982, 7993, 8004, 8016,  8027, 8038, 8049, 8060, 8070, 8081, 8091, 8102, 8112, 8122, 8132, 8143, 8152, 8162,  8172, 8182, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192};
static const Word16 kFactor2Aggressiveness3[257] = {7184, 7184, 7184, 7229, 7270, 7306,  7339, 7369, 7397, 7424, 7448, 7472, 7495, 7517, 7537, 7558, 7577, 7596, 7614, 7632,  7650, 7667, 7683, 7699, 7715, 7731, 7746, 7761, 7775, 7790, 7804, 7818, 7832, 7845,  7858, 7871, 7884, 7897, 7910, 7922, 7934, 7946, 7958, 7970, 7982, 7993, 8004, 8016,  8027, 8038, 8049, 8060, 8070, 8081, 8091, 8102, 8112, 8122, 8132, 8143, 8152, 8162,  8172, 8182, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192,  8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192, 8192};
// sum of log2(i) from table index to st->anaLen2 in Q5// Note that the first table value is invalid, since log2(0) = -infinity
static const Word16 kSumLogIndex[66] = {0,  22917,  22917,  22885,  22834,  22770,  22696,  22613,  22524,  22428,  22326,  22220,  22109,  21994,  21876,  21754,  21629,  21501,  21370,  21237,  21101,  20963,  20822,  20679,  20535,  20388,  20239,  20089,  19937,  19783,  19628,  19470,  19312,  19152,  18991,  18828,  18664,  18498,  18331,  18164,  17994,  17824,  17653,  17480,  17306,  17132,  16956,  16779,  16602,  16423,  16243,  16063,  15881,  15699,  15515,  15331,  15146,  14960,  14774,  14586,  14398,  14209,  14019,  13829,  13637,  13445};
// sum of log2(i)^2 from table index to st->anaLen2 in Q2// Note that the first table value is invalid, since log2(0) = -infinity
static const Word16 kSumSquareLogIndex[66] = {0,  16959,  16959,  16955,  16945,  16929,  16908,  16881,  16850,  16814,  16773,  16729,  16681,  16630,  16575,  16517,  16456,  16392,  16325,  16256,  16184,  16109,  16032,  15952,  15870,  15786,  15700,  15612,  15521,  15429,  15334,  15238,  15140,  15040,  14938,  14834,  14729,  14622,  14514,  14404,  14292,  14179,  14064,  13947,  13830,  13710,  13590,  13468,  13344,  13220,  13094,  12966,  12837,  12707,  12576,  12444,  12310,  12175,  12039,  11902,  11763,  11624,  11483,  11341,  11198,  11054};
// log2(table index) in Q12// Note that the first table value is invalid, since log2(0) = -infinity
static const Word16 kLogIndex[129] = {0,  0,  4096,   6492,   8192,   9511,  10588,  11499,  12288,  12984,  13607,  14170,  14684,  15157,  15595,  16003,  16384,  16742,  17080,  17400,  17703,  17991,  18266,  18529,  18780,  19021,  19253,  19476,  19691,  19898,  20099,  20292,  20480,  20662,  20838,  21010,  21176,  21338,  21496,  21649,  21799,  21945,  22087,  22226,  22362,  22495,  22625,  22752,  22876,  22998,  23117,  23234,  23349,  23462,  23572,  23680,  23787,  23892,  23994,  24095,  24195,  24292,  24388,  24483,  24576,  24668,  24758,  24847,  24934,  25021,  25106,  25189,  25272,  25354,  25434,  25513,  25592,  25669,  25745,  25820,  25895,  25968,  26041,  26112,  26183,  26253,  26322,  26390,  26458,  26525,  26591,  26656,  26721,  26784,  26848,  26910,  26972,  27033,  27094,  27154,  27213,  27272,  27330,  27388,  27445,  27502,  27558,  27613,  27668,  27722,  27776,  27830,  27883,  27935,  27988,  28039,  28090,  28141,  28191,  28241,  28291,  28340,  28388,  28437,  28484,  28532,  28579,  28626,  28672};
// determinant of estimation matrix in Q0 corresponding to the log2 tables above// Note that the first table value is invalid, since log2(0) = -infinity
static const Word16 kDeterminantEstMatrix[66] = {0,  29814,  25574,  22640,  20351,  18469,  16873,  15491,  14277,  13199,  12233,  11362,  10571,   9851,   9192,   8587,  8030,   7515,   7038,   6596,   6186,   5804,   5448,   5115,  4805,   4514,   4242,   3988,   3749,   3524,   3314,   3116,  2930,   2755,   2590,   2435,   2289,   2152,   2022,   1900,  1785,   1677,   1575,   1478,   1388,   1302,   1221,   1145,  1073,   1005,    942,    881,    825,    771,    721,    674,  629,    587,    547,    510,    475,    442,    411,    382,  355,    330};
static const Word16 rot_tab[ANA_LEN]={32767, 0, 32757, 804, 32727, 1607, 32678, 2410, 32609, 3211, 32520, 4011, 32412, 4807, 32284, 5601, 32137, 6392, 31970, 7179, 31785, 7961, 31580, 8739, 31356, 9511, 31113, 10278, 30851, 11038, 30571, 11792, 30272, 12539, 29955, 13278, 29621, 14009, 29268, 14732, 28897, 15446, 28510, 16150, 28105, 16845, 27683, 17530, 27244, 18204, 26789, 18867, 26318, 19519, 25831, 20159, 25329, 20787, 24811, 21402, 24278, 22004, 23731, 22594, 23169, 23169, 22594, 23731, 22004, 24278, 21402, 24811, 20787, 25329, 20159, 25831, 19519, 26318, 18867, 26789, 18204, 27244, 17530, 27683, 16845, 28105, 16150, 28510, 15446, 28897, 14732, 29268, 14009, 29621, 13278, 29955, 12539, 30272, 11792, 30571, 11038, 30851, 10278, 31113, 9511, 31356, 8739, 31580, 7961, 31785, 7179, 31970, 6392, 32137, 5601, 32284, 4807, 32412, 4011, 32520, 3211, 32609, 2410, 32678, 1607, 32727, 804, 32757, 0, 32767, -804, 32757, -1607, 32727, -2410, 32678, -3211, 32609, -4011, 32520, -4807, 32412, -5601, 32284, -6392, 32137, -7179, 31970, -7961, 31785, -8739, 31580, -9511, 31356, -10278, 31113, -11038, 30851, -11792, 30571, -12539, 30272, -13278, 29955, -14009, 29621, -14732, 29268, -15446, 28897, -16150, 28510, -16845, 28105, -17530, 27683, -18204, 27244, -18867, 26789, -19519, 26318, -20159, 25831, -20787, 25329, -21402, 24811, -22004, 24278, -22594, 23731, -23169, 23169, -23731, 22594, -24278, 22004, -24811, 21402, -25329, 20787, -25831, 20159, -26318, 19519, -26789, 18867, -27244, 18204, -27683, 17530, -28105, 16845, -28510, 16150, -28897, 15446, -29268, 14732, -29621, 14009, -29955, 13278, -30272, 12539, -30571, 11792, -30851, 11038, -31113, 10278, -31356, 9511, -31580, 8739, -31785, 7961, -31970, 7179, -32137, 6392, -32284, 5601, -32412, 4807, -32520, 4011, -32609, 3211, -32678, 2410, -32727, 1607, -32757, 804};
// Tables for data buffer indexes that are bit reversed and thus need to be swapped. 
static const Word16 bit_rev[1<<STAGES]={0, 128, 64, 192, 32, 160, 96, 224, 16, 144, 80, 208, 48, 176, 112, 240, 8, 136, 72, 200, 40, 168, 104, 232, 24, 152, 88, 216, 56, 184, 120, 248, 4, 132, 68, 196, 36, 164, 100, 228, 20, 148, 84, 212, 52, 180, 116, 244, 12, 140, 76, 204, 44, 172, 108, 236, 28, 156, 92, 220, 60, 188, 124, 252, 2, 130, 66, 194, 34, 162, 98, 226, 18, 146, 82, 210, 50, 178, 114, 242, 10, 138, 74, 202, 42, 170, 106, 234, 26, 154, 90, 218, 58, 186, 122, 250, 6, 134, 70, 198, 38, 166, 102, 230, 22, 150, 86, 214, 54, 182, 118, 246, 14, 142, 78, 206, 46, 174, 110, 238, 30, 158, 94, 222, 62, 190, 126, 254, 1, 129, 65, 193, 33, 161, 97, 225, 17, 145, 81, 209, 49, 177, 113, 241, 9, 137, 73, 201, 41, 169, 105, 233, 25, 153, 89, 217, 57, 185, 121, 249, 5, 133, 69, 197, 37, 165, 101, 229, 21, 149, 85, 213, 53, 181, 117, 245, 13, 141, 77, 205, 45, 173, 109, 237, 29, 157, 93, 221, 61, 189, 125, 253, 3, 131, 67, 195, 35, 163, 99, 227, 19, 147, 83, 211, 51, 179, 115, 243, 11, 139, 75, 203, 43, 171, 107, 235, 27, 155, 91, 219, 59, 187, 123, 251, 7, 135, 71, 199, 39, 167, 103, 231, 23, 151, 87, 215, 55, 183, 119, 247, 15, 143, 79, 207, 47, 175, 111, 239, 31, 159, 95, 223, 63, 191, 127, 255};

static __inline int _NormW32(Word32 a) 
{
  int zeros;

  if (a == 0) return 0;
  else if (a < 0)  a = ~a;

  if (!(0xFFFF8000 & a)) zeros = 16;
  else                   zeros = 0;

  if (!(0xFF800000 & (a << zeros))) zeros += 8;
  if (!(0xF8000000 & (a << zeros))) zeros += 4;
  if (!(0xE0000000 & (a << zeros))) zeros += 2;
  if (!(0xC0000000 & (a << zeros))) zeros += 1;

  return zeros;
}

static __inline int _NormU32(UWord32 a) 
{
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

static __inline int _NormW16(Word16 a)
{
  int zeros;

  if (a == 0) return 0;
  else if (a < 0) a = ~a;

  if (!(0xFF80 & a)) zeros = 8;
  else               zeros = 0;

  if (!(0xF800 & (a << zeros))) zeros += 4;
  if (!(0xE000 & (a << zeros))) zeros += 2;
  if (!(0xC000 & (a << zeros))) zeros += 1;

  return zeros;
}

static Word32 _norm(Word32 x) {Word32 i; for(i=1;i<32;i++) if((x^(x<<i))<0) break;  return(i-1);}

// Update the noise estimation information.
static void UpdateNoiseEstimate(NS_STRUCT* st, int offset)
{
    Word32 tmp32no1, tmp32no2;
    Word16 tmp16, i;

    for (tmp16 = -32768, i = 0; i < HALF_ANAL_BLOCKL; i++) tmp16 = max(tmp16, st->noiseEstLogQuantile[i+offset]);
    // Guarantee a Q-domain as high as possible and still fit in int16
    st->qNoise = 14 - (Word16) ((11819L*tmp16+1048576L)>>21);

    for (i = 0; i < HALF_ANAL_BLOCKL; i++)
	{
        // st->quantile[i]=exp(st->lquantile[offset+i]); // in Q21
        tmp32no2 = 11819L*st->noiseEstLogQuantile[offset + i];
        tmp32no1 = 0x00200000 | (tmp32no2 & 0x001FFFFF); // 2^21 + frac
        tmp16 = (Word16)(tmp32no2>>21) + st->qNoise - 21; // shift 21 to get result in Q0//shift to get result in Q(qNoise)
        if (tmp16>=0) tmp32no1 = tmp32no1<<tmp16;
	    else          tmp32no1 = tmp32no1>>(-tmp16);
	    st->noiseEstQuantile[i] = (Word16)max(min(tmp32no1, 32767), -32768);
	}
}

// Noise Estimation
static void NoiseEstimationC(NS_STRUCT* st, UWord16* magn, UWord32* noise, Word16* q_noise)
{
    Word16 lmagn[HALF_ANAL_BLOCKL], counter, countDiv, countProd, delta, zeros, frac, log2, tabind, logval, tmp16;
    int i, s, offset;

    tabind = 8 - st->normData;
    assert(labs(tabind) < 9);

    if (tabind < 0) logval = -NS_kLogTable[-tabind];
    else            logval = NS_kLogTable[tabind];
 
    // lmagn(i)=log(magn(i))=log(2)*log2(magn(i)) // magn is in Q(-STAGES), and the real lmagn values are: real_lmagn(i)=log(magn(i)*2^STAGES)=log(magn(i))+log(2^STAGES) // lmagn in Q8
    for (i = 0; i < HALF_ANAL_BLOCKL; i++)
	{
    	zeros = (Word16)_norm((Word32)magn[i])-16;
        frac = (magn[i]<<zeros)&0x3FFF;
        // log2(magn(i))
        log2 = (Word16)(((14-zeros) << 8) + NS_kLogTableFrac[frac>>6]);
        // log2(magn(i))*log(2) + log(2^STAGES)
        lmagn[i] = (Word16)(log2*22713L>>15);
		if (magn[i]==0) lmagn[i]=0;            
		lmagn[i] += logval;//0;
	}

    // loop over simultaneous estimates
    for (s = 0; s < SIMULT; s++)
    {
        offset = s * HALF_ANAL_BLOCKL;        
        // Get counter values from state
        counter = st->noiseEstCounter[s];
        assert(st->noiseEstCounter[s] < 201);
        countDiv = NS_kCounterDiv[counter];
        countProd = (Word16)((Word32)counter*countDiv);
        
        // quant_est(...)
        for (i = 0; i < HALF_ANAL_BLOCKL; i++)
	    {
            // compute delta
            if (st->noiseEstDensity[offset+i]>512) delta = 160<<(_norm(st->noiseEstDensity[offset+i])-16); // Get the value for delta by shifting intead of dividing.
	        else                                   delta = (st->blockIndex < END_STARTUP_LONG) ? FACTOR_Q7_STARTUP : FACTOR_Q7; // Smaller step size during startup. This prevents from using unrealistic values causing overflow.
            
            // update log quantile estimate
            tmp16 = (Word16)((Word32)delta*countDiv>>14);
            if (lmagn[i] > st->noiseEstLogQuantile[offset + i]) st->noiseEstLogQuantile[offset + i] += (tmp16+2)>>2; // +=QUANTILE*delta/(st->counter[s]+1) QUANTILE=0.25, =1 in Q2   // CounterDiv=1/(st->counter[s]+1) in Q15
	        else  st->noiseEstLogQuantile[offset + i] = max(st->noiseEstLogQuantile[offset + i] - (((tmp16+1)>>1)*3L>>1), logval);  // *(1-QUANTILE), in Q2 QUANTILE=0.25, 1-0.25=0.75=3 in Q2 // This is the smallest fixed point representation we can have, hence we limit the output.
            
            // update density estimate
            if (labs(lmagn[i] - st->noiseEstLogQuantile[offset + i]) < WIDTH_Q8) st->noiseEstDensity[offset + i] = (Word16)(((Word32)st->noiseEstDensity[offset+i]*countProd+16384)>>15) + (Word16)((21845L*countDiv+16384)>>15);

        }  // end loop over magnitude spectrum
        
        if (counter >= END_STARTUP_LONG)
	    {
            st->noiseEstCounter[s] = 0;
            if (st->blockIndex >= END_STARTUP_LONG)
	        {
                UpdateNoiseEstimate(st, offset);
            }
        }
        st->noiseEstCounter[s]++;    
    }  // end loop over simultaneous estimates
    
    // Sequentially update the noise during startup
    if (st->blockIndex < END_STARTUP_LONG) 
    {
        UpdateNoiseEstimate(st, offset);
    }
    
    for (i = 0; i < HALF_ANAL_BLOCKL; i++) noise[i] = (UWord32)(st->noiseEstQuantile[i]); // Q(qNoise)
    *q_noise = (Word16)st->qNoise;
}

void NS_CalcParametricNoiseEstimate(NS_STRUCT* st, Word16 pink_noise_exp_avg, Word32 pink_noise_num_avg, int freq_index, UWord32* noise_estimate, UWord32* noise_estimate_avg)
{ 
    Word32 tmp32no1, tmp32no2;
    Word16 int_part, frac_part;
    
    // Use pink noise estimate  // noise_estimate = 2^(pinkNoiseNumerator + pinkNoiseExp * log2(j))
    tmp32no1 = pink_noise_num_avg - ((Word32)pink_noise_exp_avg*kLogIndex[freq_index]>>15) + ((Word32)(st->minNorm - 8)<<11); // Q11
    
    // Calculate output: 2^tmp32no1 Output in Q(minNorm-STAGES)
    if (tmp32no1 > 0)
    {
        int_part = (Word16)(tmp32no1>>11);
        frac_part = (Word16)(tmp32no1 & 0x000007ff); // Q11
        // Piecewise linear approximation of 'b' in 2^(int_part+frac_part) = 2^int_part * (1 + b)  // 'b' is given in Q11 and below stored in frac_part.
        if (frac_part>>10) tmp32no2 = 2048 - ((2048 - frac_part)*311L>>8); // Upper fractional part
	    else               tmp32no2 = frac_part*201L>>8;  // Lower fractional part
        // Shift fractional part to Q(minNorm-STAGES)
        if(int_part>=11) tmp32no2 = tmp32no2<<(int_part - 11);
        else             tmp32no2 = tmp32no2>>(11 - int_part);

        *noise_estimate_avg = (1UL<<int_part) + tmp32no2;
        // Scale up to initMagnEst, which is not block averaged
        *noise_estimate = (*noise_estimate_avg) * (st->blockIndex + 1);
    }
}

// Extract thresholds for feature parameters
// histograms are computed over some window_size (given by window_pars)
// thresholds and weights are extracted every window
// flag 0 means update histogram only, flag 1 means compute the thresholds/weights
// threshold and weights are returned in: st->priorModelPars
static void NS_FeatureParameterExtraction(NS_STRUCT* st, int flag)
{
    UWord32 tmpU32;
    UWord32 histIndex;
    UWord32 posPeak1SpecFlatFX, posPeak2SpecFlatFX, posPeak1SpecDiffFX, posPeak2SpecDiffFX;    
    Word32 tmp32;
    Word32 fluctLrtFX, thresFluctLrtFX, avgHistLrtFX, avgSquareHistLrtFX, avgHistLrtComplFX;
    Word16 j, numHistLrt;    
    int i;
    int useFeatureSpecFlat, useFeatureSpecDiff, featureSum;
    int maxPeak1, maxPeak2, weightPeak1SpecFlat, weightPeak2SpecFlat, weightPeak1SpecDiff, weightPeak2SpecDiff;
    
    //update histograms
    if (!flag)
    {
        // LRT
        // Type casting to UWord32 is safe since negative values will not be wrapped to larger
        // values than HIST_PAR_EST
        histIndex = (UWord32)(st->featureLogLrt);
        if (histIndex < HIST_PAR_EST) st->histLrt[histIndex]++;
        // Spectral flatness // (st->featureSpecFlat*20)>>10 = (st->featureSpecFlat*5)>>8
        histIndex = st->featureSpecFlat * 5UL>>8;
        if (histIndex < HIST_PAR_EST) st->histSpecFlat[histIndex]++;
        
        // Spectral difference
        histIndex = HIST_PAR_EST;
	    // Guard against division by zero  // If timeAvgMagnEnergy == 0 we have no normalizing statistics and therefore can't update the histogram
        if (st->timeAvgMagnEnergy > 0) histIndex = (st->featureSpecDiff*5L>>8)/st->timeAvgMagnEnergy;
        if (histIndex < HIST_PAR_EST) st->histSpecDiff[histIndex]++;
    }
    
    // extract parameters for speech/noise probability
    if (flag)
    {
        useFeatureSpecDiff = 1;
        //for LRT feature:
        // compute the average over st->featureExtractionParams.rangeAvgHistLrt
        for (avgSquareHistLrtFX=avgHistLrtFX=numHistLrt=i=0; i<BIN_SIZE_LRT; i++)
	    {
            j = 2*i+1;
            tmp32 = (Word32)st->histLrt[i]*j;
            avgHistLrtFX += tmp32;
            numHistLrt += st->histLrt[i];
            avgSquareHistLrtFX += tmp32*j;
        }
        avgHistLrtComplFX = avgHistLrtFX;
        for (; i<HIST_PAR_EST; i++)
	    {
            j = 2*i+1;
            tmp32 = (Word32)st->histLrt[i]*j;
            avgHistLrtComplFX += tmp32;
            avgSquareHistLrtFX += tmp32*j;
        }
        fluctLrtFX = (Word32)avgSquareHistLrtFX*numHistLrt - (Word32)avgHistLrtFX*avgHistLrtComplFX;
        thresFluctLrtFX = THRES_FLUCT_LRT * numHistLrt;
        // get threshold for LRT feature:
        tmpU32 = (FACTOR_1_LRT_DIFF * (UWord32)avgHistLrtFX);
        if ((fluctLrtFX < thresFluctLrtFX)||(numHistLrt==0)||(tmpU32>100UL*numHistLrt)) st->thresholdLogLrt = st->maxLrt; //very low fluctuation, so likely noise
        else                                                                            st->thresholdLogLrt = min(max((Word32)((tmpU32 << (9 + 8)) /numHistLrt/25), st->minLrt), st->maxLrt); // check if value is within min/max range      
        if (fluctLrtFX < thresFluctLrtFX) useFeatureSpecDiff = 0; // Do not use difference feature if fluctuation of LRT feature is very low: most likely just noise state
        
        // for spectral flatness and spectral difference: compute the main peaks of histogram
        maxPeak2 = maxPeak1 = 0;
        posPeak2SpecFlatFX = posPeak1SpecFlatFX = 0;
        weightPeak2SpecFlat = weightPeak1SpecFlat = 0;
        
        // peaks for flatness
        for (i = 0; i < HIST_PAR_EST; i++)
	    {
            if (st->histSpecFlat[i] > maxPeak1) // Found new "first" peak
	        {
                maxPeak2 = maxPeak1;
                weightPeak2SpecFlat = weightPeak1SpecFlat;
                posPeak2SpecFlatFX = posPeak1SpecFlatFX;
                
                maxPeak1 = st->histSpecFlat[i];
                weightPeak1SpecFlat = st->histSpecFlat[i];
                posPeak1SpecFlatFX = 2*i+1;
            }
	        else if (st->histSpecFlat[i] > maxPeak2) // Found new "second" peak
	        {
                maxPeak2 = st->histSpecFlat[i];
                weightPeak2SpecFlat = st->histSpecFlat[i];
                posPeak2SpecFlatFX = 2*i+1;
            }
        }
        
        // for spectral flatness feature
        useFeatureSpecFlat = 1;
        // merge the two peaks if they are close
        if ((posPeak1SpecFlatFX - posPeak2SpecFlatFX < LIM_PEAK_SPACE_FLAT_DIFF) && (weightPeak2SpecFlat * LIM_PEAK_WEIGHT_FLAT_DIFF > weightPeak1SpecFlat))
	    {
            weightPeak1SpecFlat += weightPeak2SpecFlat;
            posPeak1SpecFlatFX = (posPeak1SpecFlatFX + posPeak2SpecFlatFX) >> 1;
        }
        //reject if weight of peaks is not large enough, or peak value too small
        if (weightPeak1SpecFlat < THRES_WEIGHT_FLAT_DIFF || posPeak1SpecFlatFX < THRES_PEAK_FLAT) useFeatureSpecFlat = 0;
	    else  st->thresholdSpecFlat = min(max(FACTOR_2_FLAT_Q10 * posPeak1SpecFlatFX, MIN_FLAT_Q10), MAX_FLAT_Q10); // if selected, get the threshold  // compute the threshold and check if value is within min/max range  //Q10

        // done with flatness feature        
        if (useFeatureSpecDiff)
	    {
          //compute two peaks for spectral difference
          maxPeak2 = maxPeak1 = 0;
          posPeak2SpecDiffFX = posPeak1SpecDiffFX = 0;
          weightPeak2SpecDiff = weightPeak1SpecDiff = 0;
          // peaks for spectral difference
          for (i = 0; i < HIST_PAR_EST; i++)
	      {
              if (st->histSpecDiff[i] > maxPeak1)
			  {
                  // Found new "first" peak
                  maxPeak2 = maxPeak1;
                  weightPeak2SpecDiff = weightPeak1SpecDiff;
                  posPeak2SpecDiffFX = posPeak1SpecDiffFX;
        
                  maxPeak1 = st->histSpecDiff[i];
                  weightPeak1SpecDiff = st->histSpecDiff[i];
                  posPeak1SpecDiffFX = 2*i+1;
			  }
	    	  else if (st->histSpecDiff[i] > maxPeak2)
			  {
                  // Found new "second" peak
                  maxPeak2 = st->histSpecDiff[i];
                  weightPeak2SpecDiff = st->histSpecDiff[i];
                  posPeak2SpecDiffFX = 2*i+1;
			  }
          }
        
          // merge the two peaks if they are close
          if ((posPeak1SpecDiffFX - posPeak2SpecDiffFX < LIM_PEAK_SPACE_FLAT_DIFF) && (weightPeak2SpecDiff * LIM_PEAK_WEIGHT_FLAT_DIFF > weightPeak1SpecDiff)) 
	      {
              weightPeak1SpecDiff += weightPeak2SpecDiff;
              posPeak1SpecDiffFX = (posPeak1SpecDiffFX + posPeak2SpecDiffFX) >> 1;
          }
          // get the threshold value and check if value is within min/max range //5x bigger
	      st->thresholdSpecDiff = min(max(FACTOR_1_LRT_DIFF * posPeak1SpecDiffFX, MIN_DIFF), MAX_DIFF);
          //reject if weight of peaks is not large enough
          if (weightPeak1SpecDiff < THRES_WEIGHT_FLAT_DIFF) useFeatureSpecDiff = 0;
        } // done with spectral difference feature
        
        // select the weights between the features // st->priorModelPars[4] is weight for LRT: always selected
        featureSum = 6 / (1 + useFeatureSpecFlat + useFeatureSpecDiff);
        st->weightLogLrt = featureSum;
        st->weightSpecFlat = useFeatureSpecFlat * featureSum;
        st->weightSpecDiff = useFeatureSpecDiff * featureSum;
        
        // set histograms to zero for next update
		for(i=0; i<HIST_PAR_EST; i++) st->histLrt[i]=st->histSpecDiff[i]=st->histSpecFlat[i]=0;
    }  // end of flag == 1
}

// Compute spectral flatness on input spectrum. //magn is the magnitude spectrum. spectral flatness is returned in st->featureSpecFlat
void NS_ComputeSpectralFlatness(NS_STRUCT* st, UWord16* magn)
{
    UWord32 tmpU32, avgSpectralFlatnessNum, avgSpectralFlatnessDen;
    Word32 tmp32, currentSpectralFlatness, logCurSpectralFlatness;
    Word16 zeros, frac, intPart;
    int i;
    
    // for flatness
    avgSpectralFlatnessNum = 0;
    avgSpectralFlatnessDen = st->sumMagn - magn[0]; // Q(normData-STAGES)
    
    // compute log of ratio of the geometric to arithmetic mean: check for log(0) case
    // flatness = exp( sum(log(magn[i]))/N - log(sum(magn[i])/N) ) = exp( sum(log(magn[i]))/N ) * N / sum(magn[i]) = 2^( sum(log2(magn[i]))/N - (log2(sum(magn[i])) - log2(N)) )
    for (i=1; i<HALF_ANAL_BLOCKL; i++)
    {
        // First bin is excluded from spectrum measures. Number of bins is now a power of 2
        if (magn[i])
	    {
	        zeros = (Word16)_norm((Word32)magn[i])-16;
            frac = (Word16)((magn[i]<<zeros)&0x3FFF);
            // log2(magn(i))
            tmpU32 = (UWord32)(((14 - zeros) << 8) + NS_kLogTableFrac[frac>>6]); // Q8
            avgSpectralFlatnessNum += tmpU32; // Q8
        } 
	    else
	    {
            //if at least one frequency component is zero, treat separately
            tmpU32 = (UWord32)st->featureSpecFlat*SPECT_FLAT_TAVG_Q14; // Q24
            st->featureSpecFlat -= ((UWord32)tmpU32>>14); // Q10
            return;
        }
    }
    //ratio and inverse log: check for case of log(0)
    zeros = _NormU32(avgSpectralFlatnessDen);
    frac = (Word16)(((avgSpectralFlatnessDen << zeros) & 0x7FFFFFFF) >> 23);
    // log2(avgSpectralFlatnessDen)
    tmp32 = (Word32)(((31 - zeros) << 8) + NS_kLogTableFrac[frac]); // Q8
    logCurSpectralFlatness = (Word32)avgSpectralFlatnessNum + (7L<<15) - (tmp32 << 7); // Q(8+STAGES-1)
    logCurSpectralFlatness = logCurSpectralFlatness<<2; // Q17
    tmp32 = (Word32)(0x00020000 | (labs(logCurSpectralFlatness) & 0x0001FFFF)); //Q17
    intPart = 7 - (Word16)(logCurSpectralFlatness>>17); // Shift 7 to get the output in Q10 (from Q17 = -17+10)
    if (intPart > 0) currentSpectralFlatness = tmp32>>intPart;
    else             currentSpectralFlatness = tmp32<<(-intPart);
    
    //time average update of spectral flatness feature
    tmp32 = SPECT_FLAT_TAVG_Q14*(currentSpectralFlatness - (Word32)st->featureSpecFlat); // Q24
    st->featureSpecFlat = st->featureSpecFlat + (tmp32>>14); // Q10
}

// Compute the difference measure between input spectrum and a template/learned noise spectrum
// magn_tmp is the input spectrum
// the reference/template spectrum is  st->magn_avg_pause[i]
// returns (normalized) spectral difference in st->featureSpecDiff
void NS_ComputeSpectralDifference(NS_STRUCT* st, UWord16* magnIn)
{
    // This is to be calculated: avgDiffNormMagn = var(magnIn) - cov(magnIn, magnAvgPause)^2 / var(magnAvgPause)
    UWord32 tmpU32no1, tmpU32no2;
    UWord32 varMagnUFX, varPauseUFX, avgDiffNormMagnUFX;
    Word32 tmp32no1, tmp32no2;
    Word32 avgPauseFX, avgMagnFX, covMagnPauseFX;
    Word32 maxPause, minPause;
    Word16 tmp16no1;
    int i, norm32, nShifts;
    
    // compute average quantities
    for (maxPause=minPause=st->avgMagnPause[0], avgPauseFX=i=0; i < HALF_ANAL_BLOCKL; i++)
    {
        avgPauseFX += st->avgMagnPause[i]; // Compute mean of magn_pause in Q(prevQMagn)
        maxPause = max(maxPause, st->avgMagnPause[i]);
        minPause = min(minPause, st->avgMagnPause[i]);
    }
    
	// normalize by replacing div of "HALF_ANAL_BLOCKL" with "STAGES-1" shifts
    avgPauseFX = avgPauseFX>>7;
    avgMagnFX = (Word32)st->sumMagn>>7;
    
	// Largest possible deviation in magnPause for (co)var calculations
    tmp32no1 = max(maxPause - avgPauseFX, avgPauseFX - minPause);
    // Get number of shifts to make sure we don't get wrap around in varPause
    nShifts = max(0, 18 - _NormW32(tmp32no1));
    
    varMagnUFX = 0;
    varPauseUFX = 0;
    covMagnPauseFX = 0;
    for (i = 0; i < HALF_ANAL_BLOCKL; i++)
    {
        // Compute var and cov of magn and magn_pause
        tmp16no1 = (Word16)((Word32)magnIn[i] - avgMagnFX);
        tmp32no2 = st->avgMagnPause[i] - avgPauseFX;
        varMagnUFX += (UWord32)tmp16no1*tmp16no1; // Q(2*qMagn)
        tmp32no1 = tmp32no2*tmp16no1; // Q(prevQMagn+qMagn)
        covMagnPauseFX += tmp32no1; // Q(prevQMagn+qMagn)
        tmp32no1 = tmp32no2>>nShifts; // Q(prevQMagn-minPause)
        varPauseUFX += tmp32no1*tmp32no1; // Q(2*(prevQMagn-minPause))
    }
    //update of average magnitude spectrum: Q(-2*STAGES) and averaging replaced by shifts
    st->curAvgMagnEnergy += st->magnEnergy>>(2 * st->normData + 8 - 1);
    
    avgDiffNormMagnUFX = varMagnUFX; // Q(2*qMagn)
    if ((varPauseUFX) && (covMagnPauseFX))
    {
        tmpU32no1 = (UWord32)labs(covMagnPauseFX); // Q(prevQMagn+qMagn)
        norm32 = _NormU32(tmpU32no1) - 16;
        if (norm32 > 0) tmpU32no1 = (UWord32)tmpU32no1<<norm32; // Q(prevQMagn+qMagn+norm32)
	    else            tmpU32no1 = tmpU32no1>>(-norm32); // Q(prevQMagn+qMagn+norm32)
        
        tmpU32no2 = tmpU32no1*tmpU32no1; // Q(2*(prevQMagn+qMagn-norm32))
        
        nShifts += norm32;
        nShifts <<= 1;
        if (nShifts < 0)
	    {
            varPauseUFX >>= (-nShifts); // Q(2*(qMagn+norm32+minPause))
            nShifts = 0;
        }
        if (varPauseUFX > 0) avgDiffNormMagnUFX -= min(avgDiffNormMagnUFX, (tmpU32no2/varPauseUFX)>>nShifts); // Q(2*qMagn)
	    else                 avgDiffNormMagnUFX = 0;
    }
    //normalize and compute time average update of difference feature
    tmpU32no1 = avgDiffNormMagnUFX>>(2 * st->normData);
    if (st->featureSpecDiff > tmpU32no1) st->featureSpecDiff -= (UWord32)(st->featureSpecDiff - tmpU32no1)*SPECT_DIFF_TAVG_Q8>>8; // Q(-2*STAGES)
    else                                 st->featureSpecDiff += (UWord32)(tmpU32no1 - st->featureSpecDiff)*SPECT_DIFF_TAVG_Q8>>8; // Q(-2*STAGES)
}

// Compute speech/noise probability
// speech/noise probability is returned in: probSpeechFinal  //snrLocPrior is the prior SNR for each frequency (in Q11) //snrLocPost is the post SNR for each frequency (in Q11)
void NS_SpeechNoiseProb(NS_STRUCT* st, UWord16* nonSpeechProbFinal, UWord32* priorLocSnr, UWord32* postLocSnr)
{
    UWord32 zeros, num, den, tmpU32no1, tmpU32no2, tmpU32no3;
    
    Word32 invLrtFX, indPriorFX, tmp32, tmp32no1, tmp32no2, besselTmpFX32;
    Word32 frac32, logTmp;
    Word32 logLrtTimeAvgKsumFX;
    
    Word16 indPriorFX16;
    Word16 tmp16, tmp16no1, tmp16no2, tmpIndFX, tableIndex, frac, intPart;
    
    int i, normTmp, normTmp2, nShifts;
    
    // compute feature based on average LR factor
    // this is the average over all frequencies of the smooth log LRT
    logLrtTimeAvgKsumFX = 0;
    for (i = 0; i < HALF_ANAL_BLOCKL; i++)
    {
        besselTmpFX32 = (Word32)postLocSnr[i]; // Q11
        normTmp = _NormU32(postLocSnr[i]);
        num = (UWord32)postLocSnr[i]<<normTmp; // Q(11+normTmp)
        if (normTmp > 10) den = (UWord32)priorLocSnr[i]<<(normTmp - 11); // Q(normTmp)
	    else              den = priorLocSnr[i]>>(11 - normTmp); // Q(normTmp)
        
        if (den > 0) besselTmpFX32 -= num/den; // Q11
	    else         besselTmpFX32 -= num; // Q11
        
        // st->logLrtTimeAvg[i] += LRT_TAVG * (besselTmp - log(snrLocPrior) - st->logLrtTimeAvg[i]);
        // Here, LRT_TAVG = 0.5
        zeros = _NormU32(priorLocSnr[i]);
        frac32 = (Word32)(((priorLocSnr[i] << zeros) & 0x7FFFFFFF) >> 19);
        tmp32 = (frac32*frac32*-43L)>>19;
        tmp32 += (frac32*5412L>>12);
        frac32 = tmp32 + 37;
        // tmp32 = log2(priorLocSnr[i])
        tmp32 = (Word32)(((31 - zeros) << 12) + frac32) - (11 << 12); // Q12
        logTmp = (tmp32*178)>>8; // log2(priorLocSnr[i])*log(2)
        tmp32no1 = (logTmp + st->logLrtTimeAvgW32[i])>>1; // Q12
        st->logLrtTimeAvgW32[i] += (besselTmpFX32 - tmp32no1); // Q12
        
        logLrtTimeAvgKsumFX += st->logLrtTimeAvgW32[i]; // Q12
    }
    st->featureLogLrt = logLrtTimeAvgKsumFX * 5L>>(8 + 10); // 5 = BIN_SIZE_LRT / 2
    // done with computation of LR factor
    
    //compute the indicator functions
    
    // average LRT feature
    // FLOAT code // indicator0 = 0.5 * (tanh(widthPrior * (logLrtTimeAvgKsum - threshPrior0)) + 1.0);
    tmpIndFX = 16384; // Q14(1.0)
    tmp32no1 = logLrtTimeAvgKsumFX - st->thresholdLogLrt; // Q12
    nShifts = 7 - 8; // WIDTH_PR_MAP_SHIFT - STAGES + 5;
    //use larger width in tanh map for pause regions
    if (tmp32no1 < 0)
    {
        tmpIndFX = 0;
        tmp32no1 = -tmp32no1;    
        nShifts++; //widthPrior = widthPrior * 2.0;
    }
    if(nShifts>=0) tmp32no1 = tmp32no1<<nShifts; // Q14
    else           tmp32no1 = tmp32no1>>(-nShifts); // Q14
    
    // compute indicator function: sigmoid map
    tableIndex = (Word16)(tmp32no1>>14);
    if ((tableIndex < 16) && (tableIndex >= 0)) 
    {
        tmp16no2 = kIndicatorTable[tableIndex];
        tmp16no1 = kIndicatorTable[tableIndex + 1] - kIndicatorTable[tableIndex];
        frac = (Word16)(tmp32no1 & 0x00003fff); // Q14
        tmp16no2 += (Word16)((Word32)tmp16no1*frac>>14);
        if (tmpIndFX == 0) tmpIndFX = 8192 - tmp16no2; // Q14
	    else               tmpIndFX = 8192 + tmp16no2; // Q14
    }
    indPriorFX = (Word32)st->weightLogLrt*tmpIndFX; // 6*Q14
    
    //spectral flatness feature
    if (st->weightSpecFlat) 
    {
        tmpU32no1 = (UWord32)st->featureSpecFlat*400; // Q10
        tmpIndFX = 16384; // Q14(1.0)
        //use larger width in tanh map for pause regions
        tmpU32no2 = st->thresholdSpecFlat - tmpU32no1; //Q10
        nShifts = 4;
        if (st->thresholdSpecFlat < tmpU32no1)
	    {
            tmpIndFX = 0;
            tmpU32no2 = tmpU32no1 - st->thresholdSpecFlat;
            nShifts++; //widthPrior = widthPrior * 2.0;
        }
        tmp32no1 = (Word32)((UWord32)tmpU32no2<<nShifts)/25; //Q14
        tmpU32no1 = (tmpU32no2<<nShifts)/25; //Q14
        // compute indicator function: sigmoid map // indicator1 = 0.5 * (tanh(sgnMap * widthPrior * (threshPrior1 - tmpFloat1)) + 1.0);
        tableIndex = (Word16)((UWord32)tmpU32no1>>14);
        if (tableIndex < 16)
	    {
            tmp16no2 = kIndicatorTable[tableIndex];
            tmp16no1 = kIndicatorTable[tableIndex + 1] - kIndicatorTable[tableIndex];
            frac = (Word16)(tmpU32no1 & 0x00003fff); // Q14
            tmp16no2 += (Word16)((Word32)tmp16no1*frac>>14);
            if (tmpIndFX) tmpIndFX = 8192 + tmp16no2; // Q14
	        else          tmpIndFX = 8192 - tmp16no2; // Q14
        }
        indPriorFX += (Word32)st->weightSpecFlat*tmpIndFX; // 6*Q14
    }
    
    //for template spectral-difference
    if (st->weightSpecDiff)
    {
        tmpU32no1 = 0;
        if (st->featureSpecDiff) 
	    {
            normTmp = min(20 - 8, _NormU32(st->featureSpecDiff));
            tmpU32no1 = (UWord32)st->featureSpecDiff<<normTmp; // Q(normTmp-2*STAGES)
            tmpU32no2 = st->timeAvgMagnEnergy>>(20 - 8 - normTmp);
            if (tmpU32no2 > 0) tmpU32no1 = tmpU32no1/tmpU32no2; // Q(20 - STAGES)
	        else               tmpU32no1 = (UWord32)(0x7fffffff);
        }
        tmpU32no3 = ((UWord32)st->thresholdSpecDiff<<17)/25;
        tmpU32no2 = tmpU32no1 - tmpU32no3;
        nShifts = 1;
        tmpIndFX = 16384; // Q14(1.0)
        //use larger width in tanh map for pause regions
        if (tmpU32no2 & 0x80000000) 
	    {
            tmpIndFX = 0;
            tmpU32no2 = tmpU32no3 - tmpU32no1;
            nShifts--; //widthPrior = widthPrior * 2.0;
        }
        tmpU32no1 = tmpU32no2>>nShifts;
        // compute indicator function: sigmoid map  //  indicator2 = 0.5 * (tanh(widthPrior * (tmpFloat1 - threshPrior2)) + 1.0);
        tableIndex = (Word16)((UWord32)tmpU32no1>>14);
        if (tableIndex < 16)
	    {
            tmp16no2 = kIndicatorTable[tableIndex];
            tmp16no1 = kIndicatorTable[tableIndex + 1] - kIndicatorTable[tableIndex];
            frac = (Word16)(tmpU32no1 & 0x00003fff); // Q14
            tmp16no2 += (Word16)(((Word32)tmp16no1*frac+8192)>>14);
            if (tmpIndFX) tmpIndFX = 8192 + tmp16no2;
	        else          tmpIndFX = 8192 - tmp16no2;
        }
        indPriorFX += (Word32)st->weightSpecDiff*tmpIndFX; // 6*Q14
    }

    //combine the indicator function with the feature weights
    // indPrior = 1 - (weightIndPrior0 * indicator0 + weightIndPrior1 * indicator1 + weightIndPrior2 * indicator2); // Q14
    indPriorFX16 = (Word16)((98307L-indPriorFX)/6);

    // done with computing indicator function
    //compute the prior probability // st->priorNonSpeechProb += PRIOR_UPDATE * (indPriorNonSpeech - st->priorNonSpeechProb);
    tmp16 = indPriorFX16 - st->priorNonSpeechProb; // Q14
    st->priorNonSpeechProb += (Word16)((Word32)PRIOR_UPDATE_Q14*tmp16>>14); // Q14

    //final speech probability: combine prior model with LR factor:
	for(i=0; i<HALF_ANAL_BLOCKL; i++) nonSpeechProbFinal[i]=0;

    if (st->priorNonSpeechProb > 0)
    {
        for (i = 0; i < HALF_ANAL_BLOCKL; i++)
	    {
            // invLrt = exp(st->logLrtTimeAvg[i]);
            // invLrt = st->priorSpeechProb * invLrt;
            // nonSpeechProbFinal[i] = (1.0 - st->priorSpeechProb) / (1.0 - st->priorSpeechProb + invLrt);
            // invLrt = (1.0 - st->priorNonSpeechProb) * invLrt;
            // nonSpeechProbFinal[i] = st->priorNonSpeechProb / (st->priorNonSpeechProb + invLrt);
            if (st->logLrtTimeAvgW32[i] < 65300)
	        {
                tmp32no1 = st->logLrtTimeAvgW32[i]*23637L>>14; // Q12
                intPart = (Word16)(tmp32no1>>12);
                if (intPart < -8) intPart = -8;
                frac = (Word16)(tmp32no1 & 0x00000fff); // Q12
                
                // Quadratic approximation of 2^frac
                tmp32no2 = (11L * frac * frac >>17) + (21L * frac>>5); // Q12
	    	    if(intPart >= 4) tmp32no2 = tmp32no2<<(intPart - 4);
	    	    else             tmp32no2 = tmp32no2>>(4 - intPart);
                invLrtFX = (1L<<(8 + intPart)) + tmp32no2; // Q8
                
                normTmp = _NormW32(invLrtFX);
                normTmp2 = _NormW16((Word16)(16384 - st->priorNonSpeechProb));
                if (normTmp + normTmp2 >= 7)
	    	    {
                    if (normTmp + normTmp2 < 15) invLrtFX = (invLrtFX>>(15 - normTmp2 - normTmp))*(16384 - st->priorNonSpeechProb)>>(normTmp + normTmp2 - 7); // Q14
	    	        else                         invLrtFX = invLrtFX*(16384 - st->priorNonSpeechProb)>>8; // Q14
                    nonSpeechProbFinal[i] = (UWord16)(((Word32)st->priorNonSpeechProb<<8)/((Word32)st->priorNonSpeechProb + invLrtFX)); // Q8
                }
            }
        }
    }
}

// Transform input (speechFrame) to frequency domain magnitude (magnU16)
void NS_DataAnalysis(NS_STRUCT* st, short* speechFrame, UWord16* magnU16)
{
    UWord32 tmpU32no1;

    Word32   tmp_2_w32 = 0;
    Word32   sum_log_magn = 0;
    Word32   sum_log_i_log_magn = 0;
    
    UWord16  sum_log_magn_u16 = 0;
    UWord16  tmp_u16 = 0;
    
    Word16   sum_log_i = 0;
    Word16   sum_log_i_square = 0;
    Word16   frac = 0;
    Word16   log2 = 0;
    Word16   matrix_determinant = 0;
    Word16   maxWinData;
    
    int i, j, zeros;
    int net_norm = 0;
    int right_shifts_in_magnU16 = 0;
    int right_shifts_in_initMagnEst;
    Word16 t, smax, winData_buff[ANA_LEN*2];
	Word32 low, up, vv, v2;
	Word16 k, m, n, istep, wr, wi;
    Word32 tr32, ti32, qr32, qi32;
    
    // Update analysis buffer for lower band, and window data before FFT.
	for (i=0; i<ANA_LEN-FRM_LEN; i++) st->analysisBuffer[i]=st->analysisBuffer[i+FRM_LEN];
	for (i=0; i<FRM_LEN; i++) st->analysisBuffer[i+ANA_LEN-FRM_LEN]=speechFrame[i];
    // Window data before FFT.
    for (smax=i=0; i<ANA_LEN; i++)
	{
		winData_buff[i] = (Word16)((kBlocks160w256x[i]*st->analysisBuffer[i]+8192L)>>14); // Q0
		smax = max(smax, (Word16)labs(winData_buff[i]));
	}
    // Get input energy
    t = _NormW32((Word32)smax*smax);
    if (smax == 0) st->scaleEnergyIn = 0; // Since norm(0) returns 0
	else           st->scaleEnergyIn = (t > 9) ? 0 : 9 - t;

    for (st->energyIn=i=0; i<ANA_LEN; i++) st->energyIn += (Word32)winData_buff[i]*winData_buff[i]>>st->scaleEnergyIn;

	// Acquire norm for winData_buff
	maxWinData = (short)min(smax, 32767); // Guard the case for abs(-32768).
    st->normData = _NormW16(maxWinData);
	
	// Reset zero input flag
	st->zeroInputSignal = 0;
    if (maxWinData == 0) {st->zeroInputSignal = 1;   return;} // Treat zero input separately.
    
    // Determine the net normalization in the frequency domain
    net_norm = 8 - st->normData;
    // Track lowest normalization factor and use it to prevent wrap around in shifting
    right_shifts_in_magnU16 = st->normData - st->minNorm;
    right_shifts_in_initMagnEst = max(-right_shifts_in_magnU16, 0);
    st->minNorm -= right_shifts_in_initMagnEst;
    right_shifts_in_magnU16 = max(right_shifts_in_magnU16, 0);
    
    // create real_buff as winData_buff with imag. part = zeros, normalize it // Decimation in time. Swap the elements with bit-reversed indexes.
    for (i=0; i<ANA_LEN; i++)
	{
		st->real[i] = winData_buff[bit_rev[i]]<<st->normData; // Q(normData)
	    st->imag[i] = 0; //set zeros to the imaginary parts for complex forward FFT input
	}
    // FFT output will be in st->real[], st->imag[].
	for(k=8, istep=2, n=1; n<ANA_LEN; n*=2, istep*=2, k--)
    {
        for (m = 0; m < n; m++)
        {
            j = m << k;
            wr = rot_tab[j];
            wi = -rot_tab[j+1];
            for (i = m; i < ANA_LEN; i += istep)
            {
                j = i + n;
                tr32 = (wr*st->real[j] - wi*st->imag[j] + 1L)>>1;
                ti32 = (wr*st->imag[j] + wi*st->real[j] + 1L)>>1;
                qr32 = (Word32)st->real[i]<<14;
                qi32 = (Word32)st->imag[i]<<14;
                st->real[j] = (Word16)((qr32 - tr32 + 16384L)>>15);
                st->imag[j] = (Word16)((qi32 - ti32 + 16384L)>>15);
                st->real[i] = (Word16)((qr32 + tr32 + 16384L)>>15);
                st->imag[i] = (Word16)((qi32 + ti32 + 16384L)>>15);
            }
        }
    }
    
    for (st->magnEnergy=st->sumMagn=i=0; i<HALF_ANAL_BLOCKL; i++)
	{
        tmpU32no1 = (Word32)st->real[i]*st->real[i] + (Word32)st->imag[i]*st->imag[i]; // magnitude spectrum in Q(2*(normData-STAGES))
        st->magnEnergy += tmpU32no1; // Q(2*(normData-STAGES))
        
        //magnU16[i] = SqrtFloor(tmpU32no1); // Q(normData-STAGES)
        for(low=0, up=32767, j=0; j<16; j++)
	    {
	    	vv=(low+up)>>1;
        	v2=vv*vv;		
        	if(tmpU32no1>(UWord32)v2) low=vv;
		    else if(tmpU32no1<(UWord32)v2) up=vv;
		    else {low=vv; break;}
	    }
        magnU16[i] = (UWord16)low;
        st->sumMagn += magnU16[i]; // Q(normData-STAGES)
    }

	// Gather information during startup for noise parameter estimation
    if (st->blockIndex < END_STARTUP_SHORT)
    {
		// Update initMagnEst
        for (i=0; i<HALF_ANAL_BLOCKL; i++) st->initMagnEst[i] = (st->initMagnEst[i]>>right_shifts_in_initMagnEst) + (magnU16[i]>>right_shifts_in_magnU16); // Q(minNorm-STAGES)
        
		// For pink noise estimation. Collect data neglecting lower frequency band // log2(magnU16(i)) in Q8
		for (sum_log_i_log_magn=sum_log_magn=0, i=KSTARTBAND; i<HALF_ANAL_BLOCKL; i++)
	    {
			zeros = _norm((Word32)magnU16[i])-16;
            frac = (magnU16[i]<<zeros)&0x3FFF;                
            log2 = (Word16)(((14-zeros) << 8) + NS_kLogTableFrac[frac>>6]);
			if (magnU16[i]==0) log2 = 0;
            sum_log_magn += log2; // Q8                
            sum_log_i_log_magn += (Word32)kLogIndex[i]*log2>>3;  // Q17
        }
        
        //compute simplified noise model during startup
        // Estimate White noise
		// Update the average magnitude spectrum, used as noise estimate.
        // Switch whiteNoiseLevel to Q(minNorm-STAGES) and Shift to same Q-domain as whiteNoiseLevel.
        st->whiteNoiseLevel = (st->whiteNoiseLevel>>right_shifts_in_initMagnEst) + ((UWord32)st->sumMagn*st->overdrive>>(8+8)>>right_shifts_in_magnU16); // Q(minNorm-STAGES)
        
        // Estimate Pink noise parameters
        // Denominator used in both parameter estimates.
        // The value is only dependent on the size of the frequency band (KSTARTBAND) and to reduce computational complexity stored in a table (kDeterminantEstMatrix[])
        matrix_determinant = kDeterminantEstMatrix[KSTARTBAND]; // Q0
        sum_log_i = kSumLogIndex[KSTARTBAND]; // Q5
        sum_log_i_square = kSumSquareLogIndex[KSTARTBAND]; // Q2
        // Necessary number of shifts to fit sum_log_magn in a word16
        zeros = 16 - _NormW32(sum_log_magn);
        if (zeros < 0) zeros = 0;
        
        sum_log_magn_u16 = (UWord16)(sum_log_magn<<1>>zeros);//Q(9-zeros)
        
        // Calculate and update pinkNoiseNumerator. Result in Q11.
        tmpU32no1 = (UWord32)sum_log_i_log_magn>>12; // Q5
        
        // Shift the largest value of sum_log_i and tmp32no3 before multiplication
        tmp_u16 = (UWord16)sum_log_i<<1; // Q6
        if (tmpU32no1 < (UWord32)sum_log_i) tmp_u16 = tmp_u16>>zeros;
	    else                                tmpU32no1 = tmpU32no1>>zeros;
        
        matrix_determinant = matrix_determinant>>zeros; // Q(-zeros)
        tmp_2_w32 = ((Word32)sum_log_i_square*sum_log_magn_u16 - tmpU32no1*tmp_u16)/matrix_determinant + ((Word32)net_norm<<11); // Q11
        if (tmp_2_w32 < 0) tmp_2_w32 = 0;
        
        st->pinkNoiseNumerator += tmp_2_w32; // Q11
        
        // Calculate and update pinkNoiseExp. Result in Q14.
        tmp_2_w32 = (Word32)sum_log_i*sum_log_magn_u16 - (Word32)(HALF_ANAL_BLOCKL - KSTARTBAND)*(sum_log_i_log_magn>>(3 + zeros)); // Q(14-zeros)
        if (tmp_2_w32 > 0) st->pinkNoiseExp += min(max(tmp_2_w32/matrix_determinant, 0), 16384); // If the exponential parameter is negative force it to zero, which means a flat spectrum. // Q14
    }
}

void NS_DataSynthesis(NS_STRUCT* st, short* outFrame)
{
    Word32 energyOut, tmp32;
    int scaleEnergyOut = 0;
	Word16 t, smax, energyRatio, gainFactor;
    Word16 i, j, k, m, n, istep, wr, wi, shift, outCIFFT, re[ANA_LEN], im[ANA_LEN];
    Word32 tr32, ti32, qr32, qi32, round2; 

    if (st->zeroInputSignal==0)
	{            
		// Filter the data in the frequency domain, and create spectrum.
	    for(i=0; i<HALF_ANAL_BLOCKL; i++)
	    {
			re[i]= (Word16)((Word32)st->real[i]*st->noiseSupFilter[i]>>14); // Q(normData-STAGES)
		    im[i]= -(Word16)(-(Word32)st->imag[i]*st->noiseSupFilter[i]>>14); // Q(normData-STAGES)
	    }
		// For n-point FFT, first copy the first n + 2 elements into complex FFT, then construct the remaining n - 2 elements by real FFT's conjugate-symmetric properties.
        for (i=HALF_ANAL_BLOCKL; i<ANA_LEN; i++)
	    {
            re[i] = re[ANA_LEN-i];
            im[i] = -im[ANA_LEN-i];
	    }
	    // Decimation in time. Swap the elements with bit-reversed indexes.
        for (i=0; i<ANA_LEN; i++)
	    {
            st->real[i] = re[bit_rev[i]];
            st->imag[i] = im[bit_rev[i]];
	    }   
        
        // Inverse FFT output will be in rfft_out[].
	    for(outCIFFT=0, k=8, istep=2, n=1; n<ANA_LEN; n*=2, istep*=2, k--)
        {
            // variable scaling
            for (tmp32=i=0; i < ANA_LEN; i++) tmp32 = max(tmp32, max(labs(st->real[i]), labs(st->imag[i])));
	    	shift=tmp32*19777>>28;
	    	outCIFFT += shift;
	    	round2 = 8192L<<shift;
	    	shift += 14;
        
            for (m=0; m<n; m++)
            {
                j = m << k;
                wr = rot_tab[j];
                wi = rot_tab[j+1];
                for (i=m; i<ANA_LEN; i += istep)
                {
                    j = i + n;
                    tr32 = ((Word32)wr*st->real[j] - wi*st->imag[j] + 1)>>1;
                    ti32 = ((Word32)wr*st->imag[j] + wi*st->real[j] + 1)>>1;
                    qr32 = (Word32)st->real[i]<<14;
                    qi32 = (Word32)st->imag[i]<<14;
        
                    st->real[j] = (Word16)((qr32 - tr32+round2)>>shift);
                    st->imag[j] = (Word16)((qi32 - ti32 + round2)>>shift);
                    st->real[i] = (Word16)((qr32 + tr32 + round2)>>shift);
                    st->imag[i] = (Word16)((qi32 + ti32 + round2)>>shift);
                }
            }
        }
        
        //Denormalize(st, rfft_out, outCIFFT); // Strip out the imaginary parts of the complex inverse FFT output.
		for (i=0; i<ANA_LEN; i++) st->real[i] = (Word16)max(min((Word32)st->real[i]<<outCIFFT>>st->normData, 32767), -32768); // Q0
        
        //scale factor: only do it after END_STARTUP_LONG time
        gainFactor = 8192; // 8192 = Q13(1.0)
        if ((st->gainMap==1)&&(st->blockIndex>END_STARTUP_LONG)&&(st->energyIn>0))
        {
            for (smax=0, i=0; i<ANA_LEN; i++) smax = max(smax, (Word16)labs(st->real[i]));
            t = _NormW32((Word32)smax*smax);
            if (smax == 0) scaleEnergyOut = 0; // Since norm(0) returns 0
	        else           scaleEnergyOut = (t > 9) ? 0 : 9 - t;
            for (energyOut=i=0; i < ANA_LEN; i++) energyOut += (Word32)st->real[i]*st->real[i]>>scaleEnergyOut; // Q(-scaleEnergyOut)

            if ((scaleEnergyOut==0) && !(energyOut & 0x7f800000)) 
	        {
	        	if(8 + scaleEnergyOut >= st->scaleEnergyIn) energyOut = energyOut<<(8 + scaleEnergyOut - st->scaleEnergyIn);
	        	else                                        energyOut = energyOut>>(st->scaleEnergyIn - 8 - scaleEnergyOut);
	        }
	        else
	        {
	        	st->energyIn = st->energyIn>>(8 + scaleEnergyOut - st->scaleEnergyIn); // Q(-8-scaleEnergyOut)
	        }
            
            // Limit the ratio to [0, 1] in Q8, i.e., [0, 256]  //assert(energyRatio < 257);
	        energyRatio = min(max((energyOut + (st->energyIn>>1))/st->energyIn, 0), 256);
            
            // all done in lookup tables now    
            //combine both scales with speech/noise prob: note prior (priorSpeechProb) is not frequency dependent
            // factor = st->priorSpeechProb*factor1 + (1.0-st->priorSpeechProb)*factor2; // original code
            gainFactor =  (Word16)((Word32)(16384 - st->priorNonSpeechProb)*kFactor1Table[energyRatio]>>14)
	      	            + (Word16)((Word32)st->priorNonSpeechProb*st->factor2Table[energyRatio]>>14); // Q13
        }  // out of flag_gain_map==1
        
        // Synthesis, read out fully processed segment, and update synthesis buffer.
        // synthesis
        for (i=0; i<ANA_LEN; i++) 
        {
            tmp32 = st->synthesisBuffer[i] + ((((kBlocks160w256x[i]*st->real[i]+8192L)>>14)*gainFactor+4096L)>>13); // Down shift with rounding  // Q0
	        st->synthesisBuffer[i] = (Word16)max(min(tmp32, 32767), -32768);
        }    
	}

    // read out fully processed segment
    for (i=0; i <FRM_LEN; i++) outFrame[i] = st->synthesisBuffer[i]; // Q0
    // update synthesis buffer
	for (i=0; i<ANA_LEN-FRM_LEN; i++) st->synthesisBuffer[i]=st->synthesisBuffer[i+FRM_LEN];
	for (i=0; i<FRM_LEN; i++) st->synthesisBuffer[i+ANA_LEN-FRM_LEN]=0;
}

// Initialize state // Initialization of struct //fs == 16000
void NS_Init(NS_STRUCT* st, Word32 fs, Word32 mode)
{
    int i;    
    
	st->fs = fs;
    st->thresholdLogLrt = 212644; //default threshold for LRT feature
    st->maxLrt = 524288;
    st->minLrt = 104858;

	for(i=0; i<ANA_LEN; i++) st->dataBufHBFX[i]=st->synthesisBuffer[i]=st->analysisBuffer[i]=0;
	for(i=0; i<HALF_ANAL_BLOCKL; i++) st->noiseEstQuantile[i]=0;

    for (i=0; i < SIMULT * HALF_ANAL_BLOCKL; i++) st->noiseEstLogQuantile[i] = 2048; // Q8
	for (i=0; i < SIMULT * HALF_ANAL_BLOCKL; i++) st->noiseEstDensity[i] = 153; // Q9
    
	st->noiseEstCounter[0] = 66;
	st->noiseEstCounter[1] = 133;
	st->noiseEstCounter[2] = END_STARTUP_LONG;
    // Initialize suppression filter with ones
    for (i = 0; i<HALF_ANAL_BLOCKL; i++) st->noiseSupFilter[i] = 16384;
    
    //initialize variables for new method
    st->priorNonSpeechProb = 8192; // Q14(0.5) prior probability for speech/noise
    for (i = 0; i < HALF_ANAL_BLOCKL; i++)
    {
        st->prevMagnU16[i] = 0;
        st->prevNoiseU32[i] = 0; //previous noise-spectrum
        st->logLrtTimeAvgW32[i] = 0; //smooth LR ratio
        st->avgMagnPause[i] = 0; //conservative noise spectrum estimate
        st->initMagnEst[i] = 0; //initial average magnitude spectrum
    }
    
    //feature quantities
    st->thresholdSpecDiff = 50; //threshold for difference feature: determined on-line
    st->thresholdSpecFlat = 20480; //threshold for flatness: determined on-line
    st->featureLogLrt = st->thresholdLogLrt; //average LRT factor (= threshold)
    st->featureSpecFlat = st->thresholdSpecFlat; //spectral flatness (= threshold)
    st->featureSpecDiff = st->thresholdSpecDiff; //spectral difference (= threshold)
    st->weightLogLrt = 6; //default weighting par for LRT feature
    st->weightSpecFlat = 0; //default weighting par for spectral flatness feature
    st->weightSpecDiff = 0; //default weighting par for spectral difference feature
    
    st->curAvgMagnEnergy = 0; //window time-average of input magnitude spectrum
    st->timeAvgMagnEnergy = 0; //normalization for spectral difference
    st->timeAvgMagnEnergyTmp = 0; //normalization for spectral difference
    
    //histogram quantities: used to estimate/update thresholds for features

	for(i=0; i<HIST_PAR_EST; i++) st->histLrt[i]=st->histSpecDiff[i]=st->histSpecFlat[i]=0;
    
    st->blockIndex = -1; //frame counter
    
    //st->modelUpdate    = 500;   //window for update
    st->modelUpdate = (1 << STAT_UPDATES); //window for update
    st->cntThresUpdate = 0; //counter feature thresholds updates
    
    st->sumMagn = 0;
    st->magnEnergy = 0;
    st->prevQMagn = 0;
    st->qNoise = 0;
    st->prevQNoise = 0;
    
    st->energyIn = 0;
    st->scaleEnergyIn = 0;
    
    st->whiteNoiseLevel = 0;
    st->pinkNoiseNumerator = 0;
    st->pinkNoiseExp = 0;
    st->minNorm = 15; // Start with full scale
    st->zeroInputSignal = 0;

    st->aggrMode = mode;
    if (mode == 0)
    {
        st->overdrive = 256; // Q8(1.0)
        st->denoiseBound = 8192; // Q14(0.5)
        st->gainMap = 0; // No gain compensation
    }
    else if (mode == 1)
    {
        st->overdrive = 256; // Q8(1.0)
        st->denoiseBound = 4096; // Q14(0.25)
        st->factor2Table = kFactor2Aggressiveness1;
        st->gainMap = 1;
    }
    else if (mode == 2)
    {
        st->overdrive = 282; // ~= Q8(1.1)
        st->denoiseBound = 2048; // Q14(0.125)
        st->factor2Table = kFactor2Aggressiveness2;
        st->gainMap = 1;
    } 
    else if (mode == 3)
    {
        st->overdrive = 320; // Q8(1.25)
        st->denoiseBound = 1475; // ~= Q14(0.09)
        st->factor2Table = kFactor2Aggressiveness3;
        st->gainMap = 1;
    }
}

void NS_run(NS_STRUCT* st, short* speechFrame, short* outFrame) 
{
    // main routine for noise suppression    
    UWord32 tmpU32no1, tmpU32no2, tmpU32no3;
    UWord32 satMax, maxNoiseU32;
    UWord32 tmpMagnU32, tmpNoiseU32;
    UWord32 nearMagnEst;
    UWord32 noiseUpdateU32;
    UWord32 noiseU32[HALF_ANAL_BLOCKL], postLocSnr[HALF_ANAL_BLOCKL], priorLocSnr[HALF_ANAL_BLOCKL], prevNearSnr[HALF_ANAL_BLOCKL];
    UWord32 curNearSnr;
    UWord32 priorSnr;
    UWord32 noise_estimate = 0;
    UWord32 noise_estimate_avg = 0;
    UWord32 numerator = 0;
    
    Word32 tmp32no1, tmp32no2;
    Word32 pink_noise_num_avg = 0;
    
    UWord16 tmpU16no1;
    UWord16 magnU16[HALF_ANAL_BLOCKL], prevNoiseU16[HALF_ANAL_BLOCKL], nonSpeechProbFinal[HALF_ANAL_BLOCKL];
    UWord16 gammaNoise, prevGammaNoise;
    UWord16 noiseSupFilterTmp[HALF_ANAL_BLOCKL];
    
    Word16 qMagn, qNoise;
    Word16 pink_noise_exp_avg = 0;
    
    int i;
    int nShifts, postShifts, norm32no1, norm32no2;
    int flag, sign;
    int q_domain_to_use = 0;
    
    // Store speechFrame and transform to frequency domain
    NS_DataAnalysis(st, speechFrame, magnU16);
    
    if (st->zeroInputSignal==0)
    {
        // Update block index when we have something to process
        st->blockIndex++;
        
        // Norm of magn
        qMagn = st->normData - 8;
        
        // Compute spectral flatness on input spectrum
        NS_ComputeSpectralFlatness(st, magnU16);
        
        // quantile noise estimate
        NoiseEstimationC(st, magnU16, noiseU32, &qNoise);
        
        //noise estimate from previous frame
        for (i = 0; i < HALF_ANAL_BLOCKL; i++) prevNoiseU16[i] = (UWord16)(st->prevNoiseU32[i]>>11); // Q(prevQNoise)
        
        if (st->blockIndex < END_STARTUP_SHORT)
        {
            // Noise Q-domain to be used later; see description at end of section.
            q_domain_to_use = min(qNoise, st->minNorm - 8);
            
            // Calculate frequency independent parts in parametric noise estimate and calculate the estimate for the lower frequency band (same values for all frequency bins)
            if (st->pinkNoiseExp)
	        {
                pink_noise_exp_avg = (Word16)(st->pinkNoiseExp/(st->blockIndex + 1)); // Q14
                pink_noise_num_avg = st->pinkNoiseNumerator/(st->blockIndex + 1); // Q11
                NS_CalcParametricNoiseEstimate(st, pink_noise_exp_avg, pink_noise_num_avg, KSTARTBAND, &noise_estimate, &noise_estimate_avg);
            }
	        else 
	        {
                // Use white noise estimate if we have poor pink noise parameter estimates
                noise_estimate = st->whiteNoiseLevel; // Q(minNorm-STAGES)
                noise_estimate_avg = noise_estimate / (st->blockIndex + 1); // Q(minNorm-STAGES)
            }
            for (i = 0; i < HALF_ANAL_BLOCKL; i++)
	        {
                // Estimate the background noise using the pink noise parameters if permitted
                if ((st->pinkNoiseExp) && (i >= KSTARTBAND)) 
	            {
                    // Reset noise_estimate
                    noise_estimate = 0;
                    noise_estimate_avg = 0;
                    // Calculate the parametric noise estimate for current frequency bin
                    NS_CalcParametricNoiseEstimate(st, pink_noise_exp_avg, pink_noise_num_avg, i, &noise_estimate, &noise_estimate_avg);
                }
                // Calculate parametric Wiener filter
                noiseSupFilterTmp[i] = st->denoiseBound;
                if (st->initMagnEst[i])
	            {
                    // numerator = (initMagnEst - noise_estimate * overdrive) // Result in Q(8+minNorm-STAGES)
                    tmpU32no1 = (UWord32)noise_estimate*st->overdrive;
                    numerator = (UWord32)st->initMagnEst[i]<<8;
                    if (numerator > tmpU32no1)
	        	    {
                        // Suppression filter coefficient larger than zero, so calculate.
                        numerator -= tmpU32no1;
                        
                        // Determine number of left shifts in numerator for best accuracy after division
	        	        nShifts = min(max(_NormU32(numerator), 0), 6);
                        
                        // Shift numerator to Q(nShifts+8+minNorm-STAGES)
                        numerator = (UWord32)numerator<<nShifts;
                        
                        // Shift denominator to Q(nShifts-6+minNorm-STAGES)
                        tmpU32no1 = st->initMagnEst[i]>>(6 - nShifts);          
	        	        if (tmpU32no1 == 0) tmpU32no1 = 1; // This is only possible if numerator = 0, in which case we don't need any division.
                        tmpU32no2 = numerator/tmpU32no1; // Q14
                    
	        	        noiseSupFilterTmp[i] = (UWord16)min(max(tmpU32no2, st->denoiseBound), 16384); // Q14
                    }
                }
                // Weight quantile noise 'noiseU32' with modeled noise 'noise_estimate_avg'
                // 'noiseU32 is in Q(qNoise) and 'noise_estimate' in Q(minNorm-STAGES)
                // To guarantee that we do not get wrap around when shifting to the same domain we use the lowest one. Furthermore, we need to save 6 bits for the weighting.
                // 'noise_estimate_avg' can handle this operation by construction, but 'noiseU32' may not.
                
                // Shift 'noiseU32' to 'q_domain_to_use'
                tmpU32no1 = noiseU32[i]>>(qNoise - q_domain_to_use);
                // Shift 'noise_estimate_avg' to 'q_domain_to_use'
                tmpU32no2 = noise_estimate_avg>>(st->minNorm - 8 - q_domain_to_use);
                // Make a simple check to see if we have enough room for weighting 'tmpU32no1' without wrap around
                nShifts = 0;
                if (tmpU32no1 & 0xfc000000)
	            {
                    tmpU32no1 = tmpU32no1>>6;
                    tmpU32no2 = tmpU32no2>>6;
                    nShifts = 6;
                }
                tmpU32no1 *= st->blockIndex;
                tmpU32no2 *= (END_STARTUP_SHORT - st->blockIndex);
                // Add them together and divide by startup length
                noiseU32[i] = (tmpU32no1 + tmpU32no2)/END_STARTUP_SHORT;
                // Shift back if necessary
                noiseU32[i] = (UWord32)noiseU32[i]<<nShifts;
            }
            // Update new Q-domain for 'noiseU32'
            qNoise = q_domain_to_use;
        }
        // compute average signal during END_STARTUP_LONG time:
        // used to normalize spectral difference measure
        if (st->blockIndex < END_STARTUP_LONG)
        {
            // substituting division with shift ending up in Q(-2*STAGES)
            st->timeAvgMagnEnergyTmp += st->magnEnergy>>(2*st->normData + 8 - 1);
            st->timeAvgMagnEnergy = st->timeAvgMagnEnergyTmp/(st->blockIndex + 1);
        }
        
        //start processing at frames == converged+1
        // STEP 1: compute prior and post SNR based on quantile noise estimates        
        // compute direct decision (DD) estimate of prior SNR: needed for new method
        satMax = 1048575;// Largest possible value without getting overflow despite shifting 12 steps
        
        postShifts = 6 + qMagn - qNoise;
        nShifts = 5 - st->prevQMagn + st->prevQNoise;
        for (i = 0; i < HALF_ANAL_BLOCKL; i++)
        {
            // FLOAT:
            // post SNR
            // postLocSnr[i] = 0.0;
            // if (magn[i] > noise[i]) postLocSnr[i] = magn[i] / (noise[i] + 0.0001);
            // // previous post SNR
            // // previous estimate: based on previous frame with gain filter (smooth is previous filter)
            //
            // prevNearSnr[i] = st->prevMagnU16[i] / (st->noisePrev[i] + 0.0001) * (st->smooth[i]);
            //
            // // DD estimate is sum of two terms: current estimate and previous estimate
            // // directed decision update of priorSnr (or we actually store [2*priorSnr+1])
            // priorLocSnr[i] = DD_PR_SNR * prevNearSnr[i] + (1.0 - DD_PR_SNR) * (postLocSnr[i] - 1.0);
            
            // calculate post SNR: output in Q11
            postLocSnr[i] = 2048; // 1.0 in Q11
            tmpU32no1 = (UWord32)magnU16[i]<<6; // Q(6+qMagn)
            if (postShifts < 0) tmpU32no2 = noiseU32[i]>>(-postShifts); // Q(6+qMagn)
	        else                tmpU32no2 = (UWord32)noiseU32[i]<<postShifts; // Q(6+qMagn)
            
            if (tmpU32no1 > tmpU32no2)
	        {
                // Current magnitude larger than noise
                if (tmpU32no2 > 0) postLocSnr[i] = min(satMax, ((UWord32)tmpU32no1<<11)/tmpU32no2); // Q11
	            else               postLocSnr[i] = satMax;
            }
            
            // calculate prevNearSnr[i] and save for later stead of recalculating it later
            nearMagnEst = (UWord32)(UWord16)st->prevMagnU16[i]*(UWord16)st->noiseSupFilter[i]; // Q(prevQMagn+14)
            tmpU32no1 = (UWord32)nearMagnEst<<3; // Q(prevQMagn+17)
            tmpU32no2 = st->prevNoiseU32[i]>>nShifts; // Q(prevQMagn+6)
            
            if (tmpU32no2 > 0) tmpU32no1 = min(satMax, tmpU32no1/tmpU32no2); // Q11
	        else               tmpU32no1 = satMax; // Q11
            
            prevNearSnr[i] = tmpU32no1; // Q11
            
            //directed decision update of priorSnr
            priorSnr = (UWord32)prevNearSnr[i]*DD_PR_SNR_Q11 + (UWord32)(postLocSnr[i] - 2048)*ONE_MINUS_DD_PR_SNR_Q11 + 512; // Q22 (added 512 for rounding)
            // priorLocSnr = 1 + 2*priorSnr
            priorLocSnr[i] = 2048 + (priorSnr>>10); // Q11
        }  // end of loop over frequencies
        // done with step 1: DD computation of prior and post SNR
        
        // STEP 2: compute speech/noise likelihood
        
        //compute difference of input spectrum with learned/estimated noise spectrum
        NS_ComputeSpectralDifference(st, magnU16);
        //compute histograms for determination of parameters (thresholds and weights for features)
        //parameters are extracted once every window time (=st->modelUpdate)
        //counter update
        st->cntThresUpdate++;
        flag = (int)(st->cntThresUpdate == st->modelUpdate);
        //update histogram
        NS_FeatureParameterExtraction(st, flag);
        //compute model parameters
        if (flag)
        {
            st->cntThresUpdate = 0; // Reset counter
            //update every window: get normalization for spectral difference for next window estimate
            
            // Shift to Q(-2*STAGES)
            st->curAvgMagnEnergy = st->curAvgMagnEnergy>>STAT_UPDATES;
            
            tmpU32no1 = (st->curAvgMagnEnergy + st->timeAvgMagnEnergy + 1) >> 1; //Q(-2*STAGES)
            // Update featureSpecDiff
            if ((tmpU32no1 != st->timeAvgMagnEnergy) && (st->featureSpecDiff) && (st->timeAvgMagnEnergy > 0))
	        {
                norm32no1 = 0;
                tmpU32no3 = tmpU32no1;
                while (0xFFFF0000 & tmpU32no3)
	            {
                    tmpU32no3 >>= 1;
                    norm32no1++;
                }
                tmpU32no2 = st->featureSpecDiff;
                while (0xFFFF0000 & tmpU32no2)
	            {
                    tmpU32no2 >>= 1;
                    norm32no1++;
                }
                tmpU32no3 = (UWord32)tmpU32no3*tmpU32no2/st->timeAvgMagnEnergy;
                if (_NormU32(tmpU32no3) < norm32no1) st->featureSpecDiff = 0x007FFFFF;
	            else                                 st->featureSpecDiff = min(0x007FFFFF, (UWord32)tmpU32no3<<norm32no1);
            }
            
            st->timeAvgMagnEnergy = tmpU32no1; // Q(-2*STAGES)
            st->curAvgMagnEnergy = 0;
        }
        
        //compute speech/noise probability
        NS_SpeechNoiseProb(st, nonSpeechProbFinal, priorLocSnr, postLocSnr);
        
        //time-avg parameter for noise update
        gammaNoise = NOISE_UPDATE_Q8; // Q8
        
        maxNoiseU32 = 0;
        postShifts = st->prevQNoise - qMagn;
        nShifts = st->prevQMagn - qMagn;
        for (i = 0; i < HALF_ANAL_BLOCKL; i++)
        {
          // temporary noise update: use it for speech frames if update value is less than previous
          // the formula has been rewritten into: noiseUpdate = noisePrev[i] + (1 - gammaNoise) * nonSpeechProb * (magn[i] - noisePrev[i])
          if (postShifts < 0) tmpU32no2 = magnU16[i]>>(-postShifts); // Q(prevQNoise)
	      else                tmpU32no2 = (UWord32)magnU16[i]<<postShifts; // Q(prevQNoise)
          if (prevNoiseU16[i] > tmpU32no2) {sign = -1;   tmpU32no1 = prevNoiseU16[i] - tmpU32no2;}
	      else                             {sign = 1;    tmpU32no1 = tmpU32no2 - prevNoiseU16[i];}
          noiseUpdateU32 = st->prevNoiseU32[i]; // Q(prevQNoise+11)
          tmpU32no3 = 0;
          if ((tmpU32no1) && (nonSpeechProbFinal[i]))
	      {
              // This value will be used later, if gammaNoise changes
              tmpU32no3 = tmpU32no1*nonSpeechProbFinal[i]; // Q(prevQNoise+8)
              if (0x7c000000 & tmpU32no3) tmpU32no2 = (tmpU32no3>>5)*gammaNoise; // Q(prevQNoise+11)
	          else                        tmpU32no2 = tmpU32no3*gammaNoise>>5;   // Q(prevQNoise+11)
              if (sign > 0) noiseUpdateU32 += tmpU32no2; // Q(prevQNoise+11)
	          else          noiseUpdateU32 -= tmpU32no2; // Q(prevQNoise+11)
          }
        
          //increase gamma (i.e., less noise update) for frame likely to be speech
          prevGammaNoise = gammaNoise;
          //time-constant based on speech/noise state
          //increase gamma (i.e., less noise update) for frames likely to be speech
          gammaNoise = (nonSpeechProbFinal[i] < ONE_MINUS_PROB_RANGE_Q8) ? GAMMA_NOISE_TRANS_AND_SPEECH_Q8 : NOISE_UPDATE_Q8;
        
          if (prevGammaNoise != gammaNoise)
	      {
              // new noise update
              // this line is the same as above, only that the result is stored in a different variable and the gammaNoise has changed
              // noiseUpdate = noisePrev[i] + (1 - gammaNoise) * nonSpeechProb * (magn[i] - noisePrev[i])
              if (0x7c000000 & tmpU32no3) tmpU32no2 = (tmpU32no3>>5)*gammaNoise; // Q(prevQNoise+11)
	          else                        tmpU32no2 = (UWord32)tmpU32no3*gammaNoise>>5; // Q(prevQNoise+11)
              if (sign > 0) tmpU32no1 = st->prevNoiseU32[i] + tmpU32no2; // Q(prevQNoise+11)
	          else          tmpU32no1 = st->prevNoiseU32[i] - tmpU32no2; // Q(prevQNoise+11)
              if (noiseUpdateU32 > tmpU32no1) noiseUpdateU32 = tmpU32no1; // Q(prevQNoise+11)
          }
          noiseU32[i] = noiseUpdateU32; // Q(prevQNoise+11)
          if (maxNoiseU32 < noiseUpdateU32) maxNoiseU32 = noiseUpdateU32;
        
          // conservative noise update // if (prob_speech < PROB_RANGE) st->avgMagnPause[i] = st->avgMagnPause[i] + (1.0 - gamma_pause)*(magn[i] - st->avgMagnPause[i]);
          if(nShifts>0) tmp32no2 = st->avgMagnPause[i]>>nShifts;
	      else          tmp32no2 = st->avgMagnPause[i]<<(-nShifts);           
        
          if (nonSpeechProbFinal[i] > ONE_MINUS_PROB_RANGE_Q8)
	      {
              if (nShifts < 0) tmp32no1 = (((Word32)magnU16[i] - tmp32no2)*ONE_MINUS_GAMMA_PAUSE_Q8 + 128)>>8; // Q(qMagn)
	          else             tmp32no1 = ((((Word32)magnU16[i]<<nShifts) - st->avgMagnPause[i])*ONE_MINUS_GAMMA_PAUSE_Q8 + (128 << nShifts))>>(8 + nShifts); // Q(qMagn)
              tmp32no2 += tmp32no1; // Q(qMagn)
          }
          st->avgMagnPause[i] = tmp32no2;
        }  // end of frequency loop
        
        norm32no1 = _NormU32(maxNoiseU32);
        qNoise = st->prevQNoise + norm32no1 - 5;
        // done with step 2: noise update
        
        // STEP 3: compute dd update of prior snr and post snr based on new noise estimate
        nShifts = st->prevQNoise + 11 - qMagn;
        for (i = 0; i < HALF_ANAL_BLOCKL; i++)
        {
            // // post and prior SNR
            // curNearSnr = 0.0;
            // if (magn[i] > noise[i]) curNearSnr = magn[i] / (noise[i] + 0.0001) - 1.0;
            // // DD estimate is sum of two terms: current estimate and previous estimate
            // // directed decision update of snrPrior
            // snrPrior = DD_PR_SNR * prevNearSnr[i] + (1.0 - DD_PR_SNR) * curNearSnr;
            // // gain filter
            // tmpFloat1 = st->overdrive + snrPrior;
            // theFilter[i] = tmpFloat2 = snrPrior / tmpFloat1;
        
            // calculate curNearSnr again, this is necessary because a new noise estimate has been made since then. for the original
            curNearSnr = 0; // Q11
            if (nShifts < 0)
	        {
                // This case is equivalent with magn < noise which implies curNearSnr = 0;
                tmpMagnU32 = (UWord32)magnU16[i]; // Q(qMagn)
                tmpNoiseU32 = (UWord32)noiseU32[i]<<(-nShifts); // Q(qMagn)
            } 
	        else if (nShifts > 17)
	        {
                tmpMagnU32 = (UWord32)magnU16[i]<<17; // Q(qMagn+17)
                tmpNoiseU32 = noiseU32[i]>>(nShifts-17); // Q(qMagn+17)
            } 
	        else
	        {
                tmpMagnU32 = (UWord32)magnU16[i]<<nShifts; // Q(qNoise_prev+11)
                tmpNoiseU32 = noiseU32[i]; // Q(qNoise_prev+11)
            }
            if (tmpMagnU32 > tmpNoiseU32)
	        {
                tmpU32no1 = tmpMagnU32 - tmpNoiseU32; // Q(qCur)
                norm32no2 = min(11, _NormU32(tmpU32no1));
                tmpU32no1 = (UWord32)tmpU32no1<<norm32no2; // Q(qCur+norm32no2)
                tmpU32no2 = tmpNoiseU32>>(11 - norm32no2); // Q(qCur+norm32no2-11)
                if (tmpU32no2 > 0) tmpU32no1 = tmpU32no1/tmpU32no2; // Q11
                curNearSnr = min(satMax, tmpU32no1); // Q11
            }
        
            //directed decision update of priorSnr    //priorSnr = DD_PR_SNR * prevNearSnr + (1.0-DD_PR_SNR) * curNearSnr;
            priorSnr = (UWord32)prevNearSnr[i]*DD_PR_SNR_Q11 + (UWord32)curNearSnr*ONE_MINUS_DD_PR_SNR_Q11; // Q22
            
            //gain filter
            tmpU32no1 = (UWord32)(st->overdrive) + ((priorSnr + 8192)>>14); // Q8
            
	        assert(st->overdrive > 0);
            tmpU16no1 = (UWord16)((priorSnr + (tmpU32no1 >> 1))/tmpU32no1); // Q14
            st->noiseSupFilter[i] = min(max(tmpU16no1, st->denoiseBound), 16384);
            
            // Weight in the parametric Wiener filter during startup
            if (st->blockIndex < END_STARTUP_SHORT)
	        {
                // Weight the two suppression filters
                tmpU32no1 = (UWord32)(UWord16)st->noiseSupFilter[i]*(UWord16)st->blockIndex;
                tmpU32no2 = (UWord32)(UWord16)noiseSupFilterTmp[i]*(UWord16)(END_STARTUP_SHORT - st->blockIndex);
                tmpU32no1 += tmpU32no2;
                st->noiseSupFilter[i] = (UWord16)(tmpU32no1/END_STARTUP_SHORT);
            }
        }  // end of loop over frequencies
        //done with step3
        
        // save noise and magnitude spectrum for next frame
        st->prevQNoise = qNoise;
        st->prevQMagn = qMagn;
        if (norm32no1 > 5) {for (i = 0; i < HALF_ANAL_BLOCKL; i++) st->prevNoiseU32[i] = noiseU32[i]<<(norm32no1 - 5);} // Q(qNoise+11)
        else               {for (i = 0; i < HALF_ANAL_BLOCKL; i++) st->prevNoiseU32[i] = noiseU32[i]>>(5 - norm32no1);} // Q(qNoise+11)
	    for (i=0; i<HALF_ANAL_BLOCKL; i++) st->prevMagnU16[i] = magnU16[i]; // Q(qMagn)
    }
	NS_DataSynthesis(st, outFrame);
}//line=1619, 1491, 1472, 