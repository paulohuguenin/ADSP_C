#ifndef DEFINES_H
#define DEFINES_H

//#define WIN

#define NUM_MAX_STEP 1000

// Debugs each decomposition step
//#define DBG_WATCH_DCMP_STEP

#define _MAX_PATH 260

// Load complex dictionary in memory
#define LOAD_DIC 0

// For frequency quantization multiplies Fs/Ffund
#define RF_COEF 1

// Chooses squared amp quantization
// if not defined chooses linear amp quantization
//#define SQRD_AMP_QUANT

// Enable correct compression ratio calculation
// if not defined enable real compression ratio calculation
//#define CORRECT_COMP_RATIO

#define PI acos(-1.0)

// structure book Header Bit Allocation
#define NORM_INT_NBITS 25
#define NORM_DEC_NBITS 7
#define NUM_STRUCT_NBITS 10
#define AMP_INT_NBITS 1
#define AMP_DEC_NBITS 10
#define RHO_INT_NBITS 3
#define RHO_DEC_NBITS 10
#define PHASE_INT_NBITS 3
#define PHASE_DEC_NBITS 10
	// The number of bits of the number of bits of amp, rho and phase
	// Quantization Vector
#define AMP_2NBITS 5
#define RHO_2NBITS 5
#define PHASE_2NBITS 5

//#define TOTAL_NBITS_SB_HEADER 10

// oscillogram Header Bit Allocation

#define F_NBITS 1
#define FS_NBITS 20
#define SIG_SIZE_NBITS 12
#define NUM_CHANNEL_NBITS 9

//#define TOTAL_NBITS_OSC_HEADER 10

#define NUM_QUANT 292

//#define GABEXP
//#define NOISE_LENGTH 128
//#define NOISE
//#define COMTRADE_DECOMP
//#define COMTRADE_QUANT_FILE
//#define RD_STUDY
//#define COMTRADE_QUANT_CHANNEL
//#define MSE_X_NUMSTRUCT
//#define JASP
//#define GERA_ARQ_ANAUTO
//#define CUT_COMTRADE


#endif
