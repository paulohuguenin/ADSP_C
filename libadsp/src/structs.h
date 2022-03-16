#ifndef STRUCTS_H
#define STRUCTS_H

/**
 * \brief Struct of exponential atom parameters
 */

typedef struct strtContinuousExp_tag
{
	double innerProduct;
	double rho;
	double xi;
	double phase;
	int a;
	int b;
    int nextAtom;
    int prevAtom;
    int origAtomIndex;
} strtContinuousExp;


/**
 * \brief Struct of atom parameters
 */

typedef struct strtParameter_tag
{
    double innerProduct;
    double s; // Gabor
    double rho;
    double xi;
    double phase;
    int u; // Gabor
    int a;
    int b;
    int nextAtom;
    int prevAtom;
    int origAtomIndex;
    double beta;
    int dicType;
} strtParameter;

/**
 * \brief Struct of main header parameters
 */

typedef struct strtSBBHeader_tag
{
    int     type;
    // int     dicType;
    int     numSignal;
    int     signalSize;
    int     subBlockSize;
    int     blockSize;
    int     blockHop;
    int     numBlock;
    double* norm;
    double  Fs;
    int     initBlock;
    int     finalBlock;
} strtSBBHeader;

/**
 *
 */

typedef struct strtStructBookQ_tag
{
	int innerProduct;
	int rhoSignal;
	int rho;
	//int u;
	int xi;
	int phase;
	int a;
	int b;
} strtStructBookQ;

/**
 *
 */

typedef struct strtOscHeader_tag
{
	int	F;
	int	Fs;
	int	signalSize;
	int	nChannel;
} strtOscHeader;

/**
 *
 */

typedef struct strtStructBookQHeader_tag
{
	double		norm;
	int		numberStruct;
	double		max_amp;
	int		nbits_amp;
	double		max_rho;
	double		min_rho;
	int		nbits_rho;
	double		max_phase;
	double		min_phase;
	int		nbits_phase;

} strtStructBookQHeader;

/**
 *
 */

typedef struct strtSBQHeadExp_tag
{
	double		norm;
	int		    numberStruct;
    /////////////////////////
    double		min_amp;
	double		max_amp;
	int		    nbits_amp;
    /////////////////////////
    double		min_rho;
	double		max_rho;
	int		    nbits_rho;
    /////////////////////////
    double		min_xi;
	double		max_xi;
	int		    nbits_xi;
    /////////////////////////
	double		min_phase;
	double		max_phase;
    int		    nbits_phase;
	/////////////////////////
	int		    nbits_a;
    /////////////////////////
    int		    nbits_b;

} strtSBQHeadExp;


typedef struct strtLinkSBQHeadExp_tag
{
    double      norm;
    int         numberStruct;
    /////////////////////////
    double      min_amp;
    double      max_amp;
    int         nbits_amp;
    /////////////////////////
    double      min_rho;
    double      max_rho;
    int         nbits_rho;
    /////////////////////////
    int         nbits_xi;
    /////////////////////////
    double      min_phase;
    double      max_phase;
    int         nbits_phase;
    /////////////////////////
    int         nbits_sample;
    //////////////////////////
    int         nbits_block;
} strtLinkSBQHeadExp;


//edit

typedef struct strtSBQHeadParm_tag
{
    double      norm;
    int         numberStruct;
    /////////////////////////
    double      min_amp;
    double      max_amp;
    int         nbits_amp;
    /////////////////////////
    double      min_rho;
    double      max_rho;
    int         nbits_rho;
    /////////////////////////
    double      min_xi;
    double      max_xi;
    int         nbits_xi;
    /////////////////////////
    double      min_phase;
    double      max_phase;
    int         nbits_phase;
    /////////////////////////
    int         nbits_a;
    /////////////////////////
    int         nbits_b;

} strtSBQHeadParm;


typedef struct strtLinkSBQHeadParm_tag
{
    double      norm;
    int         numberStruct;
    /////////////////////////
    double      min_amp;
    double      max_amp;
    int         nbits_amp;
    /////////////////////////
    double      min_rho;
    double      max_rho;
    int         nbits_rho;
    /////////////////////////
    int         nbits_xi;
    /////////////////////////
    double      min_phase;
    double      max_phase;
    int         nbits_phase;
    /////////////////////////
    int         nbits_sample;
    //////////////////////////
    int         nbits_block;
} strtLinkSBQHeadParm;


/**
 *
 */

// Structs used for RD optimization

typedef struct strtQuantExp_tag
{
    int nb_amp;
    int nb_rho;
    int nb_phase;
    int nb_xi;
    int nb_sample;
    int nb_deltasup;
    int nb_block;
    int i_qrho;
    double qrho;
    double R_rho;
    double H_rho;
    int i_qphi;
    double qphi;
    double R_phi;
    double H_phi;
    double R_xi;
    double R_sample;
    double R_amp;
    double H_amp;
    double R_deltasup;
    double H_deltasup;
} strtQuantExp;


//edit

typedef struct strtQuantParm_tag
{
    int nb_amp;
    int nb_rho;
    int nb_phase;
    int nb_xi;
    int nb_sample;
    int nb_deltasup;
    int nb_block;
    int i_qrho;
    double qrho;
    double R_rho;
    double H_rho;
    int i_qphi;
    double qphi;
    double R_phi;
    double H_phi;
    double R_xi;
    double R_sample;
    double R_amp;
    double H_amp;
    double R_deltasup;
    double H_deltasup;
} strtQuantParm;

/**
 *
 */

typedef struct strtRD_tag
{
    double rate;
    double dist;
} strtRD;

/**
 *
 */

typedef struct strtOpCurveExp_tag
{
    int nb_amp;
    int nb_rho;
    int nb_phase;
    int nb_xi;
    int nb_sample;
    int nb_deltasup;
    int nb_block;
    double rate;
    double dist;
    double lambda;
} strtOpCurveExp;


//edit

typedef struct strtOpCurveParm_tag
{
    int nb_amp;
    int nb_rho;
    int nb_phase;
    int nb_xi;
    int nb_sample;
    int nb_deltasup;
    int nb_block;
    double rate;
    double dist;
    double lambda;
} strtOpCurveParm;

/**
 * \brief Struct of dictionary input parameters
 */

typedef struct strtFileDictionary_tag
{
    double scale;
    int delta_tau;
    int fdiscrtype;
    double freqi;
    double freqf;
    int dicType;
    double decay;
    double rise;
} strtFileDictionary;

/**
 * \brief Struct of linked exponential atom parameters
 */

typedef struct strtLinkExp_tag
{
    double coefFirst;
    double coefSum;
    double coefSqrSum;
    double* pRho;
    int* pRhoSign;
    double xi;
    double phase;
    int initSample;
    int endSample;
    int initBlock;
    int endBlock;
    int numBlock;
} strtLinkExp;


//edit

typedef struct strtLinkParm_tag
{
    double coefFirst;
    double coefSum;
    double coefSqrSum;
    double* pRho;
    int* pRhoSign;
    double xi;
    double phase;
    int initSample;
    int endSample;
    int initBlock;
    int endBlock;
    int numBlock;
} strtLinkParm;

/**
 * \brief Struct of referenced linked atom parameters
 */

typedef struct strtRefLinkSB_tag
{
    int numBlock;
    double coefSum;
    double coefSqrSum;
    int* block;
    int* sbElement;
    int* rhoSign;
} strtRefLinkSB;


/**
 * \brief Struct of rate-distortion indices
 */

typedef struct strtRDIndex_tag
{
	double*	theta_vec;
	int*	index_vec;
	int		numElement;
} strtRDIndex;

/**
 * \brief
 */

typedef struct strtRDByAmpQuant_tag
{
	int** 			numAmpRange; // [signal,block]
	strtQuantExp*** quantExp;	 // [signal,block,ampRange]
	double***		rate;        // [signal,block,ampRange]
	int***		    nAtom;       	// [signal,block,ampRange]
	double			totalRate;
	double 			lambda;
} strtRDByAmpQuant;

typedef struct strtOpCurveQuantExp_tag
{
	int numElement;
    strtQuantExp* quantExp;
    double *rate;
    double *dist;
    double *lambda;
} strtOpCurveQuantExp;


//edit

typedef struct strtOpCurveQuantParm_tag
{
    int numElement;
    strtQuantExp* quantExp;
    double *rate;
    double *dist;
    double *lambda;
} strtOpCurveQuantParm;

#endif
