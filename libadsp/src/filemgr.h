// filemgr.h
#ifndef FILEMGR_H
#define FILEMGR_H

#include "defines.h"
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "structs.h"
#include <string.h>

using namespace std;
//=======================
// CFileGenData Class
//=======================

/**
 * \brief General data input file class
 */

class CFileGenData {

    char m_fileName[_MAX_PATH];
    int m_sigType;
    int m_dicType;
    int m_procType;
    int m_flOMP;
    int m_flANN;
    int m_blockHop;
    int m_blockSize;
    int m_initBlock;
    int m_endBlock;
    int m_nMaxStep;
    double m_coefTempSup;
    int m_flSeqBlock;
    double m_coefSeqBlock;
    int m_flEvalAtomCont;
    int m_flFindSupport;
    int m_flOptDecay;
    double m_approxRatioTarget;
    double m_snrTarget;
    int m_flImpStage;

public:
    // Constructors
    CFileGenData ();
  
    // Destructor
    ~CFileGenData(){};

    void setFileName(char* file);
    void loadData();


    // Get methods
    int getSigType();
    int getDicType();
    int getProcType();
    int getFlagOMP();
    int getFlagANN();
    int getBlockHop();
    int getBlockSize();
    int getInitBlock();
    int getEndBlock();
    int getNumMaxStep();
    double getCoefTempSup();
    int getFlagSeqBlock();
    double getCoefSeqBlock();
    int getFlagEvalAtomCont();
    int getFlagFindSupport();
    int getFlagOptDecay();
    double getApproxRatioTarget();
    double getSNRTarget();
    int getPrintDecompStage();
};



/**
 * \brief Proceedings input file class
 */

class CFileProceeding {

    char m_fileName[_MAX_PATH];
    int m_flDecomp;
    int m_flGroupBlock;
    int m_flRDTableDf;
    int m_flRDOpCurve;
    int m_flRDOpCurveAmpRange;
    int m_flEncAmpRange;
    int m_flEncOpCurve;
    int m_flEncAmpRangeArithCod;
    int m_flDecAmpRangeArithCod;
    int m_flPll;

public:
    // Constructors
    CFileProceeding ();
  
    // Destructor
    ~CFileProceeding(){};

    void setFileName(char* file);
    void loadData();

    // Get methods
    int getDecomp();
    int getGroupBlock();
    int getRDTableDf();
    int getRDOpCurve();
    int getRDOpCurveAmpRange();
    int getEncAmpRange();
    int getEncOpCurve();
    int getEncAmpRangeArithCod();
    int getDecAmpRangeArithCod();
    int getPll();
};


class CFileGroupBlock {

    char m_fileName[_MAX_PATH];
    int m_numBlockPerGroup;

public:
    // Constructors
    CFileGroupBlock ();
  
    // Destructor
    ~CFileGroupBlock(){};

    void setFileName(char* file);
    void loadData();

    // Get methods
    int getNumBlockPerGroup();
};





//=======================
// CFileBlockRange Class
//=======================

/**
 * \brief Black range input file class
 */

class CFileBlockRange {

    char m_fileName[_MAX_PATH];
    int m_numRange;
    int* m_pInitBlock;
    int* m_pFinalBlock;

public:

    // Constructors
    CFileBlockRange ();
  
    // Destructor
    ~CFileBlockRange();

    void setFileName(char* file);
    void loadData();

    int getNumRange();
    int* getPInitBlock();
    int* getPFinalBlock();
};

//=======================
// CFileDictionary Class
//=======================

/**
 * \brief Dictionary parameters input file class
 */

class CFileDictionary {

    char m_fileName[_MAX_PATH];
    int m_numDicBlock;
    strtFileDictionary* m_dicBlock;
    int m_numFreq;

public:

    // Constructors
    CFileDictionary ();
  
    // Destructor
    ~CFileDictionary();

    void setFileName(char* file);
    void loadData();
    void printToScreen();

    int getDicType(int blockInd);
    int getNumDicBlock();
    strtFileDictionary* getDicBlock();
    int getNumFreq();
    double getFreqi(int blockInd);
    double getFreqf(int blockInd);
    double getFDiscrType(int blockInd);
    double getDecay(int blockInd);
    double getRise(int blockInd);



};

//=======================
// CFileRDBitRange Class
//=======================

/**
 * \brief Bit range input file class
 */

class CFileRDBitRange {

    char m_fileName[_MAX_PATH];
    int init_nbit_amp;
    int end_nbit_amp;
    int delta_nbit_amp;
    int init_nbit_rho;
    int end_nbit_rho;
    int delta_nbit_rho;
    int init_nbit_phase;
    int end_nbit_phase;
    int delta_nbit_phase;
    char m_DfTableFName[_MAX_PATH];
public:

    // Constructors
    CFileRDBitRange ();
  
    // Destructor
    ~CFileRDBitRange(){};

    void setFileName(char* file);
    void loadData();
    int getInitNbitAmp();
    int getEndNbitAmp();
    int getDeltaNbitAmp();
    int getInitNbitRho();
    int getEndNbitRho();
    int getDeltaNbitRho();
    int getInitNbitPhase();
    int getEndNbitPhase();
    int getDeltaNbitPhase();
    char * getDfTableFName();

};

class CFileDfTable {

    char m_fileName[_MAX_PATH];
    int m_Nqrho;
    int m_Nqphi;
    double* m_qrho;
    double* m_qphi;
    double**** m_DfTable;
    int m_min_nbit_amp;
    int m_max_nbit_amp;

public:
    // Constructors
    CFileDfTable ();
  
    // Destructor
    ~CFileDfTable();

    void setFileName(char* file);
    void loadData();
    
    // Get methods
    int getBaseBlockSize();
    int getNumQRho();
    int getNumQPhi();
    double* getQRho();
    double* getQPhi();
    double**** getDfTable();
    int getMinNbitAmp();
    int getMaxNbitAmp();
};


//=======================================
// Input files of decomposition program
//=======================================

/**
 * \brief Decomp input file class
 */


class CFileDecomp {

    char m_fileName[_MAX_PATH];
    int m_sigType;
    int m_dicType;
    int m_mpType;
    int m_flOMP;
    int m_flANN;
    int m_nMaxStep;
    double m_coefTempSup;
    int m_flSeqBlock;
    double m_coefSeqBlock;
    int m_flEvalAtomCont;
    int m_flFindSupport;
    int m_flOptDecay;
    double m_approxRatioTarget;
    double m_snrTarget;
    int m_flImpStage;
    int m_flRDopt;

public:
    // Constructors
    CFileDecomp ();
  
    // Destructor
    ~CFileDecomp(){};

    void setFileName(char* file);
    void loadData();
    void echo();


    // Get methods
    int getSigType();
    int getDicType();
    int getMPType();
    int getFlagOMP();
    int getFlagANN();
    int getNumMaxStep();
    double getCoefTempSup();
    int getFlagSeqBlock();
    double getCoefSeqBlock();
    int getFlagEvalAtomCont();
    int getFlagFindSupport();
    int getFlagOptDecay();
    double getApproxRatioTarget();
    double getSNRTarget();
    int getPrintDecompStage();
    int getFlagRDopt();
};

/**
 * \brief Decomp block range input file class
 */

class CFileDecompBlockRange {

    char m_fileName[_MAX_PATH];
    int m_blockHop;
    int m_blockSize;
    int m_initBlock;
    int m_endBlock;
    int m_numRange;
    int* m_pInitBlock;
    int* m_pFinalBlock;

public:

    // Constructors
    CFileDecompBlockRange ();
  
    // Destructor
    ~CFileDecompBlockRange();

    void setFileName(char* file);
    void loadData();

    int getBlockHop();
    int getBlockSize();
    int getInitBlock();
    int getEndBlock();
    //
    int getNumRange();
    int* getPInitBlock();
    int* getPFinalBlock();
    
};


//=======================================
// Input files of encode program
//=======================================


/**
 * \brief Encode input file class
 */

class CFileEncode {

    char m_fileName[_MAX_PATH];
    int m_flGroupBlock;
    int m_numBlockPerGroup;
    int m_flRDTableDf;
    int m_flRDOpCurve;
    int m_flRDOpCurveAmpRange;
    int m_flEncAmpRange;
    int m_flEncOpCurve;
    int m_flEncAmpRangeArithCod;
    int m_flDecAmpRangeArithCod;

public:
    // Constructors
    CFileEncode ();
  
    // Destructor
    ~CFileEncode(){};

    void setFileName(char* file);
    void loadData();

    // Get methods
    int getGroupBlock();
    int getNumBlockPerGroup();
    int getRDTableDf();
    int getRDOpCurve();
    int getRDOpCurveAmpRange();
    int getEncAmpRange();
    int getEncOpCurve();
    int getEncAmpRangeArithCod();
    int getDecAmpRangeArithCod();
};


//=======================
// CFileNetwork Class
//=======================

/**
 * \brief Networl parameters input file class
 */

class CFileNetwork {

    char m_fileName[_MAX_PATH];

    double* m_x1_xoffset;
    double* m_x1_gain;
    double m_x1_ymin;
    double* m_b1;
    double ** m_w1;
    double* m_b2;
    double ** m_w2;

    int m_blockSize;
    int m_chosenNet;
    int m_size_x;
    int m_size_h;
    int m_size_y;

public:
    // Constructors
    CFileNetwork();
  
    // Destructor
    ~CFileNetwork();

    void    setFileName(char* file);
    void    setBlockSize(int blockSize);
    void    setNetwork(int num);
    void    loadData();

    double*     getX1_xoffset();
    double*     getX1_gain();
    double      getX1_ymin();
    double*     getB1();
    double**    getW1();
    double*     getB2();
    double**    getW2();

};


#endif
