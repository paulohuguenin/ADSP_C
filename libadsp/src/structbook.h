//structbook.h
#ifndef STRUCTBOOK_H
#define STRUCTBOOK_H

#include <string.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sstream>

#include "defines.h"
#include "bitstream.h"
#include "parameter.h"
#include "structs.h"
#include "filemgr.h"
#include "datasignal.h"
#include "cgmatrix.h"

using namespace std;


// ============================
//  CStructBook
// ============================

/**
 * \brief Structure book class
 */

class CStructBook
{

protected:
	int		m_iNumElement;
    double  m_norm;

public:

	CStructBook();

	virtual ~CStructBook();

	// Interface

	// Set methods
	void    setNumElement(int numElement);
    void    setNorm(double norm);

	// Get methods
	int     getNumElement();
    double     getNorm();



};

// ============================
//  CStructBookExp
// ============================

/**
 * \brief Exponential structure book class
 */

class CStructBookExp : public CStructBook
{

    strtContinuousExp* m_structBook;
    strtSBQHeadExp sbQHeadExp;
    int m_numElemOpCurve;
    strtOpCurveExp* m_opCurve;
    double m_minamp;
    double m_maxamp;
    double m_minrho;
    double m_maxrho;
    double m_minphase;
    double m_maxphase;

public:

	CStructBookExp();
	virtual ~CStructBookExp();

    //===========================================
	// Intrinsic
    void addElement(strtParameter* parm);

    void addElement(strtContinuousExp cExpStructure);
    
    void addElement ( strtContinuousExp* cExpStructure, int numEl );

    void setNextAtomIndex ( int atomIndex,
                            int nextAtom);

    void setPrevAtomIndex ( int atomIndex,
                            int prevAtom);

    void setOrigAtomIndex ( int atomIndex,
                            int origAtom);

    void removeElement(int atomIndex);

	void saveElementASCII(	FILE* stream,
							double meanApproxRatio,
							double approxRatio,
							double befSupInnerP,
							double aftSupInnerP,
                            double normRatio);

    void saveElementBin( FILE* stream);
    
    void saveElementBin( FILE* stream, int iNumElement, strtContinuousExp sb);

	void loadElementASCII(FILE* stream);

	void loadElementASCII(char* str);

	void adjustStructBook(int L);

	void saveStructBook(char* file);

	void saveStructBook(FILE* stream);

	void saveStructBookASCII(FILE* stream);

	// ... respective load functions to be written
	void loadStructBook(FILE* stream);

    void printToScreen();

    void printElementToScreen(int index);
    
    //=============================================
    void convertCoefNegativeToPositive();
    
    void findParmRange( double& min_coef,
                        double& max_coef,
                        double& min_rho,
                        double& max_rho,
                        double& min_phase,
                        double& max_phase);
                        
    int findDeltaSupMax();
    
    void sepByAmp(	strtContinuousExp* pSB,
                    int numElement,
                    double lowerAmpRangeLimit,
                    double upperAmpRangeLimit);
     
    void sepBySubBlock(	strtContinuousExp* pSB,
                        int numElement,
                        int lowerSubBlockLimit,
                        int upperSubBlockLimit);

    //===========================================
    // Quantization methods
    void setQuantConfig(    strtContinuousExp* pSB,
                            int sbNumElement,
                            double norm,
                            int nbits_amp,
						    int nbits_rho,
						    int nbits_phase,
						    double freq1, // (Ffund/ Finit) (electric/audio)
						    double freq2, // (Fs/ Fend) (electric/audio)
						    int signalSize,
                            int sigType);

    void quantStructBook(   strtContinuousExp* pSB,
                            int numElement);

    int computeNumBits();

    double computeRate(int numSamples);
    
    void setOpCurve(strtOpCurveExp* opCurve,int numElement);
    
    void setMinAmp(double);
    void setMaxAmp(double);
    void setMinRho(double);
    void setMaxRho(double);
    void setMinPhase(double);
    void setMaxPhase(double);
    

    //===========================================
    // Get methods
	strtContinuousExp*	getStructBook() const;
	strtOpCurveExp*	    getOpCurve() const;
    int getNextAtomIndex ( int atomIndex);
    int getPrevAtomIndex ( int atomIndex);
    int getOrigAtomIndex ( int atomIndex);
    int getNumElemOpCurve();
    double getMinAmp();
    double getMaxAmp();
    double getMinRho();
    double getMaxRho();
    double getMinPhase();
    double getMaxPhase();

};


// ========================
// CStructBookQExp
// ========================
/*
class CStructBookQExp : public CStructBookQ
{
	strtStructBookQHeader	structBookQuantizedHeader;
	strtStructBookQ*		m_structBookQ;

public:

	CStructBookQExp();

	~CStructBookQExp();

	void setStructBookQ(	strtStructBookQHeader sBookQHeader,
							strtStructBookQ* structBookQ,
							int iNumElementQ);
	
	void setStructBookQ();

	void saveStructBookQHeader(FILE* stream);
	
	void saveStructBookQHeader(CBitStream& bitStream);

	void saveStructBookQ(	CBitStream& bitStream,
                            double Ffund,
                            double Fs,
                            int signalSize);

	// ... respective load functions to be written

	void loadStructBookQHeader(FILE* stream);
	
	void loadStructBookQHeader(CBitStream& bitStream);

	void loadStructBookQ(   CBitStream& bitStream,
                            double Ffund,
                            double Fs,
                            int signalSize);
	
	float quantizeFloat(float input, int int_nbits, int dec_nbits);
	
	double quantizeDouble(double input, int int_nbits, int dec_nbits);
	
	void quantizeStructBookQHeader();

	void quantizeStructureBook(	CStructBook* structBook,
					int nbits_amp,
					int nbits_rho,
					int nbits_phase,
					double Ffund,
					double Fs,
					int signalSize);

	void quantizeStructureBook(	CStructBook* structBook,
					int nbits_amp,
					int nbits_rho,
					int nbits_phase,
					double Ffund,
					double Fs,
					int signalSize,
					int dummy);
	
	void dequantizeStructureBook(	CStructBook* structBook, 
					double Ffund,
					double Fs,
					int signalSize); //sqrd_amp
	
	void dequantizeStructureBook(	CStructBook* structBook, 
					double Ffund,
					double Fs,
					int signalSize,
					int dummy); // linear amp

    void			calcNumBits(	double Ffund,
                                double Fs,
                                int signalSize);

	strtStructBookQHeader	getSBQHeader() const;
	strtStructBookQ*	getStructBookQ() const;
	double			getNorm() const;

	const CStructBookQ& operator=(const CStructBookQ& structBookQ);

};
*/

// Functions that manipulate the classes

void separateByAmpExp(CStructBook* sbIn,CStructBook* sbOut, int NAmpRange, double* AmpRangeLimit );


//edit

class CStructBookParm : public CStructBook
{
    strtParameter* m_structBook;
    strtSBQHeadParm sbQHeadParm;
    int m_numElemOpCurve;
    strtOpCurveParm* m_opCurve;
    double m_minamp;
    double m_maxamp;
    double m_minrho;
    double m_maxrho;
    double m_minphase;
    double m_maxphase;

public:

	CStructBookParm();
	virtual ~CStructBookParm();

    //===========================================
	// Intrinsic
    void addElement(strtParameter* parm);

    void addElement(strtParameter cParmStructure);
    
    void addElement ( strtParameter* cParmStructure, int numEl );

    void setNextAtomIndex ( int atomIndex,
                            int nextAtom);

    void setPrevAtomIndex ( int atomIndex,
                            int prevAtom);

    void setOrigAtomIndex ( int atomIndex,
                            int origAtom);

    void removeElement(int atomIndex);

    // void newFileASCII(  int initBlock,
    //                     int finalBlock,
    //                     CDataSignal* dataSignal,
    //                     CFileDictionary* dicData);

    void saveMainHeaderASCII (  char* fileName,
                                // FILE* stream,
                                int initBlock,
                                int finalBlock,
                                CDataSignal* dataSignal);

    void saveSignalHeaderASCII (char* fileName,
                                // FILE* stream,
                                int i,
                                CDataSignal* dataSignal);

    void saveBlockHeaderASCII ( char* fileName,
                                // FILE* stream, 
                                int j,
                                double initBlockNorm);

	void saveElementASCII(	char* fileName,
                            // FILE* stream,
                            double meanApproxRatio,
							double approxRatio,
							double befSupInnerP,
							double aftSupInnerP,
                            double normRatio,
                            int chosenNet);

    void saveInnerProdASCII (   char* fileName,
                                char* fileName2,
                                cgMatrix<double>& a,
                                int step,
                                int block,
                                int initBlock,
                                long& pos);

    void saveBlockNullSamplesASCII( char* fileName/*,
                                    FILE* stream*/);

    void saveBlockEndingASCII(  char* fileName/*,
                                FILE* stream*/);

    void saveSignalEndingASCII( char* fileName/*,
                                FILE* stream*/);

    void saveSignalNullSamplesASCII(char* fileName/*,
                                    FILE* stream*/);

    void saveDecompEndingASCII( char* fileName/*,
                                FILE* stream*/);


    void saveMainHeaderBin (char* sbbFName,
                            // FILE* stream,
                            int initBlock,
                            int finalBlock,
                            CDataSignal* dataSignal,
                            CFileDictionary* dicData,
                            CFileDecomp* genData);
    
    void newFileBin(char* sbbFName,
                    // FILE* stream,
                    int initBlock,
                    int finalBlock,
                    CDataSignal* dataSignal);

    void saveSignalHeaderBin (  char* sbbFName,
                                // FILE* stream,
                                int i,
                                CDataSignal* dataSignal);

    void saveBlockHeaderBin (   char* sbbFName,
                                // FILE* stream,
                                int j,
                                double initBlockNorm);

    void saveElementBin(char* sbbFName/*,
                        FILE* stream*/);
    
    void saveElementBin(FILE* stream,
                        int iNumElement,
                        strtParameter sb);

    void saveBlockEndingBin(char* sbbFName/*,
                            FILE* stream*/);

    void saveSignalEndingBin(   char* sbbFName/*,
                                FILE* stream*/);

    void saveDecompEndingBin(   char* sbbFName/*,
                                FILE* stream*/);

	void loadElementASCII(FILE* stream);

	void loadElementASCII(char* str);

	void adjustStructBook(int L);

	void saveStructBook(char* file);

	void saveStructBook(FILE* stream);

	void saveStructBookASCII(FILE* stream);

	// ... respective load functions to be written
	void loadStructBook(FILE* stream);

    void printToScreen();

    void printElementToScreen(int index);
    
    //=============================================
    void convertCoefNegativeToPositive();
    
    void findParmRange( double& min_coef,
                        double& max_coef,
                        double& min_rho,
                        double& max_rho,
                        double& min_phase,
                        double& max_phase);
                        
    int findDeltaSupMax();
    
    void sepByAmp(	strtParameter* pSB,
                    int numElement,
                    double lowerAmpRangeLimit,
                    double upperAmpRangeLimit);
     
    void sepBySubBlock(	strtParameter* pSB,
                        int numElement,
                        int lowerSubBlockLimit,
                        int upperSubBlockLimit);

    //===========================================
    // Quantization methods
    void setQuantConfig(    strtParameter* pSB,
                            int sbNumElement,
                            double norm,
                            int nbits_amp,
						    int nbits_rho,
						    int nbits_phase,
						    double freq1, // (Ffund/ Finit) (electric/audio)
						    double freq2, // (Fs/ Fend) (electric/audio)
						    int signalSize,
                            int sigType);

    void quantStructBook(   strtParameter* pSB,
                            int numElement);

    int computeNumBits();

    double computeRate(int numSamples);
    
    void setOpCurve(strtOpCurveParm* opCurve,int numElement);
    
    void setMinAmp(double);
    void setMaxAmp(double);
    void setMinRho(double);
    void setMaxRho(double);
    void setMinPhase(double);
    void setMaxPhase(double);
    

    //===========================================
    // Get methods
	strtParameter*	getStructBook() const;
	strtOpCurveParm*	    getOpCurve() const;
    int getNextAtomIndex ( int atomIndex);
    int getPrevAtomIndex ( int atomIndex);
    int getOrigAtomIndex ( int atomIndex);
    int getNumElemOpCurve();
    double getMinAmp();
    double getMaxAmp();
    double getMinRho();
    double getMaxRho();
    double getMinPhase();
    double getMaxPhase();

    void saveHeader (char* fileName);
    void saveElement ( char* fileName,
                    double innerProd_xp,
                    double innerProd_xq,
                    double innerProd_pp,
                    double innerProd_qq,
                    double innerProd_pq);
};

void separateByAmp(CStructBook* sbIn,CStructBook* sbOut, int NAmpRange, double* AmpRangeLimit );



#endif
