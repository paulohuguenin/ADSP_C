// engine.h
#ifndef ENCODE_H
#define ENCODE_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <istream>

#include "filemgr.h"
#include "structbook.h"
#include "linkstrbook.h"
#include "dictionary.h"
#include "datasignal.h"
#include "loadfile.h"
#include "metric.h"
#include "adspdefines.h"
#include "adspsort.h"
#include "qstep.h"
#include "arithmetic_codec.h"



using namespace std;

void optimizeRDExp(	CFileRDBitRange* rdBitRaC++ Reference: A site with its main focus on a complete Standard Template Library reference, the Technical Specifications, and a list of selected non-ANSI/ISO libraries. An offline archive is available.
nge,
					CDataSignal* dataSignal,
					strtSBBHeader sbbHeader,
					CStructBook** structBook,
					int Nfreq,
					char* InputFile);

void optimizeRDExpAmpRangeOpCurve(	CFileRDBitRange* rdBitRange,
									CDataSignal* dataSignal,
									strtSBBHeader sbbHeader,
									CStructBook** structBook,
									int Nfreq,
									char* InputFile);

void quantizeStructBookExp(   double min_amp, double max_amp,
							  double min_rho,double max_rho,
							  double min_phase, double max_phase,
							  strtQuantExp quantExp,
							  strtContinuousExp* pSB, int numElement,
							  strtContinuousExp* pSBQ,int &numElementQ);

double quantizeStructBookExpDistAtom( double min_amp, double max_amp,
									  double min_rho,double max_rho,
									  double min_phase, double max_phase,
									  strtQuantExp quantExp,
									  strtContinuousExp* pSB, int numElement,
									  strtContinuousExp* pSBQ,int &numElementQ,
									  int sigSize, double signorm);

void quantizeStructBookExpRecUnquantNoReset( double min_amp, double max_amp,
									  double min_rho,double max_rho,
									  double min_phase, double max_phase,
									  strtQuantExp quantExp,
									  strtContinuousExp* pSB, int numElement,
									  strtContinuousExp* pSBQ,int &numElementQ,
									  double* recSignal, double signorm, int signalSize,
									  int* ampIndex, int* rhoIndex, int* xiIndex,
									  int* phiIndex, int* aIndex, int* bIndex,
									  int fdiscrtype, double step_xi);

int computeNumElementAmpQ(double min_amp, double max_amp,
						  int nb_amp,
						  strtContinuousExp* pSB, int numElement);

void quantizeStructBookExpAmp(  double min_amp, double max_amp,
							  	int nbamp,
							  	strtContinuousExp* pSB, int numElement,
							  	strtContinuousExp* pSBQ,int &numElementQ);

void dequantizeStructBookExpArithCodAmpVec(   int nb_amp, double min_amp, double max_amp,
											  unsigned* ampIndex, int numElement,
											  double* ampvec);

void dequantizeStructBookExpArithCod( int nb_amp,double min_amp, double max_amp,
									  int fdiscrtype, double step_xi,
									  double step_rho, double step_phase,
									  unsigned* ampIndex, int* rhoIndex, unsigned* xiIndex,
									  unsigned* phiIndex, unsigned* aIndex, unsigned* bIndex,
									  strtContinuousExp* pSB, int numElement);

double sumSqrErrorStructBookExpAmp( double min_amp, double max_amp,
									int nb_amp,
									strtContinuousExp* pSB, int numElement);

void synthSignalSBExp(   double* recSignal,
						 double norm,
						 int sigSize,
						 strtContinuousExp* sb,
						 int numElement);

void synthSignalSBExpNoReset( double* recSignal,
							 double norm,
							 int sigSize,
							 strtContinuousExp* sb,
							 int numElement);


double computeRateSBExp(strtQuantExp quantExp,int numElementQ);

double computeRateParameter(double nb,int numElementQ);

double computeRateSBExpAmp(double nb,int numElementQ);

double computeRateSBExpFreq(double nb,int numElementQ);

double computeRateSBExpSample(double nb,int numElementQ);

double computeRateUnifPdf(int nQStep,
						  int numElementQ);

double computeRateSBExpPhaseUnifPdf(int nQPhiStep,
									int numElementQ);

double computeRateSBExpDecayUnifPdf(int nQRhoStep,
									int numElementQ);

double computeRatePdf(	double* vec, int numElement,
						double* qstepEdge, double* qstepProb, int numQStep);

double computeRateSBExpDecayGGDPdf(	strtContinuousExp* sb, int numElement,
									double* qstepEdge, double* qstepProb, int numQStep);

double computeRateSBExpDecayGGDPdfPrint(	double* rhovec, int numElement,
									double* qstepEdge, double* qstepProb, int numQStep,
									int flag);

/////////////////////////////////////////////////////////////

double computeEntropyUnifPdf(int nQStep);

double computeEntropySBExpPhaseUnifPdf(int nQPhiStep);

double computeEntropySBExpDecayUnifPdf(int nQRhoStep);

double computeEntropyPdf( double* qstepProb, int numQStep);

double computeEntropySBExpDecayGGDPdf( double* qstepProb, int numQStep);



////////////////////////////////////////////////////////

void computeRDOpCurve(	int numQuant,
						strtRD* rd,
						strtRDIndex& rdIndex);

////////////////////////////////////////////////////////////////////////

void optimizeRDExpByAmp(CFileRDBitRange* rdBitRange,
						CDataSignal* dataSignal,
						strtSBBHeader sbbHeader,
						CStructBook** structBook,
						int Nfreq);

void calcTableRDExpByAmp(	CFileRDBitRange* rdBitRange,
							CDataSignal* dataSignal,
							strtSBBHeader sbbHeader,
							CStructBook** structBook,
							int Nfreq,
							char* InputFile);

double calcAtomDistFunction(strtContinuousExp sb,
							strtContinuousExp sbQ,
							int sigSize);


double indefIntegralAtomSqrNorm(double rho, double xi, double phase, double t);

double indefIntegralDfa(double rho, double xi, double phi, double t);

double indefIntegralDfb(double rho, double rhoq, double xi, double phi, double phiq, double t);

double indefIntegralDfc(double rho, double rhoq, double phi, double phiq, double t);

double computeMidThreadQuantStep( double nb, double min, double max);

double computeMidRiseQuantStep( double nb, double min, double max);

strtContinuousExp quantizeAtomExp(	double step_amp,double step_rho, double step_phase,
									double min_amp,double min_rho, double min_phase,
									strtContinuousExp pSB);

void synthSignalAtomExp( double* recSignal,
						 int sigSize,
						 strtContinuousExp sb);

////////////////////////////////////////////////////////////////////////

void groupBlockExp(	strtSBBHeader sbbHeader, CStructBook** structBook,
					strtSBBHeader sbbHeaderGroup, CStructBook** structBookGroup,
					int numBlockPerGroup, int iSignal);

////////////////////////////////////////////////////////////////////////

/**
 * \brief Encoding main function
 */

void encodeSBExp(	char*InputFile,
					double rateTarget,
					strtSBBHeader sbbHeader,
                	CStructBook** structBook,
                	int flFixedLambdaVSUnifBlockRate);

 /**
 * \brief Encoding with Dftable_byamp
 */

 void encodeSBExpAmpRange(	char*InputFile,
							double rateTarget,
							strtSBBHeader sbbHeader,
							CStructBook** structBook,
							CFileDfTable* DfTableFile,
							int Nfreq,
							CDataSignal* dataSignal);

void encodeSBExpAmpRangeArithCod(char*InputFile,
								double rateTarget,
								strtSBBHeader sbbHeader,
								CStructBook** structBook,
								CFileDfTable* DfTableFile,
								CFileDictionary* dicData,
								CDataSignal* dataSignal,
								double initLambda);

void  computeFullRDByAmpQuantArithCod(  strtSBBHeader sbbHeader,CStructBook** structBook,
										double** min_amp, double** max_amp,
										double** min_rho,double** max_rho,
										double** min_phase,double** max_phase,
										double**** DfTable,
										int N_qrho,int N_qphi,
										double* q_rho,double* q_phi,
										double nb_xi,double nb_sample,
										double*** distAmp,
										strtRDByAmpQuant* rdByAmpQuant,
										double**** ampQStepEdge,
										double**** ampQStepProb,
										double*** entropyAmpTable,
										double*** entropyPhiTable,
										double* deltaSupQStepEdge,
										double* deltaSupQStepProb,
										double H_deltasup,
										int deltaSupMax);

void encodeSBExpAmpRangeBisecSearch(char*InputFile,
									double rateTarget,
									strtSBBHeader sbbHeader,
									CStructBook** structBook,
									CFileDfTable* DfTableFile,
									int Nfreq,
									CDataSignal* dataSignal,
									int flLambdaBisecSearch,
									double initLambda);

void  getRateNbAmpGivenLambda(  double lambda,
								strtSBBHeader sbbHeader,CStructBook** structBook,
								double** min_amp, double** max_amp,
								double** min_rho,double** max_rho,
								double** min_phase,double** max_phase,
								double**** DfTable,
								int N_qrho,int N_qphi,
								double* q_rho,double* q_phi,
								int nb_xi,int nb_sample,
								double*** distAmp,
								double& totalRate,
								double& chosen_nb_amp);

void computeRDByAmpQuant(   strtSBBHeader sbbHeader, CStructBook** structBook,
							double** min_amp, double** max_amp,
							double** min_rho, double** max_rho,
							double** min_phase, double** max_phase,
							double**** DfTable,
							int N_qrho, int N_qphi,
							double* q_rho, double* q_phi,
							int nb_xi, int nb_sample,
							double chosenNbAmp,
							strtRDByAmpQuant* rdByAmpQuant);

void encodeSBExpAmpRangeBisecSearchArithCod(char*InputFile,
											double rateTarget,
											strtSBBHeader sbbHeader,
											CStructBook** structBook,
											CFileDfTable* DfTableFile,
											CFileDictionary* dicData,
											CDataSignal* dataSignal,
											int flLambdaBisecSearch,
											double initLambda);


void  getRateNbAmpGivenLambdaArithCod(  double lambda,
										strtSBBHeader sbbHeader,CStructBook** structBook,
										double** min_amp, double** max_amp,
										double** min_rho,double** max_rho,
										double** min_phase,double** max_phase,
										double**** DfTable,
										int N_qrho,int N_qphi,
										double* q_rho,double* q_phi,
										double nb_xi,double nb_sample,
										double*** distAmp,
										double& totalRate,
										double** chosen_nb_amp,
										double**** ampQStepEdge,
										double**** ampQStepProb,
										double**** rhoQStepEdge,
										double**** rhoQStepProb,
										double*** entropyAmpTable,
										double*** entropyRhoTable,
										double*** entropyPhiTable);

void computeRDByAmpQuantArithCod(   strtSBBHeader sbbHeader, CStructBook** structBook,
									double** min_amp, double** max_amp,
									double** min_rho, double** max_rho,
									double** min_phase, double** max_phase,
									double**** DfTable,
									int N_qrho, int N_qphi,
									double* q_rho, double* q_phi,
									double nb_xi, double nb_sample,
									double** chosenNbAmp,
									strtRDByAmpQuant* rdByAmpQuant,
									double**** ampQStepEdge,
									double**** ampQStepProb,
									double**** rhoQStepEdge,
									double**** rhoQStepProb,
									double*** entropyAmpTable,
									double*** entropyRhoTable,
									double*** entropyPhiTable);

void decodeSBExpAmpRangeArithCod(char*InputFile,
								 CFileDfTable* DfTableFile,
								 CFileDictionary* dicData);

void getQRhoQPhiArithCodDecode( 	double lambda,
									double min_amp, double max_amp,
									double min_rho, double max_rho,
									double min_phase, double max_phase,
									double**** DfTable,
									int N_qrho, int N_qphi,
									double* q_rho, double* q_phi,
									int ampRangeNumElement,
									double ampBar,
									int nb_amp,
									int iAmpRange,
									int& chosenQRhoIndex,
									int& chosenQPhiIndex,
									double* entropyRhoTable,
									double* entropyPhiTable);

int getIndexAmpRange(double amp, int NAmpRange, double* ampRangeLimit);

// void encode(char* InputFile, CFileGenData* genData, double rateTarget);


/**
 * \brief Encoding main function
 */

//void encodeLinkedStructBook(char* InputFile, CFileGenData* genData, double rateTarget);

/**
 * \brief Rate-distortion optimization
 */

void optimizeRDLinkedStrBook(double** origSignal,
							strtSBBHeader sbbHeader,
							CLinkStrBook* linkStrBookExpAux,
							strtQuantExp* chosenQuantExp,
							double rateTarget,
							strtQuantExp fixedQuantExp,
							int rdoptFlag,
							char* InputFile);

/**
 * \brief "Traces" a optimal rate-distortion decreasing curve
 */
void chooseBestQuantSetOptDecr(	int nChannel,
								int numQuant,
								strtRD** rd,
								strtRDIndex* rdIndex);

void quantSB(strtSBBHeader sbbHeader, CStructBook** structBook);

void synthSignal(double* recSignal,
                 double norm,
                 int sigSize,
                 strtContinuousExp* sb,
                 int numElement);

void synthSignalStrBook(double* recSignal,
						CStructBook* structBook,
                        strtSBBHeader mainHeader);

/**
 * \brief Encoding based on RD closed form formulation
 */

//void encodeRDClosedForm(char* InputFile, CFileGenData* genData, double rateTarget);

void printSBLinkedInfo(CStructBook* sbLinked,  strtSBBHeader sbbHeader,
                      double L_amp,  double L_rho, double L_phase);

/**
 * \brief Compute Rd optimum Solution
 */

void calcRDOptimumSolution(	int NAmpRange, double L_amp,  double L_rho, double L_phase,
							double* c1, double* c2, double* c3,
		                   	double& optR_amp, double* optR_rho, double* optR_phase);

/**
 * \brief Compute RD constants
 */

void calcRDConstant(	CStructBook* sbLinkedSepByAmp, double* ampBar, double Fs,
						double* c1, double* c2, double* c3, int NAmpRange, int signalSize,
						double L_amp,  double L_rho, double L_phase);

/**
 * \brief Calculate the squared norm
 */
double calcClosedFormSqrNorm(CParameter* parm, double Fs);
/**
 * \brief Calculate the distortion function 1
 */
double calcClosedFormDistortionFunc1(CParameter* parm, double Fs);
/**
 * \brief Calculate the distortion function 2
 */
double calcClosedFormDistortionFunc2(CParameter* parm, double Fs);
/**
 * \brief Calculate the distortion function 3
 */
double calcClosedFormDistortionFunc3(CParameter* parm, double Fs);

/**
 * \brief Calculate the second order equation roots
 */
void calcRoot2OrderBhaskara(	double a, double b, double c,
								CComplex& root1, CComplex& root2);

/**
 * \brief Calculate the theird order equation roots
 */
void calcRoot3OrderTartaglia(	double a, double b, double c,double d,
								CComplex& root1, CComplex& root2, CComplex& root3);

double linearQuant(double x, double step);

int linearQuantIndex(double x, double step);

double linearDeQuantIndex(int x, double step);

int geometricQuantIndex(double x, double step, int numPerOctave);

double geometricDeQuantIndex(int x, double step, int numPerOctave);


#endif
