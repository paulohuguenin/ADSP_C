#ifndef MPURSUIT_H
#define MPURSUIT_H

#include "filemgr.h"
#include "structbook.h"
#include "dictionary.h"
#include "datasignal.h"
#include "parameter.h"
#include "structs.h"
#include "cgmatrix.h"
#include <unistd.h>
// #include "libannformp.h"
#include "net.h"
#include "fftw3.h"


void matchingPursuit(   cgMatrix<double>& residue,
                        strtParameter* chosenParm,
                        int dicSize,
                        CDataSignal* dataSignal,
                        CFileDictionary* dicData,
                        CFileDecomp* genData, int step, int chosenDic,
                        strtParameter** dicAtoms,
                        int& a0,
                        int& b0,
                        int flagOMP);

void MPTradicional( cgMatrix<double>& residue,
                    CComplex* complexAtom,
                    double& maxInnerProd,
                    double& chosenOptPhase,
                    int& chosenTau,
                    double& chosenXi,
                    // CDictionary* dic,
                    int N, //int dicSize
                    // int decincAsymmFlag,
                    double s,
                    // int delta_tau,
                    // double* xi_vec,
                    // int Nfreq,
                    char* fileName);

void fastMPKolasa(  cgMatrix<double>& residue,
                    double& maxInnerProd,
                    double& chosenOptPhase,
                    double& chosenXi,
                    int N, //dicSize
                    int tau,
                    double s,
                    double* realAtom,
                    char* fileName);

/*CParameter* fastMPKolasaModified( cgMatrix<double>& residue,
                                    CDataSignal* dataSignal,
                                    CFileDictionary* dicData,
                                    CDictionary* dic);*/

void fastMPKolasaModified(  cgMatrix<double>& residue,
                            double& maxInnerProd,
                            double& chosenOptPhase,
                            int& chosenTau,
                            int N, //int dicSize
                            int decincAsymmFlag,
                            int delta_tau,
                            double xi,
                            double* realAtomWin,
                            CComplex* complexAtomXi,
                            CComplex* complexAtom2Xi,
                            double* conv_zxi0_w2,
                            char* fileName,
                            double s);

// double computeOptimumPhase( cgMatrix<double>& residue,
//                             double xi,
//                             double& innerProd,
//                             int signalSize,
//                             CComplex* complexAtom);

void setParameters (strtParameter* chosenParm,
                    int dicType,
                    double innerProd,
                    double s,
                    double xi,
                    double opt_phase,
                    int tau,
                    int a,
                    int b,
                    double beta);

void computeOptimumPhase(   cgMatrix<double>& residue,
                            double& opt_phase,
                            double& innerProd,
                            int signalSize,
                            double xi,
                            CComplex* complexAtom,
                            char* fileName,
                            double s,
                            double tau,
                            int N);

//void optimizeContinuousParmsBateman(   cgMatrix<double>& residue,
  //                                              strtParameter* parm);


void adjustParameters ( cgMatrix<double>& residue,
                        strtParameter* chosenParm);

// void setRealAtom (strtParameter* parm);

// double* getRealAtom (strtParameter* parm);

void updateResidue (cgMatrix<double>& residue,
                    double& norm,
                    int dicSize,
                    strtParameter* parm);

void updateResidue (cgMatrix<double>& residue,
                    cgMatrix<double>& signal,
                    cgMatrix<double>& b,
                    cgMatrix<double>& v,
                    cgMatrix<double>& Ai,
                    cgMatrix<double>& a,
                    double& norm,
                    int dicSize,
                    int step,
                    cgMatrix<double>& prevAtoms,
                    strtParameter* parm);


// void updateResidue (cgMatrix<double>& residue,
//                     cgMatrix<double>& gamma,
//                     double& norm,
//                     int dicSize,
//                     int step,
//                     strtParameter* parm);

void ANN(   CFileDecomp* genData,
            cgMatrix<double>& residue,
            int dicSize,
            int& chosenDic,
            int& chosenNet,
            int step,
            int L,
            double* approxRatio);

void DANNO( CFileDecomp* genData,
            cgMatrix<double>& residue,
            int dicSize,
            int& chosenDic,
            int& chosenNet,
            int step,
            int L,
            double* approxRatio);

void matchingPursuitEDA(    cgMatrix<double>& residue,
                            strtParameter* chosenParm,
                            int dicSize,
                            CDataSignal* dataSignal,
                            CFileDictionary* dicData,
                            CFileDecomp* genData, int step, int chosenDic,
                            strtParameter** dicAtoms,
                            int& a0,
                            int& b0,
                            int flagOMP);

void optimizeContinuousParms(   cgMatrix<double>& residue,
                                                strtParameter* parm);


// CDictionary* getDic

#endif
