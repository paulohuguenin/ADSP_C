#ifndef ADSPDEFINES_H
#define ADSPDEFINES_H

// Print rdbyamp files (rdbyamp_report.out, rdbyamp.out)
// #define PRINT_RDBYAMP

#define DEBUG
#define DEBUG_GRPBLOCK

// Decay variation factor used to group continued atoms between frames
const double DECAYVARIATIONFACTOR = 0.01;   // 1% of decay variation

// Use Uniform distribution or Generalized Gaussian distribution 
// to encode the decay parameter using arithmetic coding
#define USEGGD
#define USEGGDDSUP // deltasup

///////////////////////////////////////////////////
// PARAMETERS FOR STRBOOK STATISTICAL MODELING
// const double DELTASUPGGDMEAN = 0.0;
// const double DELTASUPGGDSCALE = 1685.1;
// const double DELTASUPGGDSHAPE = 1.196;
const double DELTASUPGGDMEAN1 = 0.0;
const double DELTASUPGGDSCALE1 = 10;
const double DELTASUPGGDSHAPE1 = 1;
const double DELTASUPGGDMEAN2 = 2047;
const double DELTASUPGGDSCALE2 = 10;
const double DELTASUPGGDSHAPE2 = 1;
const double DELTASUPGGDMEAN3 = 1023;
const double DELTASUPGGDSCALE3 = 1024;
const double DELTASUPGGDSHAPE3 = 100;
///////
const double AMPGGDMEAN = 0.0;
const double AMPGGDSCALE = 5.6208e-5;
const double AMPGGDSHAPE = 0.379;
//const double AMPGGDSCALE = 5.985324568452911e-9;
//const double AMPGGDSHAPE = 0.161;
//////
const double DECAYGGDMEAN = 0.0;
const double DECAYGGDSCALE = 0.0235;
const double DECAYGGDSHAPE = 0.49;


#endif