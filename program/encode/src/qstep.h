#ifndef QSTEP_H
#define QSTEP_H

#include <iostream> 
#include <fstream>
#include <iomanip> 
#include <istream>
#include <math.h>

#include "adspdefines.h"

using namespace std;

double gamma(double x);

double ggdPdf(double x, double mean, double scale, double shape);

void setQStepFeatureGGD(double* qstepCenter, double* qstepEdge, 
						double* qstepProb, double* qstepCumProb,
						double qstep, int numQStep, int leftShift,
						double mean, double scale, double shape);

void setQStepFeatureGGDDeltaSup(double* qstepCenter, double* qstepEdge, 
								double* qstepProb, double* qstepCumProb,
								double qstep, int numQStep, int leftShift);						


#endif