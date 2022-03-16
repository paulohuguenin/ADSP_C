#ifndef METRIC_H
#define METRIC_H

using namespace std;

double computeMSE(double* origSignal, 
                  double* recSignal, 
                  int sigSize);

double computeSqrErrorPerSample(  double* origSignal, 
								  double* recSignal, 
								  int a,
								  int b);
								  
double computeSqrError(   double* origSignal, 
						  double* recSignal, 
						  int a,
						  int b);

double computeSqrNorm(double* origSignal, 
					  int sigSize);						  
                  
#endif