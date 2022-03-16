#include "metric.h"

double computeMSE(double* origSignal, 
                  double* recSignal, 
                  int sigSize)
{
    double temp=0;
    double calc;
    double acum = 0;
    double sqrnorm=0;

    for(int i=0;i<sigSize;i++)
    {
        calc = recSignal[i] - origSignal[i];
        temp = temp + calc*calc;
    }

    return (double)temp/(double)sigSize;
}


double computeSqrErrorPerSample(  double* origSignal, 
								  double* recSignal, 
								  int a,
								  int b)
{
	double temp=0;
    double calc;
    double acum = 0;
    double sqrnorm=0;

    for(int i=a;i<=b;i++)
    {
        calc = recSignal[i] - origSignal[i];
        temp = temp + calc*calc;
    }
    
    return temp/static_cast<double>(b-a+1);
}


double computeSqrError( double* origSignal, 
						double* recSignal, 
						int a,
						int b)
{
	double temp=0;
    double calc;
    double acum = 0;
    double sqrnorm=0;

    for(int i=a;i<=b;i++)
    {
        calc = recSignal[i] - origSignal[i];
        temp = temp + calc*calc;
    }
    
    return temp;
}

double computeSqrNorm(double* origSignal, 
					  int sigSize)
{
    double temp=0.0;
    double calc;
    double acum = 0.0;
    double sqrnorm=0.0;

    for(int i=0;i<sigSize;i++)
    {
        calc = origSignal[i];
        temp = temp + calc*calc;
    }

    return temp;
}