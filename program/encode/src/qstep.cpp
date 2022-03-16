
#include "qstep.h"

double gamma(double x)
{
    int i,k,m;
    double ga,gr,r,z;

    static double g[] = {
        1.0,
        0.5772156649015329,
       -0.6558780715202538,
       -0.420026350340952e-1,
        0.1665386113822915,
       -0.421977345555443e-1,
       -0.9621971527877e-2,
        0.7218943246663e-2,
       -0.11651675918591e-2,
       -0.2152416741149e-3,
        0.1280502823882e-3,
       -0.201348547807e-4,
       -0.12504934821e-5,
        0.1133027232e-5,
       -0.2056338417e-6,
        0.6116095e-8,
        0.50020075e-8,
       -0.11812746e-8,
        0.1043427e-9,
        0.77823e-11,
       -0.36968e-11,
        0.51e-12,
       -0.206e-13,
       -0.54e-14,
        0.14e-14};

    if (x > 171.0) return 1e308;    // This value is an overflow flag.
    if (x == (int)x) {
        if (x > 0.0) {
            ga = 1.0;               // use factorial
            for (i=2;i<x;i++) {
               ga *= i;
            }
         }
         else
            ga = 1e308;
     }
     else {
        if (fabs(x) > 1.0) {
            z = fabs(x);
            m = (int)z;
            r = 1.0;
            for (k=1;k<=m;k++) {
                r *= (z-k);
            }
            z -= m;
        }
        else
            z = x;
        gr = g[24];
        for (k=23;k>=0;k--) {
            gr = gr*z+g[k];
        }
        ga = 1.0/(gr*z);
        if (fabs(x) > 1.0) {
            ga *= r;
            if (x < 0.0) {
                ga = -M_PI/(x*ga*sin(M_PI*x));
            }
        }
    }
    return ga;
}

double ggdPdf(double x, double mean, double scale, double shape)
{
// 	double ggd = 	(DECAYGGDSHAPE/(2.0*DECAYGGDSCALE*gamma(1.0/DECAYGGDSHAPE)))*
// 					exp(-pow( (fabs(x-DECAYGGDMEAN)/ DECAYGGDSCALE),DECAYGGDSHAPE) );
	double ggd = 	(shape/(2.0*scale*gamma(1.0/shape)))*
					exp(-pow( (fabs(x-mean)/ scale),shape) );				
	return ggd;
}

void setQStepFeatureGGD(double* qstepCenter, double* qstepEdge, 
						double* qstepProb, double* qstepCumProb,
						double qstep, int numQStep, int leftShift,
						double mean, double scale, double shape)
{
	int i;
	double scalingFactor=0.0;
	for (i=0;i<numQStep;i++)
	{
		qstepCenter[i] = (i- leftShift)*qstep;
		qstepEdge[i] = qstepCenter[i] - (qstep/2.0);
		
		qstepProb[i] = ggdPdf(qstepCenter[i], mean, scale, shape)*qstep;
		scalingFactor += qstepProb[i];
	}
	qstepEdge[numQStep] = qstepCenter[numQStep-1] + (qstep/2.0);
	//////////
	double cumProb =0.0;
	double new_scaling_factor =0.0;
	for (i=0;i<numQStep;i++)
	{
		qstepProb[i] = qstepProb[i]/scalingFactor;
		if (qstepProb[i]<1e-5) qstepProb[i] =1e-5;
		new_scaling_factor += qstepProb[i];
	}
	for (i=0;i<numQStep;i++)
	{
		qstepProb[i] = qstepProb[i]/new_scaling_factor;
		cumProb += qstepProb[i];
		qstepCumProb[i] = cumProb;
	}
}

void setQStepFeatureGGDDeltaSup(double* qstepCenter, double* qstepEdge, 
								double* qstepProb, double* qstepCumProb,
								double qstep, int numQStep, int leftShift)
{
	int i;
	double scalingFactor=0.0;
	for (i=0;i<numQStep;i++)
	{
		qstepCenter[i] = (i- leftShift)*qstep;
		qstepEdge[i] = qstepCenter[i] - (qstep/2.0);
		
		qstepProb[i] = qstep*(  0.5*ggdPdf(qstepCenter[i], DELTASUPGGDMEAN1, DELTASUPGGDSCALE1, DELTASUPGGDSHAPE1) +
								ggdPdf(qstepCenter[i], DELTASUPGGDMEAN2, DELTASUPGGDSCALE2, DELTASUPGGDSHAPE2) +
								1.2 * ggdPdf(qstepCenter[i], DELTASUPGGDMEAN3, DELTASUPGGDSCALE3, DELTASUPGGDSHAPE3) );
		scalingFactor += qstepProb[i];
	}
	qstepEdge[numQStep] = qstepCenter[numQStep-1] + (qstep/2.0);
	//////////
	double cumProb =0.0;
	double new_scaling_factor =0.0;
	for (i=0;i<numQStep;i++)
	{
		qstepProb[i] = qstepProb[i]/scalingFactor;
		if (qstepProb[i]<1e-5) qstepProb[i] =1e-5;
		new_scaling_factor += qstepProb[i];
	}
	for (i=0;i<numQStep;i++)
	{
		qstepProb[i] = qstepProb[i]/new_scaling_factor;
		cumProb += qstepProb[i];
		qstepCumProb[i] = cumProb;
	}
}



