// parameter.h
#ifndef PARAMETER_H
#define PARAMETER_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "structs.h"

using namespace std;

/**
 * \brief Parameter class
 */

class CParameter 
{

protected:
	int			m_type;

public:
	CParameter(){m_type=0;};
	~CParameter(){};

	virtual void setType(int type)=0;
	virtual int getType()=0;
	virtual void printParm2Screen()=0;
};

/**
 * \brief Exponential parameter class
 */

class CExpParm: public CParameter
{
public:
	double	innerProd;
	double	rho;
	double	xi;
	double	phase;
	int		a;
	int		b;

	CExpParm(){m_type=1;};
	~CExpParm(){};

	virtual void setType(int type);
    void setParm(strtContinuousExp exp_parm);
	virtual int getType();
	virtual void printParm2Screen();
};

#endif