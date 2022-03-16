//net.h
#ifndef NET_H
#define NET_H

#include <iostream>
#include <stdlib.h>
#include <math.h>

#include <fstream>
#include <iomanip> 
#include <istream>

#include "fftw3.h"
#include "cgfunc.h"
#include "cgmatrix.h"
#include "complex.h"
#include "defines.h"
#include "datasignal.h"
#include "parameter.h"
#include "filemgr.h"
#include "structbook.h"

using namespace std;

class CNet {

protected:

    int m_signalSize;

    int m_num;

    double* m_x;

    double* m_y;

    // double* mapMinMaxFunc(  double* x,
    //                         double* x1_xoffset,
    //                         double* x1_gain,
    //                         double x1_ymin);

    // void hiddenFunc(double* x1_xoffset,
    //                 double* x1_gain,
    //                 double x1_ymin,
    //                 double* b,
    //                 double w[][16]);

    // void outputFunc(double* b,
    //                 double w[][3]);

    void coisa(double w[][16]);

public:
    // Constructors
    CNet ();

    // Destructor
    ~CNet();

    // Set methods
    void setSignalSize(int signalSize);
    void setX(double* x);
    void setNumber(int num);
    // Get methods
    double* getY();
    void net();
};

#endif
