#ifndef ADSPSORT_H
#define ADSPSORT_H

#include <iostream> 
#include <fstream>
#include <iomanip> 
#include <istream>
#include <string.h>

using namespace std;

void SWAPINT(int& a,int& b);

void SWAPDOUBLE(double& a,double& b);

void bubble_srtint( int* a, int n );

void bubble_srtdouble( double* a, int n );  

void bubble_srtdouble_noduplicate( double* a, int n, int& n_aux);

#endif
