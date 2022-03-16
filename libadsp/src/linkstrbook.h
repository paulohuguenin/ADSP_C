// linkstrbook.h
#ifndef LINKSTRBOOK_H
#define LINKSTRBOOK_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "defines.h"
#include "bitstream.h"
#include "parameter.h"
#include "structbook.h"
#include "structs.h"
#include "dictionary.h"

using namespace std;

/**
 * \brief Linked structure book class
 */

class CLinkStrBook
{

protected:
    int     m_iNumElement;
    double  m_norm;

public:

    /**
     * \brief Constructor.
     */
    CLinkStrBook();
	CLinkStrBook(const CLinkStrBook& linkStrBook);
	CLinkStrBook(const CLinkStrBook* linkStrBook);
	void copy(const CLinkStrBook* linkStrBook);
    /**
     * \brief Destructor.
     */
    virtual ~CLinkStrBook(){};

    // Interface
    virtual void load(  CStructBook** structBook,
                        strtSBBHeader mainHeader,
                        int iSignal)=0;

    double linearQuant(double x, double step);
    
    virtual void synthSignal(  double* recSignal,
                                strtSBBHeader mainHeader)=0;
 
                                
    virtual void synthAtom(   double* recSignal,
							  strtSBBHeader mainHeader,
							  int nElement)=0;

    // Set methods
    void    setNumElement(int numElement);
    void    setNorm(double norm);

    // Get methods
    int     getNumElement() const;
    double     getNorm() const;

    // operators
    virtual CLinkStrBook& operator=(CLinkStrBook& c);
};

/**
 * \brief Exponential linked structure book class
 */

class CLinkStrBookExp : public CLinkStrBook
{
    strtLinkExp* m_linkStrBook;
    strtLinkSBQHeadExp quantHeader;

public:

    /**
     * \brief Constructor.
     */
    CLinkStrBookExp();
    CLinkStrBookExp(const CLinkStrBookExp& linkStrBook);
	CLinkStrBookExp(const CLinkStrBookExp* linkStrBook);
	void copy(const CLinkStrBookExp* linkStrBook);
    /**
     * \brief Destructor.
     */
    virtual ~CLinkStrBookExp();

    /**
     * \brief load method
     */
    virtual void load(  CStructBook** structBook,
                        strtSBBHeader mainHeader,
                        int iSignal);

    void allocateElement (int numElement);

    /**
     * \brief Allocate mamory space for more element
     */
    void allocateElement ();
    /**
     * \brief Add decaying factor information
     */
    void addDecayingFactor (int index, double rho);
    /**
     * \brief Convert Linked Structure book to unlinked structure book
     */
    void linkStrBook2StrBook(CStructBook* structBook,
                             strtSBBHeader mainHeader);
    
    /**
     * \brief Print parameters to screen
     */
    void Print();
    /**
     * \brief Find the range of the parameters
     */
    void findParmRange( double& min_coef,
                        double& max_coef,
                        double& min_rho,
                        double& max_rho,
                        double& min_phase,
                        double& max_phase);

    void configQuantizer(   int nbits_amp,
                            int nbits_rho,
                            int nbits_xi,
                            int nbits_phase,
                            int nbits_sample,
                            int nbits_block);

    void quantize();
    
    int computeNumBits();
    
    virtual void synthSignal(  double* recSignal,
                                strtSBBHeader mainHeader);
 
                                
    virtual void synthAtom(   double* recSignal,
							  strtSBBHeader mainHeader,
							  int nElement);

    strtLinkExp* getLinkStrBook() const;

    virtual CLinkStrBookExp& operator=(CLinkStrBookExp& c);

};


#endif