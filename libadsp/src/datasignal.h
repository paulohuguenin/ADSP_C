// datasignal.h

#ifndef DATASIGNAL_H
#define DATASIGNAL_H

#include <typeinfo>
#include <iostream>
#include <stdlib.h>
#include <math.h>

#include "SndInFile.h"
#include "SndOutFile.h"
#include "newcomtd_99.h"
#include "wfdb.h"

//#include "dictionary.h"

#define NMAX_CANAIS_ANALOG 300
#define NMAX_CANAIS_DIGITAL 300

/* Definitions for random functions */
#define Sem_Inicial -1234567
#define IA	    16807
#define IM 	    2147483647
#define AM 	    (1.0/IM)
#define IQ 	    127773
#define IR 	    2836
#define NTAB 	    32
#define NDIV 	    (1+(IM-1)/NTAB)
#define EPS 	    1.2e-7
#define RNMX 	    (1.0-EPS)

using namespace std;


/**
 * \brief Data signal class
 */

class CDataSignal {
protected:
	char*				m_fileName;
	int                 m_numSignal;
	int					m_signalSize;
	double**			m_signal;
	int					m_blockSize;
	int					m_blockHop;
    int                 m_numBlock;
	double*				m_norm;
	int					m_type;
	double				m_Fs;

public:

    // Constructors
    CDataSignal();

    // Destructor
    virtual ~CDataSignal();

    void        setFileName(const char* file);
    void        setNumSignal(const int);
    void        setSignalSize(const int);
    void        setBlockSize(int);
    void        setBlockHop(int);
    void        setNumBlock(int);
    void        setNorm();
    void        setNorm(int signum, double norm);
    void        setType(int);
    void        setSamplingRate(double);

    char*       getFileName();
    int         getNumSignal();
    int	        getSignalSize();
    double**    getSignal();
    int         getBlockSize();
    int         getBlockHop();
    int         getNumBlock();
    double      getNorm(int );
    int         getType();
    double		getSamplingRate();

    // Interface
    virtual void    setSignal()=0;
    virtual void    setSignal(double** signal)=0;
    virtual void    saveSignal()=0;
};


/**
 * \brief Comtrade signal class
 */

class CComtradeSignal : public CDataSignal {

	CConfigComtrade*	m_config;
	int					m_numDigitalSignal;
	char**				m_digitalSignal;
	
public:

	// Constructors
	CComtradeSignal();

	// Destructor
	~CComtradeSignal();

	// Interface

    virtual void    setSignal();
    virtual void    setSignal(double** signal){};
    virtual void    saveSignal(){};

	// Intrinsic methods
	int					setConfig();
	int					setComtradeSignal();
	void				cutComtrade(char* file, 
									int initSample, 
									int endSample);

	CConfigComtrade*	getConfig() const; 
	long				getNumberOfSamples() const;
	int					getNumberOfAnalogChannels() const;
	int					getNumberOfDigitalChannels() const;
	char**				getDigitalSignal() const;
	double				getSamplingRate(int index) const;
	double				getFundamentalFrequency() const;
};

/**
 * \brief Audio signal class
 */

class CAudioSignal : public CDataSignal {
	
public:

	// Constructors
	CAudioSignal(){m_type=2;};

	// Destructor
	virtual ~CAudioSignal(){};

	// Interface
	virtual void    setSignal();
    virtual void    setSignal(double** signal);
    virtual void    saveSignal();

	// Intrinsic
	double				getSamplingRate() const;
};

/**
 * \brief Gaussian noise signal class
 */

class CNoiseSignal : public CDataSignal {
public:
	// Constructors
	CNoiseSignal(){m_type=3;};
	// Destructor
	virtual ~CNoiseSignal(){};

	// Interface
	virtual void    setSignal();
    virtual void    setSignal(double** signal){};
    virtual void    saveSignal(){};

	// Intrinsic
	double				getSamplingRate() const;

	double gausdev(long int *sem);
	double ran1(long int *sem);
};

/**
 * \brief ECG signal class
 */

class CECGSignal : public CDataSignal {
	
public:

	// Constructors
	CECGSignal(){m_type=4;};

	// Destructor
	virtual ~CECGSignal(){};

	// Interface
    virtual void    setSignal();
    virtual void    setSignal(double** signal){};
    virtual void    saveSignal(){};

	// Intrinsic
	double				getSamplingRate() const;
};

/*class CEDASignal : public CEDASignal {

public:

	// Constructors
	CEDASignal(){m_type=5;};

	// Destructor
	virtual ~CEDASignal(){};

	// Interface
    virtual void    setSignal();
    virtual void    setSignal(double** signal){};
    virtual void    saveSignal(){};

	// Intrinsic
	double				getSamplingRate() const;

};*/

#endif
