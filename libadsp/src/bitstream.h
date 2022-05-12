// bitstream.h

#ifndef BITSTREAM_H
#define BITSTREAM_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

using namespace std;

typedef struct fpBinForm_tag 
{
	int integer;
	int int_nbits;
	int decimal;
	int dec_nbits;
} strtFpBinForm;


class CBitStream {

	unsigned char*	m_bitStream;
	int		m_iNumElement;
	int		m_bitCounter;
	int		m_byteCounter;

public:

	// Constructor
	CBitStream();
	// Destructor
	~CBitStream();

	// ==============================
	
	
	
	void setBitStream(int nbits);
	
	void setBitStream(FILE* stream);
	
	void resetBitStreamCntr();

	void writeBuffer(int sbField, int nbits);

	void readBuffer(int& sbField, int nbits);
	
	strtFpBinForm float2bin(float inputint, int int_nbits, int dec_nbits);
	
	float bin2float(strtFpBinForm input);
	
	strtFpBinForm double2bin(double inputint, int int_nbits, int dec_nbits);
	
	double bin2double(strtFpBinForm input);
	
	void setNumElement(int numElement);

	unsigned char* 	getBitStream() const;
	unsigned char 	getBitStreamElement(int i) const;
	int 		getNumElement() const;

};


#endif
