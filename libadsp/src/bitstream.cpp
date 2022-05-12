// bitstream.cpp

#include "bitstream.h"

//===============================================================
// Function: 
// Goal: Contructor
// Return:
//===============================================================

CBitStream::CBitStream()
{

	m_bitStream = NULL;
	m_iNumElement = 0;
	m_bitCounter = 0;
	m_byteCounter = 0;

}

//===============================================================
// Function: 
// Goal: Destructor
// Return:
//===============================================================

CBitStream::~CBitStream()
{
	if (m_bitStream !=NULL)
	{
		delete [] m_bitStream;
	}
}

//===============================================================
// Function: setBitStream
// Goal: 
// Return:
//===============================================================

void CBitStream::setBitStream(int nbits)
{
	m_iNumElement =(int) ceil( (double)nbits / 8.0);
	if (m_bitStream == NULL)
	{
		m_bitStream = new unsigned char [m_iNumElement+1];
	}
	for(int i=0;i<m_iNumElement+1;i++)
	{
		m_bitStream[i] =0;
	}
}

//===============================================================
// Function: setBitStream
// Goal: 
// Return:
//===============================================================

void CBitStream::setBitStream(FILE* stream)
{
	if (m_bitStream==NULL)
	{
		m_bitStream = new unsigned char [m_iNumElement+1];
	}
	
	
	fread(m_bitStream, sizeof(unsigned char), m_iNumElement, stream);
	
}

//===============================================================
// Function: resetBitStreamCntr
// Goal: 
// Return:
//===============================================================

void CBitStream::resetBitStreamCntr()
{
	m_bitCounter = 0;
	m_byteCounter = 0;
}

//===============================================================
// Function: writeBuffer (int)
// Goal: 
// Return:
//===============================================================

void CBitStream::writeBuffer(int sbField, int nbits)
{
	int sbFieldAux = sbField;
	unsigned char aux;
	int nbits_lack = nbits;
	int nbits_written=0;
	int nbits_empty=0;

	while (nbits_lack!=0)
	{
		nbits_empty = 8 - m_bitCounter;

		if (nbits_lack >= nbits_empty)
		{
			aux = 0xFF;
			aux = aux >> m_bitCounter;
			aux = aux & (unsigned char) sbFieldAux;
			aux = aux << m_bitCounter;
			m_bitStream[m_byteCounter] |= aux;
			
			nbits_written = nbits_empty;
			nbits_lack -= nbits_written;
			sbFieldAux = sbFieldAux >> nbits_written;

			m_byteCounter++;
			m_bitCounter = 0;
		}
		else
		{
			aux = 0xFF;
			aux = aux >> ( m_bitCounter + (nbits_empty - nbits_lack) );
			aux = aux & (unsigned char) sbFieldAux;
			aux = aux << m_bitCounter;
			m_bitStream[m_byteCounter] |= aux;

			nbits_written = nbits_lack;
			nbits_lack -= nbits_written;

			m_bitCounter += nbits_written;
		}
	}
}

//===============================================================
// Function: readBuffer
// Goal: 
// Return:
//===============================================================

void CBitStream::readBuffer(int& sbField, int nbits)
{
	int sbFieldAux;  
	unsigned char aux;
	int nbits_lack = nbits;
	int nbits_read=0;
	int nbits_not_read=0;

	sbField = 0;

	while (nbits_lack!=0)
	{
		nbits_not_read = 8 - m_bitCounter;

		if (nbits_lack >= nbits_not_read)
		{
			aux = 0xFF;
			aux = aux << m_bitCounter;
			aux = aux & m_bitStream[m_byteCounter];
			aux = aux >> m_bitCounter;
			sbFieldAux = (int) aux;
			sbFieldAux = sbFieldAux << (nbits - nbits_lack);
			sbField |= sbFieldAux;
			
			nbits_read = nbits_not_read;
			nbits_lack -= nbits_read;
			

			m_byteCounter++;
			m_bitCounter = 0;
		}
		else
		{
			aux = 0xFF;
			aux = aux >> ( m_bitCounter + (nbits_not_read - nbits_lack) );
			aux = aux << m_bitCounter;
			aux = aux & m_bitStream[m_byteCounter];
			aux = aux >> m_bitCounter;
			sbFieldAux = (int) aux;
			sbFieldAux = sbFieldAux << (nbits - nbits_lack);
			sbField |= sbFieldAux;

			nbits_read = nbits_lack;
			nbits_lack -= nbits_read;

			m_bitCounter += nbits_read;
		}
	}
}

//===============================================================
// Function: float2bin
// Goal: 
// Return:
//===============================================================

strtFpBinForm CBitStream::float2bin(float input, int int_nbits, int dec_nbits)
{
	strtFpBinForm aux;
	
	aux.integer = (int)floor(input);
	aux.int_nbits = int_nbits;
	
	aux.decimal=0;
	aux.dec_nbits = dec_nbits;
	float x = input - aux.integer;
	for (int i = 0; i < dec_nbits; i++)
	{
		if ((x*2.0)>=1)
		{
			aux.decimal = (aux.decimal << 1) | 0x1;
			x = x*2.0 - floor(x*2.0);
		}
		else
		{
			aux.decimal = aux.decimal << 1;
			x = x * 2.0;
		}
	}
	
	return aux;
}

//===============================================================
// Function: bin2double
// Goal: 
// Return:
//===============================================================

double CBitStream::bin2double(strtFpBinForm input)
{
	double output=0;
	int x;
	double acum = 0;
	for (int i = 0; i < input.dec_nbits; i++)
	{
		x =  (input.decimal >> i) &  0x1;
		acum = acum + (double)x * pow(2, -(input.dec_nbits-i));
	}
	
	output = ( input.integer & ~(0xffffffff << input.int_nbits) ) + acum;
	
	return output;
}

//===============================================================
// Function: double2bin
// Goal: 
// Return:
//===============================================================

strtFpBinForm CBitStream::double2bin(double input, int int_nbits, int dec_nbits)
{
	strtFpBinForm aux;
	
	aux.integer = (int)floor(input);
	aux.int_nbits = int_nbits;
	
	aux.decimal=0;
	aux.dec_nbits = dec_nbits;
	double x = input - aux.integer;
	for (int i = 0; i < dec_nbits; i++)
	{
		if ((x*2.0)>=1)
		{
			aux.decimal = (aux.decimal << 1) | 0x1;
			x = x*2.0 - floor(x*2.0);
		}
		else
		{
			aux.decimal = aux.decimal << 1;
			x = x * 2.0;
		}
	}
	
	return aux;
}

//===============================================================
// Function: bin2float
// Goal: 
// Return:
//===============================================================

float CBitStream::bin2float(strtFpBinForm input)
{
	float output=0;
	int x;
	float acum = 0;
	for (int i = 0; i < input.dec_nbits; i++)
	{
		x =  (input.decimal >> i) &  0x1;
		acum = acum + (float)x * pow(2, -(input.dec_nbits-i));
	}
	
	output = ( input.integer & ~(0xffffffff << input.int_nbits) ) + acum;
	
	return output;
}

//===============================================================
// Function: setNumElement
// Goal: 
// Return:
//===============================================================

void CBitStream::setNumElement(int numElement)
{

	m_iNumElement = numElement;
}

//===============================================================
// Function: getBitStream
// Goal: 
// Return:
//===============================================================

unsigned char* CBitStream::getBitStream() const
{
	return m_bitStream;
}

//===============================================================
// Function: getBitStreamElement
// Goal: 
// Return:
//===============================================================

unsigned char CBitStream::getBitStreamElement(int i) const
{
	return m_bitStream[i];
}




//===============================================================
// Function: getNumElement
// Goal: 
// Return:
//===============================================================

int CBitStream::getNumElement() const
{
	return m_iNumElement;
}
