#ifndef LOADFILE_H
#define LOADFILE_H

#include "filemgr.h"
#include "structbook.h"
#include "linkstrbook.h"
#include "dictionary.h"
#include "datasignal.h"

using namespace std;

/**
 * \brief Load .SBB file header
 */

strtSBBHeader loadSBHeader(char* file);

strtSBBHeader loadSBHeaderGrp(char* file);

/**
 * \brief Load SBB file for exponential dictionary
 */

void loadSBExp(		char* InputFile,
				   	CFileDecompBlockRange* blockRange,
				   	strtSBBHeader sbbHeader,
				   	CStructBook** structBook);
				   	
void loadSBExpGrp(		char* InputFile,
				   	CFileDecompBlockRange* blockRange,
				   	strtSBBHeader sbbHeader,
				   	CStructBook** structBook);				   	

void saveSBHeader(char* file, strtSBBHeader sbbHeader);

void saveSBHeaderGrp(char* file, strtSBBHeader sbbHeader);
				 
void saveSBExp(	char* InputFile,
				strtSBBHeader sbbHeader,
				CStructBook** structBook);

void saveSBExpGrp(	char* InputFile,
					strtSBBHeader sbbHeader,
					CStructBook** structBook);				
				   	
/**
 * \brief Load SBB file for exponential dictionary
 */

void loadSB(		char* InputFile,
				   	CFileDecompBlockRange* blockRange,
				   	strtSBBHeader sbbHeader,
				   	CStructBook** structBook);
               
/**
 * \brief Load rate-distortion file
 */    

void loadRDFile(int numSignal,int numQuant,strtRD** rd);

#endif