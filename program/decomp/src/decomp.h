#ifndef DECOMP_H
#define DECOMP_H

#include "filemgr.h"
#include "structbook.h"
#include "dictionary.h"
#include "datasignal.h"
#include "mpursuit.h"
#include "structs.h"
#include "time.h"
// #include "libannformp.h"


void decompElectric(CFileDecomp* genData,
		    		CFileDecompBlockRange* blockRange,
	    			CFileDictionary* dicData,
	    			char* InputFile);

void decompAudio(	CFileDecomp* genData,
		    		CFileDecompBlockRange* blockRange,
	    			CFileDictionary* dicData,
	    			char* InputFile);

void decompNoise(	CFileDecomp* genData,
		    		CFileDecompBlockRange* blockRange,
	    			CFileDictionary* dicData);

//void decompECG(		CFileDecomp* genData,
	//	    		CFileDecompBlockRange* blockRange,
	  //  			CFileDictionary* dicData,
	    //			char* InputFile);
void decompEDA(		CFileDecomp* genData,
		    		CFileDecompBlockRange* blockRange,
	    			CFileDictionary* dicData,
	    			char* InputFile);

#endif
