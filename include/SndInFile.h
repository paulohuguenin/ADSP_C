#ifndef SNDINFILE_H
#define SNDINFILE_H

#define STD_SAMPLE_RATE 44100
#define STD_CHANNELS 1
#define STD_FORMAT SF_FORMAT_WAV|SF_FORMAT_PCM_16
#define STD_SEEK_STATUS false

#include <iostream>
#include <cstdlib>
#include "sndfile.h"

using namespace std;

class SndInFile
{
	bool 		file_open;
	SF_INFO 	info;
	SNDFILE 	*in_file;

	int		numberOfSamples;
	int		numberOfChannels;
	int		signalSize;
	int 		samplingRate;
	double		signalTime;

public:

	// Constructor
			SndInFile(char *filename);

	// Destructor
	virtual 	~SndInFile(void);

	// Methods for retrieving information
	int		getNumberOfSamples();
	int		getNumberOfChannels();
	int		getSignalSize();
	int 		getSamplingRate();
	double		getSignalTime();

	// Methods for audio file handling
	bool 		getFileStatus();
	void 		closeSndFile();
	void		print();
	double		*getSignal(int channel);
};


#endif
