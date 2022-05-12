#ifndef SNDOUTFILE_H
#define SNDOUTFILE_H

#define STD_SAMPLE_RATE 44100
#define STD_CHANNELS 1
#define STD_FORMAT SF_FORMAT_WAV|SF_FORMAT_PCM_16
#define STD_SEEK_STATUS false

#include <string.h>
#include <iostream>
#include <cstdlib>
//#include "ptypes.h"
#include "sndfile.h"

using namespace std;

class SndOutFile
{
	bool file_open;
	SF_INFO info;
	SNDFILE *out_file;

public:

    SndOutFile(void);
    SndOutFile(int nchan, int samplerate);
    virtual ~SndOutFile(void);

    void createFile(char* file_name);
    void writeDouble(double *sample);
    void writeDouble(double *sample, int items);
    void writeDoubleVector(double *vector, int length);
    void setFrame(sf_count_t frame);
    void setSampleRate(int sample_rate);
    int  getSampleRate();
    void setChannels(int no_channels);
    void setFormat(int format);
    void setSeekStatus(bool status);
    bool getFileStatus();
    void closeSndFile();	
};

#endif
