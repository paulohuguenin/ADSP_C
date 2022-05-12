#include "SndInFile.h"

SndInFile::SndInFile(char *file_name)
{
	// Opening audio file in READ mode
	cout << "Opening audio file: " << file_name << endl;
	in_file = sf_open(file_name, SFM_READ, &info);
	file_open = true;

	// Retrieving number of samples in the audio file 
	numberOfSamples = info.frames;
	cout << "Number of samples: " << numberOfSamples << endl;
	
	// Retrieving number of channels in the audio file 
	numberOfChannels = info.channels;
	cout << "Number of channels: " << numberOfChannels << endl;

	// Determining signal size
	signalSize = numberOfSamples * numberOfChannels;

	// Retrieving the sampling rate for the audio file 
	samplingRate = info.samplerate;
	cout << "Sampling rate: " << samplingRate << endl;

	// Calculating audio signal time
	signalTime = (double)numberOfSamples / (double)samplingRate;
	cout << "Signal time: " << signalTime << " seconds" << endl;

	cout << endl;
}

SndInFile::~SndInFile(void)
{
	if(file_open)
		sf_close(in_file);
}

int SndInFile::getSignalSize()
{
	return signalSize;
}

int SndInFile::getNumberOfSamples()
{
	return numberOfSamples;
}

int SndInFile::getNumberOfChannels()
{
	return numberOfChannels;
}

int SndInFile::getSamplingRate()
{
	return samplingRate;
}

double SndInFile::getSignalTime()
{
	return signalTime;
}


bool SndInFile::getFileStatus()
{
	return file_open;
}

void SndInFile::closeSndFile()
{
	if(file_open)
		sf_close(in_file);

	file_open = false;
}

void SndInFile::print()
{
	sf_count_t 	readSamples;
	sf_count_t	items;

	items		= signalSize;

	double		*sample;
	sample		= new double[signalSize];
	
	readSamples 	= sf_read_double(in_file, sample, items);

	for (int k=0; k<10; k++)
	{
		cout << sample[k] << endl;	
		k++;
	}
}

double* SndInFile::getSignal(int channel = 0)
{
	sf_count_t 	readSamples;
	sf_count_t	items;

	items		= signalSize;

	double		*sample;
	sample		= new double[signalSize];
	
	double		*channelSample=NULL;
	channelSample	= new double[numberOfSamples];

	sf_seek(in_file, 0, SEEK_SET);
	readSamples 	= sf_read_double(in_file, sample, items);

	for (int k=0; k<numberOfSamples; k++)
	{
		channelSample[k] = sample[numberOfChannels*k + channel];
	}
	// delete [] sample;
	return channelSample;	
}
