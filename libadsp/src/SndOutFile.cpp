#include "SndOutFile.h"

SndOutFile::SndOutFile(void)
{
    setSampleRate(STD_SAMPLE_RATE);
    setChannels(STD_CHANNELS);
    setFormat(STD_FORMAT);
    setSeekStatus(STD_SEEK_STATUS);
    file_open = false;
}

SndOutFile::SndOutFile(int nchan, int samplerate)
{
    setSampleRate(samplerate);
    setChannels(nchan);
    setFormat(STD_FORMAT);
    setSeekStatus(STD_SEEK_STATUS);
    file_open = false;
}

SndOutFile::~SndOutFile(void)
{

    if(file_open)
        sf_close(out_file);
}

void SndOutFile::setFrame(sf_count_t frame)
{
	info.frames = frame;
}

void SndOutFile::setChannels(int no_channels)
{
	info.channels = no_channels;
}

void SndOutFile::setFormat(int format)
{
	info.format = format;
}

void SndOutFile::setSeekStatus(bool status)
{
	info.seekable = status;
}

void SndOutFile::setSampleRate(int sample_rate)
{
	info.samplerate = sample_rate;
}

int SndOutFile::getSampleRate()
{
	return info.samplerate;
}


void SndOutFile::createFile(char* file_name)
{
	char aux[10]="false";

	if(sf_format_check(&info))
	{
		out_file = sf_open(file_name,SFM_WRITE,&info);
		file_open = true;
		strcpy(aux,"true");
		//aux = "true";
	}
	sf_command(out_file,SFC_SET_CLIPPING,NULL,SF_TRUE);
	sf_command(out_file,SFC_SET_NORM_DOUBLE,NULL,SF_TRUE);

}

void SndOutFile::writeDouble(double *sample)
{
	sf_write_double(out_file, sample,1);

	//sf_write_int(out_file,(int *)sample,1);
}

void SndOutFile::writeDouble(double *sample, int items)
{
    sf_write_double(out_file, sample,items);
}

void SndOutFile::writeDoubleVector(double *vector, int length)
{
    sf_writef_double(out_file,vector,length);
}

bool SndOutFile::getFileStatus()
{
	return file_open;
}

void SndOutFile::closeSndFile()
{
    if(file_open)
        sf_close(out_file);

    file_open = false;
}
