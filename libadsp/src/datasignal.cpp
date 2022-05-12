// datasignal.cpp
#include "datasignal.h"

/////////////////////////////////////////
// CDataSignal
/////////////////////////////////////////

CDataSignal::CDataSignal()
{
	m_fileName = NULL;
	m_numSignal = 0;
	m_signalSize = 0;
	m_signal = NULL;
	m_norm = NULL;

}

CDataSignal::~CDataSignal()
{
	int i;
	if (m_fileName != NULL) delete [] m_fileName;

    if (m_signal != NULL)
    {
        for (i=0; i<m_numSignal; i++)
        {
            if (m_signal[i] != NULL) delete [] m_signal[i];
        }
        delete [] m_signal;
    }

	if (m_norm != NULL) delete [] m_norm;

    m_numSignal = 0;
    m_signalSize = 0;
    m_type =0;
    m_Fs = 0.0;
}

/////
// Set methods
/////
void CDataSignal::setFileName(const char* file)
{
	if (m_fileName==NULL)  m_fileName = new char[_MAX_PATH];
	strcpy(m_fileName, file);
}

void CDataSignal::setNumSignal(const int numSignal)
{
    m_numSignal = numSignal;
}

void CDataSignal::setSignalSize(const int signalSize)
{
	m_signalSize = signalSize;
}

void CDataSignal::setBlockSize(int blockSize)
{
	m_blockSize = blockSize;
}

void CDataSignal::setBlockHop(int blockHop)
{
	m_blockHop = blockHop;
}

void CDataSignal::setNumBlock(int numBlock)
{
	m_numBlock = numBlock;
}

void CDataSignal::setNorm()
{
	if (m_norm==NULL) m_norm = new double[m_numSignal];
	int i,j;
	for(i=0;i<m_numSignal;i++)
	{
		m_norm[i] = 0.0;
		for(j=0;j<m_signalSize;j++)
		{
			m_norm[i] += m_signal[i][j]*m_signal[i][j];
		}
		m_norm[i] = sqrt(m_norm[i]);
	}
}

void CDataSignal::setNorm(int signum, double norm)
{
	if (m_norm==NULL) m_norm = new double[m_numSignal];

	m_norm[signum-1] = norm;
}

void CDataSignal::setType(int type)
{
    m_type = type;
}

void CDataSignal::setSamplingRate(double sr)
{
    m_Fs = sr;
}

////
// Get methods
////
char* CDataSignal::getFileName()
{
	return m_fileName;
}

int CDataSignal::getNumSignal()
{
    return m_numSignal;
}

int CDataSignal::getSignalSize()
{
	return m_signalSize;
}

double** CDataSignal::getSignal()
{
	return m_signal;
}

int CDataSignal::getBlockSize()
{
	return m_blockSize;
}

int CDataSignal::getBlockHop()
{
	return m_blockHop;
}

int CDataSignal::getNumBlock()
{
	return m_numBlock;
}

double CDataSignal::getNorm(int numSignal)
{
	return m_norm[numSignal];
}

int CDataSignal::getType()
{
	return m_type;
}

double CDataSignal::getSamplingRate()
{
	return m_Fs;
}

/////////////////////////////////
// CComtradeSignal
/////////////////////////////////

//===============================================================
// Function: CComtradeSignal
// Goal: Constructor
// Return:
//===============================================================

CComtradeSignal::CComtradeSignal()
{
	m_config = NULL;
	m_numDigitalSignal =0;
	m_digitalSignal = NULL;
	m_type =1;
}

//===============================================================
// Function: ~CComtradeSignal
// Goal: Destructor
// Return:
//===============================================================

CComtradeSignal::~CComtradeSignal()
{
	int i;

	if (m_digitalSignal != NULL)
	{
		for (i=0;i<m_numDigitalSignal;i++)
		{
			if ( m_digitalSignal[i] !=NULL )
				delete [] m_digitalSignal[i];
		}
		delete [] m_digitalSignal;
	}

	if (m_config!=NULL)
	{
		delete m_config;
	}

}

void CComtradeSignal::setSignal()
{
	if (!setConfig())
	{
		cout << "Cannot open .cfg file" << endl;
	}

	m_signalSize = (int)m_config->GetNumeroAmostrasPorCanal();

	if (!setComtradeSignal())
	{
		cout << "Cannot read .dat file" << endl;
	}

}


//===============================================================
// Function: setConfig
// Goal: Set datafile configuration parameters
// Return:
//===============================================================

int CComtradeSignal::setConfig()
{
	CTradutor *pSource;
	char erro[255];

	pSource = new CNewComtrade;

	if (m_config==NULL) m_config = new CConfigComtrade;

    if (!pSource->SetNomeArquivo( m_fileName, 0, erro))
    {
	    cout << erro << endl;
	    return -1;
    }

    //if (!pSource->ExtrairConfiguracao(&config, erro))
	if (!pSource->ExtrairConfiguracao(m_config, erro))
    {
	    cout << erro << endl;
	    return -1;
	}

	delete pSource;

	return 1;
}

//===============================================================
// Function: setComtradeSignal
// Goal:
// Return:
//===============================================================

int CComtradeSignal::setComtradeSignal()
{
	CTradutor*	pSource;
	void*       vCanais[__num_max_canais];
	size_t      numpts,
		        iNumCanais;
	int         iIndiceCanais[__num_max_canais];
	long        offset;
	char		erro[255];
	long		auxL;

	//auxL= config.GetNumeroAmostrasPorCanal();
	auxL= m_config->GetNumeroAmostrasPorCanal();

	//iNumCanais = config.GetNumCanaisAnalogicos() + config.GetNumCanaisDigitais();
	iNumCanais = m_config->GetNumCanaisAnalogicos() + m_config->GetNumCanaisDigitais();

	m_numSignal = m_config->GetNumCanaisAnalogicos();
	m_numDigitalSignal = m_config->GetNumCanaisDigitais();

	int i;

	for ( i=0;i<iNumCanais;i++)
	{
		iIndiceCanais[i] = i+1;
		vCanais[i] = new float[auxL];
	}

	numpts = (size_t) auxL;
	offset = 0;

	// Modifica a extensao do arquivo para ".dat"
	char *pos;
	pos = strrchr( m_fileName, '.');
	if (pos!=NULL) strcpy( &pos[1], "dat");
	else strcat( m_fileName, ".dat");

	pSource = new CNewComtrade;

	if (!pSource->SetNomeArquivo( m_fileName, 0, erro))
    {
		cout << erro << endl;
	    	return -1;
    }

	// Extrai os valores da amostras dos canais
	if (!pSource->ExtrairDados(	m_config,//&config,
					iIndiceCanais,
					iNumCanais,
					vCanais,
                    numpts,
					offset,
					0,
					erro) )
	{
		cout << erro << endl;
		return 0;
	}


	// =====================================================

	CCanalAnalogico* pChannel;

	int k=0;

	int j;

	int iNumCanaisAnalog = m_config->GetNumCanaisAnalogicos();

	if (m_signal == NULL)
	{
		m_signal = new double*[iNumCanaisAnalog];

		for ( i=0;i<iNumCanaisAnalog;i++)
		{
			pChannel = (CCanalAnalogico*) m_config->GetCanalPtr(i+1);
			float linCoef = pChannel->GetCoefLin();
			float angCoef = pChannel->GetCoefAng();
			m_signal[i] = new double [numpts];
			for (j=0;j<numpts;j++)
			{
				//if (analogSignal[i]==NULL)
				m_signal[i][j] = ( angCoef * ( (float*)vCanais[k] )[j] ) + linCoef ;
			}
			k++;
		}
	}

	// =========================================================

	int iNumCanaisDigitais = m_config->GetNumCanaisDigitais();

	if (m_digitalSignal == NULL)
	{
		m_digitalSignal = new char* [iNumCanaisDigitais];

		for (i=0;i<iNumCanaisDigitais;i++)
		{
			m_digitalSignal[i] = new char [numpts+1];

			for (j=0;j<numpts;j++)
			{
				//if (digitalSignal[i]==NULL)
				if ( ( (char*)vCanais[k] )[j] == 1 )
					m_digitalSignal[i][j] = '1';
				else
					m_digitalSignal[i][j] = '0';
			}
			m_digitalSignal[i][j] = 0;
			k++;
		}
	}

	// =============================================================

	delete pSource;
	for (i=0;i<iNumCanais;i++)
	{
		delete [] vCanais[i];
	}

	return 1;

}

//===============================================================
// Function: cutComtrade
// Goal:
// Return:
//===============================================================

void CComtradeSignal::cutComtrade(char* file, int initSample, int endSample)
{
	CConfigComtrade*	pConfig;
	CTradutor*			pSource;
	CTradutor*			pDest;
	void*				vCanais[__num_max_canais];
	void*				vCanaisMod[__num_max_canais];

	size_t				numpts,
						iNumCanais;
	int					iIndiceCanais[__num_max_canais];
	long				offset;
	char				erro[255];
	long				auxL;
	char				newFileName[_MAX_PATH];
	char*				pos;


	pConfig = new CConfigComtrade;
	pSource = new CNewComtrade;

	if (!pSource->SetNomeArquivo( file, 0, erro))
    {
		cout << erro << endl;
    }

	if (!pSource->ExtrairConfiguracao(pConfig,erro))
    {
		cout << erro << endl;
    }

	auxL= pConfig->GetNumeroAmostrasPorCanal();
	iNumCanais = pConfig->GetNumCanaisAnalogicos() + pConfig->GetNumCanaisDigitais();

	numpts = (size_t) auxL;
	offset = 0;

	int i,k;

	for ( i=0;i<iNumCanais;i++)
	{
		iIndiceCanais[i] = i+1;
		if (i<pConfig->GetNumCanaisAnalogicos())
		{
			vCanais[i] = new float[numpts];
			for (k=0;k<numpts;k++)
			{
				( (float*) vCanais[i])[k] =0.0;

			}
			vCanaisMod[i] = new float[endSample-initSample+1];
			for (k=0;k<endSample-initSample+1;k++)
			{
				( (float*) vCanaisMod[i])[k] =0.0;

			}
		}
		else
		{
			vCanais[i] = new char[numpts];
			for (k=0;k<numpts;k++)
			{
				( (char*) vCanais[i])[k] = '0';

			}
			vCanaisMod[i] = new char[endSample-initSample+1];
			for (k=0;k<endSample-initSample+1;k++)
			{
				( (char*) vCanaisMod[i])[k] = '0';

			}
		}
	}

	pos = strrchr( file, '.');
	if (pos!=NULL) strcpy( &pos[1], "dat");
	else strcat( file, ".dat");

	if (!pSource->SetNomeArquivo( file, 0, erro))
    {
		cout << erro << endl;
    }

	// Extrai os valores da amostras dos canais
	if (!pSource->ExtrairDados(	pConfig,
					iIndiceCanais,
					iNumCanais,
					vCanais,
                    numpts,
					offset,
					0,
					erro) )
	{
		cout << erro << endl;
	}

	int j;
	for ( i=0;i<iNumCanais;i++)
	{
		if (i<pConfig->GetNumCanaisAnalogicos())
		{
			j = 0;
			for (k=initSample;k<endSample+1;k++)
			{
				( (float*) vCanaisMod[i])[j] =	( (float*) vCanais[i])[k];
				j++;

			}

		}
		else
		{
			j = 0;
			for (k=initSample;k<endSample+1;k++)
			{
				( (char*) vCanaisMod[i])[j] = ( (char*) vCanais[i])[k];
				j++;

			}
		}


	}

	pDest = new CNewComtrade;

	// Modificar a freq. de amostragem e o valor da ultima amostra
	//double fs_aux;
	//fs_aux = 60 * ceil( (double)(pConfig->GetTaxa(1))->GetTaxaAmostragem()/(double)60);
	pConfig->AlterarTaxa(1,
						 (pConfig->GetTaxa(1))->GetTaxaAmostragem(),
						 endSample-initSample+1);


	// Modificar o tempo da primeira amostra
	__instante* inst1;
	__instante* inst2;
	inst1 = pConfig->GetInstantePrimeiraAmostra();
	inst2 = pConfig->GetInstanteTrigger();
	inst1->tm_sec = inst1->tm_sec + initSample*(1.0/(double)( (pConfig->GetTaxa(1))->GetTaxaAmostragem() ));
	if ( inst2->tm_sec < inst1->tm_sec)
	{
		inst2->tm_sec = inst1->tm_sec;
		pConfig->SetInstanteTrigger(inst2);
	}

	pConfig->SetInstantePrimeiraAmostra(inst1);

	strcpy(newFileName,file);
	pos = strrchr( newFileName, '.');
	if (pos!=NULL)
	{
		strcpy( &pos[0], "_mod.cfg");
	}
	else
	{
		strcat( newFileName, "_mod.cfg");
	}

	pDest->SetNomeArquivo(newFileName,1,erro);
	pConfig->SetTipoArquivo(0);
	pDest->InserirConfiguracao(pConfig,erro);

	// Insere dat ================================
	pos = strrchr( newFileName, '.');
	strcpy( &pos[1], "dat");

	pDest->SetNomeArquivo(newFileName,1,erro);

	if (!pDest->InserirDados(pConfig,
						iIndiceCanais,
						//m_iNumSBSet,
						iNumCanais,
						vCanaisMod,
						endSample-initSample+1,//(size_t)pConfig->GetNumeroAmostrasPorCanal(),
						0,
						erro))
	{
		cout << erro <<endl;
	}

	delete pSource;
	delete pDest;
	delete pConfig;
	for (i=0;i<iNumCanais;i++)
	{
		delete [] vCanais[i];
		delete [] vCanaisMod[i];
	}

}

//===============================================================
// Function: getConfig
// Goal:
// Return:
//===============================================================

CConfigComtrade* CComtradeSignal::getConfig() const
{
	return m_config;
}


//===============================================================
// Function: getNumberOfSamples
// Goal:
// Return:
//===============================================================

long CComtradeSignal::getNumberOfSamples() const
{
	return m_config->GetNumeroAmostrasPorCanal();
}

//===============================================================
// Function: getNumberOfAnalogChannels
// Goal:
// Return:
//===============================================================

int CComtradeSignal::getNumberOfAnalogChannels() const
{
	return m_config->GetNumCanaisAnalogicos();
}

//===============================================================
// Function: getNumberOfDigitalChannels
// Goal:
// Return:
//===============================================================

int CComtradeSignal::getNumberOfDigitalChannels() const
{
	return m_config->GetNumCanaisDigitais();
}

//===============================================================
// Function: getDigitalSignal
// Goal:
// Return:
//===============================================================

char** CComtradeSignal::getDigitalSignal() const
{
	return m_digitalSignal;
}


//===============================================================
// Function: getSamplingRate
// Goal:
// Return:
//===============================================================

double CComtradeSignal::getSamplingRate(int index) const
{
	return (double)( (m_config->GetTaxa(index))->GetTaxaAmostragem() );
}

//===============================================================
// Function: getFundamentalFrequency
// Goal:
// Return:
//===============================================================

double CComtradeSignal::getFundamentalFrequency() const
{
	return (double)m_config->GetFreqLinha();
}


////////////////////////////////////////////////////////////////////////////////
// CAudioSignal
////////////////////////////////////////////////////////////////////////////////

void CAudioSignal::setSignal()
{
	int i;
	SndInFile audio(m_fileName);

	m_numSignal = audio.getNumberOfChannels();

	m_signalSize = audio.getNumberOfSamples();

	m_Fs = audio.getSamplingRate();


	if (m_signal==NULL)
	{
		m_signal = new double*[m_numSignal];
		for (i=0; i<m_numSignal; i++)
		{
			m_signal[i] = new double[m_signalSize];
			memcpy(m_signal[i],audio.getSignal(i), sizeof(double)*m_signalSize);
		}

	}
	else
	{
		for (i=0; i<m_numSignal; i++)
		{
			memcpy(m_signal[i],audio.getSignal(i), sizeof(double)*m_signalSize);
		}
	}


}

void CAudioSignal::setSignal(double** signal)
{
    int i;
    if (m_signal==NULL)
    {
        m_signal = new double*[m_numSignal];
        for (i=0; i<m_numSignal; i++)
        {
            m_signal[i] = new double[m_signalSize];
            memcpy(m_signal[i],signal[i], sizeof(double)*m_signalSize);
        }
    }
    else
    {
        for (i=0; i<m_numSignal; i++)
        {
            memcpy(m_signal[i],signal[i], sizeof(double)*m_signalSize);
        }
    }
}

void CAudioSignal::saveSignal()
{
    int i,j;

    SndOutFile audio(m_numSignal,static_cast<int>(m_Fs));

    audio.createFile(m_fileName);

    //audio.setChannels(m_numSignal);

    //audio.setSampleRate((int)m_Fs);

    int length = m_signalSize*m_numSignal;
    double* sample = new double[length];
    int k=0;
    for (i=0; i<m_numSignal; i++)
    {
        for (j=0; j<m_signalSize; j++)
        {
            sample[k] = m_signal[i][j];
            k++;
        }
    }

    audio.writeDouble(sample,length);

    delete [] sample;
}


//===============================================================
// Function: getSamplingRate
// Goal:
// Return:
//===============================================================

double CAudioSignal::getSamplingRate() const
{
	return m_Fs;
}

////////////////////////////////////////////////////////////////////////////////
// CNoiseSignal
////////////////////////////////////////////////////////////////////////////////

void CNoiseSignal::setSignal()
{

	long *aux_sem;
	long zero=0;

	aux_sem = (long *)malloc(sizeof(long));

	*aux_sem = Sem_Inicial;

	int i,j;

	if (m_signal==NULL)
	{
		m_signal = new double*[m_numSignal];
		for (i=0; i<m_numSignal; i++)
		{
			m_signal[i] = new double[m_signalSize];
			for (j=0; j<m_signalSize; j++)
			{
				if ( (i==0) && (j==0))
	      			m_signal[i][j] = gausdev(aux_sem);
				else
	      			m_signal[i][j] = gausdev(&zero);
			}
		}
	}


	free(aux_sem);

}

//===============================================================
// Function: getSamplingRate
// Goal:
// Return:
//===============================================================

double CNoiseSignal::getSamplingRate() const
{
	return m_Fs;
}


//===============================================================
// Function: gausdev
// Goal: Funcao para retornar uma amostra dentro de uma distribuicao
//		 gaussiana com media zero e variancia unitaria, utiliza a
//		 funcao ran1 (acima) e as transformacoes de Box-Miller.
//		- Ver:	Numerical Recipits in C, Press, Teukolsky,
//				Vetterling, Flannery, Cambridge University Press, pag. 289
//		sem ->	ponteiro para a semente(long int) que inicia a geracao
//				da sequencia
// Return:
//===============================================================

double CNoiseSignal::gausdev(long int *sem)
{
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;

	if (*sem<0) iset = 0;

	if (iset == 0)
    {
   		do
		{
			v1 = 2.0*ran1(sem) - 1.0;
			v2 = 2.0*ran1(sem) - 1.0;
			rsq = v1*v1 + v2*v2;
        }
        while ((rsq>=1.0)||(rsq==0.0));

		fac = sqrt(-2.0*log((double)rsq)/rsq);
		gset = v1*fac;
		iset = 1;

		return v2*fac;
	}
	else
    {
   		iset = 0;
        return gset;
    }
}

//===============================================================
// Function: ran1
// Goal: Funcao para retornar um amostra de uma distribuicao
//		 uniforme entre 0.0 e 1.0.
//		- Inicializar com idum inteiro negativo.
//		- Depois nao modificar "idum" entre amostras sucessivas de
//        uma sequencia.
//		- Ver:	Numerical Recipes in C, Press, Teukolsky, Vetterling,
//				Flannery, Cambridge University Press, pag. 280
//		sem ->	ponteiro para a semente(long int) que inicia a
//				geracao da sequencia
// Return:
//===============================================================

double CNoiseSignal:: ran1(long int *sem)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*sem <= 0 || !iy)
	{
   		if (-(*sem)< 1)
			*sem =1;
		else
			*sem = -(*sem);

		for (j=NTAB+7; j>=0; j--)
		{
			k=(*sem)/IQ;
			*sem=IA*(*sem-k*IQ)-IR*k;

			if (*sem <0) *sem += IM;
			if (j<NTAB) iv[j] = *sem;
		}

		iy = iv[0];
	}
	k = (*sem)/IQ;
	*sem = IA*(*sem -k*IQ) -IR*k;

	if (*sem <0) *sem +=IM;

	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *sem;

	if ((temp=AM*iy)>RNMX)
		return RNMX;
	else
		return temp;

}

/*
Comentado em 01/21
////////////////////////////////////////////////////////////////////////////////
// CECGSignal
////////////////////////////////////////////////////////////////////////////////

void CECGSignal::setSignal()
{
	int nsig;
	WFDB_Sample* signal;
	WFDB_Siginfo* s;
	WFDB_Frequency freq = (WFDB_Frequency)0;
	char record[_MAX_PATH];
	int j = 0;

	strcpy(record,m_fileName);

	nsig = isigopen(record, NULL, 0);
	if (freq <= (WFDB_Frequency)0) freq = sampfreq(record);

	m_Fs = freq;

	if (nsig < 1) exit(1);

	s = new WFDB_Siginfo[nsig*sizeof(WFDB_Siginfo)];
	signal = new WFDB_Sample[nsig*sizeof(WFDB_Sample)];

	if(isigopen(record,s,nsig)<1) exit(1);

	nsig = isigopen(record, s, nsig);
    m_numSignal = nsig;
    m_signalSize = s->nsamp;
    if (m_signal==NULL)
    {
        m_signal = new double*[m_numSignal];
		for (int j = 0; j<m_signalSize; j++)
		{
			if (getvec(signal)<0) break;
    		for (int i = 0; i < m_numSignal; i++) {
				if (j == 0)	m_signal[i] = new double[m_signalSize];
				m_signal[i][j] = (double)signal[i];
				// cout << m_signal[i][j] << endl;
				// cout << signal[i] << endl;
			}
		}
    }

    else
    {
        for (int j = 0; j<m_signalSize; j++)
		{
			if (getvec(signal)<0) break;
    		for (int i = 0; i < m_numSignal - 1; i++) {
				if (j == 0) m_signal[i] = new double[m_signalSize];
				m_signal[i][j] = (double)signal[i];
			}
		}
    }
    delete[] s;
    delete[] signal;
}

//===============================================================
// Function: getSamplingRate
// Goal:
// Return:
//===============================================================

double CECGSignal::getSamplingRate() const
{
	return m_Fs;
}
*/
