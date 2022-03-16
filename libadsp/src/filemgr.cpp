// filemgr.cpp

#include "filemgr.h"


//=======================
// CFileGenData Class
//=======================

CFileGenData::CFileGenData ()
{
    m_sigType = 0;
    m_dicType =0;
    m_procType =0;
    m_flOMP =0;
    m_flANN =0;
    m_blockHop = 0;
    m_blockSize = 0;
    m_initBlock = 0;
    m_endBlock = 0;
    m_nMaxStep = 0;
    m_coefTempSup = 0;
    m_flSeqBlock = 0;
    m_coefSeqBlock = 0;
    m_flEvalAtomCont = 0;
    m_flFindSupport =0;
    m_flOptDecay = 0;
    m_approxRatioTarget=0.0;
    m_snrTarget=0.0;
    m_flImpStage=0;
}

void CFileGenData::setFileName(char* file)
{
    strcpy(m_fileName, file);
}

void CFileGenData::loadData()
{
    FILE* stream;
    char auxStr[_MAX_PATH];
    char* sNumber;
    
    // open file for reading
    stream = fopen(m_fileName,"r");
    
    // Jump Comment
    fgets( auxStr, _MAX_PATH, stream );
    
    // read the signal type
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_sigType = atoi(sNumber);    // <<<<<<<<<<<<<<<<

    // read the proceeding type
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_procType = atoi(sNumber);    // <<<<<<<<<<<<<<<<

    // read the dictionary type
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_dicType = atoi(sNumber);    // <<<<<<<<<<<<<<<<

    //read usage flag for Orthogonal Matching Pursuit
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_flOMP = atoi(sNumber);  // <<<<<<<<<<<<<<<<

    //read usage flag for artificial neural network
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_flANN = atoi(sNumber);    // <<<<<<<<<<<<<<<<

    // read the block hop
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_blockHop = atoi(sNumber);    // <<<<<<<<<<<<<<<<

    // read the block number size
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_blockSize = atoi(sNumber);    // <<<<<<<<<<<<<<<<

    // read the initial block number
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_initBlock = atoi(sNumber);    // <<<<<<<<<<<<<<<<

    // read the final block number
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_endBlock = atoi(sNumber);     // <<<<<<<<<<<<<<<<

    // read the max. number of iterations
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_nMaxStep = atoi(sNumber);    // <<<<<<<<<<<<<<<<

    // read the coef. for temporal support
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_coefTempSup = atof(sNumber); // <<<<<<<<<<<<<<<<

    // read usage flag for sequence blocks heuristic
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_flSeqBlock = atoi(sNumber);  // <<<<<<<<<<<<<<<<

    // read the coef. for sequence of blocks
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_coefSeqBlock = atof(sNumber);// <<<<<<<<<<<<<<<<
    
    // read usage flag for evaluating atom continuity
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_flEvalAtomCont = atoi(sNumber);  // <<<<<<<<<<<<<<<<
    
    // read usage flag for finding support
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_flFindSupport = atoi(sNumber);  // <<<<<<<<<<<<<<<<
	
	// read usage flag for optimum decaying
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_flOptDecay = atoi(sNumber);  // <<<<<<<<<<<<<<<<
    
    // read the approximation ratio target
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_approxRatioTarget = atof(sNumber);// <<<<<<<<<<<<<<<<
    
    // read the SNR target
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_snrTarget = atof(sNumber);// <<<<<<<<<<<<<<<<
    
    // read flag that print decomp stages
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_flImpStage = atoi(sNumber);  // <<<<<<<<<<<<<<<<
    
    // close file
    fclose(stream);
}

int CFileGenData::getSigType()
{
    return m_sigType;
}

int CFileGenData::getDicType()
{
    return m_dicType;
}

int CFileGenData::getProcType()
{
    return m_procType;
}

int CFileGenData::getFlagOMP()
{
   return m_flOMP;
}

int CFileGenData::getFlagANN()
{
   return m_flANN;
}

int CFileGenData::getBlockHop()
{
    return m_blockHop;
}

int CFileGenData::getBlockSize()
{
    return m_blockSize;
}

int CFileGenData::getInitBlock()
{
    return m_initBlock;
}

int CFileGenData::getEndBlock()
{
    return m_endBlock;
}

int CFileGenData::getNumMaxStep()
{
    return m_nMaxStep;
}

double CFileGenData::getCoefTempSup()
{
    return m_coefTempSup;
}
 
int CFileGenData::getFlagSeqBlock()
{
   return m_flSeqBlock;
}

double CFileGenData::getCoefSeqBlock()
{
    return m_coefSeqBlock;
}

int CFileGenData::getFlagEvalAtomCont()
{
    return m_flEvalAtomCont;   
}

int CFileGenData::getFlagFindSupport()
{
    return m_flFindSupport;   
}

int CFileGenData::getFlagOptDecay()
{
    return m_flOptDecay;   
}

double CFileGenData::getApproxRatioTarget()
{
	return m_approxRatioTarget;
}

double CFileGenData::getSNRTarget()
{
	return m_snrTarget;
}

int CFileGenData::getPrintDecompStage()
{
    return m_flImpStage;   
}

//=======================
// CFileProceeding Class
//=======================

CFileProceeding::CFileProceeding ()
{
    m_flDecomp=0;
 	m_flGroupBlock=0;
 	m_flRDTableDf=0;
 	m_flRDOpCurve=0;
 	m_flRDOpCurveAmpRange=0;
 	m_flEncAmpRange=0;
 	m_flEncOpCurve=0;
 	m_flEncAmpRangeArithCod=0;
 	m_flDecAmpRangeArithCod = 0;
}

void CFileProceeding::setFileName(char* file)
{
    strcpy(m_fileName, file);
}

void CFileProceeding::loadData()
{
    FILE* stream;
    char auxStr[_MAX_PATH];
    char* sNumber;
    
    // open file for reading
    stream = fopen(m_fileName,"r");
    
    // Jump Comment
    fgets( auxStr, _MAX_PATH, stream );
    
    // read decomp switch
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_flDecomp = atoi(sNumber);    // <<<<<<<<<<<<<<<<

    // read group blocks switch
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_flGroupBlock = atoi(sNumber);    // <<<<<<<<<<<<<<<<

    // read RD mean Atom distortion switch
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_flRDTableDf = atoi(sNumber);    // <<<<<<<<<<<<<<<<

    // read Rd Operational curve switch
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_flRDOpCurve = atoi(sNumber);    // <<<<<<<<<<<<<<<<
    
    // read Rd Operational curve switch (Amplitude Range Separation)
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_flRDOpCurveAmpRange = atoi(sNumber);    // <<<<<<<<<<<<<<<<
    
    // read encode using Df table switch
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_flEncAmpRange = atoi(sNumber);    // <<<<<<<<<<<<<<<<
    
    // read encode using operational curve switch
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_flEncOpCurve = atoi(sNumber);    // <<<<<<<<<<<<<<<<
    
    // read encode using Df table switch (with arithmetic Coding)
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_flEncAmpRangeArithCod = atoi(sNumber);    // <<<<<<<<<<<<<<<<
    
    // read DEcode using Df table switch (with arithmetic Coding)
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_flDecAmpRangeArithCod = atoi(sNumber);    // <<<<<<<<<<<<<<<<

    // read DEcode using Df table switch (with arithmetic Coding)
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_flPll = atoi(sNumber);    // <<<<<<<<<<<<<<<<
    
    // close file
    fclose(stream);
}

int CFileProceeding::getDecomp()
{
    return m_flDecomp;
}

int CFileProceeding::getGroupBlock()
{
    return m_flGroupBlock;
}

int CFileProceeding::getRDTableDf()
{
    return m_flRDTableDf;
}

int CFileProceeding::getRDOpCurve()
{
    return m_flRDOpCurve;
}

int CFileProceeding::getRDOpCurveAmpRange()
{
    return m_flRDOpCurveAmpRange;
}

int CFileProceeding::getEncAmpRange()
{
    return m_flEncAmpRange;
}

int CFileProceeding::getEncOpCurve()
{
    return m_flEncOpCurve;
}

int CFileProceeding::getEncAmpRangeArithCod()
{
    return m_flEncAmpRangeArithCod;
}

int CFileProceeding::getDecAmpRangeArithCod()
{
	return m_flDecAmpRangeArithCod;
}

int CFileProceeding::getPll()
{
	return m_flPll;
}


//=======================
// CFileGroupBlock Class
//=======================

CFileGroupBlock::CFileGroupBlock ()
{
    m_numBlockPerGroup = 1;
}

void CFileGroupBlock::setFileName(char* file)
{
    strcpy(m_fileName, file);
}

void CFileGroupBlock::loadData()
{
    FILE* stream;
    char auxStr[_MAX_PATH];
    char* sNumber;
    
    // open file for reading
    stream = fopen(m_fileName,"r");
    
    // Jump Comment
    fgets( auxStr, _MAX_PATH, stream );
    
    // read number of blocks per group for encoding
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_numBlockPerGroup = atoi(sNumber);  // <<<<<<<<<<<<<<<<
    
    // close file
    fclose(stream);
}

int CFileGroupBlock::getNumBlockPerGroup()
{
	return m_numBlockPerGroup;
}



//=======================
// CFileBlockRange Class
//=======================

CFileBlockRange::CFileBlockRange ()
{
    m_numRange = 0;
    m_pInitBlock = NULL;
    m_pFinalBlock = NULL;
}

CFileBlockRange::~CFileBlockRange()
{
    if (m_pInitBlock !=NULL)
    {
        delete [] m_pInitBlock;
    }

    if (m_pFinalBlock !=NULL)
    {
        delete [] m_pFinalBlock;
    }
}

void CFileBlockRange::setFileName(char* file)
{
    strcpy(m_fileName, file);
}

void CFileBlockRange::loadData()
{
    FILE* stream;
    char auxStr[_MAX_PATH];
    char* sNumber;
    int numAux;
    
    // open file for reading
    stream = fopen(m_fileName,"r");

    fgets( auxStr, _MAX_PATH, stream );
    fgets( auxStr, _MAX_PATH, stream );

    while(!feof(stream))
    {
        fgets( auxStr, _MAX_PATH, stream );
        auxStr[5] = NULL;
        sNumber = &auxStr[0];
        numAux= atoi(sNumber);
        if (numAux==99999) break;
        m_numRange++;
    }

    m_pInitBlock = new int[m_numRange];
    m_pFinalBlock = new int[m_numRange];

    rewind(stream);
    fgets( auxStr, _MAX_PATH, stream );
    fgets( auxStr, _MAX_PATH, stream );
    // read the block ranges
    int i;
    for (i=0;i<m_numRange;i++)
    {
        fgets( auxStr, _MAX_PATH, stream );
        sNumber = &auxStr[11];
        m_pFinalBlock[i]= atoi(sNumber);
        auxStr[5] = NULL;
        sNumber = &auxStr[0];
        m_pInitBlock[i]= atoi(sNumber);
        printf("-Block range %d; init block %d and final block %d\n",i+1,m_pInitBlock[i],m_pFinalBlock[i]);
    }

    fclose(stream);

}

int CFileBlockRange::getNumRange()
{
    return m_numRange;
}

int* CFileBlockRange::getPInitBlock()
{
    return m_pInitBlock;
}

int* CFileBlockRange::getPFinalBlock()
{
    return m_pFinalBlock;
}


//=======================
// CFileDictionary Class
//=======================


CFileDictionary::CFileDictionary ()
{
    m_numDicBlock = 0;
    m_dicBlock = NULL;
    m_numFreq =0;
}

CFileDictionary::~CFileDictionary()
{
    if (m_dicBlock!=NULL) 
    {
        delete[] m_dicBlock;
    }
}

void CFileDictionary::setFileName(char* file)
{
   strcpy(m_fileName, file);
}

void CFileDictionary::loadData()
{
    strtFileDictionary* dicBlockAux=NULL;

    FILE* dicFile;
    char auxStr[_MAX_PATH];
    char auxStr2[_MAX_PATH];
    char* sNumber;
    // Opening the "dictionary.dat" file
    dicFile = fopen(m_fileName,"r");
    // read comment line
    fgets( auxStr, _MAX_PATH, dicFile );
    fgets( auxStr, _MAX_PATH, dicFile );
    while(!feof(dicFile))
    {
        if ( m_dicBlock==NULL )
        {
            m_dicBlock = new strtFileDictionary[1];
        }
        dicBlockAux = new strtFileDictionary[m_numDicBlock + 1];

        memcpy ( dicBlockAux,
                 m_dicBlock,
                 m_numDicBlock*sizeof ( strtFileDictionary ) );

        // read line
        fgets( auxStr, _MAX_PATH, dicFile );
        // read scale
        strcpy(auxStr2,auxStr);
        auxStr2[5] = NULL;
        sNumber = &auxStr[0];
        dicBlockAux[m_numDicBlock].scale = atof(sNumber);    // <<<<<<<<<<<<<<<<
        //cout << dicBlockAux[m_numDicBlock].scale << endl;
        if (dicBlockAux[m_numDicBlock].scale==99999) break;

        // read delta_tau
        strcpy(auxStr2,auxStr);
        auxStr2[14] = NULL;
        sNumber = &auxStr[9];
        dicBlockAux[m_numDicBlock].delta_tau = atoi(sNumber);    // <<<<<<<<<<<<<<<<
        //cout << dicBlockAux[m_numDicBlock].delta_tau << endl;
        // read fdiscrtype
        strcpy(auxStr2,auxStr);
        auxStr2[24] = NULL;
        sNumber = &auxStr[22];
        dicBlockAux[m_numDicBlock].fdiscrtype = atoi(sNumber);    // <<<<<<<<<<<<<<<<
        //cout << dicBlockAux[m_numDicBlock].fdiscrtype << endl;        
        // read freqi
        strcpy(auxStr2,auxStr);
        auxStr2[36] = NULL;
        sNumber = &auxStr[26];
        dicBlockAux[m_numDicBlock].freqi = atof(sNumber);    // <<<<<<<<<<<<<<<<
        //cout << dicBlockAux[m_numDicBlock].freqi << endl;        
        // read freqf
        strcpy(auxStr2,auxStr);
        auxStr2[50] = NULL;
        sNumber = &auxStr[40];
        dicBlockAux[m_numDicBlock].freqf = atof(sNumber);    // <<<<<<<<<<<<<<<<    
        // read dicType
        strcpy(auxStr2,auxStr);
        auxStr2[59] = NULL;
        sNumber = &auxStr[57];
        dicBlockAux[m_numDicBlock].dicType = atof(sNumber);    // <<<<<<<<<<<<<<<<    
        //read Decay
        strcpy(auxStr2,auxStr);
        auxStr2[70] = NULL;
        sNumber = &auxStr[65];
        dicBlockAux[m_numDicBlock].decay = atof(sNumber);    // <<<<<<<<<<<<<<<<  
        //read Rise
        strcpy(auxStr2,auxStr);
        auxStr2[81] = NULL;
        sNumber = &auxStr[76];
        dicBlockAux[m_numDicBlock].rise = atof(sNumber);    // <<<<<<<<<<<<<<<<  
        

        //cout << dicBlockAux[m_numDicBlock].freqf << endl;
        if ( m_dicBlock!=NULL )
        {
            delete [] m_dicBlock;
        }

        m_dicBlock = new strtFileDictionary[m_numDicBlock + 1];

        memcpy ( m_dicBlock,dicBlockAux, ( m_numDicBlock+1 ) *sizeof ( strtFileDictionary ) );

        if (dicBlockAux!=NULL)
        {
            delete [] dicBlockAux;
        }

        m_numDicBlock++;
    }
    fclose(dicFile);
    
    if (m_numDicBlock!=0)
    {
		// Considering that all scales have the same freqi and freqf
		if (dicBlockAux[0].fdiscrtype==1) // linear
		{
			m_numFreq = (int)ceil(dicBlockAux[0].freqf/dicBlockAux[0].freqi);     
		}
		else if (dicBlockAux[0].fdiscrtype==2) // geometric with quarter-tone discretization
		{
			// '24' is referred to a quarter of tone, '12' is referred to half tone
			m_numFreq = (int)ceil(24 * ( log10(dicBlockAux[0].freqf/dicBlockAux[0].freqi)/log10(2.0) ) )+1; 
		}
	}
    
}

void CFileDictionary::printToScreen()
{
    printf("Loaded dictionary:\n");
    for (int i=0; i<m_numDicBlock; i++)
    {
        printf("%10.4f %5d %2d %10.4f %10.4f %2d %5.3f %5.3f \n",  m_dicBlock[i].scale,
                                    m_dicBlock[i].delta_tau,
                                    m_dicBlock[i].fdiscrtype,
                                    m_dicBlock[i].freqi,
                                    m_dicBlock[i].freqf,
                                    m_dicBlock[i].dicType,
                                    m_dicBlock[i].decay, 
                                    m_dicBlock[i].rise);
    }
}

int CFileDictionary::getNumDicBlock()
{
    return m_numDicBlock;
}


strtFileDictionary* CFileDictionary::getDicBlock()
{
    return m_dicBlock;
}

int CFileDictionary::getNumFreq()
{
	return m_numFreq;
}

double CFileDictionary::getFreqi(int blockInd)
{
	return m_dicBlock[blockInd].freqi;
}

double CFileDictionary::getFreqf(int blockInd)
{
	return m_dicBlock[blockInd].freqf;
}

double CFileDictionary::getFDiscrType(int blockInd)
{
	return m_dicBlock[blockInd].fdiscrtype;
}

int CFileDictionary::getDicType(int blockInd)
{
    return m_dicBlock[blockInd].dicType;
}

double CFileDictionary::getDecay(int blockInd)
{
    return m_dicBlock[blockInd].decay;
}

double CFileDictionary::getRise(int blockInd)
{
    return m_dicBlock[blockInd].rise;
}


//=======================
// CFileRDBitRange Class
//=======================

CFileRDBitRange::CFileRDBitRange ()
{
    init_nbit_amp=0;
    end_nbit_amp=0;
    delta_nbit_amp=0;
    init_nbit_rho=0;
    end_nbit_rho=0;
    delta_nbit_rho=0;
    init_nbit_phase=0;
    end_nbit_phase=0;
    delta_nbit_phase=0;
}


void CFileRDBitRange::setFileName(char* file)
{
   strcpy(m_fileName, file);
}

void CFileRDBitRange::loadData()
{

    FILE* file;
    char auxStr[_MAX_PATH];
    char auxStr2[_MAX_PATH];
    char* sNumber;
    // Opening the "rdbitrange.dat" file
    file = fopen(m_fileName,"r");
    if (file==NULL)
    {
    	printf("Arquivo %s inexistente!!!\n", m_fileName);
    	exit(0);
    }
    
    // read comment line
    fgets( auxStr, _MAX_PATH, file );
    fgets( auxStr, _MAX_PATH, file );
    // read line
    fgets( auxStr, _MAX_PATH, file );
    // read init nbit
    strcpy(auxStr2,auxStr);
    auxStr2[10] = NULL;
    sNumber = &auxStr[5];
    init_nbit_amp= atoi(sNumber);
    // read end nbit
    strcpy(auxStr2,auxStr);
    auxStr2[16] = NULL;
    sNumber = &auxStr[11];
    end_nbit_amp= atoi(sNumber);
    // read delta nbit
    strcpy(auxStr2,auxStr);
    auxStr2[22] = NULL;
    sNumber = &auxStr[17];
    delta_nbit_amp= atoi(sNumber);
    //////////////////////
    // read line
    fgets( auxStr, _MAX_PATH, file );
    // read init nbit
    strcpy(auxStr2,auxStr);
    auxStr2[10] = NULL;
    sNumber = &auxStr[5];
    init_nbit_rho= atoi(sNumber);
    // read end nbit
    strcpy(auxStr2,auxStr);
    auxStr2[16] = NULL;
    sNumber = &auxStr[11];
    end_nbit_rho= atoi(sNumber);
    // read delta nbit
    strcpy(auxStr2,auxStr);
    auxStr2[22] = NULL;
    sNumber = &auxStr[17];
    delta_nbit_rho= atoi(sNumber);
    ////////////////////////
    // read line
    fgets( auxStr, _MAX_PATH, file );
    // read init nbit
    strcpy(auxStr2,auxStr);
    auxStr2[10] = NULL;
    sNumber = &auxStr[5];
    init_nbit_phase= atoi(sNumber);
    // read end nbit
    strcpy(auxStr2,auxStr);
    auxStr2[16] = NULL;
    sNumber = &auxStr[11];
    end_nbit_phase= atoi(sNumber);
    // read delta nbit
    strcpy(auxStr2,auxStr);
    auxStr2[22] = NULL;
    sNumber = &auxStr[17];
    delta_nbit_phase= atoi(sNumber);
    
    // DFTable File Name
    fgets( auxStr, _MAX_PATH, file ); // comments
    fgets( m_DfTableFName, _MAX_PATH, file );
    char* ptr = strchr(m_DfTableFName,'\n');
    ptr[0] = '\0';
    
    //int strlength = strlen(m_DfTableFName); 
    //m_DfTableFName[strlength-1] = '\0';
    
    //////////////////////////
    fclose(file);
}

int CFileRDBitRange::getInitNbitAmp()
{
	return init_nbit_amp;
}

int CFileRDBitRange::getEndNbitAmp()
{
	return end_nbit_amp;
}

int CFileRDBitRange::getDeltaNbitAmp()
{
	return delta_nbit_amp;
}
int CFileRDBitRange::getInitNbitRho()
{
	return init_nbit_rho;
}

int CFileRDBitRange::getEndNbitRho()
{
	return end_nbit_rho;
}

int CFileRDBitRange::getDeltaNbitRho()
{
	return delta_nbit_rho;
}
int CFileRDBitRange::getInitNbitPhase()
{
	return init_nbit_phase;
}

int CFileRDBitRange::getEndNbitPhase()
{
	return end_nbit_phase;
}

int CFileRDBitRange::getDeltaNbitPhase()
{
	return delta_nbit_phase;
}

char* CFileRDBitRange::getDfTableFName()
{
	return m_DfTableFName;
}

//=======================
// CFileDfTable Class
//=======================

CFileDfTable::CFileDfTable ()
{
	m_Nqrho=0;
	m_Nqphi=0;
	m_qrho=NULL;
	m_qphi=NULL;
	m_DfTable=NULL;
	m_min_nbit_amp=0;
 	m_max_nbit_amp=0;
}
CFileDfTable::~CFileDfTable ()
{
	int nb_amp,i_qrho,i_qphi;
	for (nb_amp=m_min_nbit_amp;nb_amp<=m_max_nbit_amp;nb_amp++)
	{
		for (i_qrho=0;i_qrho<m_Nqrho;i_qrho++)
		{
			for (i_qphi=0;i_qphi<m_Nqphi;i_qphi++)
			{	
				delete [] m_DfTable[nb_amp-m_min_nbit_amp][i_qrho][i_qphi];
			}
			delete [] m_DfTable[nb_amp-m_min_nbit_amp][i_qrho];
		}
		delete [] m_DfTable[nb_amp-m_min_nbit_amp];
	}
	delete [] m_DfTable;
	////
	if (m_qrho!=NULL)
	{
		delete [] m_qrho;
	}
	/////
	if (m_qphi!=NULL)
	{
		delete [] m_qphi;
	}
}

void CFileDfTable::setFileName(char* file)
{
    strcpy(m_fileName, file);
}

void CFileDfTable::loadData()
{
    FILE* iobin;
    
    int dummyint;
    double dummydouble;
    int i;
    
    // open file for reading
    cout << m_fileName << endl;
    if ( (iobin = fopen(m_fileName,"rb"))==NULL)
    {
    	printf("File %s not found.\n",m_fileName);
    }
   
	fread(&dummyint, sizeof(int), 1, iobin);
	m_Nqrho = dummyint;
	
	cout << m_Nqrho << endl; 
	
	fread(&dummyint, sizeof(int), 1, iobin);
	m_Nqphi = dummyint;
	
	cout << m_Nqphi << endl;
	
	m_qrho = new double[m_Nqrho];
	for (i=0;i<m_Nqrho;i++)
	{
		fread(&dummydouble, sizeof(double), 1, iobin);
		m_qrho[i] = dummydouble;	
		//cout << m_qrho[i] << "; " ;
	}
	//cout << endl;
	
	m_qphi = new double[m_Nqphi];	
    for (i=0;i<m_Nqphi;i++)
	{
		fread(&dummydouble, sizeof(double), 1, iobin);
		m_qphi[i] = dummydouble;	
		//cout << m_qphi[i] << "; " ;
	}
	//cout << endl;
	
	m_min_nbit_amp=1;
	m_max_nbit_amp=16;
	int nb_amp, iAmpRange;
	int i_qrho,i_qphi;
	
	m_DfTable = new double***[m_max_nbit_amp-m_min_nbit_amp+1];
    for (nb_amp=m_min_nbit_amp;nb_amp<=m_max_nbit_amp;nb_amp++)
	{
		m_DfTable[nb_amp-m_min_nbit_amp] = new double**[m_Nqrho];
		for (i_qrho=0;i_qrho<m_Nqrho;i_qrho++)
		{
			m_DfTable[nb_amp-m_min_nbit_amp][i_qrho] = new double*[m_Nqphi];
			for (i_qphi=0;i_qphi<m_Nqphi;i_qphi++)
			{
				m_DfTable[nb_amp-m_min_nbit_amp][i_qrho][i_qphi]= new double[nb_amp];
			}
		}
	}
	
	for (nb_amp=m_min_nbit_amp;nb_amp<=m_max_nbit_amp;nb_amp++)
	{
		for (iAmpRange=0;iAmpRange<nb_amp;iAmpRange++)
		{
			for (i_qrho=0;i_qrho<m_Nqrho;i_qrho++)
			{
				for (i_qphi=0;i_qphi<m_Nqphi;i_qphi++)
				{	
					fread(&dummydouble, sizeof(double), 1, iobin);
					m_DfTable[nb_amp-m_min_nbit_amp][i_qrho][i_qphi][iAmpRange] = dummydouble;
					//cout <<m_DfTable[nb_amp-m_min_nbit_amp][i_qrho][i_qphi][iAmpRange] << "; " ;
				}
			}
		}
	}
	//cout<<endl;
	
    // close file
    fclose(iobin);
}

int CFileDfTable::getNumQRho()
{
	return m_Nqrho;
}

int CFileDfTable::getNumQPhi()
{
	return m_Nqphi;
}

double* CFileDfTable::getQRho()
{
	return m_qrho;
}

double* CFileDfTable::getQPhi()
{
	return m_qphi;
}

double**** CFileDfTable::getDfTable()
{
	return m_DfTable;
}

int CFileDfTable::getMinNbitAmp()
{
	return m_min_nbit_amp;
}

int CFileDfTable::getMaxNbitAmp()
{
	return m_max_nbit_amp;
}

//=======================
// CFileDecomp Class
//=======================

CFileDecomp::CFileDecomp ()
{
    m_sigType = 0;
    m_dicType =0;
    m_mpType =0;
    m_flOMP =0;
    m_flANN =0;
    m_nMaxStep = 0;
    m_coefTempSup = 0;
    m_flSeqBlock = 0;
    m_coefSeqBlock = 0;
    m_flEvalAtomCont = 0;
    m_flFindSupport =0;
    m_flOptDecay = 0;
    m_approxRatioTarget=0.0;
    m_snrTarget=0.0;
    m_flImpStage=0;
    m_flRDopt = 0;
}

void CFileDecomp::setFileName(char* file)
{
    strcpy(m_fileName, file);
}


void CFileDecomp::loadData()
{
    FILE* stream;
    char auxStr[_MAX_PATH];
    char* sNumber;
    
    // open file for reading
    stream = fopen(m_fileName,"r");
    
    // Jump Comment
    fgets( auxStr, _MAX_PATH, stream );
    
    // read the signal type
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_sigType = atoi(sNumber);    // <<<<<<<<<<<<<<<<

    // read the MP type
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_mpType = atoi(sNumber);    // <<<<<<<<<<<<<<<<

    // read the dictionary type
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_dicType = atoi(sNumber);    // <<<<<<<<<<<<<<<<

    //read usage flag for Orthogonal Matching Pursuit
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_flOMP = atoi(sNumber);  // <<<<<<<<<<<<<<<<

    // read usage flag for artificial neural network
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_flANN = atoi(sNumber);  // <<<<<<<<<<<<<<<<

    // read the max. number of iterations
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_nMaxStep = atoi(sNumber);    // <<<<<<<<<<<<<<<<

    // read the coef. for temporal support
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_coefTempSup = atof(sNumber); // <<<<<<<<<<<<<<<<

    // read usage flag for sequence blocks heuristic
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_flSeqBlock = atoi(sNumber);  // <<<<<<<<<<<<<<<<

    // read the coef. for sequence of blocks
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_coefSeqBlock = atof(sNumber);// <<<<<<<<<<<<<<<<
    
    // read usage flag for evaluating atom continuity
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_flEvalAtomCont = atoi(sNumber);  // <<<<<<<<<<<<<<<<
    
    // read usage flag for finding support
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_flFindSupport = atoi(sNumber);  // <<<<<<<<<<<<<<<<
    
    // read usage flag for optimum decaying
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_flOptDecay = atoi(sNumber);  // <<<<<<<<<<<<<<<<
    
    // read the approximation ratio target
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_approxRatioTarget = atof(sNumber);// <<<<<<<<<<<<<<<<
    
    // read the SNR target
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_snrTarget = atof(sNumber);// <<<<<<<<<<<<<<<<
    
    // read flag that print decomp stages
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_flImpStage = atoi(sNumber);  // <<<<<<<<<<<<<<<<
    
    // read flag that activates rate distortion optimization
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_flRDopt = atoi(sNumber);  // <<<<<<<<<<<<<<<<
    


    // close file
    fflush(stream);
    fclose(stream);
}

void CFileDecomp::echo()
{
    cout << "======================================" << endl;
    cout << " ECHO OF DECOMP FILE" << endl;
    cout << "======================================" << endl;
    
    cout << "File name: " << m_fileName << endl;
    cout << "Signal type: " << m_sigType << endl;
    cout << "MP type: " << m_mpType << endl;
    cout << "m_flOMP: " << m_flOMP << endl;
    cout << "m_flANN: " << m_flANN << endl;
    cout << "Dictionary type: " << m_dicType << endl;
    cout << "Max number of steps: " << m_nMaxStep << endl;
    cout << "Coef suporte temporal: " << m_coefTempSup << endl;
    cout << "m_flSeqBlock: " << m_flSeqBlock << endl;
    cout << "m_coefSeqBlock: " << m_coefSeqBlock << endl;
    cout << "m_flEvalAtomCont: " << m_flEvalAtomCont << endl;
    cout << "m_flFindSupport: " << m_flFindSupport << endl;
    cout << "m_flOptDecay: " << m_flOptDecay << endl;
    cout << "m_approxRatioTarget: " << m_approxRatioTarget << endl;
    cout << "m_snrTarget: " << m_snrTarget << endl;
    cout << "m_flImpStage: " << m_flImpStage << endl;
    cout << "m_flRDopt: " << m_flRDopt << endl;
    
    cout << "======================================" << endl;
    cout << "======================================" << endl;
}

int CFileDecomp::getSigType()
{
    return m_sigType;
}

int CFileDecomp::getDicType()
{
    return m_dicType;
}

int CFileDecomp::getMPType()
{
    return m_mpType;
}

int CFileDecomp::getFlagOMP()
{
   return m_flOMP;
}

int CFileDecomp::getFlagANN()
{
   return m_flANN;
}

int CFileDecomp::getNumMaxStep()
{
    return m_nMaxStep;
}

double CFileDecomp::getCoefTempSup()
{
    return m_coefTempSup;
}
 
int CFileDecomp::getFlagSeqBlock()
{
   return m_flSeqBlock;
}

double CFileDecomp::getCoefSeqBlock()
{
    return m_coefSeqBlock;
}

int CFileDecomp::getFlagEvalAtomCont()
{
    return m_flEvalAtomCont;   
}

int CFileDecomp::getFlagFindSupport()
{
    return m_flFindSupport;   
}

int CFileDecomp::getFlagOptDecay()
{
    return m_flOptDecay;   
}

double CFileDecomp::getApproxRatioTarget()
{
    return m_approxRatioTarget;
}

double CFileDecomp::getSNRTarget()
{
    return m_snrTarget;
}

int CFileDecomp::getPrintDecompStage()
{
    return m_flImpStage;   
}

int CFileDecomp::getFlagRDopt()
{
    return m_flRDopt;   
}

//=======================
// CFileDecompBlockRange Class
//=======================

CFileDecompBlockRange::CFileDecompBlockRange ()
{
    m_blockHop = 0;
    m_blockSize = 0;
    m_initBlock = 0;
    m_endBlock = 0;
    m_numRange = 0;
    m_pInitBlock = NULL;
    m_pFinalBlock = NULL;
}

CFileDecompBlockRange::~CFileDecompBlockRange()
{
    if (m_pInitBlock !=NULL)
    {
        delete [] m_pInitBlock;
    }

    if (m_pFinalBlock !=NULL)
    {
        delete [] m_pFinalBlock;
    }
}

void CFileDecompBlockRange::setFileName(char* file)
{
    strcpy(m_fileName, file);
}

void CFileDecompBlockRange::loadData()
{
    FILE* stream;
    char auxStr[_MAX_PATH];
    char* sNumber;
    int numAux;
    
    // open file for reading
    stream = fopen(m_fileName,"r");
    
    // Jump Comment
    fgets( auxStr, _MAX_PATH, stream );
    
     // read the block hop
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_blockHop = atoi(sNumber);    // <<<<<<<<<<<<<<<<

    // read the block number size
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_blockSize = atoi(sNumber);    // <<<<<<<<<<<<<<<<

    // read the initial block number
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_initBlock = atoi(sNumber);    // <<<<<<<<<<<<<<<<

    // read the final block number
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_endBlock = atoi(sNumber);     // <<<<<<<<<<<<<<<<

    // Jump Comment
    fgets( auxStr, _MAX_PATH, stream );
    fgets( auxStr, _MAX_PATH, stream );

    while(!feof(stream))
    {
        fgets( auxStr, _MAX_PATH, stream );
        auxStr[5] = NULL;
        sNumber = &auxStr[0];
        numAux= atoi(sNumber);
        if (numAux==99999) break;
        m_numRange++;
    }

    m_pInitBlock = new int[m_numRange];
    m_pFinalBlock = new int[m_numRange];

    rewind(stream);
    fgets( auxStr, _MAX_PATH, stream );
    fgets( auxStr, _MAX_PATH, stream );
    fgets( auxStr, _MAX_PATH, stream );
    fgets( auxStr, _MAX_PATH, stream );
    fgets( auxStr, _MAX_PATH, stream );
    fgets( auxStr, _MAX_PATH, stream );
    fgets( auxStr, _MAX_PATH, stream );
    // read the block ranges
    int i;
    for (i=0;i<m_numRange;i++)
    {
        fgets( auxStr, _MAX_PATH, stream );
        sNumber = &auxStr[11];
        m_pFinalBlock[i]= atoi(sNumber);
        auxStr[5] = NULL;
        sNumber = &auxStr[0];
        m_pInitBlock[i]= atoi(sNumber);
        printf("-Block range %d; init block %d and final block %d\n",i+1,m_pInitBlock[i],m_pFinalBlock[i]);
    }

    fclose(stream);

}

int CFileDecompBlockRange::getBlockHop()
{
    return m_blockHop;
}

int CFileDecompBlockRange::getBlockSize()
{
    return m_blockSize;
}

int CFileDecompBlockRange::getInitBlock()
{
    return m_initBlock;
}

int CFileDecompBlockRange::getEndBlock()
{
    return m_endBlock;
}

int CFileDecompBlockRange::getNumRange()
{
    return m_numRange;
}

int* CFileDecompBlockRange::getPInitBlock()
{
    return m_pInitBlock;
}

int* CFileDecompBlockRange::getPFinalBlock()
{
    return m_pFinalBlock;
}



CFileEncode::CFileEncode ()
{
    m_flGroupBlock=0;
    m_numBlockPerGroup = 1;
    m_flRDTableDf=0;
    m_flRDOpCurve=0;
    m_flRDOpCurveAmpRange=0;
    m_flEncAmpRange=0;
    m_flEncOpCurve=0;
    m_flEncAmpRangeArithCod=0;
    m_flDecAmpRangeArithCod = 0;
}

void CFileEncode::setFileName(char* file)
{
    strcpy(m_fileName, file);
}

void CFileEncode::loadData()
{
    FILE* stream;
    char auxStr[_MAX_PATH];
    char* sNumber;
    
    // open file for reading
    stream = fopen(m_fileName,"r");
    
    // Jump Comment
    fgets( auxStr, _MAX_PATH, stream );

    // read group blocks switch
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_flGroupBlock = atoi(sNumber);    // <<<<<<<<<<<<<<<<
    
    // read number of blocks per group for encoding
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_numBlockPerGroup = atoi(sNumber);  // <<<<<<<<<<<<<<<<

    // read RD mean Atom distortion switch
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_flRDTableDf = atoi(sNumber);    // <<<<<<<<<<<<<<<<

    // read Rd Operational curve switch
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_flRDOpCurve = atoi(sNumber);    // <<<<<<<<<<<<<<<<
    
    // read Rd Operational curve switch (Amplitude Range Separation)
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_flRDOpCurveAmpRange = atoi(sNumber);    // <<<<<<<<<<<<<<<<
    
    // read encode using Df table switch
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_flEncAmpRange = atoi(sNumber);    // <<<<<<<<<<<<<<<<
    
    // read encode using operational curve switch
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_flEncOpCurve = atoi(sNumber);    // <<<<<<<<<<<<<<<<
    
    // read encode using Df table switch (with arithmetic Coding)
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_flEncAmpRangeArithCod = atoi(sNumber);    // <<<<<<<<<<<<<<<<
    
    // read DEcode using Df table switch (with arithmetic Coding)
    fgets( auxStr, _MAX_PATH, stream );
    auxStr[50] = NULL;
    sNumber = &auxStr[40];
    m_flDecAmpRangeArithCod = atoi(sNumber);    // <<<<<<<<<<<<<<<<
    
    // close file
    fclose(stream);
}

int CFileEncode::getGroupBlock()
{
    return m_flGroupBlock;
}

int CFileEncode::getNumBlockPerGroup()
{
    return m_numBlockPerGroup;
}

int CFileEncode::getRDTableDf()
{
    return m_flRDTableDf;
}

int CFileEncode::getRDOpCurve()
{
    return m_flRDOpCurve;
}

int CFileEncode::getRDOpCurveAmpRange()
{
    return m_flRDOpCurveAmpRange;
}

int CFileEncode::getEncAmpRange()
{
    return m_flEncAmpRange;
}

int CFileEncode::getEncOpCurve()
{
    return m_flEncOpCurve;
}

int CFileEncode::getEncAmpRangeArithCod()
{
    return m_flEncAmpRangeArithCod;
}

int CFileEncode::getDecAmpRangeArithCod()
{
    return m_flDecAmpRangeArithCod;
}

//=======================
// CFileNetwork Class
//=======================

CFileNetwork::CFileNetwork ()
{
    m_x1_xoffset=NULL;
    m_x1_gain=NULL;
    m_x1_ymin=0;
    m_b1=NULL;
    m_w1=NULL;
    m_b2=NULL;
    m_w2=NULL;
    m_blockSize = 0;
    m_chosenNet  = 0;
    m_size_x=0;
    m_size_h=0;
    m_size_y=0;
}

CFileNetwork::~CFileNetwork()
{
    if (m_x1_xoffset != NULL) delete [] m_x1_xoffset;

    if (m_x1_gain != NULL) delete [] m_x1_gain;

    if (m_b1 != NULL) delete [] m_b1;

    for(int i = 0; i < m_size_h; ++i)
    {    
        if (m_w1[i] != NULL) delete [] m_w1[i];
    }

    if (m_w1 != NULL) delete [] m_w1;

    if (m_b2 != NULL) delete [] m_b2;

    for(int i = 0; i < m_size_y; ++i)
    {    
        if (m_w2[i] != NULL) delete [] m_w2[i];
    }

    if (m_w2 != NULL) delete [] m_w2;
}

void CFileNetwork::setFileName(char* file)
{
    strcpy(m_fileName, file);
}

void CFileNetwork::setBlockSize(int blockSize)
{
    m_blockSize = blockSize;
}

void CFileNetwork::setNetwork(int num)
{
    m_chosenNet = num;
}

void CFileNetwork::loadData()
{
    m_size_x = 64;
    m_size_h = 16;
    m_size_y = 3;

    int net = 0;

    m_x1_xoffset = new double[m_size_x];

    m_x1_gain = new double[m_size_x];

    m_b1 = new double[m_size_h];

    m_w1 = new double*[m_size_h];
    for(int i = 0; i < m_size_h; ++i)
    {
        m_w1[i] = new double[m_size_x];
    }

    m_b2 = new double[m_size_y];

    m_w2 = new double*[m_size_y];
    for(int i = 0; i < m_size_y; ++i)
    {
        m_w2[i] = new double[m_size_h];
    }

    FILE* netFile;
    char auxStr[(24*m_size_x)+1];
    char auxStr2[(24*m_size_x)+1];
    char* sNumber;
    // Opening the "dictionary.dat" file
    netFile = fopen(m_fileName,"r");

    fgets( auxStr, _MAX_PATH, netFile );
    sNumber = &auxStr[0];
    net = atoi(sNumber);    // <<<<<<<<<<<<<<<<


    while (net != m_chosenNet)
    {
        for (int i = 0; i<42; i++)
        {
            fgets( auxStr, (24*m_size_x)+1, netFile );
        }
            fgets( auxStr, _MAX_PATH, netFile );
            sNumber = &auxStr[0];
            net = atoi(sNumber);    // <<<<<<<<<<<<<<<<
    }

    // read x1_xoffset
    fgets( auxStr, (24*m_size_x)+1, netFile );
    for (int i = 0; i < m_size_x; i++)
    {
        strcpy(auxStr2,auxStr);
        auxStr2[24*(i+1)] = NULL;
        sNumber = &auxStr[24*i];
        // cout << sNumber << endl;
        m_x1_xoffset[i] = atof(sNumber);    // <<<<<<<<<<<<<<<<
        // cout << m_x1_xoffset[i] << endl;
    }
    fgets( auxStr, _MAX_PATH, netFile );
    // read x1_gain
    fgets( auxStr, (24*m_size_x)+1, netFile );
    for (int i = 0; i < m_size_x; i++)
    {
        strcpy(auxStr2,auxStr);
        auxStr2[24*(i+1)] = NULL;
        sNumber = &auxStr[24*i];
        m_x1_gain[i] = atof(sNumber);    // <<<<<<<<<<<<<<<<
    }

    fgets( auxStr, _MAX_PATH, netFile );
    // read x1_ymin
    fgets( auxStr, 24, netFile );
    // auxStr2[24] = NULL;
    sNumber = &auxStr[0];
    m_x1_ymin = atof(sNumber);    // <<<<<<<<<<<<<<<<

    fgets( auxStr, _MAX_PATH, netFile );
    // read b1
    fgets( auxStr, (24*m_size_h)+1, netFile );
    for (int i = 0; i < m_size_h; i++)
    {
        strcpy(auxStr2,auxStr);
        auxStr2[24*(i+1)] = NULL;
        sNumber = &auxStr[24*i];
        m_b1[i] = atof(sNumber);    // <<<<<<<<<<<<<<<<
    }

    fgets( auxStr, _MAX_PATH, netFile );
    // read w1
    for (int i = 0; i < m_size_h; i++)
    {
        fgets( auxStr, (24*m_size_x)+1, netFile );
        for (int j = 0; j < m_size_x; j++)
        {
            strcpy(auxStr2,auxStr);
            auxStr2[24*(j+1)] = NULL;
            sNumber = &auxStr[24*j];
            m_w1[i][j] = atof(sNumber);    // <<<<<<<<<<<<<<<<
        }
        fgets( auxStr, _MAX_PATH, netFile );
    }

    // fgets( auxStr, _MAX_PATH, netFile );
    // read b2
    fgets( auxStr, (24*m_size_y)+1, netFile );
    for (int i = 0; i < m_size_y; i++)
    {
        strcpy(auxStr2,auxStr);
        auxStr2[24*(i+1)] = NULL;
        sNumber = &auxStr[24*i];
        m_b2[i] = atof(sNumber);    // <<<<<<<<<<<<<<<<
    }

    fgets( auxStr, _MAX_PATH, netFile );
    // read w2
    for (int i = 0; i < m_size_y; i++)
    {
        fgets( auxStr, (24*m_size_h)+1, netFile );
        for (int j = 0; j < m_size_h; j++)
        {
            strcpy(auxStr2,auxStr);
            auxStr2[24*(j+1)] = NULL;
            sNumber = &auxStr[24*j];
            m_w2[i][j] = atof(sNumber);    // <<<<<<<<<<<<<<<<
        }
        fgets( auxStr, _MAX_PATH, netFile );
    }
    fclose(netFile);
}

double* CFileNetwork::getX1_xoffset()
{
    return m_x1_xoffset;
}

double* CFileNetwork::getX1_gain()
{
    return m_x1_gain;
}

double CFileNetwork::getX1_ymin()
{
    return m_x1_ymin;
}

double* CFileNetwork::getB1()
{
    return m_b1;
}

double** CFileNetwork::getW1()
{
    return m_w1;
}

double* CFileNetwork::getB2()
{
    return m_b2;
}

double** CFileNetwork::getW2()
{
    return m_w2;
}
