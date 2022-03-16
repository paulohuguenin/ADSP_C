// linkstrbook.cpp

#include "linkstrbook.h"

// ============================
//  CStructBook Class
// ============================
CLinkStrBook::CLinkStrBook()
{
    m_iNumElement = 0;
    m_norm = 0.0;
}

CLinkStrBook::CLinkStrBook(const CLinkStrBook& linkStrBook)
{
    copy( &linkStrBook);
}

CLinkStrBook::CLinkStrBook(const CLinkStrBook* linkStrBook)
{
    copy( linkStrBook);
}

void CLinkStrBook::copy(const CLinkStrBook* linkStrBook)
{
    m_iNumElement = linkStrBook->getNumElement();
    m_norm = linkStrBook->getNorm();
}

double CLinkStrBook::linearQuant(double x, double step)
{
    double y;
    if (step==0.0)
        y = 0.0;    
    else
        y = (int)( (x + step/2) / step) * step;
        
    return y;
}

void CLinkStrBook::setNumElement ( int numElement )
{
    m_iNumElement = numElement;
}

void CLinkStrBook::setNorm( double norm )
{
    m_norm = norm;
}

int CLinkStrBook::getNumElement() const
{
    return m_iNumElement;
}

double CLinkStrBook::getNorm() const
{
    return m_norm;
}

CLinkStrBook& CLinkStrBook::operator=(CLinkStrBook& c)
{
    this->m_iNumElement = c.getNumElement();
    this->m_norm = c.getNorm();
    return *this;
}

// ============================
//  CLinkStrBookExp Class
// ============================
CLinkStrBookExp::CLinkStrBookExp()
{
    m_linkStrBook = NULL;
}

CLinkStrBookExp::CLinkStrBookExp(const CLinkStrBookExp& linkStrBook)
{
    CLinkStrBookExp::copy( &linkStrBook);
}

CLinkStrBookExp::CLinkStrBookExp(const CLinkStrBookExp* linkStrBook)
{
    CLinkStrBookExp::copy( linkStrBook);
}

void CLinkStrBookExp::copy(const CLinkStrBookExp* linkStrBook)
{
    CLinkStrBook::copy( linkStrBook);
    
    if (m_linkStrBook == NULL)
        m_linkStrBook = new strtLinkExp[m_iNumElement];
    int i,j;
    strtLinkExp* linkStrBook_aux;
    linkStrBook_aux = linkStrBook->getLinkStrBook();
    for (i=0; i < m_iNumElement; i++)
    {
        m_linkStrBook[i].coefFirst = linkStrBook_aux[i].coefFirst;
        m_linkStrBook[i].coefSum = linkStrBook_aux[i].coefSum;
        m_linkStrBook[i].coefSqrSum = linkStrBook_aux[i].coefSqrSum;
        m_linkStrBook[i].xi = linkStrBook_aux[i].xi;
        m_linkStrBook[i].phase = linkStrBook_aux[i].phase;
        m_linkStrBook[i].initSample = linkStrBook_aux[i].initSample;
        m_linkStrBook[i].endSample = linkStrBook_aux[i].endSample;
        m_linkStrBook[i].initBlock = linkStrBook_aux[i].initBlock;
        m_linkStrBook[i].endBlock = linkStrBook_aux[i].endBlock;
        m_linkStrBook[i].numBlock = linkStrBook_aux[i].numBlock;
        if (this->m_linkStrBook[i].pRho != NULL)
        {
            delete [] this->m_linkStrBook[i].pRho;
            delete [] this->m_linkStrBook[i].pRhoSign;
        }
        m_linkStrBook[i].pRho = new double[m_linkStrBook[i].numBlock];
        m_linkStrBook[i].pRhoSign = new int[m_linkStrBook[i].numBlock];
        for(j=0; j<m_linkStrBook[i].numBlock; j++)
        {
            m_linkStrBook[i].pRho[j] = linkStrBook_aux[i].pRho[j];
            m_linkStrBook[i].pRhoSign[j] = linkStrBook_aux[i].pRho[j];
        }
    }
}


CLinkStrBookExp::~CLinkStrBookExp()
{
    int i,j;
    if ( m_linkStrBook!=NULL )
    {
        for (i=0; i < m_iNumElement; i++)
        {
            delete [] m_linkStrBook[i].pRho;
            m_linkStrBook[i].pRho = NULL;
            delete [] m_linkStrBook[i].pRhoSign;
            m_linkStrBook[i].pRhoSign=NULL;
        }
        delete [] m_linkStrBook;
        m_linkStrBook = NULL;       
    }
}

void CLinkStrBookExp::load( CStructBook** structBook,
                            strtSBBHeader mainHeader,
                            int iSignal)
{
    int numSignal = mainHeader.numSignal;
    int numBlock = mainHeader.numBlock;
    strtContinuousExp *sb, *sb_aux ;
    int sbNumElement;

    int i,j,ishift;
    int prevAtom, nextAtom;
    double coefFirst;
    int sbElementTotal = 0;
    for (i=0; i< numBlock; i++)
    {
        sbElementTotal+=sbNumElement = ((CStructBookExp*)structBook[iSignal])[i].getNumElement();
    }
    allocateElement(sbElementTotal);
    for (i=numBlock-1; i>=0; i--)
    {
        sb = ((CStructBookExp*)structBook[iSignal])[i].getStructBook();
        sbNumElement = ((CStructBookExp*)structBook[iSignal])[i].getNumElement();
        //printf("sbNumElement %d\n",sbNumElement);
        //cout << ((CStructBookExp*)structBook[iSignal])[i].getNumElement() << endl;
        for (j=0;j<sbNumElement;j++)
        {
            /*printf("Coef: %f | rho: %f | xi: %f | phase: %f |a: %d| b: %d|na: %d|pa: %d|origatom: %d|\n",
                                                                sb[j].innerProduct,
                                                                sb[j].rho,
                                                                sb[j].xi,
                                                                sb[j].phase,
                                                                sb[j].a,
                                                                sb[j].b,
                                                                sb[j].nextAtom,
                                                                sb[j].prevAtom,
                                                                sb[j].origAtomIndex);
            */

            prevAtom = sb[j].prevAtom;
            nextAtom = sb[j].nextAtom;
            if (nextAtom==-1)
            {
                //allocateElement();
                m_iNumElement++;
                coefFirst = sb[j].innerProduct;
                m_linkStrBook[m_iNumElement-1].coefSum = m_linkStrBook[m_iNumElement-1].coefSum + 
                                                        sb[j].innerProduct;
                m_linkStrBook[m_iNumElement-1].coefSqrSum = m_linkStrBook[m_iNumElement-1].coefSqrSum + 
                                                            sb[j].innerProduct*sb[j].innerProduct;  
                addDecayingFactor (m_iNumElement-1, sb[j].rho);
                m_linkStrBook[m_iNumElement-1].xi = sb[j].xi;
                m_linkStrBook[m_iNumElement-1].phase = sb[j].phase;
                //m_linkStrBook[m_iNumElement-1].initSample = 0;
                m_linkStrBook[m_iNumElement-1].initSample = sb[j].a ;
                m_linkStrBook[m_iNumElement-1].endSample = sb[j].b;
                m_linkStrBook[m_iNumElement-1].initBlock=i;
                m_linkStrBook[m_iNumElement-1].endBlock=i;
                ishift = 1;
                while (prevAtom!=-1)
                {
                    sb_aux = ((CStructBookExp*)structBook[iSignal])[i-ishift].getStructBook();
                    sb_aux[prevAtom].nextAtom = i-ishift+1;
                    coefFirst = sb_aux[prevAtom].innerProduct;
                    m_linkStrBook[m_iNumElement-1].coefSum = m_linkStrBook[m_iNumElement-1].coefSum + 
                                                             sb_aux[prevAtom].innerProduct;
                    m_linkStrBook[m_iNumElement-1].coefSqrSum = m_linkStrBook[m_iNumElement-1].coefSqrSum + 
                                                                sb_aux[prevAtom].innerProduct*sb_aux[prevAtom].innerProduct;
                    addDecayingFactor (m_iNumElement-1, sb_aux[prevAtom].rho);
                    m_linkStrBook[m_iNumElement-1].xi = sb_aux[prevAtom].xi;
                    m_linkStrBook[m_iNumElement-1].phase = sb_aux[prevAtom].phase;
                    //m_linkStrBook[m_iNumElement-1].initSample = 0;
                    m_linkStrBook[m_iNumElement-1].initSample = sb_aux[prevAtom].a ;
                    m_linkStrBook[m_iNumElement-1].initBlock=i-ishift;
                    ishift++;
                    prevAtom = sb_aux[prevAtom].prevAtom;
                }
                m_linkStrBook[m_iNumElement-1].coefFirst = coefFirst;
            }
        }
    }
}

void CLinkStrBookExp::allocateElement (int numElement)
{
    if ( m_linkStrBook!=NULL )
    {
        delete [] m_linkStrBook;
    }
    m_linkStrBook = new strtLinkExp[numElement];
}

void CLinkStrBookExp::allocateElement ()
{
    strtLinkExp* aux;

    if (m_linkStrBook==NULL)
    {
        m_linkStrBook = new strtLinkExp;
    }

    aux = new strtLinkExp[m_iNumElement + 1];


    memcpy ( aux,m_linkStrBook,m_iNumElement*sizeof ( strtLinkExp ) );

    aux[m_iNumElement].coefFirst = 0.0;
    aux[m_iNumElement].coefSum = 0.0;
    aux[m_iNumElement].coefSqrSum = 0.0;
    aux[m_iNumElement].pRho = NULL;
    aux[m_iNumElement].pRhoSign = NULL;
    aux[m_iNumElement].xi = 0.0;
    aux[m_iNumElement].phase = 0.0;
    aux[m_iNumElement].initSample = 0;
    aux[m_iNumElement].endSample = 0;
    aux[m_iNumElement].initBlock = 0;
    aux[m_iNumElement].endBlock = 0;
    aux[m_iNumElement].numBlock = 0;

    if ( m_linkStrBook!=NULL )
    {
        delete [] m_linkStrBook;
    }

    m_linkStrBook = new strtLinkExp[m_iNumElement + 1];

    memcpy ( m_linkStrBook,aux, ( m_iNumElement+1 ) *sizeof ( strtLinkExp ) );

    m_iNumElement++;

    delete [] aux;
}

void CLinkStrBookExp::addDecayingFactor (int index, double rho)
{

    double* auxRho;
    int* auxRhoSign;

    int numBlock = m_linkStrBook[index].numBlock;
    auxRho = new double[numBlock+1];
    auxRhoSign = new int[numBlock+1];

    if (m_linkStrBook[index].pRho!=NULL)
    {
        memcpy ( &auxRho[1],m_linkStrBook[index].pRho,numBlock*sizeof ( double ) );
        memcpy ( &auxRhoSign[1],m_linkStrBook[index].pRhoSign,numBlock*sizeof ( int ) );
    }
    auxRho[0] = fabs(rho);
    if (rho>=0)
        auxRhoSign[0] = 1;
    else
        auxRhoSign[0] = -1;

    if ( m_linkStrBook[index].pRho!=NULL )
    {
        delete [] m_linkStrBook[index].pRho;
    }
    if ( m_linkStrBook[index].pRhoSign!=NULL )
    {
        delete [] m_linkStrBook[index].pRhoSign;
    }

    m_linkStrBook[index].pRho = new double[numBlock+1];
    m_linkStrBook[index].pRhoSign = new int[numBlock+1];


    memcpy ( m_linkStrBook[index].pRho, auxRho, ( numBlock+1 ) *sizeof ( double ) );
    memcpy ( m_linkStrBook[index].pRhoSign, auxRhoSign, ( numBlock+1 ) *sizeof ( int ) );

    m_linkStrBook[index].numBlock = numBlock+1;

    delete [] auxRho;
    delete [] auxRhoSign;

}

void CLinkStrBookExp::linkStrBook2StrBook(	CStructBook* structBook,
											strtSBBHeader mainHeader)
{
	int i,j,k,ind;
	double* rAtomBlock, *rLinkAtom;
	CExpDictionary expDic, expDicLinkedAtom;
	strtParameter expParm, expParmAux;
	double rho,prevrho, nextrho;
	int numAtomSameRho;
	
	//cout << " m_iNumElement: " << m_iNumElement << endl;
    for (i=0; i < m_iNumElement; i++)
    {
    	if (m_linkStrBook[i].numBlock==1)
		{
			//cout << " Just one atom" << endl;
		    expParm.innerProduct = m_linkStrBook[i].coefFirst;
		    expParm.xi = m_linkStrBook[i].xi;
		    rho = static_cast<double>(m_linkStrBook[i].pRho[0]*
								static_cast<double>(m_linkStrBook[i].pRhoSign[0]));
		    expParm.rho = rho;
		    expParm.phase = m_linkStrBook[i].phase;
		    expParm.a = (m_linkStrBook[i].initBlock*mainHeader.blockHop) + 
		                m_linkStrBook[i].initSample;
		    expParm.b = (m_linkStrBook[i].initBlock*mainHeader.blockHop) + 
		                m_linkStrBook[i].endSample;
		    //expParm.printParm2Screen();
			((CStructBookExp*)structBook)->addElement(&expParm);			
		}
    	else
		{
		    //cout << " Linked atom" << endl;
		    cgMatrix<double> cgSignal(1,m_linkStrBook[i].numBlock*mainHeader.blockSize,0.0);
		    expDic.setSignalSize(mainHeader.blockSize);
		    //---------------------------
		    // Fill cgSignal
		    //cout << " Fill cgSignal" << endl;
			expParm.xi = m_linkStrBook[i].xi;
			expParm.phase = m_linkStrBook[i].phase;
			for (j=0; j<m_linkStrBook[i].numBlock; j++ )
			{
				expParm.rho = (double)(m_linkStrBook[i].pRho[j]*
							           static_cast<double>(m_linkStrBook[i].pRhoSign[j]));
				if (j==0)
				{
					expParm.innerProduct = m_linkStrBook[i].coefFirst;
					expParm.a = m_linkStrBook[i].initSample;
					if (m_linkStrBook[i].numBlock==1)
						expParm.b = m_linkStrBook[i].endSample;
					else
						expParm.b = expDic.getSignalSize()-1;
						
				}
				else if (j==m_linkStrBook[i].numBlock-1)
				{
					expParm.a = 0;
					expParm.b = m_linkStrBook[i].endSample;
				}
				else
				{
					expParm.a = 0;
					expParm.b = expDic.getSignalSize()-1;
				}
				
				// ----
				expDic.setRealAtom(&expParm);
				rAtomBlock = expDic.getRealAtom();
				for (k=0; k<expDic.getSignalSize();k++)
				{
					ind = (j*expDic.getSignalSize()) + k;
					cgSignal[0][ind] = expParm.innerProduct*rAtomBlock[k];	
				}
				expDic.computeNextBlockCoefPhase(&expParm);
			}
			//--------------------------------------------------
		    // Partition the Linked Atom if necessary
		    //cout << " Partition the Linked Atom if necessary" << endl;
			double phaseCurrent = m_linkStrBook[i].phase;
			expParm.xi = m_linkStrBook[i].xi;
			int initBlock = 0;
			int endBlock =0;
			int a_aux, b_aux;
			for (j=0; j<m_linkStrBook[i].numBlock-1; j++ )
			{
				expParm.rho = (double)(m_linkStrBook[i].pRho[j]*
							           static_cast<double>(m_linkStrBook[i].pRhoSign[j]));
				nextrho = (double)(m_linkStrBook[i].pRho[j+1]*
							           static_cast<double>(m_linkStrBook[i].pRhoSign[j+1]));
		        if (fabs(nextrho-expParm.rho) > 0.01*fabs(expParm.rho))
		        {
		        	numAtomSameRho = endBlock-initBlock+1;
					//----		           
					if (initBlock==0)
					{
						expParm.a = m_linkStrBook[i].initSample;
					}
					else
					{
						expParm.a = 0;
					}
					//----
					if (endBlock==m_linkStrBook[i].numBlock-1)
					{
						
						expParm.b = ((numAtomSameRho-1)*expDic.getSignalSize()) + 
						            m_linkStrBook[i].endSample;
					}
					else
					{
						expParm.b = (numAtomSameRho*expDic.getSignalSize())-1;
					}
					//---
					expDicLinkedAtom.setSignalSize(numAtomSameRho*expDic.getSignalSize());
					expParm.phase = phaseCurrent;
					expDicLinkedAtom.setRealAtom(&expParm);
					rLinkAtom = expDicLinkedAtom.getRealAtom();
					expParm.innerProduct = cgSignal.dprodInterval(	rLinkAtom,
															initBlock*mainHeader.blockSize,
															((endBlock+1)*mainHeader.blockSize) -1);
				    a_aux = expParm.a;
				    b_aux = expParm.b;
					expParm.a = (m_linkStrBook[i].initBlock*mainHeader.blockHop) +
					            initBlock*mainHeader.blockSize+expParm.a;
					expParm.b = (m_linkStrBook[i].initBlock*mainHeader.blockHop) +
					            initBlock*mainHeader.blockSize+expParm.b;
					//expParm.printParm2Screen();
					((CStructBookExp*)structBook)->addElement(&expParm);
					// update
					expParm.a = a_aux;
					expParm.b = b_aux;
					expDicLinkedAtom.computeNextBlockCoefPhase(&expParm);
					phaseCurrent = expParm.phase;
					initBlock = endBlock+1;
					endBlock = initBlock;
				}
				else
				{
					endBlock++;
				}	
				//cout << initBlock << " " << endBlock << endl;
			}
			//cout << "Last block" << endl;
			// Last block
			expParm.rho = (double)(m_linkStrBook[i].pRho[m_linkStrBook[i].numBlock-1]*
					static_cast<double>(m_linkStrBook[i].pRhoSign[m_linkStrBook[i].numBlock-1]));
			// ---
		    numAtomSameRho = endBlock-initBlock+1;
			//----		           
			if (initBlock==0)
			{
				expParm.a = m_linkStrBook[i].initSample;
			}
			else
			{
				expParm.a = 0;
			}
			//----
			if (endBlock==m_linkStrBook[i].numBlock-1)
			{
				
				expParm.b = ((numAtomSameRho-1)*expDic.getSignalSize()) + 
							m_linkStrBook[i].endSample;
			}
			else
			{
				expParm.b = (numAtomSameRho*expDic.getSignalSize())-1;
			}
			//---
			//cout << "Configure dictionary" <<endl;
			expDicLinkedAtom.setSignalSize(numAtomSameRho*expDic.getSignalSize());
			expParm.phase = phaseCurrent;
			expDicLinkedAtom.setRealAtom(&expParm);
			rLinkAtom = expDicLinkedAtom.getRealAtom();
			expParm.innerProduct = cgSignal.dprodInterval(	rLinkAtom,
									initBlock*mainHeader.blockSize,
									((endBlock+1)*mainHeader.blockSize) -1);
			expParm.a = (m_linkStrBook[i].initBlock*mainHeader.blockHop) +
						initBlock*mainHeader.blockSize+expParm.a;
			expParm.b = (m_linkStrBook[i].initBlock*mainHeader.blockHop) +
						initBlock*mainHeader.blockSize+expParm.b;
			//expParm.printParm2Screen();
			((CStructBookExp*)structBook)->addElement(&expParm);
		}
    }
}

void CLinkStrBookExp::Print()
{
    int i,j;
    for (i=0; i < m_iNumElement; i++)
    {
        printf("->  Element: %d\n", i);
        printf("        CoefFirst  = %f\n",m_linkStrBook[i].coefFirst);
        printf("        CoefSum    = %f\n",m_linkStrBook[i].coefSum);
        printf("        CoefSqrSum = %f\n",m_linkStrBook[i].coefSqrSum);
        printf("        Xi         = %f\n",m_linkStrBook[i].xi);
        printf("        Phase      = %f\n",m_linkStrBook[i].phase);
        printf("        InitSample = %d\n",m_linkStrBook[i].initSample);
        printf("        EndSample  = %d\n",m_linkStrBook[i].endSample);
        printf("        InitBlock  = %d\n",m_linkStrBook[i].initBlock);
        printf("        EndBlock   = %d\n",m_linkStrBook[i].endBlock);
        printf("        NumBlock   = %d\n",m_linkStrBook[i].numBlock);
        for(j=0; j<m_linkStrBook[i].numBlock; j++)
        {
            printf("        Rho%d       = %f\n",j,(double)(m_linkStrBook[i].pRho[j]*(double)m_linkStrBook[i].pRhoSign[j]));
        }
    }
}

void CLinkStrBookExp::findParmRange(double& min_coef,
                                    double& max_coef,
                                    double& min_rho,
                                    double& max_rho,
                                    double& min_phase,
                                    double& max_phase)
{
    int i,j;
    double mincoef_aux, minrho_aux, minphase_aux;
    double maxcoef_aux, maxrho_aux, maxphase_aux;
    mincoef_aux     = 10000;
    minrho_aux      = 10000;
    minphase_aux    = 10000;
    maxcoef_aux     = -10000;
    maxrho_aux      = -10000;
    maxphase_aux    = -10000;
    for (i=0; i < m_iNumElement; i++)
    {
        if (mincoef_aux > m_linkStrBook[i].coefFirst)
            mincoef_aux = m_linkStrBook[i].coefFirst;
        if (maxcoef_aux < m_linkStrBook[i].coefFirst)
            maxcoef_aux = m_linkStrBook[i].coefFirst;
        if (minphase_aux > m_linkStrBook[i].phase)
            minphase_aux = m_linkStrBook[i].phase;
        if (maxphase_aux < m_linkStrBook[i].phase)
            maxphase_aux = m_linkStrBook[i].phase;
        for(j=0; j<m_linkStrBook[i].numBlock; j++)
        {
            if (minrho_aux > m_linkStrBook[i].pRho[j])
                minrho_aux = m_linkStrBook[i].pRho[j];
            if (maxrho_aux < m_linkStrBook[i].pRho[j])
                maxrho_aux = m_linkStrBook[i].pRho[j];
        }
    }
    //min_coef = mincoef_aux;
    min_coef = 0.0;
    max_coef = maxcoef_aux;
    min_rho  = minrho_aux;
    max_rho  = maxrho_aux;
    min_phase= minphase_aux;
    max_phase= maxphase_aux;
}

void CLinkStrBookExp::configQuantizer(  int nbits_amp,
                                        int nbits_rho,
                                        int nbits_xi,
                                        int nbits_phase,
                                        int nbits_sample,
                                        int nbits_block)
{

    findParmRange(quantHeader.min_amp,
                  quantHeader.max_amp,
                  quantHeader.min_rho,
                  quantHeader.max_rho,
                  quantHeader.min_phase,
                  quantHeader.max_phase);

    quantHeader.nbits_amp = nbits_amp;
    quantHeader.nbits_rho = nbits_rho;
    quantHeader.nbits_xi = nbits_xi;
    quantHeader.nbits_phase = nbits_phase;
    quantHeader.nbits_sample = nbits_sample;
    quantHeader.nbits_block = nbits_block;


    quantHeader.norm = m_norm;
    // numberStruct will set after quantization
    quantHeader.numberStruct = 0;
    
    /*printf("Quantizer header:\n");
    printf(" - quantHeader.nbits_amp: %d\n",quantHeader.nbits_amp);
    printf(" - quantHeader.nbits_rho: %d\n",quantHeader.nbits_rho);
    printf(" - quantHeader.nbits_xi: %d\n",quantHeader.nbits_xi);
    printf(" - quantHeader.nbits_phase: %d\n",quantHeader.nbits_phase);
    printf(" - quantHeader.nbits_sample: %d\n",quantHeader.nbits_sample);
    printf(" - quantHeader.nbits_block: %d\n",quantHeader.nbits_block);
    printf(" - quantHeader.norm: %f\n",quantHeader.norm);
    printf(" - quantHeader.numberStruct: %d\n",quantHeader.numberStruct);
    */
}

void CLinkStrBookExp::quantize()
{
    int i,j;
    double step_amp, nlevel_amp;
    nlevel_amp = pow(2.0, (double)quantHeader.nbits_amp) - 1;
    if (nlevel_amp==0.0)
        step_amp = 0.0;
    else
        step_amp = fabs( (quantHeader.max_amp - quantHeader.min_amp) / nlevel_amp );
    
    double step_rho, nlevel_rho;
    nlevel_rho = pow(2.0, (double)quantHeader.nbits_rho) - 1;
    if (nlevel_rho==0.0)
        step_rho=0.0;
    else
        step_rho = fabs( (quantHeader.max_rho - quantHeader.min_rho) / nlevel_rho );
        
    double step_phase, nlevel_phase;
    nlevel_phase = pow(2.0, (double)quantHeader.nbits_phase) - 1;
    if (step_phase==0.0)
        step_phase =0.0;
    else
        step_phase = fabs( (quantHeader.max_phase - quantHeader.min_phase) / nlevel_phase );
    /*printf("nlevel_amp %f step_amp %f\n",nlevel_amp,step_amp);
    printf("nlevel_rho %f step_rho %f\n",nlevel_rho,step_rho);
    printf("nlevel_phase %f step_phase %f\n",nlevel_phase,step_phase);
    */
    double coefFirst;
    double phase;
    double pRho;
    for (i=0; i< m_iNumElement; i++)
    {
        //printf(" - numElement: %d\n",i);
        coefFirst = m_linkStrBook[i].coefFirst;
        m_linkStrBook[i].coefFirst = linearQuant(coefFirst,step_amp);
        //printf("coefFirst: %f = %f\n", coefFirst,m_linkStrBook[i].coefFirst);
        if (m_linkStrBook[i].coefFirst!=0.0)
        {
            quantHeader.numberStruct++;
            phase = m_linkStrBook[i].phase;
            m_linkStrBook[i].phase = linearQuant(phase,step_phase);
            //printf("phase: %f = %f\n", phase,m_linkStrBook[i].phase);
            for (j=0; j< m_linkStrBook[i].numBlock; j++)
            {
                pRho = m_linkStrBook[i].pRho[j];
                m_linkStrBook[i].pRho[j] = linearQuant(pRho,step_rho);
                //printf("rho (%d): %f = %f\n",j, pRho,m_linkStrBook[i].pRho[j]);
            }
        }
    }
}

int CLinkStrBookExp::computeNumBits()
{
    int i,j;
    int nbits = 0;
    for (i=0; i< m_iNumElement; i++)
    {
        if (m_linkStrBook[i].coefFirst != 0.0)
        {
            nbits = nbits + quantHeader.nbits_amp; // coefFirst
            nbits = nbits + quantHeader.nbits_phase; // phase
            nbits = nbits + quantHeader.nbits_xi; // xi
            nbits = nbits + quantHeader.nbits_sample; // initSample
            nbits = nbits + quantHeader.nbits_sample; // endSample
            nbits = nbits + quantHeader.nbits_block;  // initBlock
            nbits = nbits + quantHeader.nbits_block; // numBlock;
            for (j=0; j< m_linkStrBook[i].numBlock; j++)
            {
                nbits = nbits + quantHeader.nbits_rho*m_linkStrBook[i].numBlock; // pRho
                nbits = nbits + m_linkStrBook[i].numBlock; // pRhoSign
            }
        }
    }
    return nbits;
}

void CLinkStrBookExp::synthSignal(  double* recSignal,
									strtSBBHeader mainHeader)
{
    CExpDictionary expDic;
    strtParameter* parm;
    // parm = new CExpParm;
    strtLinkExp* pLinkStrExp;
    pLinkStrExp = m_linkStrBook;
    int numElement = m_iNumElement;
    double norm = m_norm;
    int signalSize = mainHeader.signalSize;
    
    double* realAtom;
    
    int i,j,k,index;

	expDic.setSignalSize(mainHeader.blockSize);

    for (i=0; i<signalSize; i++)
        recSignal[i] = 0.0;

    for (i=0; i<numElement; i++)
    { 
        //for (j=0; j<signalSize; j++)
        //    sig_aux[j] = 0.0;
        //cout << " Sintetiza  elemento: "<< i << endl;
        if (pLinkStrExp[i].coefFirst!=0.0)
        {
            ((strtParameter*)parm)->xi = pLinkStrExp[i].xi;
            ((strtParameter*)parm)->phase = pLinkStrExp[i].phase;
            for (j=0; j<pLinkStrExp[i].numBlock; j++ )
            {
                //cout << "   Sintetiza  bloco: "<< j << endl;
                ((strtParameter*)parm)->rho = (double)(pLinkStrExp[i].pRho[j]*
                                                  static_cast<double>(pLinkStrExp[i].pRhoSign[j]));
                if (j==0)
                {
                    ((strtParameter*)parm)->innerProduct = pLinkStrExp[i].coefFirst;
                    ((strtParameter*)parm)->a = pLinkStrExp[i].initSample;
                    if (pLinkStrExp[i].numBlock==1)
                        ((strtParameter*)parm)->b = pLinkStrExp[i].endSample;
                    else
                        ((strtParameter*)parm)->b = expDic.getSignalSize()-1;
                }
                else if (j==pLinkStrExp[i].numBlock-1)
                {
                    ((strtParameter*)parm)->a = 0;
                    ((strtParameter*)parm)->b = pLinkStrExp[i].endSample;
                }
                else
                {
                    ((strtParameter*)parm)->a = 0;
                    ((strtParameter*)parm)->b = expDic.getSignalSize()-1;
                }
                //printf("Element: %d; Block: %d\n",i,j);
                //((CExpParm*)parm)->printParm2Screen();
                //cout << "   SetRealAtom " << endl;
                expDic.setRealAtom(parm);
                realAtom = expDic.getRealAtom();
                
                for (k=0; k<expDic.getSignalSize();k++)
                {
                    index = ( (pLinkStrExp[i].initBlock+j)*expDic.getSignalSize() ) + k;
                    if (index >= signalSize) break;
                    recSignal[index] = recSignal[index] + norm*((strtParameter*)parm)->innerProduct*realAtom[k];
                    //sig_aux[k] = norm*((CExpParm*)parm)->innerProd*m_realAtom[k];
                }
                //cout << "   calcula fase do proximo bloco " << endl;
                expDic.computeNextBlockCoefPhase(parm);
            }
            
        }
    }
    
    delete parm;
    
}

void CLinkStrBookExp::synthAtom(  double* recSignal,
                                  strtSBBHeader mainHeader,
                                  int nElement)
{
    CExpDictionary expDic;
    strtParameter* parm;
    // parm = new CExpParm;
    strtLinkExp* pLinkStrExp;
    pLinkStrExp = m_linkStrBook;
    int numElement = m_iNumElement;
    double norm = m_norm;
    int signalSize = mainHeader.signalSize;
    
    double* realAtom;
    
    int i,j,k,index;

	expDic.setSignalSize(mainHeader.blockSize);

    for (i=0; i<signalSize; i++)
        recSignal[i] = 0.0;

    if (pLinkStrExp[nElement].coefFirst!=0.0)
    {
        ((strtParameter*)parm)->xi = pLinkStrExp[nElement].xi;
        ((strtParameter*)parm)->phase = pLinkStrExp[nElement].phase;
        for (j=0; j<pLinkStrExp[nElement].numBlock; j++ )
        {
            //cout << "   Sintetiza  bloco: "<< j << endl;
            ((strtParameter*)parm)->rho = (double)(pLinkStrExp[nElement].pRho[j]*
                                           static_cast<double>(pLinkStrExp[nElement].pRhoSign[j]));
            if (j==0)
            {
                ((strtParameter*)parm)->innerProduct = pLinkStrExp[nElement].coefFirst;
                ((strtParameter*)parm)->a = pLinkStrExp[nElement].initSample;
                if (pLinkStrExp[nElement].numBlock==1)
                    ((strtParameter*)parm)->b = pLinkStrExp[nElement].endSample;
                else
                    ((strtParameter*)parm)->b = expDic.getSignalSize()-1;
            }
            else if (j==pLinkStrExp[nElement].numBlock-1)
            {
                ((strtParameter*)parm)->a = 0;
                ((strtParameter*)parm)->b = pLinkStrExp[nElement].endSample;
            }
            else
            {
                ((strtParameter*)parm)->a = 0;
                ((strtParameter*)parm)->b = expDic.getSignalSize()-1;
            }
            //printf("Element: %d; Block: %d\n",i,j);
            //((CExpParm*)parm)->printParm2Screen();
            //cout << "   SetRealAtom " << endl;
            expDic.setRealAtom(parm);
            realAtom = expDic.getRealAtom();
            //printf("InitBlockSample %d \n",(pLinkStrExp[i].initBlock+j)*m_signalSize);
            //cout << "   Soma ao Sinal Recontruido" << endl;
            for (k=0; k<expDic.getSignalSize();k++)
            {
                index = ( (pLinkStrExp[nElement].initBlock+j)*expDic.getSignalSize() ) + k;
                if (index >= signalSize) break;
                recSignal[index] = recSignal[index] + norm*((strtParameter*)parm)->innerProduct*realAtom[k];
                //sig_aux[k] = norm*((CExpParm*)parm)->innerProd*m_realAtom[k];
            }
            //cout << "   calcula fase do proximo bloco " << endl;
            expDic.computeNextBlockCoefPhase(parm);
        }
    }
    delete parm;
}



strtLinkExp* CLinkStrBookExp::getLinkStrBook() const
{
    return m_linkStrBook;
}

CLinkStrBookExp& CLinkStrBookExp::operator=( CLinkStrBookExp& c)
{
    this->setNumElement(c.getNumElement());
    this->setNorm(c.getNorm());
    int allocateFirstTime = 0;
    if (this->m_linkStrBook == NULL)
    {
        m_linkStrBook = new strtLinkExp[m_iNumElement];
        allocateFirstTime = 1;
    }
    int i,j;
    strtLinkExp* linkStrBook_aux;
    linkStrBook_aux = c.getLinkStrBook();
    for (i=0; i < m_iNumElement; i++)
    {
        this->m_linkStrBook[i].coefFirst = linkStrBook_aux[i].coefFirst;
        this->m_linkStrBook[i].coefSum = linkStrBook_aux[i].coefSum;
        this->m_linkStrBook[i].coefSqrSum = linkStrBook_aux[i].coefSqrSum;
        this->m_linkStrBook[i].xi = linkStrBook_aux[i].xi;
        this->m_linkStrBook[i].phase = linkStrBook_aux[i].phase;
        this->m_linkStrBook[i].initSample = linkStrBook_aux[i].initSample;
        this->m_linkStrBook[i].endSample = linkStrBook_aux[i].endSample;
        this->m_linkStrBook[i].initBlock = linkStrBook_aux[i].initBlock;
        this->m_linkStrBook[i].endBlock = linkStrBook_aux[i].endBlock;
        this->m_linkStrBook[i].numBlock = linkStrBook_aux[i].numBlock;
        if ( (this->m_linkStrBook[i].pRho != NULL) && (allocateFirstTime==0) )
        {
            delete [] this->m_linkStrBook[i].pRho;
            delete [] this->m_linkStrBook[i].pRhoSign;
        }
        this->m_linkStrBook[i].pRho = new double[this->m_linkStrBook[i].numBlock];
        this->m_linkStrBook[i].pRhoSign = new int[this->m_linkStrBook[i].numBlock];
        for(j=0; j<this->m_linkStrBook[i].numBlock; j++)
        {
            this->m_linkStrBook[i].pRho[j] = linkStrBook_aux[i].pRho[j];
            this->m_linkStrBook[i].pRhoSign[j] = linkStrBook_aux[i].pRhoSign[j];
        }
    }

    return *this;
}


