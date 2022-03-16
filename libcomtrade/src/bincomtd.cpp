/*****************************************************************************

   Projeto : SISTEMA INTEGRADO DE APOIO `A ANALISE DE PERTURBACOES - SINAPE
   Nome    : BINCOMTD.CPP
   Funcao  : Rotinas para conversão entre formato interno (CConfigCOMTRADE)
            e o formato COMTRADE (IEEE standard Common Format for Transient
            Data Exchange for power systems), IEEE C37.111-1991
   C       : C++
   Autor   : MAMR
   Data    : nov/94

*****************************************************************************/

#include <math.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include "bincomtd.h"

#ifndef ASSERT
    #define ASSERT      assert
#endif

#define COMTRADE001 "Valores minimo e maximo do sinal nao estao definidos"
#define COMTRADE002 "Nao consegue gerar arquivo de cabecalho"
#define COMTRADE003 "Arquivo com numero de caracteres maior que"
#define COMTRADE004 "Nao ha' espaco para ler cabecalho do arquivo"
#define COMTRADE005 "pode estar incompleto ou vazio"
#define COMTRADE006 "Valor de tempo das amostras maior que 999999999.5"
#define COMTRADE007 "nao possui formato compativel com o COMTRADE"
#define COMTRADE008 "Tradutor trabalha com vetor de inteiros no intervalo [%d,%d]"
#define COMTRADE009 "Erro interno do tradutor !!"
#define COMTRADE010 "Nao ha' memoria para ler configuracao !!"
#define COMTRADE011 "Nao ha' memoria para ler dados !!"

#ifndef max
#define max(a,b)        (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
#define min(a,b)        (((a) < (b)) ? (a) : (b))
#endif

#define ERRODFORMATO()  {\
        if (erro!=NULL)\
        sprintf(erro,"%s %s\nPosicao:%lu",GetNomeArquivo(),COMTRADE007,posicao);\
        free(canais_); free(auxV); return 0;}

int CBinaryComtrade::ExtrairConfiguracao(CConfigComtrade* config, char* erro) {
    return CComtrade::ExtrairConfiguracao( config, erro);
}

int CBinaryComtrade::ExtrairConfiguracao_(CConfigComtrade* config, char* erro) {
    return CComtrade::ExtrairConfiguracao_( config, erro);
}

int CBinaryComtrade::ExtrairDados( CConfigComtrade* config,
								int* numero_canais, 
								size_t tam, 
								void* canais,
                                size_t numpts, 
								long offset,
                                char AcessoRandomico, 
								char* erro) {

    char    cont;
    size_t  i;
    int     j, k, NumBytesPorAmostra, NumCanais;
    long    posicao, NumBytesTotal;
    void    *canais_;
    float   **fVetor;
    int     **iVetor;
    char    **cVetor;
    short   *pI;
    int     *pInit;
    char    *auxV, *pC, *pC0, *pI0;
    CCanal *pCanal;

    // Verifica se tipo do arquivo e' ASCII
    if (config->GetTipoArquivo()==0) 
		return CComtrade::ExtrairDados( config, 
										numero_canais, 
										tam, 
										canais, 
										numpts, 
										offset, 
										AcessoRandomico, 
										erro);

    // aloca buffers necessários
    NumBytesPorAmostra = 2*sizeof(DWORD) + 2*config->GetNumCanaisAnalogicos() +
                         2*(int)ceil((float)config->GetNumCanaisDigitais()/16);
    auxV = (char*) calloc( NumBytesPorAmostra, 1);

    j = 0;
    for (i=0; i<tam; i++) if (j<numero_canais[i]) j = numero_canais[i];

    canais_ = calloc( j, sizeof(void*));
    if ( (auxV==NULL) || (canais_==NULL) ) 
	{
        if (erro!=NULL) sprintf( erro, COMTRADE011);
        free(canais_); 
		free(auxV);  
		return 0;
    }
    memcpy( canais_, canais, j*sizeof(void*));

    // Posiciona o buffer na posicao correta
    if (AcessoRandomico) 
	{
    	if (!EncherBuffer( offset*NumBytesPorAmostra, erro)) return 0;
    };
    posicao = offset * NumBytesPorAmostra;
	NumBytesTotal = NumBytesPorAmostra * config->GetNumeroAmostrasPorCanal();

    // Le arquivo de dados
    NumCanais=config->GetNumCanaisAnalogicos()+config->GetNumCanaisDigitais();
    pI0 = auxV + 2*sizeof(DWORD);
    pC0 = pI0 + config->GetNumCanaisAnalogicos()*2;

    //printf("numpts: %d; NumBytesTotal: %d; \n", numpts, NumBytesTotal);
    
    if (m_timestamp ==NULL) m_timestamp = new double[numpts];
    //double *timestamp = new double[numpts];
    double *diffTStamp = new double[numpts-1];
    double *samplingFreq = new double[numpts-1];
    int nrate=0;
    double Fs;
    int endsample;
    FILE * pFile;
    double fator = static_cast<double>(config->GetFator());
    pFile = fopen ("timestamp.out","w");
    fprintf (pFile, "TIMESTAMP     SAMPLING FREQ    ENDSAMPLE \n");
    for (i=0; (i<numpts) && (posicao<NumBytesTotal); i++) 
	{
		
        // realinha pointers e indices
        k = 0;                          
		fVetor = (float**) canais_;
        iVetor = (int**) canais_;       
		cVetor = (char**) canais_;

        // extrai amostras indice i do buffer
        if (!ExtrairBytes( auxV, NumBytesPorAmostra))	
		{
        	long linha1 = (long)i; 
        	config->AlterarTaxa(config->GetNumeroTaxas(),-1.F,linha1); 
        	free(canais_); 
			free(auxV); 
        	return 1;  
        }
        
        	
        // Verifica se o valor lido corresponde ao número da amostra
        // New - Só se o usuário quiser ler o COMTRADE correspondente ao no da amostra
        if ( m_iVerificaAmostra == 1)	
		{
        	if ((*(long*)auxV)!=offset+(long)i)	
			{
				if ((*(long*)auxV)!=offset+(long)i+1) ERRODFORMATO();           
        	}
        }	
        
        pInit = (int*)auxV; 		
        //printf("- i:%d; posicao: %d \n",i,posicao);
        //printf("ind: %d; timestamp: %d \n", *pInit, *(pInit+1) );
        
        m_timestamp[i] = fator * static_cast<double>(*(pInit+1));
        
        if (i>0) 
        {
        	diffTStamp[i-1] = m_timestamp[i] - m_timestamp[i-1];
        	Fs = (1.0/diffTStamp[i-1])*1000000;
        	samplingFreq[i-1] = Fs;
        	endsample = i+1;
        	//printf("Fs: %f; endsample: %d \n", Fs, endsample);
        	fprintf (pFile, "%12.2f %12.2f     %10d \n",m_timestamp[i],samplingFreq[i-1],endsample);
        }
        else
        {
        	endsample = i+1;
        	fprintf (pFile, "%12.2f %12.12f    %10d \n",m_timestamp[i],0.0,endsample);
        }
                
        // Distribui nos vetores destino
        pI = (short*)pI0; 
		pC = pC0;   
		cont=0;
        for (j=1; j<=NumCanais; j++) 
		{
            pCanal = config->GetCanalPtr(j);
            if (j==numero_canais[k]) 
			{
				if (pCanal->GetAnalogico()) // analogico
				{
					switch (m_eFormatoVetorExtraido) {
					case e_inteiro: 
						*(*iVetor)++ = *pI;
					break;
					case e_real: *(*fVetor)++ = (float)*pI;
					break;
					default:
						if (erro==NULL) return 0;
						sprintf( erro, COMTRADE009);
						free(canais_); free(auxV);
						return 0;
					break;
					};
                pI++;
				} 
				else //digital
				{ 
					*(*cVetor)++ = (*pC >> cont++) & 1;
					if (cont == 8) 
					{ 
						pC++; 
						cont = 0; 
					};
				};
				k++;
				fVetor++;
				iVetor++; 
				cVetor++;
			} 
			else // canal nao interessa
			{ 
				if (pCanal->GetAnalogico()) pI++;
				else 
				{ 
					cont++; 
					if (cont==8) 
					{ 
						pC++; 
						cont=0; 
					}
				};
			};

        	if (k>=(int)tam) break;

		}; // for j

        posicao += NumBytesPorAmostra;
	}; // for i
	fclose(pFile);
	
//	delete [] timestamp;
    delete [] diffTStamp;
    delete [] samplingFreq;
    
	free(canais_); 
	free(auxV);

    return 1;
};

int CBinaryComtrade::InserirConfiguracao(CConfigComtrade* config, char* erro){
    
	int ret;
    // verifica se tipo do arquivo e' ASCII
    if (config->GetTipoArquivo()==0) 
	{
		ret = CComtrade::InserirConfiguracao( config, erro);
	} 
	else 
	{
		ret = CComtrade::InserirConfiguracao_( config, erro);
	    if (!InserirLinha("BINARY")) return 0;
	};

    return ret;
};

int CBinaryComtrade::InserirDados( CConfigComtrade* config, 
								int* numero_canais, 
								size_t tam, 
								void* canais,
                                size_t numpts, 
								long numero_ultima_amostra,
                                char* erro) {
    CCanalAnalogico *pAnalg;
    const CTaxas *pTaxa;
    char    cont;
    WORD    auxW;
    long    iMudancaTaxa, auxL;
    int     j, k, NumCanais, iTaxa=0;
	short   auxSh;
    void    *canais_;
    float   **fVetor;
    int     **iVetor;
    char    **cVetor;
    CCanal *pCanal;
    float fIntervalo, ganho[__num_max_canais], OffSet[__num_max_canais];
    double auxT;
    
    // Verifica se tipo do arquivo e' ASCII
    if (config->GetTipoArquivo()==0) 
		return CComtrade::InserirDados( config, 
									numero_canais, 
									tam, 
									canais,
									numpts, 
									numero_ultima_amostra, 
									erro);

    // Aloca buffers necessarios
    canais_ = calloc( tam, sizeof(void*));
    if (canais_==NULL) 
	{
        if (erro!=NULL) sprintf( erro, COMTRADE011);
        free(canais_); return 0;
    }
    memcpy( canais_, canais, tam*sizeof(void*));

    // Calcula instante inicial
    pTaxa = config->GetTaxa(1);
    fIntervalo =  (float)1e06 / pTaxa->GetTaxaAmostragem();
    iTaxa = 1;  
	auxT = 0;   
	iMudancaTaxa = 0;

    while (iMudancaTaxa<numero_ultima_amostra) 
	{
        pTaxa = config->GetTaxa(iTaxa++);
        fIntervalo = (float)1e06 / (pTaxa->GetTaxaAmostragem());
        auxT += min( (pTaxa->GetUltimaAmostra() - iMudancaTaxa),
                     numero_ultima_amostra)
                * fIntervalo;
        iMudancaTaxa = pTaxa->GetUltimaAmostra();
    };
    iMudancaTaxa = pTaxa->GetUltimaAmostra();
    if (iTaxa==1) iTaxa++;
    if (auxT >= 9999999999.5) 
	{
        if (erro==NULL) return 0;
        sprintf( erro, COMTRADE006);
        return 0;
    };

    // Calcula os valores para quantizar cada canal
    NumCanais=config->GetNumCanaisAnalogicos()+config->GetNumCanaisDigitais();
    
	for (j=1; j<=NumCanais; j++) 
	{
        pCanal = config->GetCanalPtr(j);
        if (pCanal->GetAnalogico()) 
		{
            pAnalg = (CCanalAnalogico*) pCanal;

            ganho[j-1]=(pAnalg->GetValorMaximoCfg()-pAnalg->GetValorMinimoCfg())/
                    COMTRADELIMITEDOQUANTIZADOR;

            if ((ganho[j-1]==0) && (pAnalg->GetFormatoDados()==e_real))
			{
                if (erro==NULL) return 0;
                sprintf( erro, COMTRADE001);
                return 0;
            };
            OffSet[j-1] =-(pAnalg->GetValorMaximoCfg()+pAnalg->GetValorMinimoCfg()) /
                        (2*ganho[j-1]);
        } 
		else {ganho[j-1]=1.F; OffSet[j-1]=0.F;};
    };

    // Escreve no arquivo de dados
    for (auxL=numero_ultima_amostra; 
		auxL < min ((long)numpts+numero_ultima_amostra,config->GetNumeroAmostrasPorCanal());
    	auxL++) 
	{
        // Realinha pointers e indices
        k = 0;  cont = 0;   auxW=0;
        fVetor = (float**) canais_;     iVetor = (int**) canais_;
        cVetor = (char**) canais_;

        // Le valores dos vetores fonte
        if (!InserirLong( auxL+1)) 
		{
			free(canais_); 
			return 0;
		};
        if (!InserirLong( (long)floor(auxT))) 
		{
			free(canais_); 
			return 0;
		};

        for (j=1; j<=NumCanais; j++) 
		{
            pCanal = config->GetCanalPtr(j);

            if (j==numero_canais[k]) 
			{
				if (pCanal->GetAnalogico()) {// analogico
					switch (((CCanalAnalogico*)pCanal)->GetFormatoDados()) {
					case e_inteiro: auxSh = (short) *(*iVetor)++;
					break;
					case e_real: auxSh = (short)
						floor( ((*(*fVetor)++)/ganho[j-1] + OffSet[j-1] + 0.5));
					break;
					default:
						if (erro!=NULL) sprintf( erro, COMTRADE009);
						free(canais_); 
						return 0;
					break;
					};

					if (!InserirInteiro( auxSh)) 
					{
						free(canais_); 
						return 0;
					};

				} 
				else //digital
				{ 
					auxW |= ((*(*cVetor)++) != 0) << cont++;
					if (cont == 16) 
					{
						if (!InserirInteiro( auxW)) 
						{
							free(canais_); 
							return 0;
						};
						auxW = 0; 
						cont = 0;
                    };
				};
				k++;
				} 
				else // canal nao interessa
				{ 
					if (pCanal->GetAnalogico()) auxSh = (short)0xFFFF;
					else 
					{
						auxW |= 1 << cont++;
                        if (cont == 16) 
						{
							if (!InserirInteiro( auxW))
							{
								free(canais_);
								return 0;
							};
                            auxW = 0; 
							cont = 0;
                        };
					};
				};
        		if (k>=(int)tam) break;
				fVetor++; 
				iVetor++; 
				cVetor++;
		}; // for j

        if (cont != 0) 
		{
            if (!InserirInteiro( auxW))
			{
				free(canais_);
				return 0;
			};

            auxW = 0; 
			cont = 0;
        };
        if (iMudancaTaxa<=auxL+1) 
		{
            pTaxa = config->GetTaxa(iTaxa++);
            if (pTaxa) 
			{
            	iMudancaTaxa = pTaxa->GetUltimaAmostra();
            	fIntervalo = (float)1e06 / pTaxa->GetTaxaAmostragem();
            };
		};
        auxT += fIntervalo;
        if (auxT >= 9999999999.5) 
		{
            if (erro==NULL) return 0;
            sprintf( erro, COMTRADE006);
            return 0;
        };
	}; // for auxL

    free(canais_);
	return 1;
};



