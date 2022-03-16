/*****************************************************************************

   Projeto : SISTEMA INTEGRADO DE APOIO `A ANALISE DE PERTURBACOES - SINAPE
   Nome    : COMTRADE.CPP
   Funcao  : Rotinas para conversão entre formato interno (CConfigCOMTRADE)
            e o formato COMTRADE (IEEE standard Common Format for Transient
            Data Exchange for power systems), IEEE C37.111-1991
   C       : C++
   Autor   : MAMR
   Data    : nov/94
   Obs.	   : Modificacao para leitura de campos vazios no COMTRADE
   Data    : fev/96

*****************************************************************************/

#include <math.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <iostream>

#include "comtrade.h"   

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
#define COMTRADE011 "Tipo do arquivo deve ser \"ASCII\" ou \"BINARY\" !"
#define COMTRADE012 "Nao ha' memoria para ler dados !"
#define COMTRADE013 "\nCampo do numero total de canais deve ser preenchido !"
#define COMTRADE014 "\nInformacao dos campos de numero de canais incoerentes !"
#define COMTRADE015 "Tem que informar horario da primeira amostra ou do disparo"  

#ifndef max
#define max(a,b)        (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
#define min(a,b)        (((a) < (b)) ? (a) : (b))
#endif

#define ERROFORMATO()\
{\
    if (erro==NULL) return 0;\
    sprintf(erro,"%s %s\nLinha:%d",GetNomeArquivo(),COMTRADE007,linha);\
    return 0;\
}

 
char* ExtrairVirgula( char* string)
{
    char *ch;
    while ((ch=strchr( string, ','))!=NULL) *ch = '.';
    return string;
};

int CComtrade::ExtrairConfiguracao(CConfigComtrade* config, char* erro) {

    char string[256];
    int linha = 1;

    // Extrai o nome da subestacao
    if (!ExtrairLexema( string, 256, ',')) ERROFORMATO();
    config->SetNomeSubstacao(string);

	// Extrai o identificador do RDP
    if (!ExtrairLinha( string, 256)) ERROFORMATO();
    config->SetIdentificadorRDP(string);

	return ExtrairConfiguracao_(config,erro);
};

//////////////////////////////////////////////////////////////////////////////
int CComtrade::LerLinha2(CConfigComtrade* config, char* erro)
{
	char    string[256], *pos;
	int     linha = 1, cont, auxI, ret, NTotCan;
    
	// O q eh isso?
    // New
    int 	m_ErroData = 0;
	
	// Extrai o numero total de canais
 	if (!ExtrairLexema( string, 256, ',')) ERROFORMATO();
    if (strlen(string)<1) 
	{
    	if (erro==NULL) return 0;
        strcpy( erro, COMTRADE013);
        return 0; 
	};
    ret = sscanf( string, "%d", &NTotCan);
   	if ((ret<=0) || (ret==EOF) || (NTotCan>__num_max_canais)) ERROFORMATO();

	// Extrai o numero de canais analogicos
    if (!ExtrairLexema( string, 256, ',')) ERROFORMATO();
    if (strlen(string)>0) 
	{
    	if ((pos = strchr( string, 'A'))==NULL) ERROFORMATO();
    	pos[0] = 0;
    	ret = sscanf( string, "%d", &auxI);
    	if ((ret==0) || (ret==EOF) || (auxI>__num_max_canais)) ERROFORMATO();
    } 
	else auxI = 0;

	for (cont=1; cont<=auxI; cont++)
	{
        if (!config->AcrescentarCanalAnalogico()) 
		{
            if (erro==NULL) return 0;
            strcpy( erro, COMTRADE010);
            return 0; 
		};
	}

	// Extrai o numero de canais digitais
    if (!ExtrairLinha( string, 256)) ERROFORMATO();
    if (strlen(string)>0) 
	{
    	if ((pos = strchr( string, 'D'))==NULL) ERROFORMATO();
    	pos[0] = 0;
	    ret = sscanf( string, "%d", &auxI);
	    if ((ret==0) || (ret==EOF) ||
	        (auxI>(__num_max_canais-config->GetNumCanaisAnalogicos()))) 
	        ERROFORMATO();
    } 
	else auxI = 0;

    for (cont=1; cont<=auxI; cont++)
	{
        if (!config->AcrescentarCanalDigital()) 
		{
            if (erro==NULL) return 0;
            strcpy( erro, COMTRADE010);
            return 0; 
		};
	}

    if (NTotCan!=config->GetNumCanaisAnalogicos()+config->GetNumCanaisDigitais()) 
	{
		if (erro==NULL) return 0;
        strcpy( erro, COMTRADE014);
		return 0; 
	};

	return ret;
};

///////////////////////////////////////////////////////////////////////

int CComtrade::LerLinhaCanalAnalogico(CConfigComtrade* config, 
									  CCanal *canal, 
									  int cont, 
									  char* erro)
{
    char    string[256];
    int     versao, ret, linha = 1, auxI;	
    long    auxL;
    CCanalAnalogico *analog;
    float auxF;
	
	sret = NULL;
	
	// Extrai o indice do canal analogico (An)
	if (!ExtrairLexema( string, 256, ',')) ERROFORMATO();
    ret = sscanf( string, "%d", &auxI);
    if ((ret<=0) || (ret==EOF)) auxI = cont;
    canal->SetNumeroCanal(auxI);
	
	// Extrai o identificador do canal (ch_id)
    if (!ExtrairLexema( string, 256, ',')) ERROFORMATO();
    canal->SetNomeCanal(string);

    if (canal->GetAnalogico()) 
	{
		analog = (CCanalAnalogico*) canal;
        
		versao = config->GetVersao();
		// Para as versoes 1997 e 1999
		if (versao > 1991)
		{	
			// Extrai a fase do canal (ph)
			if (!ExtrairLexema( string, 256, ',')) ERROFORMATO();
			if (strlen(string)>__num_char_fase)  ERROFORMATO();
			canal->SetFaseCanal(string);
			// Extrai equipamento monitorado (ccbm)
			if (!ExtrairLexema( string, 256, ',')) ERROFORMATO();
			canal->SetEquipamentoCanal(string);
		}
		
		// O q eh isso?
		analog->SetFormatoDados(m_eFormatoVetorExtraido);
        
		if (versao == 1991)	
		{
			// Extrai a fase do canal (ph) 
			if (!ExtrairLexema( string, 256, ',')) ERROFORMATO();
			if (strlen(string)>__num_char_fase)  ERROFORMATO();
			analog->SetFaseCanal(string);

			// Extrai o equipamento monitorado (ccbm)
			if (!ExtrairLexema( string, 256, ',')) ERROFORMATO();
			analog->SetEquipamentoCanal(string);
        }
		
		// Extrai a unidade do canal (uu)
		if (!ExtrairLexema( string, 256, ',')) ERROFORMATO();
		if (sret==NULL) sret = (char *) calloc(256,sizeof(char));
		double fator_mult = analog->CorrigeUnidade(string, sret);
		analog->SetUnidadeCanal(sret);
        
		// Extrai o coeficiente angular (a)
		if (!ExtrairLexema( string, 256, ',')) ERROFORMATO();
        ret = sscanf( string, "%f", &auxF);
        if ((ret<=0) || (ret==EOF)) auxF = 1.F;
        analog->SetCoefAng((float)(auxF*fator_mult));
		
		// Extrai o coeficiente linear (b)
        if (!ExtrairLexema( string, 256, ',')) ERROFORMATO();
        ret = sscanf( string, "%f", &auxF);
        if ((ret<=0) || (ret==EOF)) auxF = 0.F;
        if (auxF==0.F)
			analog->SetCoefLin(0.F);
		else
			analog->SetCoefLin((float)(auxF*fator_mult));
		
		// Extrai o time skew (skew)
        if (!ExtrairLexema( string, 256, ',')) ERROFORMATO();
        ret = sscanf( string, "%f", &auxF);
        if ((ret<=0) || (ret==EOF)) auxF = 0.F;
        analog->SetTimeSkew((float)(auxF*1E-06));
		
		// Extrai o valor minimo das amostras do canal (min)
        if (!ExtrairLexema( string, 256, ',')) ERROFORMATO();
        ret = sscanf( string, "%ld", &auxL);
        if ((ret<=0) || (ret==EOF)) auxL = -99999;
        analog->SetValorMinimo((float)auxL);
		
		// Extrai valor maximo das amostras do canal (max)
		if (versao == 1991)	
		{
			if (!ExtrairLinha( string, 256)) ERROFORMATO();
		}
		if (versao > 1991)	// (versao == 1997) || (versao == 1999)
		{	
			if (!ExtrairLexema( string, 256, ',')) ERROFORMATO();
		}
		ret = sscanf( string, "%ld", &auxL);
		if ((ret<=0) || (ret==EOF)) auxL = 99999;
		analog->SetValorMaximo((float)auxL);

        analog->AcertarTipoDoCanalAnalogico();
	}

	free(sret);

	return ret;
};

///////////////////////////////////////////////////////////////////////

int CComtrade::LerLinhaCanalDigital(CConfigComtrade* config, CCanal *canal, int cont, char* erro)
{
    char    string[256];
    int     auxI, versao, ret, linha = 1;	//cont, 

	versao = config->GetVersao();
	
	// Extrai indice do canal
	if (!ExtrairLexema( string, 256, ',')) ERROFORMATO();
   	ret = sscanf( string, "%d", &auxI);
   	if ((ret<=0) || (ret==EOF)) auxI = cont;
	canal->SetNumeroCanal(auxI);

	// Extrai nome do canal
	if (!ExtrairLexema( string, 256, ',')) ERROFORMATO();
	canal->SetNomeCanal(string);
	
	// Extrai o estado normal do status do canal (y)
	if (!ExtrairLinha( string, 256)) ERROFORMATO();
	if (versao > 1991)	// (versao == 1997) || (versao == 1999)
	{	
		if (strlen(string)>1) ERROFORMATO();
    }
	if (versao==1991) ExtrairEspacosDasBordas(string);			
	canal->SetTipoDoCanal(atoi(&string[0]));
	
	return ret;
};

///////////////////////////////////////////////////////////////////////

int CComtrade::ExtrairConfiguracao_(CConfigComtrade* config, char* erro) 
{
    char    string[256], *pos, *pos2,
    		HoraPrimPtValid, DiaPrimPtValid, HoraTrigValid, DiaTrigValid;
    int     auxI, cont, ret, linha = 1;
    long    auxL;
    CCanal *canal;
    float auxF;
    __instante *tTrigger, *tPrima;
    
	// O q eh isso?
    // New
    int 	m_ErroData = 0;
	
	// O q eh isso?
    // Valor default pois Simulacao e Comprimido nao fazem parte do COMTRADE
    config->SetSimulacao();
    config->SetComprimido();

    // Numero de canais
    linha++;
	
	// Extrai a segunda linha do arquivo ".cfg" contendo:
	// * Numero total de canais
	// * Numero de canais analogicos
	// * Numero de canais digitais
	auxI = LerLinha2(config, erro);	
	if (auxI==0) return 0;
	
	// Extrai as caracteristicas de cada canal tanto analog.
	// quanto digitais
	// * Analog -> An,ch_id,ph,ccbm,uu,a,b,skew,min,max
	// * Digital -> Dn,ch_id,y
    for (cont=1;
         cont<=config->GetNumCanaisAnalogicos()+config->GetNumCanaisDigitais();
         cont++) 
	{
        canal = config->GetCanalPtr(cont);
        ASSERT(canal!=NULL);
        linha++;

		if (canal->GetAnalogico())	
			LerLinhaCanalAnalogico(config, canal, cont, erro);
		else	
			LerLinhaCanalDigital(config,canal, cont, erro);
    
	};//for cont

    linha++;
	//Extrai frequencia da linha (lf)
    if (!ExtrairLinha( string, 256)) ERROFORMATO();
    ret = sscanf( string, "%f", &auxF);
    if ((ret<=0) || (ret==EOF)) ERROFORMATO();
	if (auxF==0)	auxF=60;
    config->SetFreqLinha(auxF);

    linha++;
	// Extrai o numero taxas de amostragem (nrates)
    if (!ExtrairLinha( string, 256)) ERROFORMATO();
    ret = sscanf( string, "%d", &auxI); 
    if ((ret<=0) || (ret==EOF) || (ret>__num_max_taxas) ) ERROFORMATO();
    
	float auxF1[100];
	long auxL1 [100];
	for(cont=0; cont<100; cont++)	
	{
		auxF1[cont] = 0;
		auxL1[cont] = 0;
	}
	
	// Extrai as taxas de amostragem e o valor da ultima
	// amostra (samp,endsamp) onde auxI = nrates
	for (cont=1; cont<=auxI; cont++) 
	{
        linha++;

        if (!ExtrairLexema( string, 256, ',')) ERROFORMATO();
        ExtrairEspacosDasBordas(string);
	    ret = sscanf( string, "%f", &auxF);   
		if ((ret<=0) || (ret==EOF)) ERROFORMATO();
        auxF1[cont-1] = auxF;
		
		if (!ExtrairLinha( string, 256)) ERROFORMATO();
        ExtrairEspacosDasBordas(string);
		ret = sscanf( string, "%ld", &auxL);   
        if ((ret<=0) || (ret==EOF)) ERROFORMATO();
		auxL1[cont-1] = auxL;

		config->AcrescentarTaxa(auxF1[cont-1],auxL1[cont-1]);
	}; //for cont

    linha++;	
	// Extrai data da primeira amostra
	DiaPrimPtValid = 1;
    tPrima = config->GetInstantePrimeiraAmostra();
    if (!ExtrairLexema( string, 256, ',')) ERROFORMATO();
    ExtrairEspacosDasBordas(string);

    if (strlen(string)>0) 
	{
		// Extrai o mes
    	if ((pos = strchr( string, '/'))==NULL) ERROFORMATO();
    	if (!pos) ERROFORMATO();
    	pos[0] = 0;
	    auxI = 0;	
		ret = sscanf( string, "%d", &auxI);
	    
	    if ((ret<=0) || (ret==EOF) || (auxI<1) || (auxI>12))   
		{
			if (auxI>12)	
			{	
	    		tPrima->tm_mday = auxI; 
	    		m_ErroData = 1;   // TRUE
	    	}
			else tPrima->tm_mon = --auxI;    

	    	if ((ret<=0) || (ret==EOF) || (auxI<1))	                     
				ERROFORMATO();                      

	    }
		else tPrima->tm_mon = --auxI;
    	
		// Extrai o dia
    	if ((pos2 = strchr( ++pos, '/'))==NULL) ERROFORMATO();
    	if (!pos2) ERROFORMATO();
    	pos2[0] = 0;
	    auxI = 0;
		ret = sscanf( pos, "%d", &auxI);
	    if ((ret<=0) || (ret==EOF) || (auxI<1) || (auxI>31)) 
			ERROFORMATO();
	    
		if (m_ErroData == 1)    // TRUE 
			tPrima->tm_mon = --auxI;
	    else 		
			tPrima->tm_mday = auxI;
    	
       	// Extrai o ano
	    auxI = 0;	
		ret = sscanf( ++pos2, "%d", &auxI);
	    if ((ret<=0) || (ret==EOF) || (auxI<0) || (auxI>99)) 
		{
	    	if ((auxI<0) || (auxI>99))	
			{
		    	if ((auxI >=2000) && (auxI<2070))
		    		auxI = abs(2000-auxI);
		    	if ((auxI>1969) && (auxI<=1999))
		    		auxI = abs(1900-auxI);  
		    }
	    	if ((ret<=0) || (ret==EOF))
	    		ERROFORMATO();
	    }
	    tPrima->tm_year = auxI;  
	} 
	else DiaPrimPtValid = 0;

	// Extrai hora da primeira amostra
    if (!ExtrairLinha( string, 256)) ERROFORMATO();
    ExtrairEspacosDasBordas(string);	
	HoraPrimPtValid = 1;

    if (strlen(string)>0) 
	{
    	// Extrai a hora
    	if ((pos = strchr( string, ':'))==NULL) ERROFORMATO();
    	if (!pos) ERROFORMATO();
    	pos[0] = 0;
	    auxI = 0;	
		ret = sscanf( string, "%d", &auxI);
	    if ((ret<=0) || (ret==EOF) || (auxI<0) || (auxI>23)) ERROFORMATO();
	    tPrima->tm_hour = auxI;
    	
		// Extrai os minutos
    	if ((pos2 = strchr( ++pos, ':'))==NULL) ERROFORMATO();
    	if (!pos2) ERROFORMATO();
    	pos2[0] = 0;
	    auxI = 0;	ret = sscanf( pos, "%d", &auxI);
	    if ((ret<=0) || (ret==EOF) || (auxI<0) || (auxI>59)) ERROFORMATO();
	    tPrima->tm_min = auxI;
    	
		// Extrai os segundos
		double auxD;
		ret = sscanf( ++pos2, "%lf", &auxD);
	    if ((ret<=0) || (ret==EOF) || (auxD<0) || (auxD>=60)) ERROFORMATO();
	    tPrima->tm_sec = auxD; 
	} 
	else HoraPrimPtValid = 0;

    
    linha++;
	// Extrai a data do trigger
	DiaTrigValid = 1;
    tTrigger = config->GetInstanteTrigger();
    if (!ExtrairLexema( string, 256, ',')) ERROFORMATO();
    ExtrairEspacosDasBordas(string);

    if (strlen(string)>0) 
	{
    	//Extrai o mes
    	if ((pos = strchr( string, '/'))==NULL) ERROFORMATO();
    	if (!pos) ERROFORMATO();
    	pos[0] = 0;
	    auxI = 0;	
		ret = sscanf( string, "%d", &auxI);
	    if ((ret<=0) || (ret==EOF) || (auxI<1) || (auxI>12))   {
	    	if (auxI>12)	{	
	    		tTrigger->tm_mday = auxI; 
	    		m_ErroData = 1;   // TRUE
	    	}else 
	    		tTrigger->tm_mon = --auxI;    
	    	 if ((ret<=0) || (ret==EOF) || (auxI<1))	                     
	    		ERROFORMATO();                      
	    }
		else tTrigger->tm_mon = --auxI;
    	
    	// Extrai o dia
    	if ((pos2 = strchr( ++pos, '/'))==NULL) ERROFORMATO();
    	if (!pos2) ERROFORMATO();
    	pos2[0] = 0;
	    auxI = 0;	ret = sscanf( pos, "%d", &auxI);
	    if ((ret<=0) || (ret==EOF) || (auxI<1) || (auxI>31)) ERROFORMATO();
	    
	    if (m_ErroData == 1)     // TRUE
	       	tTrigger->tm_mon = --auxI;
	    else
	    	tTrigger->tm_mday = auxI;
	    		
    	//Extrai o ano
	    auxI = 0;	
		ret = sscanf( ++pos2, "%d", &auxI);
	    if ((ret<=0) || (ret==EOF) || (auxI<0) || (auxI>99))	
		{
	    	if ((auxI<0) || (auxI>99))	
			{
		    	if ((auxI >=2000) && (auxI<2070))
		    		auxI = abs(2000-auxI);
		    	if ((auxI>1969) && (auxI<=1999))
		    		auxI = abs(1900-auxI);  
		    }
	    	if ((ret<=0) || (ret==EOF))
	    		ERROFORMATO();
	    }
	   	tTrigger->tm_year = auxI; 
	} 
	else DiaTrigValid = 0;

	// Extrai a hora do trigger
    if (!ExtrairLinha( string, 256)) ERROFORMATO();
    ExtrairEspacosDasBordas(string);	
	HoraTrigValid = 1;

    if (strlen(string)>0) 
	{
    	// Extrai a hora
    	if ((pos = strchr( string, ':'))==NULL) ERROFORMATO();
    	if (!pos) ERROFORMATO();
    	pos[0] = 0;
	    auxI = 0;	
		ret = sscanf( string, "%d", &auxI);
	    if ((ret<=0) || (ret==EOF) || (auxI<0) || (auxI>23)) ERROFORMATO();
	    tTrigger->tm_hour = auxI;

    	// Extrai os minutos
    	if ((pos2 = strchr( ++pos, ':'))==NULL) ERROFORMATO();
    	if (!pos2) ERROFORMATO();
    	pos2[0] = 0;
	    auxI = 0;	
		ret = sscanf( pos, "%d", &auxI);
	    if ((ret<=0) || (ret==EOF) || (auxI<0) || (auxI>59)) ERROFORMATO();
	    tTrigger->tm_min = auxI;

    	// Extrai os segundos
		double auxD;
	    ret = sscanf( ++pos2, "%lf", &auxD);
	    if ((ret<=0) || (ret==EOF) || (auxD<0) || (auxD>=60)) ERROFORMATO();
	    tTrigger->tm_sec = auxD;
	} 
	else HoraTrigValid = 0;
	
	// Executa a validacao dos horarios
	if ((!HoraPrimPtValid) && (!HoraTrigValid)) 
	{
    	if (erro==NULL) return 0;
        strcpy( erro, COMTRADE015);
        return 0;
    };
	if (HoraPrimPtValid && (!HoraTrigValid)) 
	{
		tTrigger->tm_hour = tPrima->tm_hour;
		tTrigger->tm_min = tPrima->tm_min;
		tTrigger->tm_sec = tPrima->tm_sec;
	};
	if (HoraTrigValid && (!HoraPrimPtValid)) 
	{
		tPrima->tm_hour = tTrigger->tm_hour;
		tPrima->tm_min = tTrigger->tm_min;
		tPrima->tm_sec = tTrigger->tm_sec;
	};
	if ((!DiaPrimPtValid) && (!DiaTrigValid)) 
	{
    	if (erro==NULL) return 0;
        strcpy( erro, COMTRADE015);
        return 0;
	};
	if (DiaPrimPtValid && (!DiaTrigValid)) 
	{
		tTrigger->tm_mon = tPrima->tm_mon;
		tTrigger->tm_mday = tPrima->tm_mday;
		tTrigger->tm_year = tPrima->tm_year;
	};
	if (DiaTrigValid && (!DiaPrimPtValid)) 
	{
		tPrima->tm_mon = tTrigger->tm_mon;
		tPrima->tm_mday = tTrigger->tm_mday;
		tPrima->tm_year = tTrigger->tm_year;
	};
    
    linha++;
	// Extrai o tipo do arquivo
    if (!ExtrairLinha( string, 256)) ERROFORMATO();
    ExtrairEspacosDasBordas(string);

    if (!strcmp( string, "ASCII")) config->SetTipoArquivo(0);
    else if (!strcmp( string, "BINARY")) config->SetTipoArquivo(1);
    else 
	{ 
		if (erro!=NULL) 
		{
			sprintf(erro,"Linha %d: %s",linha,COMTRADE011); 
			
		}
		return 0;
	};

    return 1;
};

///////////////////////////////////////////////////////////////

int CComtrade::ExtrairDados(CConfigComtrade* config, 
							int* numero_canais, 
							size_t tam, 
							void* canais,
							size_t numpts, 
							long offset,
							char AcessoRandomico, 
							char* erro) 
{
    char    string[256], *auxC;
    size_t  i;
    int     j, k, ret, NumCanal;
    long    auxL, linha;
    float   **fVetor = (float**) canais;
    int     **iVetor = (int**) canais;
    char    **cVetor = (char**) canais;
    CCanal *pCanal; 
    
    // Posiciona o buffer na posicao correta
    if (AcessoRandomico) 
	{
        ResetBuffer();
        auxL = 20 + config->GetNumCanaisAnalogicos()*8
                + config->GetNumCanaisDigitais()*2;
        auxC = (char*) calloc( (size_t)auxL+1, 1);
        if (auxC==NULL) {
            if (erro!=NULL) sprintf( erro, COMTRADE012);
            return 0;
          };
        for ( linha=1; linha<=offset; linha++)
            if (!ExtrairLinha( auxC, (size_t)auxL)) ERROFORMATO();
        free(auxC);
    };

    linha = offset; 
    
    // Le arquivo de dados
    NumCanal=config->GetNumCanaisAnalogicos()+config->GetNumCanaisDigitais(); 
    
	for (i=0; 
		(i<numpts) && (linha<config->GetNumeroAmostrasPorCanal()); 
		i++) 
	{
		linha++;    
        // Extrai o numero da amostra (n)
        if (!ExtrairLexema( string, 255, ','))	 
		{		
        	long linha1 = (long)linha-1;   
        	config->AlterarTaxa(config->GetNumeroTaxas(),-1.F,linha1);
        	numpts = (size_t)linha1;
        	return 1;   
        }
	   	ret = sscanf( string, "%ld", &auxL);    
        if ((ret==0) || (ret==EOF)) ERROFORMATO();    
        
		//O q eh isso?
        // New
        if ( m_iVerificaAmostra == 1)	
		{
        	if ( ((long)i+offset)!=auxL)	
			{
				if ( ((long)i+1+offset)!=auxL) ERROFORMATO();           
        	}
        }
        
        // Extrai o tempo em microsegundos da amostra (timestamp)
        if (!ExtrairLexema( string, 255, ',')) ERROFORMATO();
        ret = sscanf( string, "%ld", &auxL);
		
		// Extrai as amostras dos canais
        k = 0;
        for (j=1; j<=NumCanal; j++) 
		{
            pCanal = config->GetCanalPtr(j);

			if (j==NumCanal) 
                ExtrairLinha( string, 255);
			else 
				ExtrairLexema( string, 255, ',');
			
			ret = sscanf( string, "%ld", &auxL);
            
            if (j==numero_canais[k]) 
			{
				if (pCanal->GetAnalogico()) // analogico      
				{
              		switch (m_eFormatoVetorExtraido) {    
					case e_inteiro: 
						if ((auxL>INT_MAX) || (auxL<INT_MIN)) 
						{
							if (erro==NULL) return 0;
							sprintf( erro, COMTRADE008, INT_MIN, INT_MAX);
							return 0;
						};
						(iVetor[k++])[i] = (int) auxL;
					break;
					case e_real: (fVetor[k++])[i] = (float)auxL;
					break;
					default:
						if (erro==NULL) return 0;
						sprintf( erro, COMTRADE009);
						return 0;
					break;
					};
                } 
				else (cVetor[k++])[i] = (char) auxL; //digital
            } 
			
		}; // for j
	}; // for i
    return 1;
};

///////////////////////////////////////////////////////////////////////

int CComtrade::InserirConfiguracao_(CConfigComtrade* config, char* erro)
{
    char string[256];
    int pos, cont, pos2;
    float ganho, OffSet;
    CCanal *canal;
    CCanalAnalogico *analog;
    const CTaxas *taxa;
    __instante* instante;
	
	// Insere: 
	// * nome da subestacao
	// * identificador do RDP
    strcpy( string, config->GetNomeSubstacao());
    ExtrairVirgula( string);
    strcat( string, ",");
    pos = strlen(string);
    strcat( string, config->GetIdentificadorRDP());
	if (!InserirLinha(string)) return 0;

	// Insere:
	// * numero total de canais
	// * numero de canais analogicos
	// * numero de canais digitais
    sprintf( string, 
			"%d,%dA,%dD",
            config->GetNumCanaisAnalogicos()+config->GetNumCanaisDigitais(),
            config->GetNumCanaisAnalogicos(),
			config->GetNumCanaisDigitais());
    if (!InserirLinha(string)) return 0;
	
	// Insere as informacoes dos canais analogicos e digitais
    for (cont=1;
         cont<=config->GetNumCanaisAnalogicos()+config->GetNumCanaisDigitais();
         cont++) 
	{
        canal = config->GetCanalPtr(cont);
        ASSERT(canal!=NULL);

		// Insere o indice e nome do canal (An ou Dn, ch_id)
        pos2 = pos = sprintf(string,"%d,", canal->GetNumeroCanal());
        pos+=sprintf(string+pos,"%s",canal->GetNomeCanal());
        ExtrairVirgula(&string[pos2]); 
		strcat( string,",");
        if (!InserirString(string)) return 0;
        
		if (canal->GetAnalogico())  // Para os canais analogicos
		{
            analog = (CCanalAnalogico*) canal;

			// Insere a fase do canal (ph) e o equipamento monitorado (ccbm)
			strcpy(string,canal->GetFaseCanal());
			pos2 = pos = strlen(string)>0;
            pos+=sprintf(string+pos,",%s",analog->GetEquipamentoCanal());
            ExtrairVirgula(&string[++pos2]);
            if (!InserirString(string)) return 0;
			
			// Insere a unidade do canal (uu)
            pos =sprintf(string,",%s",analog->GetUnidadeCanal());

            ExtrairVirgula(&string[1]); 
			strcat( string, ","); 
			pos++;
            if (analog->GetValorMinimoCfg()==analog->GetValorMaximoCfg()) 
			{
                if (erro==NULL) return 0;
                sprintf( erro, COMTRADE001);
                return 0;
            };

			// Insere os coeficiente lineares e angulares (a,b)
            switch (analog->GetFormatoDados()) {
              case e_inteiro:
                pos+=sprintf(string+pos,"%G,",analog->GetCoefAng(analog->GetTipoEscalamento()));
                pos+=sprintf(string+pos,"%G,",analog->GetCoefLin(analog->GetTipoEscalamento()));
                break;
              case e_real:// sinal sera' quantizado
                ganho = COMTRADELIMITEDOQUANTIZADOR / 
					(analog->GetValorMaximoCfg()-analog->GetValorMinimoCfg());
                OffSet = (analog->GetValorMaximoCfg()+analog->GetValorMinimoCfg())/2;
                pos+=sprintf(string+pos,"%G,",
					analog->GetCoefAng(analog->GetTipoEscalamento())/ganho);
                pos+=sprintf(string+pos,"%G,",
					analog->GetCoefLin(analog->GetTipoEscalamento()) + analog->GetCoefAng(analog->GetTipoEscalamento())*OffSet);
                break;
            };

			// Insere o time skew (skew)
            pos+=sprintf(string+pos,"%G,",analog->GetTimeSkew()*1E06);

			// Insere os valores minimo  e maximo (min,max)
            switch (analog->GetFormatoDados()) {
              case e_inteiro:
                pos+=sprintf(string+pos,"%.0f,",analog->GetValorMinimoCfg());
                pos+=sprintf(string+pos,"%.0f",analog->GetValorMaximoCfg());
                break;
              case e_real:// sinal sera' quantizado
                pos+=sprintf(string+pos,"%ld,%ld",-COMTRADELIMITEDOQUANTIZADOR/2,
                                COMTRADELIMITEDOQUANTIZADOR/2);
                break;
            };
        } 
		else // canal digital
		{  
			// Insere estado normal (y)
			sprintf(string,"%d",canal->GetTipoDoCanal());
		}

        if (!InserirLinha(string)) return 0;

    };//for cont

	
	// Insere a frequencia da linha (lf)
    sprintf(string,"%d",(int)config->GetFreqLinha());
    if (!InserirLinha(string)) return 0;
    
	// Insere o numero de taxas de amostragem (nrates)
    sprintf(string,"%d",config->GetNumeroTaxas());
    if (!InserirLinha(string)) return 0;
    
	// Insere cada taxa de amostragem (em Hz) e sua correspondente 
	// ultima amostra (samp e endsamp)
    for (cont=1; cont<=config->GetNumeroTaxas(); cont++) {
        taxa = config->GetTaxa(cont);
        pos =sprintf(string,"%G,", taxa->GetTaxaAmostragem());
        pos+=sprintf(string+pos,"%li",taxa->GetUltimaAmostra());
        if (!InserirLinha(string)) return 0;
      }; //for cont
	
	// Insere a data e hora da primeira amostra
    instante = config->GetInstantePrimeiraAmostra(); 
    pos =sprintf(string,"%02d/", instante->tm_mon+1);
    pos+=sprintf(string+pos,"%02d/",instante->tm_mday);
    pos+=sprintf(string+pos,"%02d,",instante->tm_year);
    pos+=sprintf(string+pos,"%02d:",instante->tm_hour);
    pos+=sprintf(string+pos,"%02d:",instante->tm_min);
    pos+=sprintf(string+pos,"%02.6f",instante->tm_sec);
    if (!InserirLinha(string)) return 0;
	
	// Insere a data e hora do TRIGGER
    instante = config->GetInstanteTrigger();  
    pos =sprintf(string,"%02d/", instante->tm_mon+1);
    pos+=sprintf(string+pos,"%02d/",instante->tm_mday);
    pos+=sprintf(string+pos,"%02d,",instante->tm_year);
    pos+=sprintf(string+pos,"%02d:",instante->tm_hour);
    pos+=sprintf(string+pos,"%02d:",instante->tm_min);
    pos+=sprintf(string+pos,"%02.6f",instante->tm_sec);
    if (!InserirLinha(string)) return 0;

    return 1;
};

int CComtrade::InserirConfiguracao(CConfigComtrade* config, char* erro){
    int ret;

    ret = CComtrade::InserirConfiguracao_( config, erro);
	
	// Insere tipo do arquivo de dados 
	if (!InserirLinha("ASCII")) return 0;
    return ret;
};

int CComtrade::InserirDados( CConfigComtrade* config, 
							int* numero_canais, 
							size_t tam, 
							void* canais,
							size_t numpts, 
							long numero_ultima_amostra,
							char* erro) {
    int j, k, iTaxa=0;         
    long iMudancaTaxa, auxL;
    const CTaxas *pTaxa;
    CCanal *pCanal;
    CCanalAnalogico *pAnalg;
    float   **fVetor = (float**) canais;
    int     **iVetor = (int**) canais;
    char    **cVetor = (char**) canais;
    char string[80];
    float fIntervalo, ganho[__num_max_canais], OffSet[__num_max_canais];
    double auxT;
    
    // calcula instante inicial
    pTaxa = config->GetTaxa(1); 
	ASSERT(pTaxa);
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

    // calcula os valores para quantizar cada canal
    for (j=1;
         j<=config->GetNumCanaisAnalogicos()+config->GetNumCanaisDigitais();
         j++) 
	{
        pCanal = config->GetCanalPtr(j);

        if (pCanal->GetAnalogico()) 
		{
            pAnalg = (CCanalAnalogico*) pCanal;

            ganho[j-1]=(pAnalg->GetValorMaximoCfg()-pAnalg->GetValorMinimoCfg())/
                    COMTRADELIMITEDOQUANTIZADOR;

            OffSet[j-1] =-(pAnalg->GetValorMaximoCfg()+pAnalg->GetValorMinimoCfg()) /
                        (2*ganho[j-1]);

            if ((ganho[j-1]==0) && (pAnalg->GetFormatoDados()==e_real))
			{
                if (erro==NULL) return 0;
                sprintf( erro, COMTRADE001);
                return 0;
            };
        } 
		else 
		{ 
			ganho[j-1]=1.F; 
			OffSet[j-1] = 0.F;
		};

	};

    // Escreve arquivo de dados
    for (auxL=0; 
		(auxL<(long)numpts) && (auxL+numero_ultima_amostra<config->GetNumeroAmostrasPorCanal());
    	auxL++) 
	{
		// Insere numero da amostra e o timestamp (n,timestamp)
        sprintf( string,"%010ld,%010.0f",auxL+numero_ultima_amostra+1,floor(auxT));
        if (!InserirString(string)) return 0;

		// Insere as amostras dos canais analogicos e digitais por linha
        k = 0;
        for (j=1;
			j<=config->GetNumCanaisAnalogicos()+config->GetNumCanaisDigitais();
			j++) 
		{
            pCanal = config->GetCanalPtr(j);

            if (j==numero_canais[k]) 
			{
				if (pCanal->GetAnalogico()) // analogico
				{
					pAnalg = (CCanalAnalogico*) pCanal;

					switch (pAnalg->GetFormatoDados()) {
					case e_inteiro: sprintf( string,",%06d",(iVetor[k++])[auxL]);
					break;
					case e_real: // quantiza os valores
						sprintf( string,",%06.0f", 
								floor( (fVetor[k++])[auxL] / ganho[j-1]+OffSet[j-1]+0.5));
					break;
					default:
						if (erro==NULL) return 0;
						sprintf( erro, COMTRADE009);
						return 0;
					break;
					};
				} 
				else sprintf( string,",%d",(cVetor[k++])[auxL]); //digital

			}
			else // vetor nao valido 
			{	
				sprintf( string,"%s", ",999999");
			};

            if (!InserirString(string)) return 0;

        	if (k>=(int)tam) break;

		}; // for j

        if (!InserirLinha("")) return 0;

        if (iMudancaTaxa<=auxL+numero_ultima_amostra+1) 
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

    return 1;

};

int CComtrade::InserirCabecalho( char* texto, size_t tam, char *erro) {
    
	if (InserirBytes( texto, tam)==0)
	{
        if (erro==NULL) return 0;
        sprintf( erro, COMTRADE002);
        return 0;
    };
    return SalvarBuffer( 0, erro);
};

int CComtrade::ExtrairCabecalho( char* texto, size_t *tam, char *erro) {
    
	strcpy( erro, "\0");
	long tamanho = filesize(erro);

    if (strlen(erro)>0) return 0;
    if (tamanho>UINT_MAX) 
	{
        if (erro==NULL) return 0;
        sprintf(erro,"%s %d",COMTRADE003,UINT_MAX);
        return 0;
    } 
	else *tam = (size_t) tamanho;

    texto = (char*) malloc(*tam);
    if (texto==NULL) 
	{
        if (erro==NULL) return 0;
        sprintf( erro, "%s %s",COMTRADE004,GetNomeArquivo());
        return 0;
    };
    if (ExtrairBytes( texto, *tam)==0) 
	{
        if (erro==NULL) return 0;
        sprintf(erro,"%s %s",GetNomeArquivo(),COMTRADE005);
        return 0;
    };
    return 1;
};

double* CComtrade::GetTimeStamp()
{
	return m_timestamp;
}
