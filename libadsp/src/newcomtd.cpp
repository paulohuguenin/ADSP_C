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
#include "newcomtd.h"

//#include <afx.h>


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
#define COMTRADE016 "Tipo do identificador de escalamento dos dados deve ser \"P\" ou \"S\" ou \"p\" ou \"s\"!"

#ifndef max
#define max(a,b)        (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
#define min(a,b)        (((a) < (b)) ? (a) : (b))
#endif

#define ERROFORMATO()   {\
    if (erro==NULL) return 0;\
    sprintf(erro,"%s %s\nLinha:%d",GetNomeArquivo(),COMTRADE007,linha);\
    return 0;}

int CComtrade97::ExtrairConfiguracao(CConfigComtrade* config, char* erro) 
{
    char string[256];
    int ret, linha = 1, auxI;

    // Extrai o nome da subestacao
    if (!ExtrairLexema( string, 256, ',')) ERROFORMATO();
    config->SetNomeSubstacao(string);

	// Extrai o identificador do RDP
    if (!ExtrairLexema( string, 256, ',')) ERROFORMATO();
    config->SetIdentificadorRDP(string);
	
	// Extrai a versao/ano do comtrade
    if (!ExtrairLinha( string, 256)) ERROFORMATO();
    if (string[0]!=0) 
	{
    	auxI = atoi(string);
    	config->SetVersao(atoi(string));
		ret = ExtrairConfiguracao_(config,erro);
    } 
	else ret = CBinaryComtrade::ExtrairConfiguracao_(config,erro);

    return ret;
};


////////////////////////////////////////////////////////////////////////////
int CComtrade97::LerLinhaCanalAnalogico(CConfigComtrade* config, 
										CCanal *canal, 
										int cont, 
										char* erro)
{
	char    string[256];
	int  ret, linha = 1, versao;
	float auxF;
	CCanalAnalogico *analog;
		
	// Extrai os seguintes itens de configuracao do canal analogico:
	// An,ch_id,ph,ccbm,uu,a,b,skew,min,max
	ret = CComtrade::LerLinhaCanalAnalogico(config, canal, cont, erro);

	analog = (CCanalAnalogico*) canal;

	// Extrai o primario (primary)
	if (!ExtrairLexema( string, 256, ',')) ERROFORMATO();
    ret = sscanf( string, "%f", &auxF);
    if ((ret<=0) || (ret==EOF)) auxF = 0.F;
    analog->SetPrimario(auxF);
        
	// Extrai o secundario (secondary)
	versao = config->GetVersao();
	if (versao <= 1997)
	{		
		if (!ExtrairLinha( string, 256)) ERROFORMATO();
		ExtrairEspacosDasBordas(string);
		
		// Pegando algum campo depois da vírgula, caso haja
		char * Dados;
		int result;
		Dados = strchr( string, ',');
		if (Dados !=NULL)	
		{
			result = Dados - string + 1;
			Dados = string + result;
			///////////////// 
			// Campo para indicar se os valores estão escalados para primário
			// ou secundário
			if ((!strcmp( Dados, "P") ) || (!strcmp( Dados, "p") )) 
				analog->SetTipoEscalamento(0);
			else 
			{
				if ((!strcmp( Dados, "S")) || (strcmp( Dados, "s"))) 
					analog->SetTipoEscalamento(1);
				else 
				{ 
					if (erro!=NULL) 
						sprintf(erro,"Linha %d: %s",linha,COMTRADE016); 
					return 0;
				}
			}
		}
	}
	if (versao > 1997)	
	{
		if (!ExtrairLexema( string, 256, ',')) ERROFORMATO();
	}
    ret = sscanf( string, "%f", &auxF);
    if ((ret<=0) || (ret==EOF)) auxF = 0.F;
    analog->SetSecundario(auxF);

	return ret;
};
//////////////////////////////////////////////////////////////////////////////

int CComtrade97::LerLinhaCanalDigital(CConfigComtrade* config, 
									  CCanal *canal, 
									  int cont, 
									  char* erro)
{
    char    string[256];
	int  	auxI, ret, linha = 1;	// cont,
	
	// Extrai o indice do canal (Dn)
	if (!ExtrairLexema( string, 256, ',')) ERROFORMATO();
    ret = sscanf( string, "%d", &auxI);
    if ((ret<=0) || (ret==EOF)) auxI = cont;
    canal->SetNumeroCanal(auxI);
	
	// Extrai o nome do canal (ch_id)
    if (!ExtrairLexema( string, 256, ',')) ERROFORMATO();
    canal->SetNomeCanal(string);

	// Extrai a fase  do canal (ph)
	if (!ExtrairLexema( string, 256, ',')) ERROFORMATO();
	if (strlen(string)>__num_char_fase)  ERROFORMATO();
	canal->SetFaseCanal(string);

	// Extrai o equipamento monitorado (ccbm)
	if (!ExtrairLexema( string, 256, ',')) ERROFORMATO();
	canal->SetEquipamentoCanal(string);

	// Extrai o estado normal do status do canal (y)
    if (!ExtrairLinha( string, 256)) ERROFORMATO();
    if (strlen(string)>1) ERROFORMATO();	
	ExtrairEspacosDasBordas(string);
	canal->SetTipoDoCanal(atoi(&string[0]));

	return ret;
};

//////////////////////////////////////////////////////////////////////////////


int CComtrade97::ExtrairConfiguracao_(CConfigComtrade* config, char* erro) 
{
    char    string[256], *pos, *pos2,
    		HoraPrimPtValid, DiaPrimPtValid, HoraTrigValid, DiaTrigValid;
    int     auxI, cont, ret, linha = 1;
    long    auxL; 
    CCanal *canal;
    float auxF;
    __instante *tTrigger, *tPrima;  
	
    // Valor default pois Simulacao e Comprimido nao fazem parte do COMTRADE
    config->SetSimulacao();
    config->SetComprimido();

    // Numero de canais
    linha++;

	// Extrai a segunda linha do arquivo ".cfg" contendo
	// * Numero total de canais
	// * Numero de canais analógicos
	// * Numero de canais digitais
	auxI = LerLinha2(config, erro);	
	if (auxI==0) return 0;
	

	// Extrai as caracteristicas de cada canal tanto analog quanto digitais
	// * Analog -> An,ch_id,ph,ccbm,uu,a,b,skew,min,max,primary,secondary
	// * Digitais -> Dn,ch_id,ph,ccbm,y
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
			LerLinhaCanalDigital(config, canal, cont, erro);
    };//for cont
    
    linha++;
	// Extrai a frequencia da linha (lf)
    if (!ExtrairLinha( string, 256)) ERROFORMATO();
    ret = sscanf( string, "%f", &auxF);
    if ((ret<=0) || (ret==EOF)) ERROFORMATO();
	if (auxF==0)	auxF=60;
    config->SetFreqLinha(auxF);

    
    linha++;
	// Extrai o numero de taxas de amostragem (nrates)
    if (!ExtrairLinha( string, 256)) ERROFORMATO();
    ret = sscanf( string, "%d", &auxI);
    printf("nrates: ret - %d; auxI - %d\n",ret,auxI);   
    if ((ret<=0) || (ret==EOF) || (ret>__num_max_taxas) ) ERROFORMATO();
    
	float auxF1[100];
	long auxL1 [100];
	for(cont=0; cont<100; cont++)	{
		auxF1[cont] = 0;
		auxL1[cont] = 0;
	}


	// Extrai as taxas de amostragem e o valor da ultima
	// amostra (samp,endsamp) onde auxI = nrates
	if (auxI==0)  auxI = 1;
	
	
	for (cont=1; cont<=auxI; cont++) 
	{
        linha++;

        if (!ExtrairLexema( string, 256, ',')) ERROFORMATO();
        ExtrairEspacosDasBordas(string);
	    ret = sscanf( string, "%f", &auxF);
	    printf("SAMP: ret - %d; auxF - %f\n",ret,auxF); 
		if ((ret<=0) || (ret==EOF)) ERROFORMATO();
        auxF1[cont-1] = auxF;

		if (!ExtrairLinha( string, 256)) ERROFORMATO();
        ExtrairEspacosDasBordas(string);
		ret = sscanf( string, "%ld", &auxL);
		printf("ENDSAMP: ret - %d; auxL - %ld\n",ret,auxL); 
        if ((ret<=0) || (ret==EOF)) ERROFORMATO();
		auxL1[cont-1] = auxL;

		config->AcrescentarTaxa(auxF1[cont-1], auxL1[cont-1]);

	}; //for cont

	
    linha++;	
	// Extrai data da primeira amostra
	DiaPrimPtValid = 1;
    tPrima = config->GetInstantePrimeiraAmostra();
    if (!ExtrairLexema( string, 256, ',')) ERROFORMATO();
    ExtrairEspacosDasBordas(string);

    if (strlen(string)>0) 
	{
    	// Extrai o dia
    	if ((pos = strchr( string, '/'))==NULL) ERROFORMATO();
    	if (!pos) ERROFORMATO();
    	pos[0] = 0;
	    auxI = 0;	
		ret = sscanf( string, "%d", &auxI);
	    if ((ret<=0) || (ret==EOF) || (auxI<1) || (auxI>31)) ERROFORMATO();
	    tPrima->tm_mday = auxI;

    	// Extrai o mes
    	if ((pos2 = strchr( ++pos, '/'))==NULL) ERROFORMATO();
    	if (!pos2) ERROFORMATO();
    	pos2[0] = 0;
	    auxI = 0;	
		ret = sscanf( pos, "%d", &auxI);
	
	    if ((ret<=0) || (ret==EOF) || (auxI<1) || (auxI>12))   
		{
	    	if (auxI>12)	
			{
	    		tPrima->tm_mon = --tPrima->tm_mday;
	    		tPrima->tm_mday = auxI; 
	     	}
			else 
	    		tPrima->tm_mon = --auxI;    
	    	if ((ret<=0) || (ret==EOF) || (auxI<1))	                     
	    		ERROFORMATO();                      
	    }
		else tPrima->tm_mon = --auxI;   
	    
    	// Extrai o ano   
    	auxI = 0;	
		ret = sscanf( ++pos2, "%d", &auxI);  
    	if ((ret<=0) || (ret==EOF) || (auxI<1900) || (auxI>9999))	
		{
    	   	if ((auxI<1900) || (auxI>9999))	
	    		tPrima->tm_year = auxI;
		    if ((ret<=0) || (ret==EOF))
	    		ERROFORMATO();	
		}	 
	    if (auxI>=1900) 
			tPrima->tm_year = auxI-1900;    
	    if (auxI>1999)
	        tPrima->tm_year = auxI-2000;
		   
	   
	} 
	else DiaPrimPtValid = 0;
	  		
	// Extrai a hora da primeira amostra
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
	    auxI = 0;	
		ret = sscanf( pos, "%d", &auxI);
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
    	// Extrai o dia
    	if ((pos = strchr( string, '/'))==NULL) ERROFORMATO();
    	if (!pos) ERROFORMATO();
    	pos[0] = 0;
	    auxI = 0;	
		ret = sscanf( string, "%d", &auxI);
	    if ((ret<=0) || (ret==EOF) || (auxI<1) || (auxI>31)) ERROFORMATO();
	    tTrigger->tm_mday = auxI;

    	// Extrai o mes
    	if ((pos2 = strchr( ++pos, '/'))==NULL) ERROFORMATO();
    	if (!pos2) ERROFORMATO();
    	pos2[0] = 0;
	    auxI = 0;	
		ret = sscanf( pos, "%d", &auxI);
    	
	    if ((ret<=0) || (ret==EOF) || (auxI<1) || (auxI>12))   
		{
	    	if (auxI>12)	
			{
	    		tTrigger->tm_mon = --tTrigger->tm_mday;
	    		tTrigger->tm_mday = auxI; 
	     	}
			else 
	    		tTrigger->tm_mon = --auxI;    
	    	if ((ret<=0) || (ret==EOF) || (auxI<1))	                     
	    		ERROFORMATO();                      
	    }
		else tTrigger->tm_mon = --auxI;
	    	   
    	// Extrai o ano
	    auxI = 0;	
		ret = sscanf( ++pos2, "%d", &auxI);
	    if ((ret<=0) || (ret==EOF) || (auxI<1900) || (auxI>9999))    
		{
	    	if ((auxI<1900) || (auxI>9999))	
		 		tTrigger->tm_year = auxI;
		    if ((ret<=0) || (ret==EOF))
	    		ERROFORMATO();	
	    }
	        
	    if (auxI>=1900)
			tTrigger->tm_year = auxI-1900;    
	    if (auxI>1999)
	        tTrigger->tm_year = auxI-2000; 
	  
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
	// Extrai o tipo de arquivo
    if (!ExtrairLinha( string, 256)) ERROFORMATO();
    ExtrairEspacosDasBordas(string);	
    if (!strcmp( string, "ASCII")) config->SetTipoArquivo(0);
    else if (!strcmp( string, "BINARY")) config->SetTipoArquivo(1);
    else 
	{ 
		if (erro!=NULL) sprintf(erro,"Linha %d: %s",linha,COMTRADE011); 
		return 0;
	};

	linha++;
	// Extrai o fator de  multiplicacao
	if (!ExtrairLinha( string, 256)) strcpy(string, "1");
    ret = sscanf( string, "%f", &auxF);
    if ((ret<=0) || (ret==EOF) || (auxF<0)) ERROFORMATO();
	config->SetFator(auxF);

    return 1;
};

int CComtrade97::InserirConfiguracao_(CConfigComtrade* config, char* erro)
{
    char string[256], auxS[256];
    int pos, cont, pos2;
    float ganho, OffSet;
    CCanal *canal;
    CCanalAnalogico *analog;
    const CTaxas *taxa;
    __instante* instante;

	// Insere: 
	// * nome da subestacao
	// * identificador do RDP
	// * versao 97
    strcpy( string, config->GetNomeSubstacao());
    ExtrairVirgula( string);
    strcat( string, ",");
    pos = strlen(string);
    strcat( string, config->GetIdentificadorRDP());
    ExtrairVirgula( &string[pos]);
    strcat( string, ",");
	strcpy(auxS, "1997");
    strcat( string, auxS);
    if (!InserirLinha(string)) return 0;
	
	// Insere:
	// * numero total de canais
	// * numero de canais analogicos
	// * numero de canais digitais
    sprintf( string, "%d,%dA,%dD",
            config->GetNumCanaisAnalogicos()+config->GetNumCanaisDigitais(),
            config->GetNumCanaisAnalogicos(),
			config->GetNumCanaisDigitais());
    if (!InserirLinha(string)) return 0;

	// Insere as caracteristicas dos canais analogicos e digitais
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

		// Insere a fase do canal (ph) e o equipamento monitorado (ccbm)
		strcpy(string,canal->GetFaseCanal());
		pos2 = pos = strlen(string);
        pos+=sprintf(string+pos,",%s",canal->GetEquipamentoCanal());
        ExtrairVirgula(&string[++pos2]);
        if (!InserirString(string)) return 0;

        if (canal->GetAnalogico()) 
		{
            analog = (CCanalAnalogico*) canal;

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
					pos+=sprintf(string+pos,"%G,",analog->GetCoefAng(0));
					pos+=sprintf(string+pos,"%G,",analog->GetCoefLin(0));
					break;
				case e_real:// sinal sera' quantizado
					ganho = COMTRADELIMITEDOQUANTIZADOR / 
						(analog->GetValorMaximoCfg()-analog->GetValorMinimoCfg());

					OffSet = (analog->GetValorMaximoCfg()+analog->GetValorMinimoCfg())/2;
					pos+=sprintf(string+pos,"%G,",analog->GetCoefAng(0)/ganho);
					pos+=sprintf(string+pos,"%G,", 
						analog->GetCoefLin(0) + analog->GetCoefAng(0)*OffSet);
                break;
			};

			// Insere o time skew (skew)
            pos+=sprintf(string+pos,"%G,",analog->GetTimeSkew()*1E06);

			// Insere os valores minimo  e maximo (min,max)
            switch (analog->GetFormatoDados()) 
			{
				case e_inteiro:
					pos+=sprintf(string+pos,"%.0f,",analog->GetValorMinimoCfg());
					pos+=sprintf(string+pos,"%.0f",analog->GetValorMaximoCfg());
					break;
				case e_real:// sinal sera' quantizado
					pos+=sprintf(string+pos,"%ld,%ld",-COMTRADELIMITEDOQUANTIZADOR/2, 
						COMTRADELIMITEDOQUANTIZADOR/2);
                break;
            };

			// Insere primario e secundario (primary,secondary)
            pos+=sprintf(string+pos,",%G,",analog->GetPrimario());
            pos+=sprintf(string+pos,"%G",analog->GetSecundario()); 
        } 
		else  // canal digital
		{
			// Insere estado normal (y)
            sprintf(string,",%d",canal->GetTipoDoCanal());
		}
        if (!InserirLinha(string)) return 0;

	};//for cont
	
	// Insere a frequencia da linha (lf)
    sprintf(string,"%d",(int)config->GetFreqLinha());
    if (!InserirLinha(string)) return 0;    

	// Insere o numero de taxas de amostragem (nrates)
    sprintf(string,"%d",config->GetNumeroTaxas());
    if (!InserirLinha(string)) return 0;    

	// Insere as taxas e a ultima amostra correspondente (samp,endsamp)
    for (cont=1; cont<=config->GetNumeroTaxas(); cont++) 
	{
        taxa = config->GetTaxa(cont);
        pos =sprintf(string,"%G,", taxa->GetTaxaAmostragem());
        pos+=sprintf(string+pos,"%li",taxa->GetUltimaAmostra());
        if (!InserirLinha(string)) return 0;
    }; //for cont

	// Insere data e hora da primeira amostra
    instante = config->GetInstantePrimeiraAmostra();
    pos =sprintf(string,"%02d/",instante->tm_mday);
    pos+=sprintf(string+pos,"%02d/",instante->tm_mon+1);
    
    if ((instante->tm_year >= 0) && (instante->tm_year < 70))
    	pos+=sprintf(string+pos,"%04d,",instante->tm_year+2000); 
    if (instante->tm_year >= 1999)
    	pos+=sprintf(string+pos,"%04d,",instante->tm_year);
    if (instante->tm_year > 69 && instante->tm_year <= 100)
    	pos+=sprintf(string+pos,"%04d,",instante->tm_year+1900);		
    
    pos+=sprintf(string+pos,"%02d:",instante->tm_hour);
    pos+=sprintf(string+pos,"%02d:",instante->tm_min);
    pos+=sprintf(string+pos,"%02.6f",instante->tm_sec);
    if (!InserirLinha(string)) return 0;

	// Insere data e hora do Trigger
    instante = config->GetInstanteTrigger();
    pos =sprintf(string,"%02d/",instante->tm_mday);
    pos+=sprintf(string+pos,"%02d/",instante->tm_mon+1);
        
    if ((instante->tm_year >= 0) && (instante->tm_year < 70))
    	pos+=sprintf(string+pos,"%04d,",instante->tm_year+2000); 
    if (instante->tm_year >= 1999)
    	pos+=sprintf(string+pos,"%04d,",instante->tm_year);
    if (instante->tm_year > 69 && instante->tm_year <= 100)
    	pos+=sprintf(string+pos,"%04d,",instante->tm_year+1900);		
    
    pos+=sprintf(string+pos,"%02d:",instante->tm_hour);
    pos+=sprintf(string+pos,"%02d:",instante->tm_min);
    pos+=sprintf(string+pos,"%02.6f",instante->tm_sec);
    if (!InserirLinha(string)) return 0;

	// Insere o tipo do arquivo
	if (config->GetTipoArquivo()) 
	{
		if (!InserirLinha("BINARY")) return 0;
	} 
	else 
	{
		if (!InserirLinha("ASCII")) return 0;
	};

	// Insere o fator de multiplicacao
    sprintf(string,"%G", config->GetFator());
    if (!InserirLinha(string)) return 0;

    return 1;
};

int CComtrade97::InserirConfiguracao(CConfigComtrade* config, char* erro)
{
    int ret;

	ret = CComtrade97::InserirConfiguracao_( config, erro);

    return ret;
};

