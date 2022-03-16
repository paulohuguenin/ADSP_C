/*****************************************************************************

   Projeto : SINAPE
   Nome    : BUFTIPOS.CPP
   Funcao  : Arquivo com implementacao de objetos contendo formatos internos
             de representacao dos dados dos RDPs, baseados no COMTRADE
   C       : C++
   Autor   : MAMR
   Data    : 14/11/94
*****************************************************************************/
//#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "buftipos.h"

// New
#define CNewComtrade CComtrade99;

#ifndef ASSERT
	#define ASSERT		assert
#endif

#ifndef max
#define max(a,b)        (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
#define min(a,b)        (((a) < (b)) ? (a) : (b))
#endif

void StringToUpper (char* &string)
{	
	char* strAux;
	int len = strlen(string);
	
	strAux = new char[len+1];
	
	for(int i=0;i<len;i++)	
	{
		strAux[i]= toupper(string[i]);
	}
	
	strAux[len] = '\0';
   	strcpy(string,strAux);
	
	delete [] strAux;	
}

void StringReversal(char* &string)
{
	char *strAux;
	int len = strlen(string);
	
	strAux = new char[len+1];
	
	for (int i=0;i<len;i++)
	{
		strAux[i]=string[len-1-i];
	}
	
	strAux[len] = '\0';
   	strcpy(string,strAux);
	
	delete [] strAux;		
}

void ExtrairEspacosDasBordas( char* string)
{
    int i = strspn( string, " \t");
    strcpy( string, &string[i]);
	StringReversal(string);
//    strrev(string);
    i = strspn( string, " \t");
    strcpy( string, &string[i]);
	StringReversal(string);
//    strrev(string);
};


/*******************             class CCanal            ********************/

CCanal::CCanal(){
    nome = (char*) calloc( 1, sizeof(char)); 
	*nome = 0;
    numero = 1;
    estado_normal = 2;
	fase = (char*) calloc( __num_char_fase+1, sizeof(char));
	*fase = 0X0020;		
	*(fase+1) = 0;
	equipamento_monitorado = (char*) calloc( 1, sizeof(char));
	*equipamento_monitorado = 0;
};

CCanal::CCanal(const CCanal& canal)
{
	nome = NULL;
	equipamento_monitorado = NULL;
	Copiar( &canal);
};

CCanal::CCanal(const CCanal* canal)
{
	nome = NULL;
	equipamento_monitorado = NULL;
	Copiar( canal);
};

void CCanal::Copiar(const CCanal* canal)
{
	SetNumeroCanal(canal->GetNumeroCanal());
	SetNomeCanal(canal->GetNomeCanal());
	SetTipoDoCanal(canal->GetTipoDoCanal());
	SetFaseCanal(canal->GetFaseCanal());//Orig de CCanalAnalogico
	SetEquipamentoCanal(canal->GetEquipamentoCanal());//Orig de CCanalAnalogico
};

CCanal::~CCanal()
{
    if (nome!=NULL) free( nome);
	if (equipamento_monitorado!=NULL) free(equipamento_monitorado);
	if (fase!=NULL) free(fase);
	nome=NULL;
	equipamento_monitorado=NULL;
	fase=NULL;
};

unsigned char   CCanal::GetAnalogico() const {return estado_normal>1;};
int     CCanal::GetNumeroCanal() const {   return numero;};
char*   CCanal::GetNomeCanal() const { return nome;};
unsigned char CCanal::GetTipoDoCanal() const { return estado_normal;};

int     CCanal::SetNumeroCanal(int numero_) {
    if (numero_>=__num_max_canais) return 0;
    else numero = numero_;
    return 1;
};

int CCanal::SetNomeCanal(const char* nome_) {
    int tam;
    char *auxC;
    tam = strlen(nome_);
    auxC = (char*) calloc( tam+1, sizeof(char));
    if (!auxC) return 0;
    strcpy( auxC, nome_);
    ExtrairEspacosDasBordas( auxC);
    tam = strlen(auxC);
    if (nome!=NULL) free( nome);
	nome = (char*) calloc( tam+1, sizeof(char));
    if (!nome) {free(auxC); return 0;};
    strcpy( nome, auxC);
	free(auxC);
    return 1;
  };

int CCanal::SetTipoDoCanal(unsigned char tipo) {
    if (tipo<= 6) 	estado_normal = tipo;
    else return(0);
    return 1;
};

char*   CCanal::GetFaseCanal() const {return fase;};
char*   CCanal::GetEquipamentoCanal() const {return equipamento_monitorado;};

int     CCanal::SetFaseCanal(const char* fase_) {
    int c=1;
    char fase__ = (char) toupper(fase_[0]);
    if (strlen(fase_)==1) 
	{
	    switch (fase__) {
	     case 'A': strcpy(fase,"A"); break;
	     case 'B': strcpy(fase,"B"); break;
	     case 'C': strcpy(fase,"C"); break;
	     case 'N': strcpy(fase,"N"); break;
	     case 'R': strcpy(fase,"R"); break;
	     case 'Y': strcpy(fase,"Y"); break;
	     case ' ': strcpy(fase," "); break;
	     default:
	                strncpy(fase,fase_,2);
					fase[2]=0;
	                c=0; 
					break;
	    };
	} 
	else 
	{ 
	  	strncpy( fase, fase_, min(__num_char_fase,strlen(fase_)) ); 
	  	fase[min(__num_char_fase,strlen(fase_))]=0;
	};
    return c;
};

int CCanal::SetEquipamentoCanal(const char* equipamento) {
    int tam;
    char *auxC;
    tam = strlen(equipamento);
    auxC = (char*) calloc( tam+1, sizeof(char));
    if (!auxC) return 0;
    strcpy( auxC, equipamento);
    ExtrairEspacosDasBordas(auxC);
    tam = strlen(auxC);
    if (equipamento_monitorado!=NULL) free(equipamento_monitorado);
    equipamento_monitorado = (char*) calloc( tam+1, sizeof(char));
    if (!equipamento_monitorado) {free(auxC);return 0;};
    strcpy( equipamento_monitorado, auxC);
	free(auxC);
    return 1;
  };

/*******************         class CCanalAnalogico       ********************/

CCanalAnalogico::CCanalAnalogico(){
    formato = e_real;
    unidade = (char*) calloc( 1, sizeof(char)); *unidade = 0;
    coef_ang = 1.F;
    coef_lin = 0.F;
    time_skew = 0.F;
    valor_minimo_range = -1.F;
    valor_maximo_range = 1.F;
	primario = 1.F;
	secundario = 1.F;
	tipo_escala = 0;
};

CCanalAnalogico::CCanalAnalogico(const CCanalAnalogico& canal)
{
    unidade = NULL;
	CCanalAnalogico::Copiar( &canal);
};

CCanalAnalogico::CCanalAnalogico(const CCanalAnalogico* canal)
{
	unidade = NULL;
	CCanalAnalogico::Copiar( canal);
};

void CCanalAnalogico::Copiar(const CCanalAnalogico* canal)
{
	sret = NULL;
	CCanal::Copiar( canal);

	SetFormatoDados(canal->GetFormatoDados());

	// Alocando memoria
	if (sret==NULL) sret = (char *) calloc(256,sizeof(char));

	double fator_mult = CorrigeUnidade(canal->GetUnidadeCanal(), sret);

	SetUnidadeCanal(sret);
	SetCoefAng((float)(canal->GetCoefAng(canal->GetTipoEscalamento())*fator_mult));
	SetCoefLin((float)(canal->GetCoefLin(canal->GetTipoEscalamento())*fator_mult));
	SetTimeSkew(canal->GetTimeSkew());
	SetValorMinimo(canal->GetValorMinimoCfg());
	SetValorMaximo(canal->GetValorMaximoCfg());
	SetPrimario(canal->GetPrimario());
	SetSecundario(canal->GetSecundario());
	SetTipoEscalamento( canal->GetTipoEscalamento());

	free(sret);
};

CCanalAnalogico::~CCanalAnalogico()
{
    if (unidade!=NULL) free(unidade);
	unidade=NULL;
};

int     CCanalAnalogico::GetFormatoDados() const {return formato;};
char*   CCanalAnalogico::GetUnidadeCanal() const { return unidade;};

float   CCanalAnalogico::TransformaPrimSec( float valor, int PrimSec, int valorPrimitivo) const
{ 

	if (PrimSec!=valorPrimitivo)	
	{
		float primario = GetPrimario();
		float secundario = GetSecundario();

		if (PrimSec!=0)	
		{
			if ((secundario==0.0F) && (primario==0.0F))
				valor = 0;
			else
				valor = valor*(secundario/primario);
		}
		else	
		{
			if ((secundario==0.0F) && (primario==0.0F))
				valor = 0;
			else
				valor = valor*(primario/secundario);
		}
	}

	return valor;
};

float   CCanalAnalogico::TransformaPrimSec( float valor, int PrimSec) const
{ 
	return TransformaPrimSec( valor, PrimSec, tipo_escala);
};

float   CCanalAnalogico::GetCoefAng(int PrimSec) const
{ 
	float valor;
	valor = TransformaPrimSec(coef_ang, PrimSec);
	return valor;

};	

float   CCanalAnalogico::GetCoefLin(int PrimSec) const
{ 
	return TransformaPrimSec(coef_lin, PrimSec);
};	

float   CCanalAnalogico::GetValorMinimo(int PrimSec) const
{
	return GetValorMinimoCfg()*GetCoefAng(PrimSec) +
										GetCoefLin(PrimSec);
};
float   CCanalAnalogico::GetValorMaximo(int PrimSec) const
{
	return GetValorMaximoCfg()*GetCoefAng(PrimSec) +
										GetCoefLin(PrimSec);
};

float   CCanalAnalogico::GetTimeSkew() const { return time_skew;};
float   CCanalAnalogico::GetValorMinimoCfg() const { return valor_minimo_range;};
float   CCanalAnalogico::GetValorMaximoCfg() const { return valor_maximo_range;};
float   CCanalAnalogico::GetPrimario() const { return primario;};
float   CCanalAnalogico::GetSecundario() const { return secundario;};
char   CCanalAnalogico::GetTipoEscalamento() const { return tipo_escala;};

int     CCanalAnalogico::SetFormatoDados(int formato_) {
    if (formato_>e_dupla) return 0;
                    else formato=formato_;
    return 1;
};

double CCanalAnalogico::CorrigeUnidade(const char* unidade_, char* unidade_corrigida)
// Retorna o multiplicador (0.01, 1, 1000, 1000000)
// Retorna em unidade_corrigida a unidade pura, sem o fator m, k ou M
// fator = CorrigeUnidade("kV",unidade_corrigida);
// fator = 1000; unidade_corrigida = "V";
// Precisa prever os casos de mv, mV, kv, kV, Kv, KV, Mv, MV, Volts, Vs ...
// Precisa prever os casos de ma, mA, ka, kA, Ka, KA, Ma, MA, Amperes, Amps ...
// Precisa prever os casos de mva, mvA, mVa, mVA, k ..., K..., M ...
// Precisa prever os casos de mw, mW, kw, kW, Kw, KW, M ..., Watts
{
	
	const int NUMOPTS = 29;
	const int NUMOPS = 6;

	const char possibilit[NUMOPTS][16] = {{"mv"},{"mV"},{"Mv"},{"MV"},{"kv"},
		{"kV"},{"Kv"},{"KV"},{"Volts"},{"Vs"},{"ma"},{"mA"},{"Ma"},{"MA"},{"ka"},
		{"kA"},{"Ka"},{"KA"},{"Ampere"},{"Amps"},{"mw"},{"mW"},{"Mw"},{"MW"},
		{"kw"},{"kW"},{"Kw"},{"KW"},{"Watts"}};

	const char possibilit_[NUMOPS] = {'V','v','A','a','W','w'};
	
	const char corrigido[NUMOPTS] = {'V','V','V','V','V',
		'V','V','V','V','V','A','A','A','A','A',
		'A','A','A','A','A','W','W','W','W',
		'W','W','W','W','W'};

	const double fator[NUMOPTS] = {0.001,0.001,1000000.0,1000000.0,
		1000.0,1000.0,1000.0,1000.0,1.0,1.0,0.001,0.001,
		1000000.0,1000000.0,1000.0,1000.0,1000.0,1000.0,1.0,1.0,
		0.001,0.001,1000000.0,1000000.0,1000.0,1000.0,1000.0,
		1000.0,1.0};
	
	int i;
	char* sAux;
	
	sAux = new char[4];
		
	strcpy(sAux, unidade_);			//char sAux = unidade_;
	StringToUpper(sAux);
	//_strupr(sAux);					// maiúscula
	ExtrairEspacosDasBordas(sAux);	//sAux.LimparBordas();

	char *result = strstr( sAux, "VA");		// procura por VA

	if (result==sAux) 
	{
		strcpy(unidade_corrigida,"VA");
		return 1.0;
	} 
	else 
	{
		if (result==NULL) 
		{
			//não achou VA, vejo outras possibilidades
			strcpy(sAux, unidade_);
			ExtrairEspacosDasBordas(sAux);		//sAux.LimparBordas();
			i=0; unsigned char achou = 0;
			//int nCount = strlen(sAux);
			while ((i<NUMOPTS) && (!achou))	
			{
				if (strncmp(sAux,possibilit[i],2)==0) 
				{
					sAux[0] = corrigido[i];
					sAux[1] = 0;
					achou = 1;
				}
				else	
				{
					if((i<NUMOPS)&&(sAux[0]==possibilit_[i])) 
					{
						sAux[0] = (char) possibilit_[i];	
						sAux[1] = 0;
						StringToUpper(sAux);
						//_strupr(sAux);
						sAux[1] = 0;
						achou = 1;
						i=28;
					}
				}
				i++;
			}
			
			strcpy(unidade_corrigida,sAux);
			return fator[i-1];
		} 
		else 
		{
			// encontrou VA em qq posição > 0
			if (result-sAux > sizeof(char))	// não está na primeira posição ? não sei que unidade é essa	
			{
				strcpy(unidade_corrigida,unidade_);
				return 1.0;
			}
			else // está na primeira posição
			{  
				strcpy(unidade_corrigida,"VA");
				switch (unidade_[0])	{
					case 'k': 
					case 'K': 
						i = 7;
						break;
					case 'm':
						i = 3;
						break;
					case 'M':
						i = 5;
						break;
					default:
						i = 1;
						strcpy(unidade_corrigida,unidade_);
						break;
				}
				return fator[i-1];
			};
		};
	};
	delete [] sAux;
}

int CCanalAnalogico::SetUnidadeCanal(const char* unidade_) 
{
    int c=0, tam;
    char *auxC;

    tam = strlen(unidade_);
    auxC = (char*) calloc( tam+1, sizeof(char));

    if (!auxC) return 0;
    strcpy( auxC, unidade_);
    ExtrairEspacosDasBordas( auxC);
    tam = strlen(auxC);
    if (unidade!=NULL) free(unidade);
    unidade = (char*) calloc( tam+1, sizeof(char));
    if (!unidade) {free(auxC);return 0;};
	strcpy(unidade,unidade_);

	free(auxC);
	return c;
};

int     CCanalAnalogico::SetCoefAng(float coeficiente) {
    coef_ang = coeficiente;
    return 1;
};

int     CCanalAnalogico::SetCoefLin(float coeficiente) {
    coef_lin = coeficiente;
    return 1;
};

int     CCanalAnalogico::SetTimeSkew(float skew) {
    time_skew = skew;
    return 1;
};

int     CCanalAnalogico::SetValorMinimo(float valor) {
    valor_minimo_range = valor;
    if (valor_minimo_range>=valor_maximo_range) return 0;
        else return 1;
};

int     CCanalAnalogico::SetValorMaximo(float valor) {
    valor_maximo_range = valor;
    if (valor_minimo_range>=valor_maximo_range) return 0;
        else return 1;
};

int     CCanalAnalogico::SetPrimario(float valor) {
    if (valor!=0) {
    	primario = valor;
    	return 1;
    } else  return 0;
};

int     CCanalAnalogico::SetSecundario(float valor) {
    if (valor!=0) {
    	secundario = valor;
    	return 1;
    } else  return 0;
};

int     CCanalAnalogico::SetTipoEscalamento( char escala) {
	tipo_escala = escala;
    return 1;
};

int CCanalAnalogico::AcertarTipoDoCanalAnalogico() {
    char *pC;
    //acertar o tipo usando a unidade
    if (GetTipoDoCanal()!=6) return 0;
    pC = strchr( unidade, 'V');
    if (pC) { SetTipoDoCanal(2); return 1;}
    pC = strchr( unidade, 'v');
    if (pC) { SetTipoDoCanal(2); return 1;}
    pC = strchr( unidade, 'A');
    if (pC) { SetTipoDoCanal(3); return 1;}
    pC = strchr( unidade, 'a');
    if (pC) { SetTipoDoCanal(3); return 1;}
    return 0;
};

/*******************             class CTaxas            ********************/

CTaxas::CTaxas() {
    taxa_amostragem = 1.F;
    ultima_amostra = 0;
};

float CTaxas::GetTaxaAmostragem() const {return taxa_amostragem;};

long CTaxas::GetUltimaAmostra() const {return ultima_amostra;};

int CTaxas::SetTaxaAmostragem(float taxa) {
    if (taxa>0) taxa_amostragem = taxa;
    else    return 0;
    return 1;
};

int CTaxas::SetUltimaAmostra(long indice) {
    if (indice>=0) ultima_amostra = indice;
    else    return 0;
    return 1;
};

/*******************        class CConfigComtrade        ********************/

CConfigComtrade::CConfigComtrade() {
    int cont;
    time_t long_time;
    tm *pH;
    simulacao = 0;
    comprimido = 0;	
    //nome_subestacao = (char*) calloc( 1, sizeof(char));
	//*nome_subestacao = 0;
	nome_subestacao = NULL;
    //identificador_RDP = (char*) calloc( 1, sizeof(char));
	//*identificador_RDP = 0;
	identificador_RDP = NULL;
    num_canais_analogicos = 0;
    num_canais_digitais = 0;
    for (cont=0; cont<__num_max_canais; cont++) canal[cont] = NULL;
    freq_linha = 60.F;
    numero_taxas = 0;
    time(&long_time);
    pH = localtime(&long_time);

    if (pH!=NULL) 
	{
	    primeira_amostra.tm_sec = (float)pH->tm_sec;
	    primeira_amostra.tm_mday=pH->tm_mday;
	    primeira_amostra.tm_mon= pH->tm_mon;
	    primeira_amostra.tm_year= pH->tm_year;
	    primeira_amostra.tm_hour= pH->tm_hour;
	    primeira_amostra.tm_min= pH->tm_min;
	    primeira_amostra.tm_sec= (float)pH->tm_sec;
	    if (pH->tm_year>99)
	    	primeira_amostra.tm_year= pH->tm_year+1900;
	};

    memcpy( &trigger, &primeira_amostra, sizeof(__instante));
    offset_inicio_dados = 0;
    tipo_arquivo = 2;
	versao = 1991;
	fator = 1;

};

CConfigComtrade::CConfigComtrade(const CConfigComtrade& config_orig)
{
	int i, can;
	CCanal *pC;
	const CTaxas *taxa;
	__instante inst;

    nome_subestacao = NULL;    
	identificador_RDP = NULL;

    SetSimulacao( config_orig.GetSimulacao());
    SetComprimido( config_orig.GetComprimido());
    SetNomeSubstacao( config_orig.GetNomeSubstacao());
    SetIdentificadorRDP( config_orig.GetIdentificadorRDP());
    SetFreqLinha( config_orig.GetFreqLinha());
	inst=config_orig.GetInstantePrimeiraAmostra();
	SetInstantePrimeiraAmostra( &inst);
	inst = config_orig.GetInstanteTrigger();	
	SetInstanteTrigger( &inst);
	SetOffsetInicioDados(config_orig.GetOffsetInicioDados());
	SetTipoArquivo(config_orig.GetTipoArquivo());
	SetVersao(config_orig.GetVersao());
	SetFator(config_orig.GetFator());
	
	// acerta canais
	num_canais_analogicos = num_canais_digitais = 0;
    for (i=0; i<__num_max_canais; i++) canal[i] = NULL;
	for (i=0; i<config_orig.GetNumCanaisAnalogicos(); i++) 
	{
		pC = config_orig.GetCanalPtr(i+1);
		ASSERT(pC->GetTipoDoCanal()>1);
		can = AcrescentarCanalAnalogico();
		if (can>0) ((CCanalAnalogico*)GetCanalPtr(can))->Copiar((CCanalAnalogico*)pC);
	};
	for (i=config_orig.GetNumCanaisAnalogicos();
		i<config_orig.GetNumCanaisAnalogicos()+config_orig.GetNumCanaisDigitais(); 
		i++) 
	{
		pC = config_orig.GetCanalPtr(i+1);
		ASSERT(pC->GetTipoDoCanal()<2);
		can = AcrescentarCanalDigital();
		if (can>0) GetCanalPtr(can)->Copiar(pC);
	};
	// acerta taxas
	numero_taxas = 0;
	for (i=0;i<config_orig.GetNumeroTaxas(); i++) 
	{
		taxa = config_orig.GetTaxa(i+1); 
		ASSERT(taxa);
		AcrescentarTaxa( taxa->GetTaxaAmostragem(), taxa->GetUltimaAmostra());
	};
};

CConfigComtrade::~CConfigComtrade() 
{
    int cont;

    if (nome_subestacao!=NULL) 
	{
    	free(nome_subestacao);
    	nome_subestacao=NULL;
   	};
    if (identificador_RDP!=NULL) 
	{
    	free(identificador_RDP);
    	identificador_RDP=NULL;
    };
    if (num_canais_analogicos<0) 
		num_canais_analogicos = 0;
    if (num_canais_digitais>__num_max_canais) 
		num_canais_digitais = 0;
    if (num_canais_analogicos+num_canais_digitais>__num_max_canais)
        num_canais_analogicos = num_canais_digitais = 0;

    for (cont=0; cont<num_canais_analogicos; cont++)
	{
		if (canal[cont]!=NULL)
		{ 
			delete((CCanalAnalogico*)canal[cont]); 
			canal[cont] = NULL;
		}
	}
    for (cont=num_canais_analogicos;
		 cont<num_canais_digitais+num_canais_analogicos;
         cont++)
	{
        if (canal[cont]!=NULL)
		{ 
			delete((CCanal*)canal[cont]); 
			canal[cont] = NULL;
		};
	}

    num_canais_analogicos = 0;
    num_canais_digitais = 0;
    numero_taxas = 0;
    offset_inicio_dados = 0;
};


char    CConfigComtrade::GetSimulacao() const 
			{return simulacao;};
char    CConfigComtrade::GetComprimido() const 
			{return comprimido;};
char*   CConfigComtrade::GetNomeSubstacao() const 
			{return nome_subestacao;};
char*   CConfigComtrade::GetIdentificadorRDP() const
			{return identificador_RDP;};
int     CConfigComtrade::GetNumCanaisAnalogicos() const
			{return num_canais_analogicos;};
int     CConfigComtrade::GetNumCanaisDigitais() const
			{return num_canais_digitais;};

CCanal* CConfigComtrade::GetCanalPtr(int num_canal) const {
    if ( (num_canal>__num_max_canais) || 
		(num_canal<1) ||
        (num_canal>num_canais_analogicos+num_canais_digitais) )
		return NULL;
    else 
		return canal[num_canal-1];
};

float  CConfigComtrade::GetFreqLinha() const 
			{return freq_linha;};
char   CConfigComtrade::GetNumeroTaxas() const 
			{return numero_taxas;};

int CConfigComtrade::SetNumeroTaxas(const char num) {
	numero_taxas = num; 
	return 1;
};

const CTaxas* CConfigComtrade::GetTaxa(int indice) const {
    if ( (indice > __num_max_taxas) || 
		(indice > numero_taxas) || 
		(indice<1) )
        return NULL;
    else  
		return &taxas[indice-1];
};

long CConfigComtrade::GetNumeroAmostrasPorCanal() const {
	if (numero_taxas<1) 
		return 0L;
	else 
		return taxas[numero_taxas-1].GetUltimaAmostra();
};

__instante* CConfigComtrade::GetInstantePrimeiraAmostra()
			{return &primeira_amostra;};
__instante* CConfigComtrade::GetInstanteTrigger() 
			{return &trigger;};
__instante CConfigComtrade::GetInstantePrimeiraAmostra() const
			{return primeira_amostra;};
__instante CConfigComtrade::GetInstanteTrigger() const 
			{return trigger;};
long CConfigComtrade::GetOffsetInicioDados() const
			{return offset_inicio_dados;};
char   CConfigComtrade::GetTipoArquivo() const 
			{return tipo_arquivo;};
int	CConfigComtrade::GetVersao() const 
			{return versao;};
float CConfigComtrade::GetFator() const 
			{return fator;};

int CConfigComtrade::SetSimulacao(char simulacao_)
			{simulacao=simulacao_; return 1;};
int CConfigComtrade::SetComprimido( char comprimido_)
			{comprimido=comprimido_; return 1;};

int CConfigComtrade::SetNomeSubstacao( const char * nome){
    int tam;
    char *auxC;

    tam = strlen(nome);
    auxC = (char*) calloc( tam+1, sizeof(char));
    if (!auxC) return 0;
    strcpy( auxC, nome);
    ExtrairEspacosDasBordas( auxC);
    tam = strlen(auxC);
    if (nome_subestacao!=NULL) free(nome_subestacao);
    nome_subestacao = (char*) calloc( tam+1, sizeof(char));
    if (!nome_subestacao) {free(auxC);return 0;};
    strcpy( nome_subestacao, auxC);

    free(auxC);
    return 1;
  };

int CConfigComtrade::SetIdentificadorRDP( const char *identificador_){
    int tam;
    char *auxC;

    tam = strlen(identificador_);
    auxC = (char*) calloc( tam+1, sizeof(char));
    if (!auxC) return 0;
    strcpy( auxC, identificador_);
    ExtrairEspacosDasBordas( auxC);
    tam = strlen(auxC);
    if (identificador_RDP!=NULL) free(identificador_RDP);
    identificador_RDP = (char*) calloc( tam+1, sizeof(char));
    if (!identificador_RDP) {free(auxC);return 0;};
    strcpy( identificador_RDP, auxC);

    free(auxC);
    return 1;
  };

int CConfigComtrade::AcrescentarCanalAnalogico(int* num) {
	
	ASSERT((num_canais_analogicos>=0)&&(num_canais_analogicos>=0));
	int i, canais=num_canais_analogicos+num_canais_digitais;

	if (canais<__num_max_canais) 
	{
		for (i=canais; i>num_canais_analogicos; i--)
			canal[i] = canal[i-1];
		canal[num_canais_analogicos] = new CCanalAnalogico;
		if (canal[num_canais_analogicos]==NULL) return 0;
		canal[num_canais_analogicos]->SetTipoDoCanal(6);
		num_canais_analogicos++;
	} 
	else return 0;
	if (num!=NULL) *num = num_canais_analogicos;

	return num_canais_analogicos+num_canais_digitais;
};

int CConfigComtrade::AcrescentarCanalDigital(int* num) {

	ASSERT((num_canais_analogicos>=0)&&(num_canais_analogicos>=0));
	int canais=num_canais_analogicos+num_canais_digitais;
	if (canais<__num_max_canais) 
	{	
        canal[canais] = new CCanal;
        if (canal[canais]==NULL) return 0;
        canal[canais]->SetTipoDoCanal(0);
        num_canais_digitais++;
  	} 
	else return 0;	
	if (num!=NULL) *num = num_canais_analogicos+num_canais_digitais;

	return num_canais_analogicos+num_canais_digitais;
};

int CConfigComtrade::ApagarCanal(int num_canal) {

    int i, can = num_canais_analogicos+num_canais_digitais;

	if ((num_canal<1)||(num_canal>can)) return 0;
	if (canal[num_canal-1]==NULL) return 0;
	if (num_canal==can) 
	{
		if (canal[num_canal-1]->GetTipoDoCanal()<2) 
		{
			num_canais_digitais--;
			delete((CCanal*)canal[num_canal-1]);
			canal[num_canal-1] = NULL;
		} 
		else 
		{
			num_canais_analogicos--;
			delete((CCanalAnalogico*)canal[num_canal-1]);
			canal[num_canal-1] = NULL;
		};
		return 1;
	};
	if (canal[num_canal-1]->GetTipoDoCanal()<2) 
	{
		num_canais_digitais--;
		delete((CCanal*)canal[num_canal-1]);
		canal[num_canal-1] = NULL;
	} 
	else 
	{
		num_canais_analogicos--;
		delete((CCanalAnalogico*)canal[num_canal-1]);
		canal[num_canal-1] = NULL;
	};
	for (i=num_canal;i<can;i++) 
	{
		if (canal[i]==NULL) return 0;
 		canal[i-1] = canal[i];
		canal[i] = NULL;
 	};
	return 1;
};

int CConfigComtrade::InicializarMinimosEMaximos(float minimo, float maximo) 
{
    int k;
    for (k=0; k<num_canais_analogicos+num_canais_digitais; k++)
	{
        if ( (canal[k]!=NULL) && (canal[k]->GetAnalogico() ) ) 
		{
            ((CCanalAnalogico*)canal[k])->SetValorMinimo( minimo);
            ((CCanalAnalogico*)canal[k])->SetValorMaximo( maximo);
        };
	}
    return 1;
};

int CConfigComtrade::SetFreqLinha(float frequencia){
    if (frequencia>=0) freq_linha = frequencia;
    else return 0;
    return 1;
};

int CConfigComtrade::AcrescentarTaxa(float taxa, long indice) {
	
	ASSERT(numero_taxas>=0);
    if (numero_taxas<__num_max_taxas) 
	{
        taxas[numero_taxas].SetTaxaAmostragem(taxa);
        taxas[numero_taxas].SetUltimaAmostra(indice);
        numero_taxas++;
    } 
	else return 0;
    return numero_taxas;
};

int CConfigComtrade::AlterarTaxa(int numTaxa, float taxa, long indice)
{
	if ( (numTaxa<=numero_taxas) && (numTaxa>0) ){
		if (taxa>=0) taxas[numTaxa-1].SetTaxaAmostragem(taxa);
		if (indice>=0) taxas[numTaxa-1].SetUltimaAmostra(indice);
	  } else return 0;
	return numTaxa;
}

int CConfigComtrade::SetInstantePrimeiraAmostra(const __instante* instante){
    memcpy( &primeira_amostra, instante, sizeof(__instante));
    return 1;
};

int CConfigComtrade::SetInstanteTrigger(const __instante* instante){
    memcpy( &trigger, instante, sizeof(__instante));
    return 1;
};

int CConfigComtrade::SetOffsetInicioDados(long offset) {
    if (offset>=0) offset_inicio_dados = offset;
    else return 0;
    return 1;
};

int CConfigComtrade::SetTipoArquivo( char tipo) {
    tipo_arquivo = tipo;
    return 1;
};


int CConfigComtrade::SetVersao( int num) {
    versao = num;
    return 1;
};

int CConfigComtrade::SetFator( float fator_){
    if (fator>=0) fator = fator_;
    else return 0;
    return 1;
};

void CConfigComtrade::Copiar( const CConfigComtrade* config_orig)
{
	int canais[__num_max_canais];
	int c, N;
	N = config_orig->GetNumCanaisAnalogicos() +
	    config_orig->GetNumCanaisDigitais();
	for (c=0; c<N; c++) canais[c] = c+1;
	Copiar(config_orig, N, canais, 0,
		   config_orig->GetNumeroAmostrasPorCanal());
};

void CConfigComtrade::Copiar(
				const CConfigComtrade *config_orig,
				int numero_canais, 
				int *canais,
	    	    unsigned long amostra_inicial, 
				unsigned long amostra_final,
	    	    char soma)
{
	int i, can, j;
	unsigned long ultima;
	CCanal *pC;
	const CTaxas *taxa;
	float tempo;
	__instante inst;

    SetSimulacao( config_orig->GetSimulacao());
    SetComprimido( config_orig->GetComprimido());
    SetNomeSubstacao( config_orig->GetNomeSubstacao());
    SetIdentificadorRDP( config_orig->GetIdentificadorRDP());
    SetFreqLinha( config_orig->GetFreqLinha());
	inst=config_orig->GetInstantePrimeiraAmostra();
	SetInstantePrimeiraAmostra( &inst);
	inst = config_orig->GetInstanteTrigger();	SetInstanteTrigger( &inst);
	SetOffsetInicioDados(config_orig->GetOffsetInicioDados());
	SetTipoArquivo(config_orig->GetTipoArquivo());
	SetVersao(config_orig->GetVersao());
	SetFator(config_orig->GetFator());
	// acerta canais
	num_canais_analogicos = num_canais_digitais = 0;
	for (i=0; i<numero_canais; i++) 
	{
	  pC = config_orig->GetCanalPtr(canais[i]);
	  if (pC->GetTipoDoCanal()>1) //analogicos
	  {
		can = AcrescentarCanalAnalogico();
		if (can>0)
		  ((CCanalAnalogico*)GetCanalPtr(can))->Copiar((CCanalAnalogico*)pC);
	  };
	};
	for (i=0; i<numero_canais; i++) 
	{
	  pC = config_orig->GetCanalPtr(canais[i]);
	  if (pC->GetTipoDoCanal()<2) //digitais
	  {
		can = AcrescentarCanalDigital();
		if (can>0) GetCanalPtr(i+1)->Copiar(pC);
	  };
	};
	// acerta taxas
	numero_taxas = 0; 
	tempo = 0.F; 
	ultima = 0; 
	j=0;
	for (i=0;(i<config_orig->GetNumeroTaxas()) && (ultima<amostra_final); i++) 
	{
		taxa = config_orig->GetTaxa(i+1);
		if ((unsigned long)taxa->GetUltimaAmostra()> amostra_inicial) 
		{
			AcrescentarTaxa( taxa->GetTaxaAmostragem(),
						     min( amostra_final,(unsigned long)taxa->GetUltimaAmostra() ) - amostra_inicial);
			if (!j) tempo += (amostra_inicial-ultima) / taxa->GetTaxaAmostragem();
			j=1;
		} 
		else 
		{
		    if (!j)
		    tempo+=(taxa->GetUltimaAmostra()-ultima)/taxa->GetTaxaAmostragem();
		};
	  	ultima=taxa->GetUltimaAmostra();
	};
    // acerta tempo inicial
    if (soma) SomarAoTempoDaPrimeiraAmostra( tempo);
    else SomarAoTempoDaPrimeiraAmostra( -tempo);
};

void CConfigComtrade::SomarAoTempoDaPrimeiraAmostra( float tempo)
{
	__instante *t;
	struct tm tm1, *tm2;
	time_t segundos;
	float tmp;
	
    t = (__instante*) calloc( 1, sizeof(__instante));
    memcpy( t, GetInstantePrimeiraAmostra(), sizeof(__instante));

    tmp = (float)t->tm_sec + tempo;
    if ( (tmp>=60) || (tmp<0) )
	{
        tm1.tm_year=t->tm_year; 
		tm1.tm_mon=t->tm_mon; 
		tm1.tm_mday=t->tm_mday;

        tm1.tm_hour=t->tm_hour;
		tm1.tm_min=t->tm_min;
		tm1.tm_sec=0;
		
        tm1.tm_isdst=0;

        if (tmp<0) tmp -= 1;

        segundos = mktime(&tm1) + (int)tmp;
        tm2 = localtime( &segundos);

        if (tm2!=NULL) 
		{
			t->tm_year=tm2->tm_year; 
			t->tm_mon=tm2->tm_mon; 
			t->tm_mday=tm2->tm_mday;

			t->tm_hour=tm2->tm_hour;
			t->tm_min=tm2->tm_min;
			t->tm_sec=(float)tm2->tm_sec;

			if (tmp<0) t->tm_sec += 1 + (tmp - (int)tmp);
			else t->tm_sec += tmp - (int)tmp;
        };
    } 
	else t->tm_sec = tmp;

	SetInstantePrimeiraAmostra(t);

	free(t);
};

#ifndef __WINDOWSMFC
/*******************            class CFile            ********************/

CFile::CFile() {
    strcpy( arquivo,"");
    aberto = 0;
};

CFile::~CFile() { Close();};

char CFile::Open( const char* arquivo_, unsigned int modo, char* erro) 
{
    if (strlen(arquivo_)>_MAX_PATH) 
	{
        if (erro==NULL) return 0;
        sprintf(erro,"Nome do arquivo possui mais de %d caracteres",_MAX_PATH);
        return 0;
    };
    strcpy( arquivo, arquivo_);
    if ((modo & 0xE000)!=0x8000) 
	{
        if (erro==NULL) return 0;
        sprintf(erro,"Arquivo deve ser aberto em modo binario",_MAX_PATH);
        return 0;
    };
    if (aberto) Close();

    switch(modo & 3) {
      case 0: // leitura
        if ((arq = fopen(arquivo, "rb")) == NULL) {
            if (erro==NULL) return 0;
            sprintf( erro, "Arquivo %s nao pode ser aberto.", arquivo);
            return 0;
          };
        break;
      case 1: // escrita
        if ((arq = fopen(arquivo, "wb")) == NULL) {
            if (erro==NULL) return 0;
            sprintf( erro, "Arquivo %s nao pode ser aberto.", arquivo);
            return 0;
          };
        break;
      case 2: // acrescenta ao conteudo ja existente (append)
        if ((arq = fopen(arquivo, "a+b")) == NULL) {
            if (erro==NULL) return 0;
            sprintf( erro, "Arquivo %s nao pode ser aberto.", arquivo);
            return 0;
          };
        break;
      };
    aberto = 1;
    return 1;
};

void CFile::Close() {
    if (aberto) fclose( arq);
    aberto = 0;
};

unsigned int CFile::Read(char* lpBuf, unsigned int nCount) {
    if (!aberto) return 0;
    return fread( lpBuf, 1, (size_t)nCount, arq);
};

void CFile::Write( const char* lpBuf, unsigned int nCount ) {
    if (!aberto) return;
    fwrite( lpBuf, 1, (size_t)nCount, arq);
};

long CFile::Seek( long lOff, unsigned int nFrom ) {
    int auxUI;

    if (strcmp(arquivo,"")==0) return 0;
    if (!aberto)
        if ((arq = fopen(arquivo, "rb")) == NULL) return -1;
    switch (nFrom) {
      case begin: auxUI = SEEK_SET; break;
      case current: auxUI = SEEK_CUR; break;
      case end: auxUI = SEEK_END; break;
      };
    if (fseek( arq, lOff, auxUI)!=0) return -1;
    aberto = 1;

    return lOff;
};

unsigned long CFile::GetLength()
{
	unsigned long curpos, length;
	if (strcmp(arquivo,"")==0) return 0;
	if (!aberto)
		if ((arq = fopen(arquivo, "rb")) == NULL) return 0;
	curpos = ftell(arq);
	fseek(arq, 0L, SEEK_END);
	length = ftell(arq);
	fseek(arq, curpos, SEEK_SET);
	if (!aberto) fclose( arq);

	return length;
};

#endif
/*******************           class CBuffer           ********************/

CBuffer::CBuffer(size_t memoria_disponivel, char* arquivo_) {
    long tambuf;
	// inicializa ponteiros e variáveis
    arqptr = bufptr = bufend = 0L;
    haDadosSalvar = 0;
    arq = NULL;
	// determina memória disponível
    memoria_disponivel = min(memoria_disponivel,(size_t)65000);
    tambuf =( ((memoria_disponivel)==(0)) ? (65000*2):
                                            ((long)memoria_disponivel*2) );
	// ALoca espaço para o buffer de leitura
    buffer = NULL;
    while (buffer==NULL) 
	{
        tambuf = tambuf / 2;
        buffer = (char*) malloc(size_t(tambuf));
    };
	buffer[0]=0;
    buftam = size_t(tambuf);
	// Prepara arquivo para leitura
    if (arquivo_!=NULL) SetNomeArquivo( arquivo_);
                else strcpy( Arquivo, "");
};

CBuffer::~CBuffer() { 
    if (haDadosSalvar) SalvarBuffer(0);
    if (arq!=NULL) 
	{
    	delete arq;
    	arq = NULL;
    };
    if (buffer!=NULL) 
	{
    	free(buffer);
    	buffer = NULL;
    };
};

void CBuffer::ResetBuffer() { 
    buffer[0]=0;
    arqptr = bufptr = bufend = 0L;
    haDadosSalvar = 0;
};

int CBuffer::SetNomeArquivo(const char* arquivo_,unsigned int modo,char* erro){
    if (haDadosSalvar) SalvarBuffer(0);
    if (strlen(arquivo_)>_MAX_PATH) 
	{
        if (erro==NULL) return 0;
        sprintf(erro,"Nome do arquivo possui mais de %d caracteres",_MAX_PATH);
        return 0;
    };
    strcpy( Arquivo, arquivo_);
    if (arq!=NULL) delete arq;
    arq = new CFile();
    if (!arq->Open( arquivo_, modo | CFile::typeBinary)) 
	{
        if (erro==NULL) return 0;
        sprintf(erro,"Nao foi possivel abrir arquivo %s", Arquivo);
        return 0;
    };
    buffer[0]=0;
    arqptr = bufptr = bufend = 0L;
    haDadosSalvar = 0;
    return 1;
};

void CBuffer::SetCFilePtr(CFile *pF) {
    if (pF==NULL) return;
    if (haDadosSalvar) SalvarBuffer(0);
    if (arq!=NULL) delete arq;
    arq = pF->Duplicate();
    buffer[0]=0;
    arqptr = bufptr = bufend = 0L;
    haDadosSalvar = 0;
};

long CBuffer::filesize( char *erro) 
			{ return (long)arq->GetLength(); }

int CBuffer::EncherBuffer( long posicao, char *erro) {
    if (arq->Seek( posicao, CFile::begin)!=posicao) 
	{
        if (erro==NULL) return 0;
        sprintf( erro, "Nao consegui localizar posicao %ld", posicao);
        return 0;
    };
    bufend = arq->Read( buffer, buftam);
    if (bufend==0)  
	{
        if (erro==NULL) return 0;
        sprintf( erro, "Nao consegui ler nenhum byte do arquivo %s", Arquivo);
        return 0;
    };
    bufptr = 0;     
	arqptr = posicao;
    
	return 1;
};

int CBuffer::EncherBuffer(char * arquivo_, long posicao, char *erro) {
    SetNomeArquivo( arquivo_);
    if (!EncherBuffer(posicao, erro)) 
	{
        if (erro==NULL) return 0;
        sprintf( erro, "Arquivo %s nao pode ser aberto.", Arquivo);
        return 0;
    };
    return 1;
};

int CBuffer::SalvarBuffer(int dummy, char *erro) {
    arq->Write( buffer, bufptr);
    bufend = bufptr = 0;
    haDadosSalvar = 0;
    dummy = 1; return dummy;    
};

int CBuffer::SalvarBuffer(char * arquivo_, char *erro) {
    SetNomeArquivo( arquivo_);
    if (erro==NULL) 
	{
        if (!SalvarBuffer(0,NULL)) return 0;
    } 
	else 
	{
        if (!SalvarBuffer(0,erro)) return 0;
    };

    return 1;
};

BYTE CBuffer::ExtrairByte(char* Fim) {
    BYTE auxB;
    if (Fim!=NULL) 
	{
        if (ExtrairBytes( &auxB, sizeof(BYTE))==0) *Fim = 1;
		else *Fim = 0;    
	};

    return auxB;
};

WORD CBuffer::ExtrairInteiro(char* Fim) {
    WORD auxW;
    if (Fim!=NULL) 
	{
        if (ExtrairBytes( &auxW, sizeof(WORD))==0) *Fim = 1;
        else *Fim = 0;    
	};
    return auxW;
};
DWORD CBuffer::ExtrairLong(char* Fim) {
    DWORD auxD;
    if (Fim!=NULL) 
	{
        if (ExtrairBytes( &auxD, sizeof(DWORD))==0) *Fim = 1;
		else *Fim = 0;    
	};
    return auxD;
};

float CBuffer::ExtrairReal(char* Fim) {
    float auxR;
    if (Fim!=NULL) 
	{
        if (ExtrairBytes( &auxR, sizeof(float))==0) *Fim = 1;
		else *Fim = 0;    
	};
    return auxR;
};

int CBuffer::ExtrairString(char* string, size_t tam) { // terminada por \0
    int ret = 1;
    size_t cont = 0;
    if (string==NULL) return 0;
    if (bufptr>=bufend) ret = EncherBuffer(bufptr+arqptr);

    while ((cont<tam) && (buffer[bufptr]!=0) && (ret!=0)) 
	{
        string[cont] = buffer[bufptr];
        cont++; bufptr++;
        if (bufptr>=bufend) ret = EncherBuffer(bufptr+arqptr);
    };
    string[cont] = '\0';    bufptr++;
    if (bufptr>=bufend) ret = EncherBuffer(bufptr+arqptr);
    return (cont<=tam) && ret;
};

int CBuffer::PularLinha() {
//terminada por <CR>, <LF> ou <CR>+<LF>
    int ret = 1, c = 0;
  //  if (string==NULL) return 0;
    if (bufptr>=bufend) ret = EncherBuffer(bufptr+arqptr);
    while ((buffer[bufptr]!=13) && (buffer[bufptr]!=10) && (ret!=0)) 
	{
        bufptr++;
		c++;
        if (bufptr>=bufend) ret = EncherBuffer(bufptr+arqptr);
    };
    bufptr++;
    if (bufptr>=bufend) ret = EncherBuffer(bufptr+arqptr);
    if (buffer[bufptr]==10) bufptr++;

    return c;
};

int CBuffer::ExtrairLinha(char* string, size_t tam) {
//terminada por <CR>, <LF> ou <CR>+<LF>
    int ret = 1;
    size_t cont = 0;
    if (string==NULL) return 0;     
    if (bufptr>=bufend) ret = EncherBuffer(bufptr+arqptr);
    while ((cont<tam) && (buffer[bufptr]!=13) &&
            (buffer[bufptr]!=10) && (ret!=0)) 
	{
        string[cont] = buffer[bufptr];
        cont++; bufptr++;
        if (bufptr>=bufend) ret = EncherBuffer(bufptr+arqptr);
    };
    string[cont] = '\0';     
	bufptr++;
    
    if (bufptr>=bufend)
        ret = EncherBuffer(bufptr+arqptr) | (strlen(string)>0);
    
    if (buffer[bufptr]==10) bufptr++;

    return (cont<=tam) && ret;
};

int CBuffer::ExtrairLexema(char* string, size_t tamMax, char separador)
{
    int ret = 1;
    size_t cont = 0;

    if (string==NULL) return 0;
    if (bufptr>=bufend) ret = EncherBuffer(bufptr+arqptr);

    while ((cont<tamMax) && (buffer[bufptr]!=separador) && 
        (ret!=0) && (buffer[bufptr]!='\r') && (buffer[bufptr]!='\n')) 
	{
        string[cont] = buffer[bufptr];
        cont++; bufptr++;
        if (bufptr>=bufend) ret = EncherBuffer(bufptr+arqptr);
    };

    if (buffer[bufptr]=='\n') bufptr--;	
        
    string[cont] = '\0';    
	bufptr++;

    if (bufptr>=bufend) ret = EncherBuffer(bufptr+arqptr);

    return (cont<=tamMax) && ret;
};

int CBuffer::ExtrairBytes( void* vetor, size_t tam) { //sequencia de tam bytes
    int ret = 1;
    size_t resto = tam, auxI;
    char* vet = (char*)vetor;

    while ((resto>0) && (ret!=0))
	{
	    if ((bufend-bufptr)<resto)
            ret = EncherBuffer(bufptr+arqptr);
        auxI = min(resto,bufend);
        if (vetor!=NULL) 
			memcpy( &vet[tam-resto], &buffer[bufptr], auxI);
        resto -= auxI;      
		bufptr += auxI;
    };
    return ret;
};
    
int CBuffer::InserirByte(BYTE valor)
			{return InserirBytes( &valor, sizeof(BYTE));};
int     CBuffer::InserirInteiro(WORD valor)
			{return InserirBytes( &valor, sizeof(WORD));};
int CBuffer::InserirLong(DWORD valor)
			{return InserirBytes( &valor, sizeof(DWORD));};
int CBuffer::InserirReal(float valor)
			{return InserirBytes( &valor, sizeof(float));};
int CBuffer::InserirString(char* string)  //terminada por \0
			{return InserirBytes( string, strlen(string));};

int CBuffer::InserirBytes( void* vetor, size_t tam) {
    int ret = 1;
    size_t resto = tam, auxI;
    char* vet = (char*)vetor;

    while ((resto>0) && (ret!=0))
	{
        if ((buftam-bufptr)<resto) ret = SalvarBuffer(0);
        auxI = min(resto,buftam);
        memcpy( &buffer[bufptr], &vet[tam-resto], auxI);
        resto -= auxI;      
		bufptr += auxI;     
		bufend = bufptr;
    };
    haDadosSalvar = 1;
    return ret;
};

int CBuffer::InserirLinha(char* string){ 
	//string e' terminada por \0
    //sera' inserida terminada por <CR>+<LF>
    InserirString( string);
    InserirByte(13);
    return InserirByte(10);
};


/*******************            class CTradutor          ********************/

int CTradutor::ExtrairConfiguracao(CConfigComtrade* config, char* erro) {
    if (erro==NULL) return 0;
    sprintf( erro, "ExtrairConfiguracao nao esta' implementada");
    return 0;
};

int CTradutor::InicializarParametrosAdicionais( CConfigComtrade* config,
                                                char* erro)
{
    int cont;
    CCanalAnalogico *pAna;

    for ( cont=1;
          cont<=config->GetNumCanaisAnalogicos()+config->GetNumCanaisDigitais();
          cont++) 
	{
        pAna = (CCanalAnalogico *) config->GetCanalPtr(cont);
        if (pAna->GetAnalogico()) 
		{
            if (pAna->GetValorMinimo()==pAna->GetValorMaximoCfg())
                pAna->SetValorMaximo(pAna->GetValorMaximoCfg() + 1);
        }
    };
	
    return 1;
};

int CTradutor::ExtrairDados( CConfigComtrade* config, 
                        int* numero_canais, 
						size_t tam, 
						void* canais,
                        size_t numpts, 
						long offset,
                        char AcessoRandomico, 
						char* erro) {
    if (erro==NULL) return 0;
    sprintf( erro, "ExtrairDados nao esta' implementada");
    return 0;
};

int CTradutor::InserirConfiguracao(CConfigComtrade* config, char* erro) {
    if (erro==NULL) return 0;
    sprintf( erro, "InserirConfiguracao nao esta' implementada");
    return 0;
};

int CTradutor::InserirDados(CConfigComtrade* config, 
							int* numero_canais,
							size_t tam, void* canais, 
							size_t numpts, 
							long numero_ultima_amostra, 
							char* erro){
    if (erro==NULL) return 0;
    sprintf( erro, "InserirDados nao esta' implementada");
    return 0;
};

size_t AlocarVetoresDeDados( void *vetor_apontadores, 
							long numero_pontos,
                            int numero_canais, 
							int* posicoes,
                            CConfigComtrade *config)
{
    size_t pontos;
    void **vet;
    char bSemMem = 1;
    vet = (void**)vetor_apontadores;

    for (pontos=0; pontos<(unsigned)numero_canais; pontos++)
        vet[pontos] = NULL;

    numero_pontos *= 2;
    while (bSemMem && (numero_pontos>2)) 
	{
        bSemMem = 0;    
		numero_pontos /= 2;

        for (pontos=0; pontos<(unsigned)numero_canais; pontos++)
            if(vet[pontos]!=NULL)free(vet[pontos]);

        for (pontos=0; pontos<(unsigned)numero_canais; pontos++)
		{
            if ((config->GetCanalPtr(posicoes[pontos]))->GetAnalogico()) 
			{
                switch (((CCanalAnalogico*)
                        (config->GetCanalPtr(posicoes[pontos])))->GetFormatoDados()){
                  case e_inteiro:
                    bSemMem |= (vet[pontos] = (int*)
                            malloc((size_t)numero_pontos*sizeof(int)))==NULL;
                  break;
                  case e_real:
                    bSemMem |= (vet[pontos] = (float*)
                            malloc((size_t)numero_pontos*sizeof(float)))==NULL;
                  break;
                  case e_dupla: // double
                    bSemMem |= (vet[pontos] = (double*)
                            malloc((size_t)numero_pontos*sizeof(double)))==NULL;
                  break;
                }; // switch
            } 
			else 
			{
				bSemMem |= (vet[pontos]=(char*)malloc((size_t)numero_pontos*sizeof(char)))==NULL;
            }; //if Analogico
		} //for
    }; //while

    return (size_t)numero_pontos;
};
