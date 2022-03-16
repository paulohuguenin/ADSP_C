#ifndef BUFTIPOS_H
#define BUFTIPOS_H
/*****************************************************************************

   Projeto : SINAPE
   Nome    : BUFTIPOS.H
   Funcao  : Arquivo com declaracao dos formatos internos de representacao
            dos dados dos RDPs, baseados no COMTRADE
   C       : C++
   Autor   : MAMR
   Data    : 14/11/94
*****************************************************************************/

#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>

/*  estruturas usuais de arquivos */
typedef unsigned char       BYTE;
typedef unsigned short      WORD;
typedef unsigned int       DWORD;

#define _MAX_PATH 260

/*  estruturas baseadas nos formatos de arquivos do COMTRADE */
#ifndef _DEMO
	#define __num_max_canais 500
	#define __num_max_taxas  5
#endif

#ifdef _DEMO
	#define __num_max_canais 12
	#define __num_max_taxas  1
#endif

#define __erro_tam  255  // tamanho maximo da string de erro
#define __num_char_fase	 11

enum e_formato_buffer    {e_inteiro=0, e_real=1, e_dupla=2};
enum e_tipo_canal        {e_tensao=2,e_corrente,e_impedancia,e_DFT,e_outros};
enum e_tipo_arquivo      {e_ascii=0, e_binary=1};

/*  rotinas de uso generico */
void ExtrairEspacosDasBordas( char* string);

void StringReversal(char* &string);

void StringToUpper (char* &string);

/*  classes de traducao */


class CCanal {

    int numero;
    char *nome;
    unsigned char estado_normal; /* vide Obs.*/
	char *fase;				// eg.: A,B,C,N,R,Y,B
	char *equipamento_monitorado;

public:
	CCanal();
	CCanal(const CCanal& canal);
	CCanal(const CCanal* canal);
	void Copiar(const CCanal* canal);	
	~CCanal();

	// Funcoes Get
	int            GetNumeroCanal() const;
	char*          GetNomeCanal() const;
	unsigned char  GetTipoDoCanal() const;
	unsigned char  GetAnalogico() const;
	char*	       GetFaseCanal() const;
	char*          GetEquipamentoCanal() const;

	// Funcoes Set
	int            SetNumeroCanal(int numero_);
	int            SetNomeCanal(const char* nome_60);
	int            SetTipoDoCanal(unsigned char tipo);
	// Obs.: 0 -> canal digital normalmente aberto (0=fechado; 1=aberto);
	//       1 -> canal digital normalmente fechado(0=aberto; 1=fechado);
	//       2 -> tensao; 
	//       3 -> corrente; 
	//       4 -> impedancia; 
	//       5 -> DFT; 
	//       6 -> outros.
	int            SetFaseCanal(const char* fase_);
	int            SetEquipamentoCanal(const char* equipamento);
};

class CCanalAnalogico : public CCanal{

    int		formato;			// {e_inteiro, e_real, e_dupla}
    char	*unidade;			// eg.: kV, kA, Ohms
    float	coef_ang;			// Obs.1
    float	coef_lin;
    float	time_skew;		// tempo do inicio da amostragem `a amostragem do canal
    float	valor_minimo_range; // vide Obs.2
    float	valor_maximo_range;
	float	primario;
	float	secundario;
	char	tipo_escala;		// 0-P e 1-S
	char	*sret;
		
public: 
	CCanalAnalogico();
	CCanalAnalogico(const CCanalAnalogico& canal);
	CCanalAnalogico(const CCanalAnalogico* canal);
	void Copiar(const CCanalAnalogico* canal);
	~CCanalAnalogico();
	
	// Funcoes Get
	int     GetFormatoDados() const;
	char*   GetUnidadeCanal() const;
	float   GetCoefAng(int PrimSec=0) const;
	float   GetCoefLin(int PrimSec=0) const;
	float   GetTimeSkew() const;
	float   GetValorMinimo(int PrimSec=0) const;
	float   GetValorMaximo(int PrimSec=0) const;
	float   GetValorMinimoCfg() const;
	float   GetValorMaximoCfg() const;
	float   GetPrimario() const;
	float   GetSecundario() const;
	char	GetTipoEscalamento() const;

	// Funcoes Set
	int     SetFormatoDados(int formato_=e_real);
	int     SetUnidadeCanal(const char* unidade_);
	int     SetCoefAng(float coeficiente=1);
	int     SetCoefLin(float coeficiente=0);
	int     SetTimeSkew(float skew);
	int     SetValorMinimo(float valor);
	int     SetValorMaximo(float valor);
	int     SetPrimario(float valor);
	int     SetSecundario(float valor);
	int     SetTipoEscalamento( char escala); // 0 é P (Primário) e 1 é S (Secundário)

	// Demais Funcoes
	float   TransformaPrimSec( float valor, int PrimSec, int valorPrimitivo) const;
	float   TransformaPrimSec( float valor, int PrimSec=0) const;
	// Corrige a unidade do arquivo  de acordo com sua identificação - V, A ou W.
	static double CorrigeUnidade(const char* unidade_, char* unidade_corrigida=NULL);
	// Usado para corrigir o campo de estado_normal; retorna 0 se nada alterou
	int     AcertarTipoDoCanalAnalogico();

// Obs.1:   valor = coef_ang*amostra + coef_lin;
// Obs.2:   se iguais e' por que informacao nao esta disponivel.
//          Se o formato for real contem os valores minimo e maximo,
//          se inteiro, contem o valor minimo e maximo do range.
//          coef_ang & coef_lin devem atuar sobre estes valores
// Obs.3:   time_skew e' guardado em us no arquivo comtrade porem em 
//          segundos no objeto CConfigComtrade
};


class CTaxas    {

	float	taxa_amostragem;  /* Hz */
    long	ultima_amostra; /* vide Obs1 */

public:
	CTaxas();

	// Funcoes Get
	float   GetTaxaAmostragem() const;
	long    GetUltimaAmostra() const;

	// Funcoes Set
	int     SetTaxaAmostragem(float taxa=1);
	int     SetUltimaAmostra(long indice);

// Obs.1:   ultima amostra do canal, i.e., contando apenas as
//          amostras de um canal, com a referida taxa de amostragem
};

// modificacao da estrutura tm para maior precisao nos segundos
typedef struct __instante_tag {
	double	tm_sec; // Seconds after the minute
	int		tm_min;    // Minutes after the hour - [0,59]
	int		tm_hour;   // Hours since midnight - [0,23]
	int		tm_mday;   // Day of the month - [1,31]
	int		tm_mon;    // Months since January - [0,11]
	int		tm_year;   // Years since 1900 
	int		tm_year1;  // Novo
	int		tm_wday;   // Days since Sunday - [0,6]
	int		tm_yday;   // Days since January 1 - [0,365]
	int		tm_isdst;  // Daylight-saving-time flag
  } __instante;

class CConfigComtrade   {
    
	unsigned char simulacao;   /* false=0 true>0 */
    unsigned char comprimido;   /* false=0 true>0 */
    char		*nome_subestacao;
    char		*identificador_RDP;
    int			num_canais_analogicos;
    int			num_canais_digitais;
    CCanal		*canal[__num_max_canais];
    float		freq_linha;    /* em Hz */
    char		numero_taxas;/*num. de difer. taxas de amost.*/
    CTaxas		taxas[__num_max_taxas];
    __instante	primeira_amostra;
    __instante	trigger;
    long		offset_inicio_dados;
    char		tipo_arquivo;/* valor definido pelo tradutor*/
	int			versao; //ano da versao do Padrao Comtrade
	float		fator; //f.multiplic. da coluna tempo .DAT

public:
	CConfigComtrade();
    CConfigComtrade(const CConfigComtrade& config_orig);
    ~CConfigComtrade();

	//Funcoes Get
	char			GetSimulacao() const;
	char			GetComprimido() const;
	char*			GetNomeSubstacao() const;
	char*			GetIdentificadorRDP() const;
	int				GetNumCanaisAnalogicos() const;
	int				GetNumCanaisDigitais() const;
	CCanal*			GetCanalPtr(int num_canal) const;// num_canal >=1
	float			GetFreqLinha() const;
	char			GetNumeroTaxas() const;
	const CTaxas*	GetTaxa(int indice) const;
	long			GetNumeroAmostrasPorCanal() const;
	__instante*		GetInstantePrimeiraAmostra();
	__instante*		GetInstanteTrigger();
	__instante		GetInstantePrimeiraAmostra() const;
	__instante		GetInstanteTrigger() const;
	long			GetOffsetInicioDados() const;
	char			GetTipoArquivo() const;
	int				GetVersao() const;
	float			GetFator() const;

	//Funcoes Set
	int				SetSimulacao(char simulacao_=0);
	int				SetComprimido( char comprimido_=0);
	int				SetNomeSubstacao( const char * nome);
	int				SetIdentificadorRDP( const char *identificador);
	int				SetFreqLinha( float frequencia);
	int			    SetInstantePrimeiraAmostra(const __instante* instante);
	int				SetInstanteTrigger(const __instante* instante);
	int				SetOffsetInicioDados(long offset);
	int				SetTipoArquivo( char tipo);
	int				SetNumeroTaxas(const char num);
	int 			SetVersao(int num);
	int 			SetFator(float fator_);

	// Demais Funcoes
	int     AcrescentarCanalAnalogico(int* num=NULL);//num=canal acrescentado >=1
	int     AcrescentarCanalDigital(int* num=NULL); //num=canal acrescentado,>=1
	int     InicializarMinimosEMaximos(float minimo=+3.4E+38,float maximo=-3.4E+38);
	int     AcrescentarTaxa(float taxa, long indice=0);//retorna num.taxa,>=1
	int     AlterarTaxa(int numTaxa, float taxa=1, long indice=0); // se taxa<0 sera' ignorado; idem para indice
	int		ApagarCanal(int num_canal);//apaga o canal num_canal >=1

		// inicializa objeto CConfigComtrade com um ja' existente
	void    Copiar( const		CConfigComtrade* config_orig,
				int				numero_canais, 
				int				*canais,
				unsigned long	amostra_inicial, 
				unsigned long	amostra_final,
				char			soma=0);
	void    Copiar( const CConfigComtrade* config_orig);
	void    SomarAoTempoDaPrimeiraAmostra( float tempo);
};

#ifndef __WINDOWSMFC
// Adaptacao para interfacear com MFC do Visual C++ do WINDOWS
class CFile {
    char    arquivo[_MAX_PATH];//nome do arquivo com configuracao/dados
    char    aberto;
    FILE    *arq;
    
public:
    enum OpenFlags {    // combinar com |
                modeRead =          0x0000,
                modeWrite =         0x0001,
                modeReadWrite =     0x0002,
                modeCreate=         0x1000,
                typeText =          0x4000,
                typeBinary =   (int)0x8000
                };
    enum SeekPosition { begin = 0x0, current = 0x1, end = 0x2 };

    CFile();
    ~CFile();

	char    Open(const char* arquivo_, unsigned int modo, char*erro=NULL);
	void    Close();
	unsigned int Read(char* lpBuf, unsigned int nCount);
	void    Write(const char* lpBuf, unsigned int nCount );
	long    Seek( long lOff, unsigned int nFrom );
	unsigned long GetLength();
	CFile*  Duplicate() { return NULL;};
};

#endif

#ifdef __WINDOWSMFC
    #include "stdafx.h"
#endif

//typedef unsigned int size_t
class CBuffer   {
    char*       buffer; // buffer para conter copia de arquivo original
    size_t      bufptr; // proxima posicao de acesso ao buffer
    size_t      bufend; // tamanho do conteudo do buffer
    size_t      buftam; // tamanho fisico do buffer
    CFile       *arq;   // objeto de arquivo
    char        Arquivo[_MAX_PATH];//nome do arquivo com configuracao/dados
    long        arqptr; // posicao do inicio do buffer em relacao ao arquivo
    char        haDadosSalvar;// indica se buffer deve ser salvo ou nao

public:
        
	CBuffer(size_t memoria_disponivel=0, char* arquivo_=NULL);
	virtual ~CBuffer();

	void ResetBuffer();
	size_t  posicao() {return bufptr;};
	size_t  tamanho() {return buftam;};

	char*	GetNomeArquivo(){ return Arquivo;};
	int		SetNomeArquivo(const char* arquivo_, unsigned int modo=0, char* erro=NULL);
		// modo ==modeRead : abre arquivo para leitura
		// modo ==(modeCreate|modeWrite) : cria arquivo para escrita,
		//        apagando versao anterior
		// modo ==modeReadWrite :
		//              no WindowsMFC: abre arquivo para escrita e leitura
		//              no DOS abre arquivo para escrita, adicionando ao final
		// no WindowsMFC usar CFile::modeXXXX
	void	SetApenasNomeArquivo(const char* arquivo){ strcpy( Arquivo, arquivo);};
		// usado no Visual C/MFC apenas para tornar mensagens de erro melhores
	void	SetCFilePtr(CFile *pF);
		// usado no Visual C/MFC no lugar de SetNomeArquivo()
		// o arquivo ja' deve ter sido previamente aberto, via File::Open()
	CFile*	GetCFilePtr(){return arq;};
	long	GetFilePos() {return bufptr+arqptr;};
	long	filesize( char *erro=NULL);//retorna tamanho do arquivo corrente

	virtual int	EncherBuffer( long posicao=0,char *erro=NULL);
	virtual int SalvarBuffer(int dummy, char *erro=NULL);
	int			EncherBuffer(char * arquivo_, long posicao=0, char *erro=NULL);
	int			SalvarBuffer(char * arquivo_, char *erro=NULL);

	BYTE    ExtrairByte(char* Fim=NULL);// Fim=1 se for fim do arquivo
	WORD    ExtrairInteiro(char* Fim=NULL);
	DWORD   ExtrairLong(char* Fim=NULL);
	float   ExtrairReal(char* Fim=NULL);
		// 4 rotinas a seguir retornam 0 se for fim do arquivo ou string for pequena
		// string dever estar alocada c/ tam+1 bytes em ExtrairString,
		//    ExtrairLinha e ExtrairLexema
	int     ExtrairBytes( void* vetor, size_t tam);//sequencia de tam bytes
	int     ExtrairString(char* string, size_t tam); //terminada por \0
	int     ExtrairLinha(char* string, size_t tam); 
	int		PularLinha(); // retorna numero de caracteres na linha ou 0=erro
                          //terminada por <CR>, <LF> ou <CR>+<LF>
	int     ExtrairLexema(char* string, size_t tamMax, char separador);

	int     InserirByte(BYTE valor);
	int     InserirInteiro(WORD valor);
	int     InserirLong(DWORD valor);
	int     InserirReal(float valor);
	int     InserirBytes( void* vetor, size_t  tam);
	int     InserirString(char* string); //terminada por \0
	int     InserirLinha(char* string);  //string e' terminada por \0
										 //sera' inserida terminada por <CR>+<LF>
};

class CTradutor: public CBuffer {

protected:
    long    arqtam;//tamanho do arquivo de oscilogramas
    char    m_bDummyOp;//operacao dummy (*)
    //(*) Utilizado quando necessario por exemplo ler o arquivo de dados
    //    apenas para determinar os limites maximos e minimos de cada canal
    //    e/ou o numero de pontos. A rotina ExtrairDados() pode utilizar
    //    esta informacao para tomar acoes diferentes em cada caso.

public:

    CTradutor(size_t memoria_disponivel=0, char* arquivo_=NULL)
        :CBuffer( memoria_disponivel, arquivo_){m_bDummyOp=0;};
    ~CTradutor() {};

	// preenche config a partir do arquivo registrado em SetNomeArquivo
	virtual int ExtrairConfiguracao(CConfigComtrade* config, char* erro=NULL);

	// Le dados para extrair, se necessario parametros do tipo minimos e maximos
	// numero de pontos etc, quando necessario.
	// Implementacao default apenas testa se valores maximos e minimos sao iguais
	// tornando-os difierentes se for o caso
	virtual int InicializarParametrosAdicionais( CConfigComtrade* config, char* erro=NULL);

	// * numero_canais ->	possui os numeros dos canais a serem extraidos, >=1.
	//   Os numeros dos canais nao devem ser repetidos e devem estar em ordem
	//   crescente.
	// * tam -> e' o tamanho dos vetores numero_canais e canais
	// * canais -> e' um vetor de apontadores para vetores do tipo int ou float ou char
	//   que conterao os valores extraidos
	//   canais[n] -> recebe o conteudo do canal "numero_canais[n]"
	// * numpts -> é o numero de pontos, em cada canal a ser lido
	// * offset (deve ser > 0) -> indica o offset em amostras da primeira amostra,
	//   em cada canal separadamente. Se "==0" assume acesso sequencial
	// * AcessoRandomico -> "==1" obriga o arquivo a ser buscado desde o inicio ate'
	virtual int ExtrairDados( CConfigComtrade*	config, 
		                    int*				numero_canais, 
							size_t				tam, 
							void*				canais,
							size_t				numpts, 
							long				offset=0,
							char				AcessoRandomico=0, 
							char*				erro=NULL);

	// Cria arquivo a partir de config no formato do tradutor
	virtual int InserirConfiguracao(CConfigComtrade* config, char* erro=NULL);

	// Acrescenta dados ao arquivo sequencialmente
	// * canais -> e' um vetor de apontadores para vetores do tipo float, char ou int
	//   que conterao os valores. Cada 1 dos 3 casos precisa ser tratado.
	// * canais[n] -> e' inserido no arquivo na posicao indicada em "numero_canais[n]"
	//   e as posicoes intermediarias sao preenchidas com valores invalidos
	// * numero_ultima_amostra e' o valor da ultima amostra acrescentada ao arquivo
	//   de dados do COMTRADE. E' usado no caso de multiplas chamadas a rotina para
	//   o mesmo arquivo.
	virtual int InserirDados( CConfigComtrade*	config, 
							int*				numero_canais, 
							size_t				tam, 
							void*				canais,
							size_t				numpts, 
							long				numero_ultima_amostra=0,
							char*				erro=NULL);


	// Obs.: CBuffer::SetNomeArquivo() deve ser usado para inicializar o nome do
	//      arquivo correspondente

	char SetEstadoDummyOp(char estado) {m_bDummyOp = estado;return m_bDummyOp;};

};

size_t AlocarVetoresDeDados( void *				vetor_apontadores, 
							long				numero_pontos,
                            int					numero_canais, 
							int*				posicoes,
                            CConfigComtrade		*config);

#endif
