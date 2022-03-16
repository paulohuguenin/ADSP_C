#ifndef COMTRADE_H
#define COMTRADE_H

/*****************************************************************************

   Projeto : SISTEMA INTEGRADO DE APOIO `A ANALISE DE PERTURBACOES - SINAPE
   Nome    : COMTRADE.H
   Funcao  : Rotinas para conversão entre formato interno (CConfigCOMTRADE)
            e o formato COMTRADE (IEEE standard Common Format for Transient
            Data Exchange for power systems), IEEE C37.111-1991
   C       : C++
   Autor   : MAMR
   Data    : nov/94

   Data    : data da primeira modificacao
   Modific.: descricao dos itens modificados
   ...
   Data    : data da ultima modificacao
   Modific.: descricao dos itens modificados

*****************************************************************************/

#include "buftipos.h"

#define COMTRADELIMITEDOQUANTIZADOR 65534  // quantizacao para 16 bits

char* ExtrairVirgula( char* string);

class CComtrade: public CTradutor  {

protected:
    e_formato_buffer m_eFormatoVetorExtraido;  
    int m_iVerificaAmostra;
    
public:
    CComtrade(size_t memoria_disponivel=0, char* arquivo_=NULL)
               :CTradutor( memoria_disponivel, arquivo_) {
		m_eFormatoVetorExtraido = e_real; 
		m_iVerificaAmostra = 0;
		m_timestamp = NULL;
	};

    ~CComtrade() {if (m_timestamp!=NULL) delete [] m_timestamp;};

	// preenche config a partir do arquivo registrado em SetNomeArquivo
	virtual int ExtrairConfiguracao(CConfigComtrade* config, char* erro=NULL);

	virtual int ExtrairConfiguracao_(CConfigComtrade* config, char* erro=NULL);

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
	virtual int ExtrairDados(CConfigComtrade* config, 
							int* numero_canais, 
							size_t tam, 
							void* canais,
							size_t numpts, 
							long offset=0,
							char AcessoRandomico=0, 
							char* erro=NULL);

	// criar arquivo a partir de config no formato do tradutor
	virtual int InserirConfiguracao(CConfigComtrade* config, char* erro=NULL);

	// Acrescenta dados ao arquivo sequencialmente
	// canais e' um vetor de apontadores para vetores do tipo float, char ou int
	//   que conterao os valores. Cada 1 dos 3 casos precisa ser tratado.
	// "canais[n]" e' inserido no arquivo na posicao indicada em "numero_canais[n]"
	//   e as posicoes intermediarias sao preenchidas com valores invalidos
	// numero_ultima_amostra e' o valor da ultima amostra acrescentada ao arquivo
	//   de dados do COMTRADE. E' usado no caso de multiplas chamadas a rotina para
	//   o mesmo arquivo.
	virtual int InserirDados( CConfigComtrade* config, 
							int* numero_canais, 
							size_t tam, 
							void* canais,
							size_t numpts, 
							long numero_ultima_amostra=0,
							char* erro=NULL);

	// Copia o conteudo do texto para o arquivo
	// Nao desaloca espaco ocupado por texto
	int InserirCabecalho( char* texto, size_t tam, char *erro=NULL);

	// Aloca o vetor texto com tam bytes contento o conteudo do arquivo
	// texto nao deve estar alocado
	// O usuario explicitamente desaloca texto usando "free(texto)"
	int ExtrairCabecalho( char* texto, size_t *tam, char *erro=NULL);

	// Usado para forcar o tradutor a transformar os valores lidos do arquivo
	// COMTRADE para reais. Neste caso pode ser interessante calcular o valor
	// exato da amostra, i.e., usando os coeficientes angular e linear,
	// especificados na norma.
	void SetVetorExtraido( e_formato_buffer tipo=e_inteiro)
				{m_eFormatoVetorExtraido = tipo;};                       
    

	void SetVerificaAmostra( int Verifica)
				{m_iVerificaAmostra = Verifica;};                 

// Obs.: CBuffer::SetNomeArquivo() deve ser usado para inicializar o nome do
//      arquivo correspondente

	int LerLinha2(CConfigComtrade* config, char* erro=NULL);

	int LerLinhaCanalAnalogico(CConfigComtrade* config, CCanal *canal,int cont, char* erro=NULL);

	int LerLinhaCanalDigital(CConfigComtrade* config, CCanal *canal, int cont, char* erro=NULL);
	
	double* GetTimeStamp();

	char* sret;
	double* m_timestamp;

protected:
	// Auxiliar de InserirConfiguracao()
	int InserirConfiguracao_(CConfigComtrade* config, char* erro=NULL);
};

#endif
