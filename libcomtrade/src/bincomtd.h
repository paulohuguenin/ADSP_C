#ifndef BINCOMTD_H
#define BINCOMTD_H
/*****************************************************************************

   Projeto : SISTEMA INTEGRADO DE APOIO `A ANALISE DE PERTURBACOES - SINAPE
   Nome    : BINCOMTD.H
   Funcao  : Rotinas para conversão entre formato interno (CConfigCOMTRADE)
             e o formato COMTRADE (IEEE standard Common Format for Transient
             Data Exchange for power systems), IEEE C37.111-1991.
             Este tradutor gera COMTRADE em formato binario
   C       : C++
   Autor   : MAMR
   Data    : mar/95

   Data    : data da primeira modificacao
   Modific.: descricao dos itens modificados
   ...
   Data    : data da ultima modificacao
   Modific.: descricao dos itens modificados

*****************************************************************************/

#include "comtrade.h"

class CBinaryComtrade: public CComtrade  {

public:           

    CBinaryComtrade(size_t memoria_disponivel=0, char* arquivo_=NULL)
               :CComtrade( memoria_disponivel, arquivo_) {};

    ~CBinaryComtrade() {};

	// preenche config a partir do arquivo registrado em SetNomeArquivo
	int ExtrairConfiguracao(CConfigComtrade* config, char* erro=NULL);

	int ExtrairConfiguracao_(CConfigComtrade* config, char* erro=NULL);

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
	int ExtrairDados( CConfigComtrade* config, 
					int* numero_canais, 
					size_t tam, 
					void* canais,
                    size_t numpts, 
					long offset=0,    
					char AcessoRandomico=0, 
					char* erro=NULL);

	// criar arquivo a partir de config no formato do tradutor
	int InserirConfiguracao(CConfigComtrade* config, char* erro=NULL);

	// Acrescenta dados ao arquivo sequencialmente
	// canais e' um vetor de apontadores para vetores do tipo float, char ou int
	//   que conterao os valores. Cada 1 dos 3 casos precisa ser tratado.
	// "canais[n]" e' inserido no arquivo na posicao indicada em "numero_canais[n]"
	//   e as posicoes intermediarias sao preenchidas com valores invalidos
	// numero_ultima_amostra e' o valor da ultima amostra acrescentada ao arquivo
	//   de dados do COMTRADE. E' usado no caso de multiplas chamadas a rotina para
	//   o mesmo arquivo.
	int InserirDados( CConfigComtrade* config, 
					int* numero_canais, 
					size_t tam, 
					void* canais,
					size_t numpts, 
					long numero_ultima_amostra=0,
                    char* erro=NULL);

// Obs.: CBuffer::SetNomeArquivo() deve ser usado para inicializar o nome do
//      arquivo correspondente

};

#endif
