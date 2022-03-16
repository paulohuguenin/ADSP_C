#ifndef NEWCOMTD_H
#define NEWCOMTD_H

/*****************************************************************************

   Projeto : SISTEMA INTEGRADO DE APOIO `A ANALISE DE PERTURBACOES - SINAPE
   Nome    : NEWCOMTD.H
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

#include "bincomtd.h"
//#include "newcomtd_99.h"

class CComtrade97: public CBinaryComtrade  {

public:           

    CComtrade97(size_t memoria_disponivel=0, char* arquivo_=NULL)
               :CBinaryComtrade( memoria_disponivel, arquivo_) {};

    ~CComtrade97() {};

	// preenche config a partir do arquivo registrado em SetNomeArquivo
	int ExtrairConfiguracao(CConfigComtrade* config, char* erro=NULL);

	int ExtrairConfiguracao_(CConfigComtrade* config, char* erro=NULL);

	// criar arquivo a partir de config no formato do tradutor
	int InserirConfiguracao(CConfigComtrade* config, char* erro=NULL);

	int InserirConfiguracao_(CConfigComtrade* config, char* erro=NULL);

	// New

	int LerLinhaCanalAnalogico(CConfigComtrade* config, CCanal *canal, int cont, char* erro=NULL);

	int LerLinhaCanalDigital(CConfigComtrade* config, CCanal *canal, int cont, char* erro=NULL);


};


#endif
