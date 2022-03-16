#ifndef NEWCOMTD99_H
#define NEWCOMTD99_H

/*****************************************************************************

   Projeto : SISTEMA INTEGRADO DE APOIO `A ANALISE DE PERTURBACOES - SINAPE
   Nome    : NEWCOMTD.H
   Funcao  : Rotinas para conversão entre formato interno (CConfigCOMTRADE)
             e o formato COMTRADE (IEEE standard Common Format for Transient
             Data Exchange for power systems), IEEE Std C37.111-1999.
             Este tradutor gera COMTRADE em formato binario
   C       : C++
   Autora  : Suelaine
   Data    : jun/2002

   Data    : data da primeira modificacao
   Modific.: descricao dos itens modificados
   ...
   Data    : data da ultima modificacao
   Modific.: descricao dos itens modificados

*****************************************************************************/

#include "newcomtd.h"

class CComtrade99: public CComtrade97  {

public:           

    CComtrade99(size_t memoria_disponivel=0, char* arquivo_=NULL)
               :CComtrade97( memoria_disponivel, arquivo_) {};

    ~CComtrade99() {};

// preenche config a partir do arquivo registrado em SetNomeArquivo
int ExtrairConfiguracao(CConfigComtrade* config, char* erro=NULL);

int ExtrairConfiguracao_(CConfigComtrade* config, char* erro=NULL);

// criar arquivo a partir de config no formato do tradutor
int InserirConfiguracao(CConfigComtrade* config, char* erro=NULL);

int InserirConfiguracao_(CConfigComtrade* config, char* erro=NULL);

// New
int LerLinhaCanalAnalogico(CConfigComtrade* config, CCanal *canal, int cont, char* erro);

// New
int SetConversao (char prim_sec);

};


#define UltimoComtrade CComtrade99

class CNewComtrade: public UltimoComtrade  {

public:           

    CNewComtrade(size_t memoria_disponivel=0, char* arquivo_=NULL)
               :UltimoComtrade( memoria_disponivel, arquivo_) {};

    ~CNewComtrade() {};
};

#endif
