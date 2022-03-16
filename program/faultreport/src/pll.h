#ifndef PLL_H
#define PLL_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "filemgr.h"
#include "structbook.h"
#include "dictionary.h"
#include "datasignal.h"

/************************************************************************/
/*****************************LOSONHO************************************/

/*
Correcoes feitas:

	-> A funcao abs() foi substituida por fabs(), que é a contra parte que aceita doubles
	como argumento. Em alguns compiladores C++, abs() pode ser usada também, mas a 
	mudança foi mantida por	portabilidade
	-> A funcao nao era contida, ou seja, devia ser chamada varias vezes, uma para cada
	elemento, e se fosse necessario usa-la novamente com um vetor diferente, nao era 
	possivel...  (0_o)  Todas as variaveis globais se tornaram locais e os parametros mudaram

Possíveis bugs futuros:

	-> A variável <periodo_interval_estim> é declarada como int, mas parece ser
	um double, pois é inicializada depois de uma divisão. Isso acontece em outras partes 
	do código, como se o objetivo fosse truncar o resultado da divisão. Nesse caso, é
	melhor usar round() ou floor(), pois o resultado de uma conversão deste tipo é
	inesperado */

//#define	FreqAmostragem												((double)(24390.0))								// 12000, 18000	, 24000	
#define	Ff_inicial										((double)(60.0))
//#define	comprimento_mem_sinal						200												// Fs/Ff_inicial*1.5
//#define	comprimento_mem_sinal_extra				    500												// 2000										
#define	ganho_laco										((double)(1.0))									// (estabiliza + lento)1 2 3 4 5 6 7 8 10 (esbaliza + rapido)
#define	pipi											((double)(6.283185305179586476925286766559))
//#define	pi												((double)(3.1415926535897932384626433832795))
#define	duracao_processamento						    0.0000075
#define	down_sampling_factor							1.0												// Fator de redução do custo computacional 20 pra 24000


void pll_lovo__(double* Vin, double* V_estim, int V_len, double* freq_estim, 
                double* amp_estim, double* wt_estim, double FreqAmostragem);

             //( Entrada Vetor , Fundamental estimada , Tamanho amostra , Frequencia da fundamental , 
             //  Amplitude da fundamental , Fase da fundamental )
             
#endif             
