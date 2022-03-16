#include "pll.h"


/***********************************************************************************/
/*************************************	LOVO	************************************/

// V_len é o tamanho dos vetores Vin, V_estim, freq_estim, amp_estim e wt_estim, que deve ser igual (CUIDADO COM MEMORY LEAK!)
//       A função não checa se o espaço foi corretamente alocado
// OBS.: Usar argc e argv, para que a funcaoo nao gaste tanto espaco no HEAP

void pll_lovo__(double* Vin, double* V_estim, int V_len, double* freq_estim, double* amp_estim, double* wt_estim, double FreqAmostragem)
            // (  entrada  ,     saida      ,  entrada ,      saida        ,        saida     ,      saida      )          
            // ( Entrada Vetor , Fundamental estimada , Tamanho amostra , Frequencia da fundamental , Amplitude da fundamental , Fase da fundamental )
{
	double fase_estimada[2];
	double diff_phi[2], diff_phi_corrigida[2];
	double saida_real[2], saida_imag[2];
	double r_filt[2],/*[2]*/ i_filt[2]; //[2];		//double r_filt1; //double i_filt1;
	double real_volta=0, imag_volta=0;
    //double entrada=0;
	double senoide_sinc=0;
	double freq_estimada=0;
	double total_prod_real=0, total_prod_imag=0;
	//double total_prod_real=0, total_prod_imag=0;
	//double total_prod_real=0, total_prod_imag=0;
	double produto_interno=0;
	double Ff_var=0;
	double temp_c=0, temp_s=0;
	double aux, aux2, aux3, aux_real=0, aux_imag=0, aux_real1=0, aux_imag1=0;
   
	int periodo_interval_estim=0;
	//int contador_iteracoes=1;
	int i=0, index=0, pos=0, atual=0, ant=1;

	for(register int contador_iteracoes = 1; contador_iteracoes <= V_len; contador_iteracoes++)
	{
		//cout << "contador_iteracoes: " << contador_iteracoes << endl;
		/* Leitura da Entrada */	
        // entrada = Vin[contador_iteracoes]; 
			
		// Na primeira iteracao faz-se tudo em funcao da fundamental e um sistema relaxado
		if(contador_iteracoes == 1)
		{
			freq_estimada = Ff_inicial;
			periodo_interval_estim = (int)(FreqAmostragem/(freq_estimada));
			saida_real[atual] =0;
			saida_imag[atual] =0;
			//r_filt[atual]=1.0;
			//i_filt[atual]=0.0;
		}

		/* Fim das condicoes de atualizacao das variaveis */

		//armazena a nova amostra da entrada
		//sinal[pos] = entrada; 
			
		/* Achando a Fase */	
		/* Calculo do produto interno */
		
		aux = (pipi*freq_estimada)/FreqAmostragem;
		total_prod_real = 0;
		total_prod_imag = 0;
		
		for(i=0; i<periodo_interval_estim;i=i+down_sampling_factor)
		{
			aux2 = aux*((double)(i));
			temp_c = (double)(cos(aux2));
			temp_s = (double)(sin(aux2));
 			index = pos - periodo_interval_estim + i;
			//index = index % (comprimento_mem_sinal);
			index = (index < 0) ? 0 : index;
			total_prod_real = total_prod_real + (Vin[index])*temp_c;
			total_prod_imag = total_prod_imag + (Vin[index])*temp_s;
		}
		
		total_prod_real = ((double)down_sampling_factor)*total_prod_real/((double)periodo_interval_estim/2.0);
		total_prod_imag = ((double)down_sampling_factor)*total_prod_imag/((double)periodo_interval_estim/2.0);

    	//if(contador_iteracoes>1)
		//{
		//Transformação DQ 	
		//aux3= fase_estimada[ant]+freq_estimada*pipi/Fs;
		//saida_real[atual] = cos(aux3)*total_prod_real - sin(aux3)*total_prod_imag;
		//saida_imag[atual] = sin(aux3)*total_prod_real + cos(aux3)*total_prod_imag;
		//saida_real_atual = cos(fase_estimada_anterior + pipi*freq_estimada/Fs)*total_prod_real - sin(fase_estimada_anterior+ pipi*freq_estimada/Fs)*total_prod_imag;
		//saida_imag_atual = sin(fase_estimada_anterior + pipi*freq_estimada/Fs)*total_prod_real + cos(fase_estimada_anterior+ pipi*freq_estimada/Fs)*total_prod_imag;

		//Filtragem das componentes DQ 

		//aux_real = (double)(0.0063)*(saida_real[atual] + saida_real[ant]);
		//aux_imag = (double)(0.0063)*(saida_imag[atual] + saida_imag[ant]);
		//aux_real1 =(((double)0.9877)*((double)(r_filt[ant])));
		//aux_imag1 = (((double)0.9877)*((double)(i_filt[ant])));
		//r_filt[atual] = (double)((double)(aux_real) + (double)(aux_real1));
		//i_filt[atual] = (double)((double)(aux_imag) + (double)(aux_imag1));
		//aux_real = (((double)0.9877)*((double)(r_filt))) + (((double)0.0063)*((double)(aux_real)));
	    //aux_imag = (((double)0.9877)*((double)(i_filt))) + (((double)0.0063)*((double)(aux_imag)));

	   		/******* a Volta *******/

		//total_prod_real = (double)(r_filt[atual])*cos(aux3) + (double)(i_filt[atual])*sin(aux3);
		//total_prod_imag = -((double)(r_filt[atual])*sin(aux3))  + (double)(i_filt[atual])*cos(aux3);

		//fase_estimada[atual] = atan2(-imag_volta,real_volta);

		fase_estimada[atual] = atan2(-total_prod_imag,total_prod_real);
		
			/****************** derivada *********************/

		//}
		//else
		//{
			//fase_estimada[atual] = atan2(-total_prod_imag,total_prod_real);
		//}
		
		diff_phi[atual] = (fase_estimada[atual]-fase_estimada[ant]);
		
		//Este calculo de derivada pode ser melhorado
		if ( ((double) fabs(diff_phi[atual])) > pi )	//Correcao do problema decorrente da descontinuidade da fase
		{
			diff_phi[atual] = diff_phi[ant];
		}
		diff_phi_corrigida[atual] = diff_phi[atual] - aux;
		
			/***** Gerando a senoide sincronizada *****/

		produto_interno = (total_prod_real*cos(fase_estimada[atual]) - total_prod_imag*sin(fase_estimada[atual]));
		//produto_interno = (real_volta*cos(fase_estimada[atual]) - imag_volta*sin(fase_estimada[atual]));
		senoide_sinc = cos(fase_estimada[atual]);	
		
		 if(contador_iteracoes > periodo_interval_estim)
		{
			/***** Variacao da fundamental: Integral da derivada da fase menos 2*pi*Ff *****/
			
			Ff_var = (diff_phi_corrigida[atual]+diff_phi_corrigida[ant])/2.0;
			//     		Ff_var = (diff_phi_corrigida_atual+diff_phi_corrigida_anterior+diff_phi_corrigida_anterior2)/3.0;   // nova
			freq_estimada = (freq_estimada+ganho_laco*Ff_var);
		    periodo_interval_estim = (int)(FreqAmostragem/(freq_estimada));
	  
		   // printf(" %i \n",periodo_interval_estim);
           
			/* Atualizacao do periodo de integracao do correlator: Quantidade de amostras existentes em 1 ciclo da fundamental com a freq. Estimada */
		}	
		
			/***** Saida dos parametros estimados *****/
		
		V_estim   [contador_iteracoes - 1] = (double)(senoide_sinc*produto_interno);	/* Fundamental estimada */
		freq_estim[contador_iteracoes - 1] = (double)(freq_estimada);					/* Frequencia da fundamental */
		amp_estim [contador_iteracoes - 1] = (double)(produto_interno);					/* Amplitude da fundamental */
        wt_estim  [contador_iteracoes - 1] = (double)(fase_estimada[atual]);			/* Fase da fundamental */
			
			/***** Atualizacao do contador *****/	
		
		//if(contador_iteracoes < (comprimento_mem_sinal+1))
		//{
		//	contador_iteracoes++;
		//}

		pos++; // 0 ate 7318
		//pos = pos % comprimento_mem_sinal;
		atual++;
		ant++;
		atual = atual % 2;
		ant = ant % 2;
		//if(pos==(comprimento_mem_sinal))
		//pos = 0;
	}
	return;
}