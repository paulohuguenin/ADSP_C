/*****************************************************************************

   Projeto : SINAPE - Interface do Banco de oscilografia
   Nome    : CComplex.h
   Funcao  : Implementa operações com números complexos
   C       : C++
   Autores : André Luiz Lins Miranda
   Data    : 05/2002
   Obs.	   : Modificação da classe Complex para melhorar a
			 sua funcionalidade

*****************************************************************************/

#if !defined(AFX_CCOMPLEX_H__A26F9429_A41C_46F8_AB84_59F01174339C__INCLUDED_)
#define AFX_CCOMPLEX_H__A26F9429_A41C_46F8_AB84_59F01174339C__INCLUDED_

#include "defines.h"

/**  \brief Complex numeric manipulation class
 *
 *   Complex numeric manipulation class
 */
class CComplex
{
// Implementation
private:
        double re;
		double im;
        double mod; 
		double ang;

public:
    // constructors
    CComplex();
	// Se ret>0 val1 e val2 estão em coordenadas retangulares. senão eles
	// estão em coordenadas polares. val2 deve estar em radianos
    CComplex(double val1, double val2=0.0, int ret=1);

    // manipulações complexas
    double Real(void) const;   // Parte real
    double Imag(void) const;   // Parte imaginária
	void setReal(double);   // Parte real
    void setImag(double) ;   // Parte imaginária
    CComplex Conj(void) const;  // Complexo conjugado
	double Abs(void) const;	// Magnitude
    double Arg(int rad=1) const;    // Ângulo no plano. Se rad>0 retorna em radiano
									// Senão, retorna em graus
	double ArgWrap(int rad=1) const;	// Ângulo no intervalo entre 0 e pi (rad>0) ou 0 e 180 graus

		// Funções Matemáticas
	int Finito(void) const;		// >0 se infinito
	CComplex Sqrt(void) const;	// raiz quadrada
	CComplex Log(void) const;	// Logaritmo neperiano
	CComplex LogQ(double Q) const;	// Logaritmo na base Q
	CComplex Exp(void) const;		// exponencial

		// Funcoes trigonométricas
	#define CPLXvezesJ(complexo) CComplex(-complexo.Imag(),complexo.Real()) // produto de um complexo por j
	CComplex Sin(void) const;	// seno
	CComplex Cos(void) const;	// coseno
	CComplex Tan(void) const;	// tangente
	CComplex Cotan(void) const;	// cotangente

		// Funções hiperbólicas
	CComplex Cosh(void) const;	// Cosseno Hiperbólico
	CComplex aCosh(void) const;
	CComplex Sinh(void) const;	// Seno Hiperbólico
	CComplex aSinh(void) const;
	CComplex Tanh(void) const;	// Tangente Hiperbólica
	CComplex aTanh(void) const;	// arc-Tangente Hiperbólica
	CComplex CoTanh(void) const;// Cotangente Hiperbólica
	CComplex aCoTanh(void) const;//arc-Cotangente Hiperbólica
	CComplex CoSec(void) const;// Cossecante Hiperbólica
	CComplex aCoSec(void) const;//arc-Cossecante Hiperbólica
	CComplex Sec(void) const;// Secante Hiperbólica
	CComplex aSec(void) const;//arc-secante Hiperbólica

    // Operadores Binários
	CComplex operator+() const;
    CComplex operator-() const;
         
	CComplex operator+(const double real) const;
	CComplex operator-(const double real) const;
	CComplex operator+(const CComplex& c) const;
	CComplex operator-(const CComplex& c) const;
	
	const CComplex& operator+=(const double real);
	const CComplex& operator-=(const double real);
	const CComplex& operator+=(const CComplex& c);
	const CComplex& operator-=(const CComplex& c);

	CComplex operator*(const double real) const;
	CComplex operator/(const double real) const;
	CComplex operator*(const CComplex& c) const;
	CComplex operator/(const CComplex& c) const;
	CComplex operator^(const double& c) const;		// exponenciação

	const CComplex& operator*=(const double real);
	const CComplex& operator/=(const double real);
	const CComplex& operator*=(const CComplex& c);
	const CComplex& operator/=(const CComplex& c);

	const CComplex& operator=(const CComplex& c); // Operador de atribuição

	int operator==(CComplex c) const;	// Verdadeiro se returna 1, falso se 0
	int operator!=(CComplex c) const;

	// Faz a rotação do complexos em graus
	CComplex operator>>(const double gr) const;
	CComplex operator<<(const double gr) const;
	const CComplex& operator>>=(const double gr);
	const CComplex& operator<<=(const double gr);


};

#endif // !defined(AFX_CCOMPLEX_H__A26F9429_A41C_46F8_AB84_59F01174339C__INCLUDED_)
