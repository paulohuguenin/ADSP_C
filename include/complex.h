/*****************************************************************************

   Projeto : SINAPE - Interface do Banco de oscilografia
   Nome    : CComplex.h
   Funcao  : Implementa opera��es com n�meros complexos
   C       : C++
   Autores : Andr� Luiz Lins Miranda
   Data    : 05/2002
   Obs.	   : Modifica��o da classe Complex para melhorar a
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
	// Se ret>0 val1 e val2 est�o em coordenadas retangulares. sen�o eles
	// est�o em coordenadas polares. val2 deve estar em radianos
    CComplex(double val1, double val2=0.0, int ret=1);

    // manipula��es complexas
    double Real(void) const;   // Parte real
    double Imag(void) const;   // Parte imagin�ria
	void setReal(double);   // Parte real
    void setImag(double) ;   // Parte imagin�ria
    CComplex Conj(void) const;  // Complexo conjugado
	double Abs(void) const;	// Magnitude
    double Arg(int rad=1) const;    // �ngulo no plano. Se rad>0 retorna em radiano
									// Sen�o, retorna em graus
	double ArgWrap(int rad=1) const;	// �ngulo no intervalo entre 0 e pi (rad>0) ou 0 e 180 graus

		// Fun��es Matem�ticas
	int Finito(void) const;		// >0 se infinito
	CComplex Sqrt(void) const;	// raiz quadrada
	CComplex Log(void) const;	// Logaritmo neperiano
	CComplex LogQ(double Q) const;	// Logaritmo na base Q
	CComplex Exp(void) const;		// exponencial

		// Funcoes trigonom�tricas
	#define CPLXvezesJ(complexo) CComplex(-complexo.Imag(),complexo.Real()) // produto de um complexo por j
	CComplex Sin(void) const;	// seno
	CComplex Cos(void) const;	// coseno
	CComplex Tan(void) const;	// tangente
	CComplex Cotan(void) const;	// cotangente

		// Fun��es hiperb�licas
	CComplex Cosh(void) const;	// Cosseno Hiperb�lico
	CComplex aCosh(void) const;
	CComplex Sinh(void) const;	// Seno Hiperb�lico
	CComplex aSinh(void) const;
	CComplex Tanh(void) const;	// Tangente Hiperb�lica
	CComplex aTanh(void) const;	// arc-Tangente Hiperb�lica
	CComplex CoTanh(void) const;// Cotangente Hiperb�lica
	CComplex aCoTanh(void) const;//arc-Cotangente Hiperb�lica
	CComplex CoSec(void) const;// Cossecante Hiperb�lica
	CComplex aCoSec(void) const;//arc-Cossecante Hiperb�lica
	CComplex Sec(void) const;// Secante Hiperb�lica
	CComplex aSec(void) const;//arc-secante Hiperb�lica

    // Operadores Bin�rios
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
	CComplex operator^(const double& c) const;		// exponencia��o

	const CComplex& operator*=(const double real);
	const CComplex& operator/=(const double real);
	const CComplex& operator*=(const CComplex& c);
	const CComplex& operator/=(const CComplex& c);

	const CComplex& operator=(const CComplex& c); // Operador de atribui��o

	int operator==(CComplex c) const;	// Verdadeiro se returna 1, falso se 0
	int operator!=(CComplex c) const;

	// Faz a rota��o do complexos em graus
	CComplex operator>>(const double gr) const;
	CComplex operator<<(const double gr) const;
	const CComplex& operator>>=(const double gr);
	const CComplex& operator<<=(const double gr);


};

#endif // !defined(AFX_CCOMPLEX_H__A26F9429_A41C_46F8_AB84_59F01174339C__INCLUDED_)
