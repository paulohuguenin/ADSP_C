/*****************************************************************************

   Projeto : SINAPE - Interface do Banco de oscilografia
   Nome    : CComplex.cpp
   Funcao  : Implementa operações com números complexos
   C       : C++
   Autores : André Luiz Lins Miranda
   Data    : 05/2002
   Obs.	   : Modificação da classe Complex para melhorar a
			 sua funcionalidade

*****************************************************************************/

#include "complex.h"
#include <math.h>
#include <float.h>

CComplex::CComplex()
{
    re = im = 0.0;
    mod = ang = 0.0;
}

CComplex::CComplex(double val1, double val2, int ret)
{
//	if(ret>0) {
//		re = val1;
//		im = val2;
//	} else {
//		re = val1*cos(val2);
//		im = val1*sin(val2);
//	}
	if(ret>0) 
	{
		re = val1;
		im = val2;
		mod = sqrt(val1*val1 + val2*val2);
		ang = atan2(val2, val1);
	} 
	else 
	{
		mod = val1;
		ang = val2;
		re = val1*cos(val2);
		im = val1*sin(val2);
	}

}

double CComplex::Real(void) const
{
	return re;
//	return mod*cos(ang);
}

double CComplex::Imag(void) const
{
	return im;
//	return mod*sin(ang);
}

void CComplex::setReal(double real_part)
{
	re = real_part;
}

void CComplex::setImag(double imag_part)
{
	im = imag_part;
} 

CComplex CComplex::Conj(void) const
{
//  return CComplex(re, -im);
    return CComplex(this->Real(),-this->Imag());
}

double CComplex::Abs(void) const
{
//	return sqrt(re*re + im*im);
	return mod;
}

double CComplex::Arg(int rad) const
{
	double angulo = (rad>0) ? ang : 180.0*ang/atan2(0.0,-1.0);
//	if(rad>0) angulo = atan2(im, re);
//	else angulo = (180.0*atan2(im, re))/atan2(0,-1);
	return angulo;
}

double CComplex::ArgWrap(int rad) const
{
	double X = this->Real();
	double Y = this->Imag();
	double angulo;
	if(rad>0) angulo = atan2(Y, X);
	else angulo = (180.0*atan2(Y, X))/atan2(0.0,-1.0);
	return angulo;
}

int CComplex::Finito(void) const
{
#ifndef WIN
	return (finite(this->Real())==1) & (finite(this->Imag())==1);
#else
	return (_finite(this->Real())==1) & (_finite(this->Imag())==1);
#endif

}

CComplex CComplex::Sqrt(void) const
{
	double m1,a1;
	m1 = this->Abs();
	a1 = this->Arg();
    return CComplex(sqrt(m1), a1/2.0,0);
}

CComplex CComplex::Log(void) const
{
	double a, b;
	a = this->Abs();
	a = a<=0 ? DBL_EPSILON*100 : a;
	b = this->Arg();
	return CComplex(log(a),b);
}

CComplex CComplex::LogQ(double Q) const
{
	Q = Q<=0 ? DBL_EPSILON*100 : Q;
	return this->Log() / log(Q);
}

CComplex CComplex::Exp(void) const
{
	return CComplex(exp(this->Real()),this->Imag(),0);
}

		// Funcoes trigonométricas
CComplex CComplex::Sin(void) const
{
	CComplex x = *this;
	CComplex y = - CPLXvezesJ(x);
	x = CPLXvezesJ(x);
	CComplex doisjota = CComplex(0,2);
	return (x.Exp() - y.Exp())/doisjota;
}

CComplex CComplex::Cos(void) const
{
	CComplex x = *this;
	CComplex y = - CPLXvezesJ(x);
	x = CPLXvezesJ(x);
	return (x.Exp() + y.Exp())/2.0;
}

CComplex CComplex::Tan(void) const
{
	CComplex c = this->Cos();
	c = c.Abs() < DBL_EPSILON ?  DBL_EPSILON : c;
	return this->Sin() / c;
}

CComplex CComplex::Cotan(void) const
{
	CComplex s = this->Sin();
	s = s.Abs() < DBL_EPSILON ?  DBL_EPSILON : s;
	return this->Cos() / s;
}

		// Funções hiperbólicas
CComplex CComplex::Cosh(void) const
{
	double a, b;
	a = cosh(this->Real())*cos(this->Imag());
	b = sinh(this->Real())*sin(this->Imag());
	return CComplex(a,b);
}

CComplex CComplex::aCosh(void) const
{
	CComplex x = *this;
	x = x*x - 1;
	x = x.Sqrt();
	x += *this;
	return x.Log();
}

CComplex CComplex::Sinh(void) const
{
	double a, b;
	a = sinh(this->Real())*cos(this->Imag());
	b = cosh(this->Real())*sin(this->Imag());
	return CComplex(a,b);
}

CComplex CComplex::aSinh(void) const
{
	CComplex x = *this;
	x = x*x + 1;
	x = x.Sqrt();
	x += *this;
	return x.Log();
}

CComplex CComplex::Tanh(void) const
{
	CComplex a, b;
	a = this->Sinh();
	b = this->Cosh();
	if (!b.Finito())
		if (!a.Finito()) return 1;
			else	return 0;
	return a/b;
}

CComplex CComplex::aTanh(void) const
{
	CComplex a = (*this + 1) / (CComplex(1,0) - *this);
	return a.Log()/2;
}

CComplex CComplex::CoTanh(void) const
{
	CComplex a, b;
	a = this->Cosh();
	b = this->Sinh();
	if (!b.Finito())
		if (!a.Finito()) return 1;
			else	return 0;
	return a/b;
}

CComplex CComplex::aCoTanh(void) const
{
	CComplex a = (*this + 1) / (*this - 1);
	return a.Log()/2;
}

CComplex CComplex::CoSec(void) const{
	CComplex b;
	b = this->Sinh();
	if (!b.Finito()) return 0;
	return CComplex(1)/b;

};

CComplex CComplex::aCoSec(void) const{
	CComplex x = *this;
	x = (CComplex(1) / (x*x)) + 1;
	x = (CComplex(1) / *this) + x.Sqrt();
	return x.Log();
};

CComplex CComplex::Sec(void) const{
	CComplex b;
	b = this->Cosh();
	if (!b.Finito()) return 0;
	return CComplex(1)/b;

};

CComplex CComplex::aSec(void) const{
	CComplex x = *this;
	x = (CComplex(1) / (x*x)) - 1;
	x = (CComplex(1) / *this) + x.Sqrt();
	return x.Log();
};

CComplex  CComplex::operator+() const
{return *this;}

CComplex  CComplex::operator-() const
{return CComplex(-this->Real(), -this->Imag());}

CComplex CComplex::operator+(const double real) const
{return CComplex(this->Real()+real,this->Imag());}

CComplex CComplex::operator-(const double real) const
{return CComplex(this->Real()-real,this->Imag());}

CComplex CComplex::operator+(const CComplex& c) const
{return	CComplex(this->Real()+c.Real(),this->Imag()+c.Imag());}

CComplex CComplex::operator-(const CComplex& c) const
{return	CComplex(this->Real()-c.Real(),this->Imag()-c.Imag());}

const CComplex& CComplex::operator+=(const double real)
{
//	re += real;
	double X = this->Real() + real;
	double Y = this->Imag();
	mod = sqrt(X*X + Y*Y);
	ang = atan2(Y, X);
	return *this;
}

const CComplex& CComplex::operator-=(const double real)
{
//	re -= real;
	double X = this->Real() - real;
	double Y = this->Imag();
	mod = sqrt(X*X + Y*Y);
	ang = atan2(Y, X);
	return *this;
}

const CComplex& CComplex::operator+=(const CComplex& c)
{
//	re += c.Real();
//	im += c.Imag();
	double X = this->Real() + c.Real();
	double Y = this->Imag() + c.Imag();
	mod = sqrt(X*X + Y*Y);
	ang = atan2(Y, X);
	return *this;
}

const CComplex& CComplex::operator-=(const CComplex& c)
{
//	re -= c.Real();
//	im -= c.Imag();
	double X = this->Real() - c.Real();
	double Y = this->Imag() - c.Imag();
	mod = sqrt(X*X + Y*Y);
	ang = atan2(Y, X);
	return *this;
}

//CComplex CComplex::operator*(const double real) const
//{return CComplex(re*real,im*real);}
CComplex CComplex::operator*(const double real) const
{
	double R = ( this->Real() )*real;
	double I = ( this->Imag() )*real;
	return CComplex(R,I);
}

CComplex CComplex::operator/(const double real) const
{
	double d1,d2;
	if(fabs(real)<1e-12) {
		d1 = d2 = 0.0;
	} else {
		d1 = this->Real()/real;
		d2 = this->Imag()/real;
	}
	return CComplex(d1,d2);
}

CComplex CComplex::operator*(const CComplex& c) const
{
	double d1,d2;
	d1 = this->Real()*c.Real() - this->Imag()*c.Imag();
	d2 = this->Real()*c.Imag() + this->Imag()*c.Real();
    return CComplex(d1,d2);
}
	
CComplex CComplex::operator/(const CComplex& c) const
{
	double aux_re,aux_im, aux_den;
	aux_den = c.Real()*c.Real() + c.Imag()*c.Imag();
	if(fabs(aux_den)<1e-12) {
		aux_re = aux_im = 0.0;aux_den=1;
	} else {
		aux_re = (this->Real()*c.Real() + this->Imag()*c.Imag());
		aux_im = (this->Imag()*c.Real() - this->Real()*c.Imag());
	}
#ifndef WIN
	if (!finite(aux_den)) {
		if (!finite(aux_re) || !finite(aux_im))
					return 1;
			else	return 0;
	};
#else
	if (!_finite(aux_den)) {
		if (!_finite(aux_re) || !_finite(aux_im))
					return 1;
			else	return 0;
	};
#endif	

	return CComplex(aux_re/aux_den, aux_im/aux_den);
}

CComplex CComplex::operator^(const double& c) const
{
	double a,b;
	a = this->Abs();
	b = this->Arg();
	return CComplex(exp(c*log(a)),b*c,0);
}

const CComplex& CComplex::operator*=(const double real)
{
//	re *= real;
//	im *= real;
	mod *= real;
	return *this;
}
	
const CComplex& CComplex::operator/=(const double real)
{
//	if(fabs(real)<1e-12) {
//		re = im = 0.0;
//	} else {
//		re /= real;
//		im /= real;
//	}

	if(fabs(real)<1e-12) {
		mod = ang = 0.0;
	} else {
		mod /= real;
	}
	return *this;
}
	
const CComplex& CComplex::operator*=(const CComplex& c)
{
//	double d1,d2;
//	d1 = re*c.Real() - im*c.Imag();
//	d2 = re*c.Imag() + im*c.Real();
//	re = d1;
//	im = d2;
	mod *= c.Abs();
	ang += c.Arg();
	return *this;
}

const CComplex& CComplex::operator/=(const CComplex& c)
{
//	double aux_re,aux_im, aux_den;
//	aux_den = c.Real()*c.Real() + c.Imag()*c.Imag();
//	if(fabs(aux_den)<1e-12) {
//		aux_re = aux_im = 0.0;
//	} else {
//		aux_re = (re*c.Real() + im*c.Imag())/aux_den;
//		aux_im = (im*c.Real() - re*c.Imag())/aux_den;
//	}
//	re = aux_re;
//	im = aux_im;
	mod /= c.Abs();
	ang -= c.Arg();
	return *this;
}

const CComplex& CComplex::operator=(const CComplex& c) // Operador de atribuição
{
	this->re = c.Real();
	this->im = c.Imag();
	this->mod = c.Abs();
	this->ang = c.Arg();
	return *this;
}

int CComplex::operator==(const CComplex c) const
//{return (re == c.Real())&&(im == c.Imag());}
{return (this->Real() == c.Real())&&(this->Imag() == c.Imag());}

int CComplex::operator!=(const CComplex c) const
//{return (re != c.Real())||(im != c.Imag());}
{return (this->Real() != c.Real())||(this->Imag() != c.Imag());}

CComplex CComplex::operator>>(const double gr) const
{
//	double th1,th2;
//	th1 = atan2(im, re);
//	th2 = atan2(0,-1)*gr/180.0;
//	return CComplex(re*cos(th1+th2)/cos(th1),im*sin(th1+th2)/sin(th1));
	double R = this->Abs();
	double F = this->Arg() + atan2(0.0,-1.0)*gr/180.0;
	return CComplex(R,F,0);
}

CComplex CComplex::operator<<(const double gr) const
{
//	double th1,th2;
//	th1 = atan2(im, re);
//	th2 = atan2(0,-1)*gr/180.0;
//	return CComplex(re*cos(th1-th2)/cos(th1),im*sin(th1-th2)/sin(th1));
	double R = this->Abs();
	double F = this->Arg() - atan2(0.0,-1.0)*gr/180.0;
	return CComplex(R,F,0);
}

const CComplex& CComplex::operator>>=(const double gr)
{
//	double th1,th2;
//	th1 = this->Arg();
//	th2 = atan2(0,-1)*gr/180.0;
//	re *= cos(th1+th2)/cos(th1);
//	im *= sin(th1+th2)/sin(th1);
	ang += atan2(0.0,-1.0)*gr/180.0;
	return *this;
}

const CComplex& CComplex::operator<<=(const double gr)
{
//	double th1,th2;
//	th1 = this->Arg();
//	th2 = atan2(0,-1)*gr/180.0;
//	re *= cos(th1-th2)/cos(th1);
//	im *= sin(th1-th2)/sin(th1);
	ang -= atan2(0.0,-1.0)*gr/180.0;
	return *this;
}

