// parameter.cpp

#include "parameter.h"


void CExpParm::setType(int type)
{
	m_type = type;
}

void CExpParm::setParm(strtContinuousExp exp_parm)
{
    innerProd = exp_parm.innerProduct;
	rho = exp_parm.rho;
	xi = exp_parm.xi;
	phase = exp_parm.phase;
	a = exp_parm.a;
	b = exp_parm.b;

}


int CExpParm::getType()
{
	return m_type;
}

void CExpParm::printParm2Screen()
{
//    printf("Coef: %10.8f | Rho: %10.8f| Xi: %10.8f| Phase: %10.8f| a: %8i| b: %8i\n", innerProd,
 //                                                                                   rho,
 //                                                                                   xi,
 //                                                                                   phase,
  //                                                                                  a,
   //                                                                                 b);
    cout << "AQUI" << endl;
    cout << "Coef: "<<innerProd<<"| Rho: "<<rho<<"| Xi: "<<xi<<"| Phase: "<<phase<<"| a: "<<a<<"| b: "<<b<< endl; 

/*    cout<< "Inner Product: " << innerProd << endl;
	cout<< "Rho: " << rho << endl;
	cout<< "Xi: " << xi << endl;
	cout<< "Phase: " << phase << endl;
	cout<< "a: " << a << endl;
	cout<< "b: " << b << endl;*/
}
