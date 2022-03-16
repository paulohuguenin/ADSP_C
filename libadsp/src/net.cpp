//net.cpp

#include "net.h"
// ========================================
//  CNet Class
// ========================================

//===============================================================
// Function: Contructors
// Goal:
// Return:
//===============================================================
CNet::CNet()
{
	int m_signalSize = 0;
	int m_num = 0;
	m_x = NULL;
	m_y = NULL;
}


//===============================================================
// Function: Destructor
// Goal:
// Return:
//===============================================================

CNet::~CNet()
{
	if (m_x) delete[] m_x;
	if (m_y) delete[] m_y;
}

void CNet::setSignalSize(int signalSize)
{
	m_signalSize = signalSize;
    if (m_x) delete [] m_x;
    if (m_y) delete [] m_y;
    m_x = new double[m_signalSize/8];
    m_y = new double[3];
}

void CNet::setX(double* x)
{
	m_x = x;
}

void CNet::setNumber(int num)
{
	m_num = num;
}

double* CNet::getY()
{
	return m_y;
}

// double* CNet::mapMinMaxFunc(double* x,
// 							double* x1_xoffset,
// 							double* x1_gain,
// 							double x1_ymin)
// {
// 	for (int i=0; i<m_signalSize; i++)
// 	{
// 		x[i] = m_x[i] - x1_xoffset[i];
// 		x[i] *= x1_gain[i];
// 		x[i] += x1_ymin;
// 	}	

// 	return x;
// }

// void CNet::hiddenFunc(	double* x1_xoffset,
// 						double* x1_gain,
// 						double x1_ymin,
// 						double* b,
// 						double w[][16])
// {
// 	double* x = new double[m_signalSize/8];

// 	x = mapMinMaxFunc(m_x, x1_xoffset, x1_gain, x1_ymin);


// 	for (int i=0; i<m_signalSize; i++)
// 	{
// 		m_h[i] = ( 2 / ( 1 + exp( -2*x[i]*w[i][0] ) ) )- 1;
// 	}
// }

// void CNet::outputFunc(	double b[],
// 						double w[][3])
// {
// 	double max_h = -1e18;
// 	double* num = new double[m_signalSize/32];
// 	double den = 1;

// 	for (int i=0; i<m_signalSize; i++)
// 	{
// 		num[i] = 0.0;
// 	}

// 	for (int i=0; i<m_signalSize; i++)
// 	{
// 		if (max_h > m_h[i]) max_h = m_h[i];
// 	}

// 	for (int i=0; i<m_signalSize; i++)
// 	{
// 		num[i] = exp( m_h[i] );
// 	}

// 	for (int i=0; i<m_signalSize; i++)
// 	{
// 		den += num[i];
// 	}

// 	if (den == 0.0) den = 1.0;

// 	for (int i=0; i<m_signalSize; i++)
// 	{
// 		m_y[i] = num[i]/den;
// 	}
// }

void CNet::net()
{
	// Loading the network file
    CFileNetwork* ann;
    ann = new CFileNetwork;
    ann->setFileName("panelNetwork.dat");
    ann->setNetwork(m_num);
    ann->loadData();

	int size_x = 64;
	int size_h = 16;
	int size_y = 3;

    double* h;
    h = new double[size_h];

	double* x;
	x = new double[size_x];

	double max_h;

	double* num;
	num = new double[size_h];

	double den = 0.0;

	double* w1x;
	w1x = new double[size_h];

	double* w2h;
	w2h = new double[size_y];

	double* w1xb1;
	w1xb1 = new double[size_h];

	double* w2hb2;
	w2hb2 = new double[size_y];

	double* x1_xoffset;
	// x1_xoffset = new double[size_x];

    double* x1_gain;
    // x1_gain = new double[size_x];

    double x1_ymin = 0;

    double* b1;
    // b1 = new double[size_h];

    double ** w1;
    // w1 = new double*[size_h];

    // for(int i = 0; i < size_h; ++i)
    // {
    //     w1[i] = new double[size_x];
    // }

    double* b2;
    // b2 = new double[size_y];

    double ** w2;
    // w2 = new double*[size_y];

    // for(int i = 0; i < size_y; ++i)
    // {
    //     w2[i] = new double[size_h];
    // }

    // cout << ann->getX1_xoffset()<< endl;

    // cout << x1_xoffset[0] << endl;
    

 //    memcpy(x1_xoffset, ann->getX1_xoffset(),size_x);
	// memcpy(x1_gain, ann->getX1_gain(), size_x);
	// memcpy(b1, ann->getB1(), size_h);
	// memcpy(w1, ann->getW1(), size_h);
	// memcpy(&b2, ann->getB2(), size_y);
	// memcpy(&w2, ann->getW2(), size_y);
	x1_xoffset = ann->getX1_xoffset();
	x1_gain = ann->getX1_gain();
	x1_ymin = ann->getX1_ymin();

	b1 = ann->getB1();
	w1 = ann->getW1();
	b2 = ann->getB2();
	w2 = ann->getW2();


	// cout << x1_gain[0] << endl;
	// cout << x1_ymin << endl;

	for (int i=0; i<size_h; i++)
	{
		w1x[i] = 0.0;
	}

	for (int i=0; i<size_y; i++)
	{
		w2h[i] = 0.0;
	}

	// minmax
	for (int i=0; i<size_x; i++)
	{
		// cout << m_x[i] << endl;
		x[i] = m_x[i] - x1_xoffset[i];
		x[i] *= x1_gain[i];
		x[i] += x1_ymin;
		// cout << x[i] << endl;
	}

	// innerprod1
	for (int i=0; i<size_h; i++)
	{
		for (int j=0; j<size_x; j++)
		{
			w1x[i] += w1[i][j]*x[j];
			// cout << w1[i][j] << endl;
			// cout << x[j] << endl;
		}
		// cout << w1x[i] << endl;
	}

	for (int i=0; i<size_h; i++)
	{
		w1xb1[i] = w1x[i] + b1[i];
		// cout << b1[i] << endl;			
	}

	// cout << w1x[0] << endl;

	// sigmoid
	for (int i=0; i<size_h; i++)
	{
		h[i] = ( 2 / ( 1 + exp( -2 * (w1xb1[i]) ) ) )- 1;
		// cout << exp( -2 * (w1xb1[i]) )  << endl;
		// cout << h[i] << endl;
	}

	// innerprod2
	for (int i=0; i<size_y; i++)
	{
		for (int j=0; j<size_h; j++)
		{
			w2h[i] += w2[i][j]*h[j];
		}
		// cout << w2h[i] << endl;
	}

	for (int i=0; i<size_y; i++)
	{
		w2hb2[i] = w2h[i] + b2[i];
		// cout << w2h[i] << endl;
	}

	// softmax
	for (int i=0; i<size_y; i++)
	{
		num[i] = 0.0;
	}

	max_h = w2hb2[0];

	for (int i=1; i<size_y; i++)
	{
		if ( w2hb2[i] > max_h ) max_h = w2hb2[i];
	}


	// cout << max_h << endl;

	for (int i=0; i<size_y; i++)
	{
		num[i] = exp( w2hb2[i] - max_h);
	}

	for (int i=0; i<size_y; i++)
	{
		den += num[i];
	}
	if (den == 0.0) den = 1.0;
// 		

	for (int i=0; i<size_y; i++)
	{
		m_y[i] = num[i]/den;
		// cout << m_y[i] << endl;
	}
	
	if (h != NULL) delete [] h;
	if (x != NULL) delete [] x;
	if (num != NULL) delete [] num;
	if (w1x != NULL) delete [] w1x;
	if (w2h != NULL) delete [] w2h;
	if (w1xb1 != NULL) delete [] w1xb1;
	if (w2hb2 != NULL) delete [] w2hb2;

	// if (x1_xoffset != NULL) delete [] x1_xoffset;

 //    if (x1_gain != NULL) delete [] x1_gain;

 //    if (b1 != NULL) delete [] b1;

 //    for(int i = 0; i < size_h; ++i)
 //    {    
 //        if (w1[i] != NULL) delete [] w1[i];
 //    }

 //    if (w1 != NULL) delete [] w1;

 //    if (b2 != NULL) delete [] b2;

 //    for(int i = 0; i < size_y; ++i)
 //    {    
 //        if (w2[i] != NULL) delete [] w2[i];
 //    }

 //    if (w2 != NULL) delete [] w2;

	if (ann != NULL) delete ann;

}