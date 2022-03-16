
#include "adspsort.h"

void SWAPINT(int& a,int& b)   
{ 
	int t; 
	t=a; 
	a=b; 
	b=t; 
}  // Macro for swapping


void bubble_srtint( int* a, int n )  
{   
    int i, j;
       
    for(i = 0; i < n; i++)         // Make a pass through the array for each element
    {              
        for(j = 1; j < (n-i); j++) // Go through the array beginning to end
        {              
           if(a[j-1] > a[j])       // If the the first number is greater, swap it 
              SWAPINT(a[j-1],a[j]);   
        }
    }
}

void SWAPDOUBLE(double& a,double& b)   
{ 
	double t; 
	t=a; 
	a=b; 
	b=t; 
}  // Macro for swapping

void bubble_srtdouble( double* a, int n )  
{   
    int i, j;
       
    for(i = 0; i < n; i++)         // Make a pass through the array for each element
    {              
        for(j = 1; j < (n-i); j++) // Go through the array beginning to end
        {              
           if(a[j-1] > a[j])       // If the the first number is greater, swap it 
              SWAPDOUBLE(a[j-1],a[j]);   
        }
    }
}

void bubble_srtdouble_noduplicate( double* a, int n, int& n_aux )  
{   
    int i, j;
    double* aux;
    
    aux =new double[n];
       
    for(i = 0; i < n; i++)         // Make a pass through the array for each element
    {              
        for(j = 1; j < (n-i); j++) // Go through the array beginning to end
        {              
           if(a[j-1] > a[j])       // If the the first number is greater, swap it 
              SWAPDOUBLE(a[j-1],a[j]);   
        }
    }
    
    
    aux[0] = a[0];
    n_aux=1;
    for(i = 1; i < n; i++)         // Make a pass through the array for each element
    {
    	if(a[i-1] != a[i])
    	{
    		aux[n_aux] = a[i];
    		n_aux++;
    	}
    }
    
    memcpy(a,aux,sizeof(double)*n_aux);   
    
    delete [] aux;
    
}
