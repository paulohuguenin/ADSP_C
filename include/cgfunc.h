// CGFUNC.H
#ifndef CGFUNC_H
#define CGFUNC_H

#include <iostream>
#include "cgmatrix.h"

using namespace std;

/*_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_*/
/*			FUNCTIONS PROTOTYPES				     */
/*_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_*/

// Point-to-point matrix multiplication
template <class MatrixType>
cgMatrix<MatrixType>		ppMult(cgMatrix<MatrixType>, cgMatrix<MatrixType>);

// Point-to-point matrix division
template <class MatrixType>
cgMatrix<MatrixType>		ppDiv(cgMatrix<MatrixType>, cgMatrix<MatrixType>);

// Sequences Convolution
template <class MatrixType>
cgMatrix<MatrixType>		Conv(cgMatrix<MatrixType>, cgMatrix<MatrixType>);

// Vertical concatenation
template <class MatrixType>
cgMatrix<MatrixType>		VertConc(cgMatrix<MatrixType>, cgMatrix<MatrixType>);

// Horizontal concatenation
template <class MatrixType>
cgMatrix<MatrixType>		HorizConc(cgMatrix<MatrixType>, cgMatrix<MatrixType>);

// LU decomposition
template <class MatrixType>
void				LUDecomp(cgMatrix<MatrixType>, cgMatrix<MatrixType> &, cgMatrix<MatrixType> &);

// LDU decomposition
template <class MatrixType>
void				LDUDecomp(cgMatrix<MatrixType>, cgMatrix<MatrixType> &, cgMatrix<MatrixType> &, cgMatrix<MatrixType> &);

// Matrix Inversion
template <class MatrixType>
cgMatrix<MatrixType>		Inv(cgMatrix<MatrixType> );

// Sequence FFT
template <class MatrixType>
cgMatrix<MatrixType>		FFT(cgMatrix<MatrixType>);


/*_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_*/
/*			FUNCTIONS DEFINITIONS				     */
/*_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_*/

// Point-to-point matrix multiplication
template <class MatrixType>
cgMatrix<MatrixType>		ppMult(cgMatrix<MatrixType> A, cgMatrix<MatrixType> B)
{
	cgMatrix<MatrixType> C(A.getLines(), A.getColumns(), 0);

	for(int i=0; i < A.getLines(); i++)
		for(int j=0; j < A.getColumns(); j++)
			C[i][j] = A[i][j] * B[i][j];
			
	return C;			
}

// Point-to-point matrix division
template <class MatrixType>
cgMatrix<MatrixType>		ppDiv(cgMatrix<MatrixType> A, cgMatrix<MatrixType> B)
{
	cgMatrix<MatrixType> C(A.getLines(), A.getColumns(), 0);

	for(int i=0; i < A.getLines(); i++)
		for(int j=0; j < A.getColumns(); j++)
			C[i][j] = A[i][j] / B[i][j];
			
	return C;
}

// Sequences Convolution
template <class MatrixType>
cgMatrix<MatrixType>		Conv(cgMatrix<MatrixType> A, cgMatrix<MatrixType> B)
{
	cgMatrix<MatrixType> C(1, A.getColumns() + B.getColumns() - 1, 0);
	MatrixType Accum	=  0, 	Prod = 0;
	int B_Columns 	 	=  B.getColumns();
	int A_Columns		=  A.getColumns();
	int C_Columns		=  C.getColumns();	
	
  	// Local variables
  	int i, j;

  	// Multiplying the polynomials
  	for (j = 0; j <= (A_Columns + B_Columns); j++) {
      		C[0][j] = 0;
      		for (i = 0; i <= A_Columns; i++) {
			if ((j - i) >= 0 && (j - i) <= B_Columns)
	    			C[0][j] += A[0][i] * B[0][j - i];
		}
    	}
	
	return C;
}

// Vertical concatenation
template <class MatrixType>
cgMatrix<MatrixType>		VertConc(cgMatrix<MatrixType> A, cgMatrix<MatrixType> B)
{
	cgMatrix<MatrixType> C(A.getLines() + B.getLines(), A.getColumns());
					
	for(int i = 0; i < C.getLines(); i++)
		for(int j = 0; j < C.getColumns(); j++)
			if (i < A.getLines()) {
				C[i][j] = A[i][j];
			}
			else {
				C[i][j] = B[i - A.getLines()][j];
			}
				
	return C;
}

// Horizontal concatenation
template <class MatrixType>
cgMatrix<MatrixType>		HorizConc(cgMatrix<MatrixType> A, cgMatrix<MatrixType> B)
{
	cgMatrix<MatrixType> C(	A.getLines(), A.getColumns() + B.getColumns());
				
	for(int i = 0; i < C.getLines(); i++)
		for(int j = 0; j < C.getColumns(); j++)
			if (j < A.getColumns())
				C[i][j] = A[i][j];
			else
				C[i][j] = B[i][j - A.getColumns()];
	
	return C;
}

// LU Decomposition
template <class MatrixType>
void				LUDecomp(cgMatrix<MatrixType> A, cgMatrix<MatrixType> &L, cgMatrix<MatrixType> &U)
{
	L.setEye();
	for(int i=0; i < A.getLines(); i++) {
		for(int j=i+1; j < A.getLines(); j++) {
			L[j][i]	= A[j][i] / A[i][i];
			for(int k=0; k<A.getColumns(); k++)
				A[j][k] = A[j][k] - L[j][i]*A[i][k];
		}
	}
	U = A;
}

// LDU Decomposition
template <class MatrixType>
void				LDUDecomp(cgMatrix<MatrixType> A, cgMatrix<MatrixType> &L, cgMatrix<MatrixType> &D, cgMatrix<MatrixType> &U)
{
	A.Print();	
	L.setEye();
	for(int i=0; i < A.getLines(); i++) {
		for(int j=i+1; j < A.getLines(); j++) {
			L[j][i]	= A[j][i] / A[i][i];
			for(int k=0; k<A.getColumns(); k++)
				A[j][k] = A[j][k] - L[j][i]*A[i][k];
		}
	}
	
	for(int i=0; i < A.getLines(); i++)
		for(int j=0; j < A.getColumns(); j++)
			if (i == j)
				D[i][j] = A[i][j];
			
	for(int i=0; i < A.getLines(); i++)
		for(int j=0; j < A.getColumns(); j++)
			if (i == j)
				U[i][j] = 1;
			else
				U[i][j] = A[i][j];
}

// Matrix Inversion
template <class MatrixType>
cgMatrix<MatrixType>		Inv(cgMatrix<MatrixType> A)
{
	cgMatrix<MatrixType> I(A.getLines(), A.getColumns());
	I.setEye();
	
	cgMatrix<MatrixType> S(A.getLines(), 2*A.getColumns());
	S = HorizConc(A, I);
	
	cgMatrix<MatrixType> L(S.getLines(), S.getColumns());
	cgMatrix<MatrixType> U(S.getLines(), S.getColumns());
	
	LUDecomp(S, L, U);	
	
	cout << endl;
	U.Print();
	cout << endl;
	
	cgMatrix<MatrixType> J(A.getLines(), A.getColumns());
	for(int i=A.getLines()-1; i >= 0; i--) { 
		for(int j=i-1; j >= 0; j--) {
			J[j][i]	= U[i][i] / U[j][i];
			for(int k=0; k<U.getColumns(); k++)
				U[j][k] = U[i][k] - J[j][i]*U[j][k];
		}
	}
	
	for(int i=0; i < A.getLines(); i++) 
		for(int j=U.getColumns()-1; j >=0 ; j--) 
			U[i][j]	= U[i][j] / U[i][i];
	
	cout << endl;
	U.Print();
	cout << endl;
	
	cgMatrix<MatrixType> temp(A.getLines(), A.getColumns());
	
	temp = U(0, A.getLines()-1, A.getColumns(), 2*A.getColumns()-1);
	
	return(temp);
}

template <class MatrixType>
cgMatrix<MatrixType>		FFT(cgMatrix<MatrixType> X)
{
	return X;
}

#endif
