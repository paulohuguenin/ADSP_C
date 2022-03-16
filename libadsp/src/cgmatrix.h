// CGMATRIX  Generalized Matlab-like matrix manipulation library
//
// Set of functions created in order to make easier the matrix
// manipulation inside C++ codes.
//
// cgtotal March 2004.

// Using IFNDEF in order to prevent defining twice the parameters
// included in this header file.
#ifndef CGMATRIX_H
#define CGMATRIX_H

#include <iostream>
#include <fstream>
#include <new>
#include <stdio.h>
#include <string.h>
#include "defines.h"

using namespace std;

// Creating a generic stack class
// for manipulating matrices
template <class MatrixType> class cgMatrix {
private:
  			// Information data stored in the object
			MatrixType 		**Data;

			// Number of Lines and Columns of the matrix object
  			int			NumberOfLines;
  			int			NumberOfColumns;

public:
			// Default Constructor
			cgMatrix();

			// Constructor
  			cgMatrix(int, int, MatrixType);

			// Constructor with external vector loading
			cgMatrix(int, int, MatrixType*);

			// Constructor with external matrix loading
			cgMatrix(int, int, MatrixType**);

			// Copy Constructor
			cgMatrix(const cgMatrix<MatrixType> &cgCopy);

			// Destructor
			~cgMatrix();

			// Manipulating Private Data
			void					setCGMatrix(int, int, MatrixType);
			void					setCGMatrix(int, int, MatrixType*);
			void					setCGMatrix(int, int, MatrixType**);
			MatrixType**			getDataPtr() const;
			MatrixType*				getDataLinePtr(int) const;
  			MatrixType				getData(int, int) const;
  			void					setData(int, int, MatrixType);
			int						getLines() const;
			int						getColumns() const;
			void					setNumberOfLines(int);
			void					setNumberOfColumns(int);
			cgMatrix<MatrixType>	size() const;
			void					setEye();
			cgMatrix<MatrixType>	transpose() const;
			double					norm() const;
			void					fillVector(MatrixType*);
			void					fillVector(MatrixType*,int num);
			void					fillVector(MatrixType*,int a, int b);
			void					fillMatrix(MatrixType**);
			void					zeros();
			void					rotate();

			cgMatrix<MatrixType>	fast_dprod(const MatrixType* vector, int a, int b);
			double               	dprodInterval(const MatrixType* vector, int a, int b);

			// Printing the Information data on the screen
			void					Print();
			// Printing the Information data to file
			void					PrintToFile(char* fileName);

			// Creating and Overloading common operators
			cgMatrix<MatrixType>	operator+(const cgMatrix<MatrixType>& operand) const;
			const cgMatrix<MatrixType>	operator+=(const cgMatrix<MatrixType>& operand);
			cgMatrix<MatrixType>	operator-(const cgMatrix<MatrixType>& operand) const;

			cgMatrix<MatrixType>				operator*(const MatrixType) const;
			const cgMatrix<MatrixType>&			operator*=(const MatrixType);
			cgMatrix<MatrixType>				operator*(const cgMatrix<MatrixType>& operand) const;
			cgMatrix<MatrixType>				operator*(const MatrixType* vector) const;
			// For Windows
#ifdef WIN
			friend	cgMatrix<MatrixType> operator*(const MatrixType, const cgMatrix<MatrixType>&);

#else
			// The element <> is needed to run the following overload, as
			// explained in http://www.cygwin.com/ml/cygwin/1999-03/msg00477.html
			//friend	cgMatrix<MatrixType>		operator*<>(const MatrixType, const cgMatrix<MatrixType>&);
#endif

			cgMatrix<MatrixType>			operator/(const MatrixType) const;
			const cgMatrix<MatrixType>&		operator/=(const MatrixType);

			const cgMatrix<MatrixType>&		operator=(const cgMatrix<MatrixType>& operand);

			MatrixType*						operator[](const int);

			MatrixType						operator()(const int,const int) const;
			cgMatrix<MatrixType>			operator()(const int, const int, const int, const int) const;
			cgMatrix<MatrixType>			operator()(const int, const int, const char) const;
			cgMatrix<MatrixType>			operator()(const char, const int, const int) const;

			bool							operator==( const cgMatrix<MatrixType>& operand) const;
			bool							operator!=( const cgMatrix<MatrixType>& operand) const;
};


// Default Constructor method
template <class MatrixType> cgMatrix<MatrixType>::cgMatrix()
{
	Data = NULL;
	NumberOfLines = 0;
	NumberOfColumns = 0;
}

// Constructor method. Used for allocating memory
// for the objects.
template <class MatrixType> cgMatrix<MatrixType>::cgMatrix(int NumLines, int NumCol, MatrixType Init)
{
	// Allocating memory for a matrix, interpreted
	// as a vector of lines.
	try {
		this->Data = new MatrixType*[NumLines];
	} catch (bad_alloc xa) {
		cout << "Allocation Failure!" << endl;
	}

	// Allocating memory for the columns along the
	// lines.
	try {
		for(int i = 0; i < NumLines; i++)
			this->Data[i] = new MatrixType[NumCol];
	} catch (bad_alloc xa) {
		cout << "Allocation Failure!" << endl;
	}

	// Defining the Number of lines and columns, using parameters of the function.
	NumberOfLines	=	NumLines;
	NumberOfColumns	=	NumCol;

	// Initializing the information data using the parameter Init.
	for(int i = 0; i < NumberOfLines; i++)
		for(int j = 0; j < NumberOfColumns; j++)
			this->Data[i][j] = Init;
}

// Constructor method. Used for allocating memory
// for the objects. Load the external VECTOR content.
template <class MatrixType> cgMatrix<MatrixType>::cgMatrix(int NumLines, int NumCol, MatrixType* vector)
{
	// Allocating memory for a matrix, interpreted
	// as a vector of lines.
	try {
		this->Data = new MatrixType*[NumLines];
	} catch (bad_alloc xa) {
		cout << "Allocation Failure!" << endl;
	}

	// Allocating memory for the columns along the
	// lines.
	try {
		for(int i = 0; i < NumLines; i++)
			this->Data[i] = new MatrixType[NumCol];
	} catch (bad_alloc xa) {
		cout << "Allocation Failure!" << endl;
	}

	// Defining the Number of lines and columns, using parameters of the function.
	NumberOfLines	=	NumLines;
	NumberOfColumns	=	NumCol;

	// Initializing the information data using the parameter Init.
	for(int i = 0; i < NumberOfLines; i++)
		memcpy(Data[i],vector, NumberOfColumns * sizeof(MatrixType));
	//	for(int j = 0; j < NumberOfColumns; j++)
	//		this->Data[i][j] = vector[j];
}

// Constructor method. Used for allocating memory
// for the objects. Load the external MATRIX content.
template <class MatrixType> cgMatrix<MatrixType>::cgMatrix(int NumLines, int NumCol, MatrixType** matrix)
{
	// Allocating memory for a matrix, interpreted
	// as a vector of lines.
	try {
		this->Data = new MatrixType*[NumLines];
	} catch (bad_alloc xa) {
		cout << "Allocation Failure!" << endl;
	}

	// Allocating memory for the columns along the
	// lines.
	try {
		for(int i = 0; i < NumLines; i++)
			this->Data[i] = new MatrixType[NumCol];
	} catch (bad_alloc xa) {
		cout << "Allocation Failure!" << endl;
	}

	// Defining the Number of lines and columns, using parameters of the function.
	NumberOfLines	=	NumLines;
	NumberOfColumns	=	NumCol;

	// Initializing the information data using the parameter Init.
	for(int i = 0; i < NumberOfLines; i++)
		memcpy(Data[i],matrix[i], NumberOfColumns * sizeof(MatrixType));
	//	for(int j = 0; j < NumberOfColumns; j++)
	//		this->Data[i][j] = matrix[i][j];
}


// Copy constructor. Used in operation overloading.
// The same as the constructor, but the function
// parameter is passed by reference.
template <class MatrixType> cgMatrix<MatrixType>::cgMatrix(const cgMatrix<MatrixType> &cgCopy)
{

	// Allocating memory for a matrix, interpreted
	// as a vector of lines.
	try {
		this->Data = new MatrixType*[cgCopy.NumberOfLines];
	} catch (bad_alloc xa) {
		cout << "Allocation Failure!" << endl;
	}

	// Allocating memory for the columns along the
	// lines.
	try {
		for(int i = 0; i < cgCopy.NumberOfLines; i++)
			this->Data[i] = new MatrixType[cgCopy.NumberOfColumns];
	} catch (bad_alloc xa) {
		cout << "Allocation Failure!" << endl;
	}

	// Defining the Number of lines and columns, using parameters of the function.
	NumberOfLines	=	cgCopy.NumberOfLines;
	NumberOfColumns	=	cgCopy.NumberOfColumns;

	MatrixType** pCopy;
	pCopy = cgCopy.getDataPtr();
	// Initializing the information data using the parameter Init.
	for(int i = 0; i < cgCopy.NumberOfLines; i++)
	{
		memcpy(	Data[i],
				pCopy[i],
				NumberOfColumns * sizeof(MatrixType));
	}
	//		for(int j = 0; j < cgCopy.NumberOfColumns; j++)
	//		this->Data[i][j] = cgCopy.Data[i][j];
}

// Destructor. Used for freeing memory.
template <class MatrixType> cgMatrix<MatrixType>::~cgMatrix()
{
	if (this->Data != NULL)
	{
		// Freeing the column members along a line.
		for(int i = 0; i < NumberOfLines; i++)
			delete [] this->Data[i];

		// Freeing the vector of lines.
		delete [] this->Data;
	}
}



template <class MatrixType> void cgMatrix<MatrixType>::setCGMatrix(int NumLines, int NumCol, MatrixType Init)
{
	// Allocating memory for a matrix, interpreted
	// as a vector of lines.
	try {
		this->Data = new MatrixType*[NumLines];
	} catch (bad_alloc xa) {
		cout << "Allocation Failure!" << endl;
	}

	// Allocating memory for the columns along the
	// lines.
	try {
		for(int i = 0; i < NumLines; i++)
			this->Data[i] = new MatrixType[NumCol];
	} catch (bad_alloc xa) {
		cout << "Allocation Failure!" << endl;
	}

	// Defining the Number of lines and columns, using parameters of the function.
	NumberOfLines	=	NumLines;
	NumberOfColumns	=	NumCol;

	// Initializing the information data using the parameter Init.
	for(int i = 0; i < NumberOfLines; i++)
		for(int j = 0; j < NumberOfColumns; j++)
			this->Data[i][j] = Init;
}


template <class MatrixType> void cgMatrix<MatrixType>::setCGMatrix(int NumLines, int NumCol, MatrixType* vector)
{
	// Allocating memory for a matrix, interpreted
	// as a vector of lines.
	try {
		this->Data = new MatrixType*[NumLines];
	} catch (bad_alloc xa) {
		cout << "Allocation Failure!" << endl;
	}

	// Allocating memory for the columns along the
	// lines.
	try {
		for(int i = 0; i < NumLines; i++)
			this->Data[i] = new MatrixType[NumCol];
	} catch (bad_alloc xa) {
		cout << "Allocation Failure!" << endl;
	}

	// Defining the Number of lines and columns, using parameters of the function.
	NumberOfLines	=	NumLines;
	NumberOfColumns	=	NumCol;

	// Initializing the information data using the parameter Init.
	for(int i = 0; i < NumberOfLines; i++)
		memcpy(Data[i],vector, NumberOfColumns * sizeof(MatrixType));
	//	for(int j = 0; j < NumberOfColumns; j++)
	//		this->Data[i][j] = vector[j];
}


template <class MatrixType> void cgMatrix<MatrixType>::setCGMatrix(int NumLines, int NumCol, MatrixType** matrix)
{
	// Allocating memory for a matrix, interpreted
	// as a vector of lines.
	try {
		this->Data = new MatrixType*[NumLines];
	} catch (bad_alloc xa) {
		cout << "Allocation Failure!" << endl;
	}

	// Allocating memory for the columns along the
	// lines.
	try {
		for(int i = 0; i < NumLines; i++)
			this->Data[i] = new MatrixType[NumCol];
	} catch (bad_alloc xa) {
		cout << "Allocation Failure!" << endl;
	}

	// Defining the Number of lines and columns, using parameters of the function.
	NumberOfLines	=	NumLines;
	NumberOfColumns	=	NumCol;

	// Initializing the information data using the parameter Init.
	for(int i = 0; i < NumberOfLines; i++)
		memcpy(Data[i],matrix[i], NumberOfColumns * sizeof(MatrixType));
	//	for(int j = 0; j < NumberOfColumns; j++)
	//		this->Data[i][j] = matrix[i][j];
}



// Getting Data Pointer
template <class MatrixType> MatrixType** cgMatrix<MatrixType>::getDataPtr() const
{
	return this->Data;
}

// Getting Data Line Pointer
template <class MatrixType> MatrixType* cgMatrix<MatrixType>::getDataLinePtr(int NumLine) const
{
	return this->Data[NumLine];
}

// Getting the private attribute named Data.
template <class MatrixType> MatrixType cgMatrix<MatrixType>::getData(int NumLine, int NumCol) const
{
	return(Data[NumLine][NumCol]);
}

// Manipulating private data (attribute NumberOfLines).
template <class MatrixType> int cgMatrix<MatrixType>::getLines() const
{
	return(NumberOfLines);
}

// Manipulating private data (attribute NumberOfColumns).
template <class MatrixType> int cgMatrix<MatrixType>::getColumns() const
{
	return(NumberOfColumns);
}

template <class MatrixType> void cgMatrix<MatrixType>::setNumberOfLines(int val)
{
	NumberOfLines = val;
}
template <class MatrixType> void cgMatrix<MatrixType>::setNumberOfColumns(int val)
{
	NumberOfColumns = val;
}

// Setting the private attribute named Data.
template <class MatrixType> void cgMatrix<MatrixType>::setData(int NumLine, int NumCol, MatrixType Val)
{
	Data[NumLine][NumCol] = Val;
}

// Printing the attributes NumberOfLines and NumberOfColumns on the Screen
template <class MatrixType> cgMatrix<MatrixType> cgMatrix<MatrixType>::size() const
{
	// cout << "Dimensions: " << this->getLines() << "x" << this->getColumns() << endl;

	cgMatrix<MatrixType> temp(1, 2, 0.0);
	temp[0][0]	=	this->getLines();
	temp[0][1]  =	this->getColumns();

	return temp;
}

// Setting matrix eye
template <class MatrixType> void cgMatrix<MatrixType>::setEye()
{
	for(int i=0; i < getLines(); i++)
		for(int j=0; j < getColumns(); j++)
			if (i == j)
				Data[i][j] = 1;
			else
				Data[i][j] = 0;
}

// Transpose matrix
template <class MatrixType> cgMatrix<MatrixType> cgMatrix<MatrixType>::transpose() const
{
	// Creating an auxiliary object
	cgMatrix<MatrixType> temp(NumberOfColumns, NumberOfLines,0.0);

	for(int i=0; i<NumberOfLines; i++)
		for(int j=0; j<NumberOfColumns; j++)
			temp[j][i] = getData(i, j);

	return temp;
}

// Calculate the matrix norm
template <class MatrixType> double cgMatrix<MatrixType>::norm() const
{
	double acumulator = 0;

	for(int i = 0; i < NumberOfLines; i++)
	{
		for(int j = 0; j < NumberOfColumns; j++)
		{
			acumulator += Data[i][j]*Data[i][j];
		}
	}

	return sqrt(acumulator);
}

// Filling Vector
template <class MatrixType> void cgMatrix<MatrixType>::fillVector(MatrixType* vector)
{
	// Initializing the information data using the parameter Init.
	memcpy( Data[0],vector,NumberOfColumns * sizeof(MatrixType) );
	/*for(int j = 0; j < NumberOfColumns; j++)
	{
		this->Data[0][j] = vector[j];
	}*/
}

// Filling Vector
template <class MatrixType> void cgMatrix<MatrixType>::fillVector(MatrixType* vector, int num)
{
	// Initializing the information data using the parameter Init.
	memcpy( Data[0],vector, num * sizeof(MatrixType) );
	/*for(int j = 0; j < num; j++)
	{
		this->Data[0][j] = vector[j];
	}
	*/
}

// Filling Vector
template <class MatrixType> void cgMatrix<MatrixType>::fillVector(MatrixType* vector, int a, int b)
{
	// Initializing the information data using the parameter Init.
	memcpy( Data[a],vector, (b-a+1) * sizeof(MatrixType) );
	/*for(int j = 0; j < num; j++)
	{
		this->Data[0][j] = vector[j];
	}
	*/
}


// Filling Matrix
template <class MatrixType> void cgMatrix<MatrixType>::fillMatrix(MatrixType** matrix)
{
	// Initializing the information data using the parameter Init.
	for(int i = 0; i < NumberOfLines; i++)
	{
		memcpy( Data[i],
				matrix[i],
				NumberOfColumns * sizeof(MatrixType) );
	//	for(int j = 0; j < NumberOfColumns; j++)
	//	{
	//		this->Data[i][j] = matrix[i][j];
	//	}
	}
}

template <class MatrixType> void	cgMatrix<MatrixType>::zeros()
{
	for(int i = 0; i < NumberOfLines; i++)
	{
		for(int j = 0; j < NumberOfColumns; j++)
		{
			this->Data[i][j] = 0;
		}
	}
}

template <class MatrixType> void	cgMatrix<MatrixType>::rotate()
{
	// Creating an auxiliary object
	cgMatrix<MatrixType> temp(NumberOfLines, NumberOfColumns, 0.0);
	int i;
	int j;
	for( i=0; i<NumberOfLines; i++)
	{
		for( j=0; j<NumberOfColumns; j++)
		{
			temp[i][NumberOfColumns-1-j] = this->Data[i][j];
		}
	}

	for( i = 0; i < NumberOfLines; i++)
	{
		for( j = 0; j < NumberOfColumns; j++)
		{
			this->Data[i][j] = temp[i][j];
		}
	}
}



template <class MatrixType> cgMatrix<MatrixType> cgMatrix<MatrixType>::fast_dprod(const MatrixType* vector, int a, int b)

{

	MatrixType acumulator;
	acumulator = 0;
	for(int i = a; i <= b; i++)
	{
		acumulator += this->getData(0,i)* vector[i];
	}
	cgMatrix<MatrixType> temp(1, 1, acumulator);
	return temp;

}

template <class MatrixType> double cgMatrix<MatrixType>::dprodInterval(const MatrixType* vector, int a, int b)

{

	MatrixType acumulator;
	acumulator = 0;
    int j=0;
	for(int i = a; i <= b; i++)
	{
		acumulator += this->getData(0,i)* vector[j];
		j++;
	}

	return acumulator;

}


// Printing the attribute Data on the Screen
template <class MatrixType> void cgMatrix<MatrixType>::Print()
{
// 	int i, j;
//
// 	// Showing the information data from the
// 	// caller of the method.
// 	for(i = 0; i < NumberOfLines; i++){
// 		for(j = 0; j < NumberOfColumns; j++) {
// 			cout << getData(i, j) << " ";
// 		}
// 		cout << endl;
// 	}

	if ((NumberOfLines > 0) && (NumberOfColumns > 0)) {
		int i, j;
		for (i = 0; i < NumberOfLines; i++) {
			for (j = 0; j < NumberOfColumns; j++) {
				cout.fill(' ');
				cout.width(10);
				cout.precision(4);
				cout << Data[i][j];
			}
			cout << ";" << endl;
		}
	}
}

// Printing the attribute Data to File
template <class MatrixType> void cgMatrix<MatrixType>::PrintToFile(char* fileName)
{

	if ((NumberOfLines > 0) && (NumberOfColumns > 0))
	{
		ofstream outMatrixFile(fileName,ios::out);

		if (!outMatrixFile)
		{
			cerr << "Cannot open file" << endl;
		}

		int i, j;
		for (i = 0; i < NumberOfLines; i++)
		{
			for (j = 0; j < NumberOfColumns; j++)
			{
				outMatrixFile << this->Data[i][j] << '\n' ;
			}
			//outMatrixFile << '\n';
		}
	}
}


// ==============================================
// OPERATOR OVERLOADING SECTION
// ==============================================

// Operator overloading (operator +).
template <class MatrixType> cgMatrix<MatrixType> cgMatrix<MatrixType>::operator+(const cgMatrix& operand) const
{
	// Creating an auxiliary object
	cgMatrix<MatrixType> temp(NumberOfLines, NumberOfColumns, 0);

	// Summing the information data from the parameter operand and the
	// caller of the method. The result is stored in the auxiliary object.
	for(int i = 0; i < this->NumberOfLines; i++)
		for(int j = 0; j < this->NumberOfColumns; j++) {
			temp.setData(i, j, this->getData(i, j) + operand.getData(i, j));
		}

	return temp;
}

// Operator overloading (operator +=).
template <class MatrixType> const cgMatrix<MatrixType> cgMatrix<MatrixType>::operator+=(const cgMatrix& operand)
{
	// Summing the information data from the parameter operand and the
	// caller of the method. The result is stored in the auxiliary object.
	for(int i = 0; i < this->NumberOfLines; i++)
		for(int j = 0; j < this->NumberOfColumns; j++) {
			this->setData(i, j, this->getData(i, j) + operand.getData(i, j));
		}

	return *this;
}

// Operator overloading (operator -).
template <class MatrixType> cgMatrix<MatrixType> cgMatrix<MatrixType>::operator-(const cgMatrix& operand) const
{
	// Creating an auxiliary object
	cgMatrix<MatrixType> temp(NumberOfLines, NumberOfColumns, 0.0);

	// Subtracting the information data from the parameter operand from the
	// caller of the method. The result is stored in the auxiliary object.
	for(int i = 0; i < this->NumberOfLines; i++)
		for(int j = 0; j < this->NumberOfColumns; j++) {
			temp.setData(i, j, this->getData(i, j) - operand.getData(i, j));
		}

	return temp;
}

// Operator overloading (operator * - right multiplying by scalars).
template <class MatrixType> cgMatrix<MatrixType> cgMatrix<MatrixType>::operator*(const MatrixType scalar) const
{
	// Creating an auxiliary object
	cgMatrix<MatrixType> temp(NumberOfLines, NumberOfColumns, 0.0);

	// Multiplying the parameter scalar and the
	// caller of the method. The result is stored in the auxiliary object.
	for(int i = 0; i < this->NumberOfLines; i++)
		for(int j = 0; j < this->NumberOfColumns; j++) {
			temp.setData(i, j, this->getData(i, j) * scalar);
		}

	return temp;
}

// Operator overloading (operator *=).
template <class MatrixType> const cgMatrix<MatrixType>& cgMatrix<MatrixType>::operator*=(const MatrixType scalar)
{
	for(int i = 0; i < this->NumberOfLines; i++)
		for(int j = 0; j < this->NumberOfColumns; j++) {
			this->Data[i][j] *= scalar;
		}
	return *this;
}

// Operator overloading (operator * - multiplying by cgmatrices).
template <class MatrixType> cgMatrix<MatrixType> cgMatrix<MatrixType>::operator*(const cgMatrix& operand) const
{
	// Creating auxiliary objects.
	cgMatrix<MatrixType> temp(NumberOfLines, NumberOfColumns, 0.0);
	MatrixType acumulator;
	int LinesResult, ColumnsResult;

	// Testing if the multiplication is possible.
	if (this->NumberOfColumns != operand.getLines()) {
		cout << " Invalid Operation " << endl;
		cout << " Number of Columns of the left operand must equals the Number of ";
		cout << " Lines of the right operand."<< endl;
		exit(1);
	}

	// Defining the number of lines and
	// columns of the result.
	LinesResult = this->NumberOfLines;
	temp.setNumberOfLines(LinesResult);
	ColumnsResult = operand.getColumns();
	temp.setNumberOfColumns(ColumnsResult);

	// Subtracting the information data from the parameter operand from the
	// caller of the method. The result is stored in the auxiliary object.
	for(int i = 0; i < LinesResult; i++)
	{
		for(int j = 0; j < ColumnsResult; j++)
		{
			acumulator = 0;
			for (int k = 0; k < this->NumberOfColumns; k++)
			{
				acumulator += this->getData(i,k)*operand.getData(k,j);
			}
			temp.setData(i,j,(double)acumulator);
		}
	}

	return temp;
}

// Operator overloading (operator * - multiplying by vector).
template <class MatrixType> cgMatrix<MatrixType> cgMatrix<MatrixType>::operator*(const MatrixType* vector) const
{
	MatrixType acumulator;

	acumulator = 0;
	for(int i = 0; i < NumberOfColumns; i++)
	{
		acumulator += this->getData(0,i)* vector[i];
	}

	cgMatrix<MatrixType> temp(1, 1, acumulator);

	return temp;
}


// Operator overloading (operator * - left multiplying by scalars).
// template <class MatrixType> cgMatrix<MatrixType> operator*(const MatrixType scalar, const cgMatrix<MatrixType>& operand)
// {
// 	// Creating an auxiliary object
// 	cgMatrix<MatrixType> temp(operand.NumberOfLines, operand.NumberOfColumns, 0.0);
//
// 	// Multiplying the parameter scalar and the
// 	// caller of the method. The result is stored in the auxiliary object.
// 	for(int i = 0; i < operand.NumberOfLines; i++)
// 		for(int j = 0; j < operand.NumberOfColumns; j++) {
// 			temp.setData(i, j, operand.getData(i, j) * scalar);
// 		}
//
// 	return temp;
// }

// Operator overloading (operator =).
template <class MatrixType> const cgMatrix<MatrixType>& cgMatrix<MatrixType>::operator=(const cgMatrix<MatrixType>& operand)
{
	if (this->Data == NULL)
	{
		this->NumberOfLines = operand.NumberOfLines;
		this->NumberOfColumns = operand.NumberOfColumns;

		// Allocating memory for a matrix, interpreted
		// as a vector of lines.
		try {
			this->Data = new MatrixType*[NumberOfLines];
		} catch (bad_alloc xa) {
			cout << "Allocation Failure!" << endl;
		}

		// Allocating memory for the columns along the
		// lines.
		try {
			for(int i = 0; i < NumberOfLines; i++)
				this->Data[i] = new MatrixType[NumberOfColumns];
		} catch (bad_alloc xa) {
			cout << "Allocation Failure!" << endl;
		}
	}

	// Testing if the equality is possible.
	if ( (this->NumberOfLines != operand.getLines()) ||
		(this->NumberOfColumns != operand.getColumns()) )
	{
		cout << " Invalid Operation " << endl;
		cout << " Number of lines and columns of the left operand ";
		cout << " must equals the number of lines and columns of the right operand."<< endl;
		exit(1);
	}


	// Passing the information data from the parameter operand
	// to the caller of the method
	MatrixType** opDataPtr;
	opDataPtr = operand.getDataPtr();
	for(int i = 0; i < operand.getLines(); i++)
	{
		memcpy(	Data[i],
				opDataPtr[i],
				operand.getColumns() * sizeof(MatrixType));
	}

	/*	for(int j = 0; j < operand.getColumns(); j++)
		{
			this->setData(i, j, operand.getData(i, j));
		}
	*/

	// An auxiliary is not used in order
	// to allow multiple assignments.
	return *this;
}

// Operator overloading (operator /).
template <class MatrixType> cgMatrix<MatrixType> cgMatrix<MatrixType>::operator/(const MatrixType scalar) const
{
	// Creating an auxiliary object
	cgMatrix<MatrixType> temp(NumberOfLines, NumberOfColumns, 0.0);

	// Avoiding division by zero
	if (scalar == 0) {
		cout << " Division by zero " << endl;
		exit(1);
	}

	// Multiplying the parameter scalar and the
	// caller of the method. The result is stored in the auxiliary object.
	for(int i = 0; i < this->NumberOfLines; i++)
		for(int j = 0; j < this->NumberOfColumns; j++) {
			temp.setData(i, j, this->getData(i, j) / scalar);
		}

	return temp;
}

// Operator overloading (operator /=).
template <class MatrixType> const cgMatrix<MatrixType>& cgMatrix<MatrixType>::operator/=(const MatrixType scalar)
{
	// Avoiding division by zero
	if (scalar == 0) {
		cout << " Division by zero " << endl;
		exit(1);
	}

	// Multiplying the parameter scalar and the
	// caller of the method. The result is stored in the auxiliary object.
	for(int i = 0; i < this->NumberOfLines; i++)
		for(int j = 0; j < this->NumberOfColumns; j++) {
			this->Data[i][j] /= scalar;
			//temp.setData(i, j, this->getData(i, j) / scalar);
		}
	return *this;
}


// Operator overloading (operator []).
template <class MatrixType> MatrixType* cgMatrix<MatrixType>::operator[](const int i)
{
	return(Data[i]);
}

// Operator overloading (operator ()).
template <class MatrixType> MatrixType cgMatrix<MatrixType>::operator()(const int i, const int j) const
{
	return(Data[i][j]);
}

// Operator overloading (operator ()).
template <class MatrixType> cgMatrix<MatrixType> cgMatrix<MatrixType>::operator()(const int a, const int b,  const int c, const int d) const
{
	// Creating an auxiliary object
	cgMatrix<MatrixType> temp(b-a+1, d-c+1,0.0);

	for(int i=a; i<=b; i++)
		for(int j=c; j<=d; j++)
			temp[i - a][j - c] = getData(i, j);

	return temp;
}

// Operator overloading (operator ()).
template <class MatrixType> cgMatrix<MatrixType> cgMatrix<MatrixType>::operator()(const int a, const int b, const char c) const
{
	// Creating an auxiliary object
	cgMatrix<MatrixType> temp(b-a+1, NumberOfColumns,0.0);

	for(int i=a; i<=b; i++)
		for(int j=0; j<NumberOfColumns; j++)
			temp[i - a][j] = getData(i, j);

	return temp;
}

// Operator overloading (operator ()).
template <class MatrixType> cgMatrix<MatrixType> cgMatrix<MatrixType>::operator()(const char c, const int a, const int b) const
{
	// Creating an auxiliary object
	cgMatrix<MatrixType> temp(NumberOfLines, a-b+1,0.0);

	for(int i=0; i<NumberOfLines; i++)
		for(int j=a; j<=b; j++)
			temp[i][j - a] = getData(i, j);

	return temp;
}

// Operator overloading (operator == - comparison).
template <class MatrixType> bool cgMatrix<MatrixType>::operator==(const cgMatrix& operand) const
{
	bool result = 1;

	// Seeing if the dimensions of the matrices are diferent
	if ((this->NumberOfLines != operand.getLines()) || (this->NumberOfColumns != operand.getColumns())) {
		result = 0;
		return result;
	}
	else{

		// Comparing the correspondent elements of each matrix
		for(int i = 0; i < this->NumberOfLines; i++)
			for(int j = 0; j < this->NumberOfColumns; j++) {
				if (this->Data[i][j] != operand.getData(i,j))
					result = 0;
			}

		return result;
	}
}

// Operator overloading (operator != - right multiplying by scalars).
template <class MatrixType> bool cgMatrix<MatrixType>::operator!=(const cgMatrix& operand) const
{
	bool result = 0;

	// Seeing if the dimensions of the matrices are diferent
	if ((this->NumberOfLines != operand.getLines()) || (this->NumberOfColumns != operand.getColumns())) {
		result = 1;
		return result;
	}
	else
	{

		// Comparing the correspondent elements of each matrix
		for(int i = 0; i < this->NumberOfLines; i++)
			for(int j = 0; j < this->NumberOfColumns; j++) {
				if (this->Data[i][j] != operand.getData(i,j))
					result = 1;
			}

		return result;
	}
}


#endif
