//define fourmatrix class

#ifndef FOURMATRIX_H
#define FOURMATRIX_H

#include "fourvector.h"

class fourmatrix
{
private:
	double matrix[16]; //array of matrix elements (column by column)

public:

//default constructor
	fourmatrix()
	{
		for (int i=0;i<16;i++){
			matrix[i]=0.0;
		}
	}

//constructor
	fourmatrix(double a1, double a2, double a3, double a4, 
		       double b1, double b2, double b3, double b4, 
			   double c1, double c2, double c3, double c4, 
			   double d1, double d2, double d3, double d4)
	{
		matrix[0]=a1; matrix[4]=b1; matrix[8]=c1; matrix[12]=d1; 
		matrix[1]=a2; matrix[5]=b2; matrix[9]=c2; matrix[13]=d2; 
		matrix[2]=a3; matrix[6]=b3; matrix[10]=c3; matrix[14]=d3; 
		matrix[3]=a4; matrix[7]=b4; matrix[11]=c4; matrix[15]=d4; 
	}


//Access method for matrix elements
	double get (int i)
	{
		return matrix[i];
	}

}; 

#endif

