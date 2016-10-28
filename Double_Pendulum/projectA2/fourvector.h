//define fourvector class

#ifndef FOURVECTOR_H
#define FOURVECTOR_H

#include "fourmatrix.h"

class fourvector
{
private:
	double vector[4]; //array of vector elements

public:
//default constructor
	fourvector()
	{
		for (int i=0;i<4;i++){
			vector[i]=0.0;
		}
	}

//constructor
	fourvector(double a1, double a2, double a3, double a4) //a1,a2,a3,a4: 1st,2nd,3rd,4th row
	{
		vector[0]=a1;
		vector[1]=a2; 
		vector[2]=a3;
		vector[3]=a4;
	}


//Access method for vector elements
	double get (int i)
	{
		return vector[i];
	}


//Modifier method for vector elements
	void set(double a1, double a2, double a3, double a4)
	{
		vector[0]=a1;
		vector[1]=a2; 
		vector[2]=a3;
		vector[3]=a4;
	}


//Overloading * operator as scalar multiplication
	fourvector operator*(double factor)
	{
		fourvector product(vector[0]*factor,vector[1]*factor,vector[2]*factor,vector[3]*factor);
		return product;
	}

//Overloading * operator as multiplying fourmatrix
	fourvector operator*(fourmatrix T)
	{
		double a1 = vector[0]*T.get(0) + vector[1]*T.get(4) + vector[2]*T.get(8) + vector[3]*T.get(12);
		double a2 = vector[0]*T.get(1) + vector[1]*T.get(5) + vector[2]*T.get(9) + vector[3]*T.get(13);
		double a3 = vector[0]*T.get(2) + vector[1]*T.get(6) + vector[2]*T.get(10) + vector[3]*T.get(14);
		double a4 = vector[0]*T.get(3) + vector[1]*T.get(7) + vector[2]*T.get(11) + vector[3]*T.get(15);
		fourvector product (a1, a2, a3, a4);
		return product;
	}


//Overloading + operator
	fourvector operator+(fourvector v_other)
	{
		fourvector v_sum(v_other.get(0) + vector[0],
			             v_other.get(1) + vector[1],
						 v_other.get(2) + vector[2],
						 v_other.get(3) + vector[3]);
		return v_sum;
	}

}; 

#endif

