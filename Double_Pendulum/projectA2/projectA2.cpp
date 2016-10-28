#include <iostream> 
#include <fstream>
#include "fourvector.h"
#include <cstring>
#include <cstdlib>

#include "fourmatrix.h"
using namespace std;

//Program to evaluate damped/undamped double pendulum using RK4 method

int main(void)
{

	double R; //mass ratio of lower to upper pendulum
	cout << "Mass ratio R = ";
	cin >> R;
 
	double G; //damping constant
	cout << "Damping constant G = ";
	cin >> G;

	double h; //time step-size
	cout << "Time step-size dt = ";
	cin >> h;

	double t=0.0; //initialise time
	double g=10.0; //gravitational acceleration

    //set matrix operator
	fourmatrix M(0.0, 0.0, -R-1.0, R+1.0, 0.0, 0.0, R, -R-1.0, 1.0, 0.0, -G, G*(1.0-1.0/R), 0.0, 1.0, 0.0, -G/R); 


//RK4 method

	ofstream rk4results("rk4results.txt"); //output file for results

    //initalise state vector (upper and lower angles)(upper and lower angular velocities)
	fourvector v0(0.1, 0.0, 0.0, 0.0); 
 
	for (t=0.0; t<=100.01; t+=h) //time-range of computation
	{
		double a1 = v0.get(0); //upper angle
		double a2 = v0.get(1); //lower angle
		double a3 = v0.get(2); //upper angular velocity
		double a4 = v0.get(3); //lower angular velocity

		//calculate total energy of system
		double E = g*(1.0+R)*a3*a3/2.0 + g*R*a4*a4/2.0 + g*R*a3*a4 + g*(1.0+R)*a1*a1/2.0 + R*g*a2*a2/2.0;

		//print results
		rk4results << t << '\t' << //time
		v0.get(0) << '\t' << v0.get(1) << '\t' << //angles
		v0.get(2) << '\t' << v0.get(3) << '\t' << //angular velocities
		E << endl; //total energy

		fourvector k1 = (v0*M)*h;
		fourvector v1 = v0 + k1*0.5;
		fourvector k2 = (v1*M)*h;
		fourvector v2 = v0 + k2*0.5;
		fourvector k3 = (v2*M)*h;
		fourvector v3 = v0 + k3;
		fourvector k4 = (v3*M)*h;

		v0 = v0 + (k1 + k2*2.0 + k3*2.0 + k4)*(1.0/6.0); //RK4 iteration
	}

	return 0;

}
