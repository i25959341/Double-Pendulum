//Program to model a single pendulum for small angular displacements
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstring>
#include <cstdlib>

using namespace std;

#define I_TYPE double, double, double, double, double, double, double, int
#define CALL_I h,gamma_eff,theta_Initial,u_Initial,l,m,g,imax
#define I_FCN double h,double gamma_eff,double theta_Initial,double u_Initial,double l,double m,double g,int imax


void euler(I_FCN);
void leapfrog(I_FCN);
void rk4(I_FCN);

double calculate_energy(double, double, double , double , double );


int main()
{
	double h =0.01; //step size
	int imax = 10000; //Define max number of steps

	double m=1.0; //mass
	double gamma=0.2; //dampening constant, zero for undamped
	double l=1.0; //length
	double g=1.0; //gravitational constant


	double gamma_eff=(gamma/m)*sqrt(1.0/l*g);

	// Initiate initial conditions

	double theta_Initial = 0.1; //Initial value of theta 
	double u_Initial = 0.0;		//Initial value for derivative of theta

	//Using Switch function to let user choose which method to use

	int choice;
	cout<<"Select Method:\n\n1 Euler\n2 Leapfrog\n3 Runge-Kutta 4\n4 All\n\n";
	cin>>choice;

	switch (choice) {
		case 1:
			euler(CALL_I);
			break;

		case 2:
			leapfrog(CALL_I);
			break;

		case 3:
			rk4(CALL_I);
			break;
		case 4:
			euler(CALL_I);
			leapfrog(CALL_I);
			rk4(CALL_I);
			break;

		default:
			cout<<"\nPlease enter an integer 1-3\n\n";
			exit (99);
			break;
	};


	return 0;
}


void euler(I_FCN)
{
	ofstream outEul("Sing_Eul_Data.txt");
	outEul<< "t\t"<<"Theta Euler\t"<<"Energy\n" ;

	double theta = theta_Initial ;
	double u= u_Initial ;

	double theta1;
	double u1;

	double t=0.0; //Dimensionless time
	double energy;

	//EULER BEGIN******************************************

	for (int i=1; i<=imax; i++){

		t=i*h;

		theta1= theta+h*u;
		u1= -h*theta+(1.0-h*gamma_eff)*u;

		energy = calculate_energy(l, theta1, m, g, u1);

		outEul << t << "\t" << theta1 << "\t"<< energy << "\n";

		theta= theta1;
		u= u1;

	}
	return;
}


void leapfrog(I_FCN)	{

		ofstream outLF("Sing_LF_Data.txt"); //Leapfrog Stream
		outLF<< "t\t" << "Theta Leapfrog\t" << "Energy \n" ;

		double theta = theta_Initial ; //Set initial conditions
		double u= u_Initial ;

		// Variables used for LF Update
		double theta1;
		double u1;
		double theta2;
		double u2;

		double t=0.0; //Dimensionless time
		double energy;

		// LF BEGIN***************

		theta1= theta+h*u;
		u1= -h*theta+(1.0-h*gamma_eff)*u;

		outLF<<0<<"\t"<<theta<<"\n";//output initial conditions
		outLF<<h<<"\t"<<theta1<<"\n"; //output 1st Euler step

		//LF iteration
		for (int i = 2; i<=imax; i++){
			t=i*h; //Dimensionless time

			theta2= theta+2.0*h*u1;
			u2=u+2.0*h*((-theta1)-gamma_eff*u1);

			energy = calculate_energy(l, theta1, m, g, u1);

			outLF << t << "\t" << theta1 << "\t"<< energy << "\n";
			//Update theta and u

			theta=theta1;
			theta1=theta2;

			u=u1;
			u1=u2;

		}
		return;

	}


void rk4 (I_FCN){

		ofstream outRK("Sing_RK4_Data.txt"); //rk4 Stream
		outRK<< "t\t" << "Theta RK\t" << "Energy \n" ;

		double theta = theta_Initial ; //Set initial conditions
		double u= u_Initial ;
		double t=0.0; //Dimensionless time

		double theta1;
		double u1;

		vector<double>k1(2);
		vector<double>k2(2);
		vector<double>k3(2);
		vector<double>k4(2);

		double energy = calculate_energy(l,theta,m,g,u); //variable for calculating energy at each loop

		for (int i = 1; i<=imax; i++){
			
		t=i*h;			

		k1[0]=h*u;
		k1[1]=h*(-theta-gamma_eff*u);

		k2[0]=h*(u+0.5*k1[1]);
		k2[1]=h*(-(theta+0.5*k1[0])-gamma_eff*(u+0.5*k1[1]));

		k3[0]=h*(u+0.5*k2[1]);
		k3[1]=h*(-(theta+0.5*k2[0])-gamma_eff*(u+0.5*k2[1]));

		k4[0]=h*(u+k3[1]);
		k4[1]=h*(-(theta+k3[0])-gamma_eff*(u+k3[1]));

		theta1=theta+(1.0/6.0)*(k1[0]+2.0*k2[0]+2.0*k3[0]+k4[0]);
		u1= u + (1.0/6.0)*(k1[1]+2.0*k2[1]+2.0*k3[1]+k4[1]);

		energy = calculate_energy(l, theta1, m, g, u1);

		outRK << t << "\t" << theta1 << "\t"<< energy << "\n";

				
		u=u1;
		theta=theta1;	

		}
		
		return;
}


//Function monitoring energy
double calculate_energy(double l, double theta, double mass, double g, double u)
{

	double energy= 0.5*(theta*theta+u*u); //Total Energy

	return energy;

}
