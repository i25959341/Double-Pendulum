//Program to model a double pendulum for small angular displacements
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstring>
#include <cstdlib>

using namespace std;

#define I_TYPE double, double, double, double, double, double, double, double, double, double, int
#define CALL_I h,gamma_eff,theta_Initial,u_Initial, omega_Initial, v_Initial, l,m, r,g,imax
#define I_FCN double h,double gamma_eff,double theta_Initial,double u_Initial,double omega_Initial,double v_Initial,double l,double m,double r, double g,int imax

void rk4(I_FCN);
double calculate_energy(double, double, double , double , double );

int main()
{

	double h =0.001; //step size
	int imax = 100000; //Define max number of steps

	double m=1; //mass
	double r=100;
	double gamma=0.0; //dampening constant, zero for undamped
	double l=1.0; //length
	double g=1.0; //gravitational constant

	double gamma_eff=(gamma/m)*sqrt(1.0/l*g);

	// Initiate initial conditions

	double theta_Initial = 0.0; //Initial value of theta
	double u_Initial = 0.1;		//Initial value for derivative of theta
	double omega_Initial= 0.0;
	double v_Initial = 0.0;

	//Using Switch function to let user choose which method to use

	int choice;
	cout<<"Type one to run the Double Pendulum:\n\n1 Run\n\n";
	cin>>choice;

	switch (choice) {

		case 1:
			rk4(CALL_I);
			break;

		default:
			cout<<"\nPlease enter an integer 1\n\n";
			exit (99);
			break;
	};


	return 0;
}


void rk4 (I_FCN){

		ofstream outRK("RK4_Double_Data.txt"); //rk4 Stream
		outRK<< "t\t" << "Theta RK\t" << "Omega RK\t" << "u\t"<<"v\t" << "Energy\n";

		double theta = theta_Initial ; //Set initial conditions
		double u= u_Initial ;
		double omega= omega_Initial ;
		double v=v_Initial ;
		double t=0.0; //Dimensionless time

		double theta1;
		double u1;
		double omega1;
		double v1;

		vector<double>k1(4);
		vector<double>k2(4);
		vector<double>k3(4);
		vector<double>k4(4);

		double energy = calculate_energy(r, theta, omega, u, v); //variable for calculating energy at each loop

		for (int i = 1; i<=imax; i++){

		t=i*h;


		k1[0]=h*(u);
		k1[1]=h*(v);
		k1[2]=h*(-(r+1.0)*theta+r*omega-gamma_eff*u);
		k1[3]=h*((r+1.0)*theta-(r+1)*omega+gamma_eff*(1-(1/r)*u-(gamma_eff/r)*v));


		k2[0]=h*(u+0.5*k1[2]);
		k2[1]=h*(v+0.5*k1[3]);
		k2[2]=h*(-(r+1.0)*(theta+0.5*k1[0])+r*(omega+0.5*k1[1])-gamma_eff*(u+0.5*k1[2]));
		k2[3]=h*((r+1.0)*(theta+0.5*k1[0])-(r+1)*(omega+0.5*k1[1])+gamma_eff*(1-(1/r)*(u+0.5*k1[2])-(gamma_eff/r)*(v+0.5*k1[3])));


		k3[0]=h*(u+0.5*k2[2]);
		k3[1]=h*(v+0.5*k2[3]);
		k3[2]=h*(-(r+1.0)*(theta+0.5*k2[0])+r*(omega+0.5*k2[1])-gamma_eff*(u+0.5*k2[2]));
		k3[3]=h*((r+1.0)*(theta+0.5*k2[0])-(r+1)*(omega+0.5*k2[1])+gamma_eff*(1-(1/r)*(u+0.5*k2[2])-(gamma_eff/r)*(v+0.5*k2[3])));


		k4[0]=h*(u+k3[2]);
		k4[1]=h*(v+k3[3]);
		k4[2]=h*(-(r+1.0)*(theta+k3[0])+r*(omega+k3[1])-gamma_eff*(u+k3[2]));
		k4[3]=h*((r+1.0)*(theta+k3[0])-(r+1)*(omega+k3[1])+gamma_eff*(1-(1/r)*(u+k3[2])-(gamma_eff/r)*(v+k3[3])));


		theta1 = theta + (1.0/6.0)*(k1[0]+2.0*k2[0]+2.0*k3[0]+k4[0]);
		omega1 = omega + (1.0/6.0)*(k1[1]+2.0*k2[1]+2.0*k3[1]+k4[1]);
		u1 = u + (1.0/6.0)*(k1[2]+2.0*k2[2]+2.0*k3[2]+k4[2]);
		v1 = v+ (1.0/6.0)*(k1[3]+2.0*k2[3]+2.0*k3[3]+k4[3]);


		energy = calculate_energy(r, theta1, omega1, u1, v1);

		outRK << t << "\t" << theta1 << "\t"<< omega1 << "\t" << u1 << "\t" << v1 << "\t" << energy << "\n" ;


		u=u1;
		theta=theta1;
		v=v1;
		omega=omega1;

		}

		return;
}


double calculate_energy(double r, double theta, double omega, double u, double v)
{

	double energy= 0.5*(r*(v+u)*(v+u)+u*u+(1+r)*theta*theta+r*omega*omega);

	return energy;

}

