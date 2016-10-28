//Algorithm in optimizing the critical step size for stability for RK4 Method using energy consideration
// Start with an unstable case and iterate it to a stable solution

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstring>
#include <cstdlib>

using namespace std;

double energy_check (double, double);

int main(){

// Initiate variables

	int tmax= 1e6; // max time

	double m= 1.0;//mass
	const double gamma=0.2; //damping constant
	double l=1.0; // length
	double g=1.0; // gravitational constant

	double gamma_eff=(gamma/m)*sqrt(1.0/(l*g)); // effective gamma

	ofstream outfile ("RK4Stability.txt"); // initiate outfile stream

	double h=3.0; // start with a big step size

	double theta_Inital=0.2;
	double u_Initial=0.0;
	double t=0.0;

	double energy_Initial= energy_check(theta_Inital,u_Initial); // check initial energy

	double tolerance=1.1; // Energy tolerance
	double reducing_factor=1-1e-7; // Reducing stepsize

	double energy;
	double theta;
	double theta1;
	double u;
	double u1;


		vector<double>k1(2);
		vector<double>k2(2);
		vector<double>k3(2);
		vector<double>k4(2);

	do{
		theta= theta_Inital;
		u=u_Initial;

		int imax= static_cast<int>(tmax/h); // using tmax to find max number of steps

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


				energy = energy_check(theta1, u1);

					if (energy >= tolerance*energy_Initial){
						h=h*reducing_factor;//reduce step size
						cout << "Break"<<"\t"<< i*h << "\t" << h<<endl; // useful for debugging

						break;
						}

					if (i==imax){

					outfile<<"RK4 Critical Step Alogrithm \n \n \n ";
					outfile<< "The optimal step size is " <<h<<endl;
					outfile << "Tolerence:"<< tolerance<< "\t" << endl ;
					outfile << "Initial_Energy:" << energy_Initial<< endl;
					outfile << "imax ="<< imax<< endl;
					outfile << "Time=" << imax*h <<endl;
					exit(99);
				}

				theta=theta1;
				u=u1;
		}
	}while (h>0.1);

	return 0;
}


double energy_check (double theta, double u ){

	double energy= (1.0/2.0)*(u*u+theta*theta);

	return energy;

}
