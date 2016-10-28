//Algorithm in optimizing the critical step size for stability for LF Method using energy consideration
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

	const double gamma=0.0; //damping constant

	double l=1.0; // length
	double g=1.0; // gravitational constant

	double gamma_eff=(gamma/m)*sqrt(1.0/(l*g)); // effective gamma

	ofstream outfile ("LFStability.txt"); // initiate outfile stream

	double h=1.1; // start with a big step size

	double theta_Inital=0.2;
	double u_Initial=0.0;
	double t=0.0;

	double energy_Initial= energy_check(theta_Inital,u_Initial); // check initial energy

	double tolerance=30.0; // Energy tolerance
	double reducing_factor=1-1e-6; // Reducing stepsize

	double energy;

	double theta;
	double theta1;
	double theta2;

	double u;
	double u1;
	double u2;

	do{

		theta= theta_Inital;
		u=u_Initial;

		theta1= theta+h*u;
                u1= -h*theta+(1.0-h*gamma_eff)*u;


		int imax= static_cast<int>(tmax/h); // using tmax to find max number of steps

			for (int i = 2; i<=imax; i++){

                        t=i*h; //Dimensionless time

                        theta2= theta+2.0*h*u1;
                        u2=u+2.0*h*((-theta1)-gamma_eff*u1);

                        energy = energy_check( theta1, u1);




					if (energy >= tolerance*energy_Initial){
						h=h*reducing_factor;//reduce step size
						cout << "Break"<<"\t"<< i*h << "\t" << h<<endl; // useful for debugging

						break;
						}

					if (i==imax){

					outfile<<"LF Critical Step Alogrithm \n \n \n ";
					outfile<< "The optimal step size is " <<h<<endl;
					outfile << "Tolerence:"<< tolerance<< "\t" << endl ;
					outfile << "Initial_Energy:" << energy_Initial<< endl;
					outfile << "imax ="<< imax<< endl;
					outfile << "Time=" << imax*h <<endl;
					exit(99);
				}

                        theta=theta1;
                        theta1=theta2;

                        u=u1;
                        u1=u2;


		}
	}while (h>0.1);

	return 0;
}


double energy_check (double theta, double u ){

	double energy= (1.0/2.0)*(u*u+theta*theta);

	return energy;

}
