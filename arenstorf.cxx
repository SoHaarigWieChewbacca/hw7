#include <iostream>
#include <cmath>

using namespace std;

void calc(double x, double x_dot, double y, double y_dot, double* k, const double mu);
void step_size_control(double* xy_RK4, double* xy_RK5, const double tol, double& dt);


int main() {
    const double mu = 0.012277471;
    double t = 0;
    const double T = 17.065216560157;
    double dt = 1e-3;
    const double tol = 1e-5;

    double xy_RK5[4];
    double xy_RK4[4];
    xy_RK4[0] = 0.994; // x
    xy_RK4[1] = 0; // x_dot
    xy_RK4[2] = 0; // y
    xy_RK4[3] = -2.00158510637908; // y_dot

    double k1[4];
    double k2[4];
    double k3[4];
    double k4[4];
    double k5[4];
    double k6[4];
    double k7[4];


    while(t <= T){
	calc(xy_RK4[0],
	xy_RK4[1],
	xy_RK4[2],
	xy_RK4[3], k1, mu);
	calc(xy_RK4[0] + dt*0.2*k1[0],
	xy_RK4[1] + dt*0.2*k1[1],
	xy_RK4[2] + dt*0.2*k1[2],
	xy_RK4[3] + dt*0.2*k1[3], k2, mu);
	calc(xy_RK4[0] + dt*(3.0/40.0*k1[0] + 9.0/40.0*k2[0]),
	xy_RK4[1] + dt*(3.0/40.0*k1[1] + 9.0/40.0*k2[1]),
	xy_RK4[2] + dt*(3.0/40.0*k1[2] + 9.0/40.0*k2[2]),
	xy_RK4[3] + dt*(3.0/40.0*k1[3] + 9.0/40.0*k2[3]), k3, mu);
	calc(xy_RK4[0] + dt*(44.0/45.0*k1[0] - 56.0/15.0*k2[0] + 32.0/9.0*k3[0]),
	xy_RK4[1] + dt*(44.0/45.0*k1[1] - 56.0/15.0*k2[1] + 32.0/9.0*k3[1]),
	xy_RK4[2] + dt*(44.0/45.0*k1[2] - 56.0/15.0*k2[2] + 32.0/9.0*k3[2]),
	xy_RK4[3] + dt*(44.0/45.0*k1[3] - 56.0/15.0*k2[3] + 32.0/9.0*k3[3]), k4, mu);
	calc(xy_RK4[0] + dt*(19372.0/6561.0*k1[0] - 25360.0/2187.0*k2[0] + 64448.0/6561.0*k3[0] - 212.0/729.0*k4[0]),
	xy_RK4[1] + dt*(19372.0/6561.0*k1[1] - 25360.0/2187.0*k2[1] + 64448.0/6561.0*k3[1] - 212.0/729.0*k4[1]),
	xy_RK4[2] + dt*(19372.0/6561.0*k1[2] - 25360.0/2187.0*k2[2] + 64448.0/6561.0*k3[2] - 212.0/729.0*k4[2]),
	xy_RK4[3] + dt*(19372.0/6561.0*k1[3] - 25360.0/2187.0*k2[3] + 64448.0/6561.0*k3[3] - 212.0/729.0*k4[3]), k5, mu);
	calc(xy_RK4[0] + dt*(9017.0/3168.0*k1[0] - 355.0/33.0*k2[0] + 46732.0/5247.0*k3[0] + 49.0/176.0*k4[0] - 5103.0/18656.0*k5[0]),
	xy_RK4[1] + dt*(9017.0/3168.0*k1[1] - 355.0/33.0*k2[1] + 46732.0/5247.0*k3[1] + 49.0/176.0*k4[1] - 5103.0/18656.0*k5[1]),
	xy_RK4[2] + dt*(9017.0/3168.0*k1[2] - 355.0/33.0*k2[2] + 46732.0/5247.0*k3[2] + 49.0/176.0*k4[2] - 5103.0/18656.0*k5[2]),
	xy_RK4[3] + dt*(9017.0/3168.0*k1[3] - 355.0/33.0*k2[3] + 46732.0/5247.0*k3[3] + 49.0/176.0*k4[3] - 5103.0/18656.0*k5[3]), k6, mu);
	calc(xy_RK4[0] + dt*(35.0/384.0*k1[0] + 500.0/1113.0*k3[0] + 125.0/192.0*k4[0] - 2187.0/6784.0*k5[0] + 11.0/84.0*k6[0]),
	xy_RK4[1] + dt*(35.0/384.0*k1[1] + 500.0/1113.0*k3[1] + 125.0/192.0*k4[1] - 2187.0/6784.0*k5[1] + 11.0/84.0*k6[1]),
	xy_RK4[2] + dt*(35.0/384.0*k1[2] + 500.0/1113.0*k3[2] + 125.0/192.0*k4[2] - 2187.0/6784.0*k5[2] + 11.0/84.0*k6[2]),
	xy_RK4[3] + dt*(35.0/384.0*k1[3] + 500.0/1113.0*k3[3] + 125.0/192.0*k4[3] - 2187.0/6784.0*k5[3] + 11.0/84.0*k6[3]), k7, mu);

	for(int j = 0; j < 4; j++) {
	    xy_RK5[j] = xy_RK4[j] + dt*(35.0/384.0*k1[j] + 500.0/1113.0*k3[j] + 125.0/192.0*k4[j] - 2187.0/6784.0*k5[j] + 11.0/84.0*k6[j]);
	    xy_RK4[j] += dt*(5179.0/57600.0*k1[j] + 7571.0/16695.0*k3[j] + 393.0/640.0*k4[j] - 92097.0/339200.0*k5[j] + 187.0/2100.0*k6[j] + 1.0/40.0*k7[j]);
	}

	step_size_control(xy_RK4, xy_RK5, tol, dt);
	t += dt;

	cout << t << "\t" << xy_RK4[0] << "\t" << xy_RK4[1] << "\t" << xy_RK4[2] << "\t" << xy_RK4[3] << "\t" << dt << endl;
    }


    return 0;
}


void calc(double x, double x_dot, double y, double y_dot, double* k, const double mu) {
    double r = sqrt((x + mu)*(x + mu) + y*y);
    double s = sqrt((x - 1 + mu)*(x - 1 + mu) + y*y);

    k[0] = x_dot;
    k[1] = x + 2*y_dot - ((1-mu)*(x+mu))/(r*r*r) - mu*(x-1+mu)/(s*s*s);
    k[2] = y_dot;
    k[3] = y - 2*x_dot - (1-mu)*y/(r*r*r) - mu*y/(s*s*s);

}

void step_size_control(double* xy_RK4, double* xy_RK5, const double tol, double& dt) {
    double max_error = 0;
    double error;

    for(int i = 0; i < 4; i++) {
	error = abs(xy_RK4[i] - xy_RK5[i]);

	if(error > max_error)
	max_error = error;
    }

    dt *= pow(tol/max_error, 1.0/5.0);
}