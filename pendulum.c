// (Loosely) based on code by Anton Hawthorne, 16/03/2016, for PHYC20013

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cpgplot.h>

//Define G as +9.81m/s, the differential equations assume g is a +ve value.
#define G 9.81
#define PI 3.14159

// Define constants of the problem
const double M1 = 1;//in kg
const double M2 = 1;
const double L1 = 0.5;//in m
const double L2 = 0.5;

// Calculates an array of values for theta1, theta2, omega1, and omega2, given initial conditions. Also takes in a step size h and a number of iterations N.
void RungeKutta(double *theta1vals, double *theta2vals, double *omega1vals, double *omega2vals,
        double theta1_0, double theta2_0, double omega1_0, double omega2_0,
        double h, double N);

//Throughout: the variables labelled 1 belong to the first pendulum, and those labelled 2 correspond to the pendulum that is attached to the first.
//so total extended length is L1+L2. the angle theta 2 is the angle from vertical, does not depend on theta1 inherently (although its dynamics do)
double dtheta1_dt(double theta1, double theta2, double omega1, double omega2) {
    return omega1;
}

double dtheta2_dt(double theta1, double theta2, double omega1, double omega2) {
    return omega2;
}

//The right hand side of the second order differential equation: returns d(omega1)/dt = d2(theta1)/dt2 given the input masses, lengths, angles and angular velocities
double domega1_dt(double theta1, double theta2, double omega1, double omega2){
    double a = M2*L2*pow(omega2,2)*sin(theta2-theta1);
    double b = -(M1+M2)*G*sin(theta1);
    double c = M2*cos(theta2-theta1)*(L1*pow(omega1,2)*sin(theta2-theta1) -G*sin(theta2));
    double d = M1 + M2*pow(sin(theta2-theta1),2);
    return (a + b + c) / d;
}

//returns d(omega2)/dt = d2(theta2)/dt2 given the input masses, lengths, angles and angular velocities
double domega2_dt(double theta1, double theta2, double omega1, double omega2){
    double a = -(M1+M2)*(L1*pow(omega2,2)*sin(theta2-theta1) + G*sin(theta2));
    double b = cos(theta2-theta1)*((M1+M2)*G*sin(theta1) - M2*L2*pow(omega2,2)*sin(theta2-theta1));
    double c = L2*(M1+M2*pow(sin(theta2-theta1),2));
    return (a + b) / c;
}

int main() {
    //the initial angles and angular velocities
    double perturbation = 0;
    double theta1Initial = PI/2 + perturbation;//in radians
    double theta2Initial = PI/2 + perturbation;
    double omega1Initial = 0;//in radians per sec
    double omega2Initial = 0;
    
    // Notify whether the initial values satisfy the flip condition.
    char *result = ((3 * cos(theta1Initial) + cos(theta2Initial)) > 2) ? "not " : "";
    printf("It is %spossible for the pendula to flip.\n", result);

    double timeInterval = 10;//in seconds
    int stepNo = 10000;
    double stepSize = timeInterval / stepNo;

    // allocate arrays
    double *theta1vals = calloc(stepNo + 1, sizeof(double));
    double *theta2vals = calloc(stepNo + 1, sizeof(double));
    double *omega1vals = calloc(stepNo + 1, sizeof(double));
    double *omega2vals = calloc(stepNo + 1, sizeof(double));

    float *x1vals = calloc(stepNo + 1, sizeof(float));
    float *y1vals = calloc(stepNo + 1, sizeof(float));
    float *x2vals = calloc(stepNo + 1, sizeof(float));
    float *y2vals = calloc(stepNo + 1, sizeof(float));
    float *x1trailvals = calloc(stepNo + 1, sizeof(float));
    float *y1trailvals = calloc(stepNo + 1, sizeof(float));
    float *x2trailvals = calloc(stepNo + 1, sizeof(float));
    float *y2trailvals = calloc(stepNo + 1, sizeof(float));
    

    const int stepsPerIteration = 20;
    int count = 0;
    char fname[256];
    for (int N = 0; N < stepNo; N += stepsPerIteration) {
        RungeKutta(theta1vals, theta2vals, omega1vals, omega2vals, theta1Initial, theta2Initial, omega1Initial, omega2Initial, stepSize, N);
        // calculate origin points
        for (int i = 0; i <= N; ++i) {
            x1vals[i] = 0;
            y1vals[i] = 0;
            x2vals[i] = L1 * sin(theta1vals[N]);
            y2vals[i] = -L1 * cos(theta1vals[N]);
        }
        // calculate end points
        for (int i = N - stepsPerIteration; i <= N; ++i) {
            x1vals[i] = L1 * sin(theta1vals[N]);
            y1vals[i] = -L1 * cos(theta1vals[N]);
            x2vals[i] = L1 * sin(theta1vals[N]) + L2 * sin(theta2vals[N]);
            y2vals[i] = -L1 * cos(theta1vals[N]) - L2 * cos(theta2vals[N]);
        }
        // calculate trail points
        for (int i = 0; i <= N; ++i) {
            x1trailvals[i] = L1 * sin(theta1vals[i]);   
            y1trailvals[i] = -L1 * cos(theta1vals[i]);
            x2trailvals[i] = L1 * sin(theta1vals[i]) + L2 * sin(theta2vals[i]);
            y2trailvals[i] = -L1 * cos(theta1vals[i]) - L2 * cos(theta2vals[i]);
        }

        // plot results -- sorry for the disgusting but it's a FORTRAN API, what am I gonna do
        // create the window
        cpgbeg(0,"/PNG",1,1);
        cpgsls(1);
        cpgsch(1.);
        cpgswin(-1,1,-1,1);
        cpgbox("BCNST", 0.0, 0, "BCNST", 0.0, 0);
        cpgmtxt("B",3,.5,.5,"x axis");
        cpgmtxt("L",3,.5,.5,"y axis");
        cpgsch(2.);
        cpgmtxt("T",1,.5,.5,"Double Pendulum");
        // draw the first pendulum's trajectory
        cpgscr(4, 0, 0.5, 0);
        cpgsls(4);
        cpgsci(4);
        cpgline(N, x1trailvals, y1trailvals);
        // draw the second pendulum's trajectory
        cpgscr(15, 0.5, 0.5, 0.5);
        cpgsci(15);
        cpgline(N, x2trailvals, y2trailvals);
        // draw the first pendulum
        cpgsls(1);
        cpgsci(3);
        cpgline(N, x1vals, y1vals);
        // draw the second pendulum
        cpgsci(1);
        cpgline(N, x2vals, y2vals);
        cpgend();
        sprintf(fname, "output/pgplot%d.png", count++);
        rename("pgplot.png", fname);
    }
    return 0;
}
void RungeKutta(double *theta1vals, double *theta2vals, double *omega1vals, double *omega2vals,
        double theta1_0, double theta2_0, double omega1_0, double omega2_0,
        double h, double N) {
    theta1vals[0] = theta1_0;
    theta2vals[0] = theta2_0;
    omega1vals[0] = omega1_0;
    omega2vals[0] = omega2_0;

    for (int n = 0; n < N; ++n) {
        double j1 = dtheta1_dt(theta1vals[n], theta2vals[n], omega1vals[n], omega2vals[n]);
        double k1 = dtheta2_dt(theta1vals[n], theta2vals[n], omega1vals[n], omega2vals[n]);
        double l1 = domega1_dt(theta1vals[n], theta2vals[n], omega1vals[n], omega2vals[n]);
        double m1 = domega2_dt(theta1vals[n], theta2vals[n], omega1vals[n], omega2vals[n]);

        double j2 = dtheta1_dt(theta1vals[n] + h/2 * j1, theta2vals[n] + h/2 * k1, omega1vals[n] + h/2 * l1, omega2vals[n] + h/2 * m1);
        double k2 = dtheta2_dt(theta1vals[n] + h/2 * j1, theta2vals[n] + h/2 * k1, omega1vals[n] + h/2 * l1, omega2vals[n] + h/2 * m1);
        double l2 = domega1_dt(theta1vals[n] + h/2 * j1, theta2vals[n] + h/2 * k1, omega1vals[n] + h/2 * l1, omega2vals[n] + h/2 * m1);
        double m2 = domega2_dt(theta1vals[n] + h/2 * j1, theta2vals[n] + h/2 * k1, omega1vals[n] + h/2 * l1, omega2vals[n] + h/2 * m1);

        double j3 = dtheta1_dt(theta1vals[n] + h/2 * j2, theta2vals[n] + h/2 * k2, omega1vals[n] + h/2 * l2, omega2vals[n] + h/2 * m2);
        double k3 = dtheta2_dt(theta1vals[n] + h/2 * j2, theta2vals[n] + h/2 * k2, omega1vals[n] + h/2 * l2, omega2vals[n] + h/2 * m2);
        double l3 = domega1_dt(theta1vals[n] + h/2 * j2, theta2vals[n] + h/2 * k2, omega1vals[n] + h/2 * l2, omega2vals[n] + h/2 * m2);
        double m3 = domega2_dt(theta1vals[n] + h/2 * j2, theta2vals[n] + h/2 * k2, omega1vals[n] + h/2 * l2, omega2vals[n] + h/2 * m2);

        double j4 = dtheta1_dt(theta1vals[n] + h * j3, theta2vals[n] + h * k3, omega1vals[n] + h * l3, omega2vals[n] + h * m3);
        double k4 = dtheta2_dt(theta1vals[n] + h * j3, theta2vals[n] + h * k3, omega1vals[n] + h * l3, omega2vals[n] + h * m3);
        double l4 = domega1_dt(theta1vals[n] + h * j3, theta2vals[n] + h * k3, omega1vals[n] + h * l3, omega2vals[n] + h * m3);
        double m4 = domega2_dt(theta1vals[n] + h * j3, theta2vals[n] + h * k3, omega1vals[n] + h * l3, omega2vals[n] + h * m3);

        theta1vals[n + 1] = theta1vals[n] + h/6 * (j1 + 2*j2 + 2*j3 + j4);
        theta2vals[n + 1] = theta2vals[n] + h/6 * (k1 + 2*k2 + 2*k3 + k4);
        omega1vals[n + 1] = omega1vals[n] + h/6 * (l1 + 2*l2 + 2*l3 + l4);
        omega2vals[n + 1] = omega2vals[n] + h/6 * (m1 + 2*m2 + 2*m3 + m4);
    }
}
