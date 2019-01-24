// 2.5 GHz Intel Core i5, clang-800.0.42.1, macOS High Sierra version 10.13.6

#include <stdio.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>

using namespace std;

// number of equations/objects
// 3 for three_body, 2 for harmonic_osc
int N=3;

// Physical constants for 3-body problem (in MKS units)
const double Me = 5.976e24;
const double Mm = 0.0123*Me;
const double Mm2 = 0.20*Mm;
const double Tm = 648.; // hrs
const double G = 6.674e-11;
const double Re = 6378.0e3;
const double Rm = 3476.0e3;
const double Rm2 = 0.5*Rm;
const double mass[3] = {Me,Mm,Mm2};
const double radius[3] = {Re,Rm,Rm2};

// Dormand-Prince Coefficients
const double c2 = 1./5;
const double c3 = 3./10;
const double c4 = 4./5;
const double c5 = 8./9;
const double c6 = 1.;
const double a21 = 1./5;
const double a31 = 3./40;
const double a32 = 9./40;
const double a41 = 44./45;
const double a42 = -56./15;
const double a43 = 32./9;
const double a51 = 19372./6561;
const double a52 = -25360./2187;
const double a53 = 64448./6561;
const double a54 = -212./729;
const double a61 = 9017./3168;
const double a62 = -355./33;
const double a63 = 46732./5247;
const double a64 = 49./176;
const double a65 = -5103./18656;

// coefficients for 5th order solution
const double b1 = 35./384;
const double b2 = 0.;
const double b3 = 500./1113;
const double b4 = 125./192;
const double b5 = -2187./6784;
const double b6 = 11./84;

// coefficients for 4th order solution
const double b1_prime = 5179./57600;
const double b2_prime = 0.;
const double b3_prime = 7571./16695;
const double b4_prime = 393./640;
const double b5_prime = -92097./339200;
const double b6_prime = 187./2100;

/* yin[0] = position, yin[1] = derivative
 * fout[0] = dy0/dt, fout[1] = dy1/dt
 */
void harmonic_osc(double t, double yin[], double fout[]) {
    fout[0] = yin[1];
    fout[1] = -yin[0];
    return;
}

/*
    0th object = Earth, 1st = moon, 2nd = asteroid
    y[i+xo] = x pos of ith object         xo = 0
    y[i+yo] = y pos of ith object         yo = N
    y[i+vxo] = vx velocity of ith object  vxo = 2*N
    y[i+vyo] = vy velocity of ith object  vyo = 3*N
 
    fout[i+dxo] = dx/dt of ith object     dxo = 0
    fout[i+dyo] = dy/dt of ith object     dyo = N
    fout[i+dvxo] = dvx/dt of ith object   dvxo = 2*N
    fout[i+dvyo] = dvy/dt of ith object   dvyo = 3*N
 */
void three_body(double t, double y[], double fout[]) {
    
    // dXi/dt eqs
    fout[0] = y[2*N];
    fout[1] = y[2*N+1];
    fout[2] = y[2*N+2];

    // dYi/dt eqs
    fout[N] = y[3*N];
    fout[N+1] = y[3*N+1];
    fout[N+2] = y[3*N+2];
    double dx, dy, dist;

    /* dVxi/dt and dVyi/dt
       This nested for loop executes Cowell's method
       computing the x and y accelerations of each ith object
    */
    for (int i=0; i < N; i++) {
        fout[2*N+i] = 0.0;
        fout[3*N+i] = 0.0;
        for (int j=0; j < N; j++) {
            if (i != j) {
                // compute euclidean distance
                dx = y[j] - y[i];
                dy = y[N+j] - y[N+i];
                dist = sqrt(dx*dx + dy*dy);
                // collision detection
                if (dist < 0.25*(radius[i] + radius[j])) {
                    cout << "COLLISION!!! A catastrophic \
			     interstellar collision wiped \
			     all of humanity from \
			     existence. Nice job! \
			     Simulation terminating..." << endl;
                    exit(1);
                }

                // dVxi/dt
                fout[2*N+i] += G*mass[j]*dx/(dist*dist*dist);

                // dVyi/dt
                fout[3*N+i] += G*mass[j]*dy/(dist*dist*dist);
            }
        }
    } // end Cowell's method
       
    return;

} // end three_body

/*
 Adaptive stepsize algorithm
 1. calculate k1,k2,...k6 using dormand prince coefficients
 2. calculate y_n+1 and y*_n+1 from k's
 3. estimate the current relative error
    Ein+1= |((yin+1)-(y*in+1))/scalei|
    scalei ~ |yni|+|h*(dyni/dt)| + |Emin|
    (last term is a safety error term)
 4. find the max error delta_max(|Ein+1|) for i=1,2,...,N
 5. if delta_max << Emax increase
    h=h_opt=h_current*(Emax/delta_max)^(1/5)
 6. if delta_max > Emax then h = h/2 or h/5
    and repeat current step
 7. solve y_n+1 and set t_n+1=t_n+h_used repeat step 1
 */
void adaptive_rk(double yold[], double ynew[], double ynew_prime[], 
		 double h, double Emax, double tinit, double tmax, 
		 int n_eqs, void frhs (double t, double[], double[]), 
		 ofstream& output_file) {
    
    int i, nstep, istep;
    double *k1, *k2, *k3, *k4, *k5, *k6, *k7, *temp, *Error, Ecurrent, tn, h_used;
    
    // pointers to array
    k1 = new double[9*n_eqs];
    k2 = k1 + n_eqs;
    k3 = k2 + n_eqs;
    k4 = k3 + n_eqs;
    k5 = k4 + n_eqs;
    k6 = k5 + n_eqs;
    k7 = k6 + n_eqs;
    temp = k7 + n_eqs;
    Error = temp + n_eqs;
    
    h_used = h;
    nstep = (int) ((tmax-tinit)/h) + 1;
    
    for (istep=0; istep<nstep;  istep++) {
        tn = tinit + istep*h;
        
        // first step
        frhs(tn,yold,k1);
        for (i=0; i< n_eqs; i++) temp[i] = yold[i] + h_used*a21*k1[i];
        
        // second step
        frhs(tn+c2*h_used,temp,k2);
        for (i=0; i< n_eqs; i++) temp[i] = yold[i] + h_used*(a31*k1[i] 
				           + a32*k2[i]);
        
        // third step
        frhs(tn+c3*h_used,temp,k3);
        for (i=0; i< n_eqs; i++) temp[i] = yold[i] + h_used*(a41*k1[i] 
				           + a42*k2[i] + a43*k3[i]);
        
        // fourth step
        frhs(tn+c4*h_used,temp,k4);
        for (i=0; i< n_eqs; i++) temp[i] = yold[i] + h_used*(a51*k1[i] 
					   + a52*k2[i] + a53*k3[i] 
					   + a54*k4[i]);
        
        // fifth step
        frhs(tn+c5*h_used,temp,k5);
        for (i=0; i< n_eqs; i++) temp[i] = yold[i] + h_used*(a61*k1[i] 
					   + a62*k2[i] + a63*k3[i] 
					   + a64*k4[i] + a65*k5[i]);
        
        // sixth step
        frhs(tn+c6*h_used,temp,k6);
        
        // calculate 5th order accurate solution ynew
        for (i=0; i< n_eqs; i++) ynew[i] = yold[i] + h_used*(b1*k1[i] 
					   + b2*k2[i] + b3*k3[i] 
					   + b4*k4[i] + b5*k5[i] 
					   + b6*k6[i]);
        
        /* compute error using 4th and 5th order 
           dormand prince coefficients
        */
        for (i=0; i< n_eqs; i++) Error[i] = h_used*(
				 (b1-b1_prime)*k1[i] 
				+ (b2-b2_prime)*k2[i] 
			        + (b3-b3_prime)*k3[i] 
				+ (b4-b4_prime)*k4[i] 
	                        + (b5-b5_prime)*k5[i] 
			        + (b6-b6_prime)*k6[i]);
                
        // Scaling and computing max computation error
        double epsilon = 1.0e-30; // ensure no divide by zero
        double scale[6];
        for (i=0;i<2*N;i++) {
            scale[i] = (abs(ynew[i]) + abs(ynew[i+2*N]) + epsilon);
        }

        /* Ecurrent cannot be zero when scaling h_used
         (see end of function)
        */
        Ecurrent = abs(Error[0]/scale[0]); 
        
        for (i=1; i< 2*N; i++) {
            if(abs(Error[i]/scale[i]) > Ecurrent) {
                Ecurrent = abs(Error[i]/scale[i]);
            }
        }
        
        // if our error is too high recompute with smaller step size
        while(Ecurrent > Emax) {
            h_used = .5*h_used;

        // first step
        frhs(tn,yold,k1);
        for (i=0; i< n_eqs; i++) temp[i] = yold[i] + h_used*a21*k1[i];
        
        // second step
        frhs(tn+c2*h_used,temp,k2);
        for (i=0; i< n_eqs; i++) temp[i] = yold[i] + h_used*(a31*k1[i] 
				           + a32*k2[i]);
        
        // third step
        frhs(tn+c3*h_used,temp,k3);
        for (i=0; i< n_eqs; i++) temp[i] = yold[i] + h_used*(a41*k1[i] 
				           + a42*k2[i] + a43*k3[i]);
        
        // fourth step
        frhs(tn+c4*h_used,temp,k4);
        for (i=0; i< n_eqs; i++) temp[i] = yold[i] + h_used*(a51*k1[i] 
					   + a52*k2[i] + a53*k3[i] 
					   + a54*k4[i]);
        
        // fifth step
        frhs(tn+c5*h_used,temp,k5);
        for (i=0; i< n_eqs; i++) temp[i] = yold[i] + h_used*(a61*k1[i] 
					   + a62*k2[i] + a63*k3[i] 
					   + a64*k4[i] + a65*k5[i]);
        
        // sixth step
        frhs(tn+c6*h_used,temp,k6);
        
        // calculate 5th order accurate solution ynew
        for (i=0; i< n_eqs; i++) ynew[i] = yold[i] + h_used*(b1*k1[i] 
					   + b2*k2[i] + b3*k3[i] 
					   + b4*k4[i] + b5*k5[i] 
					   + b6*k6[i]);
        
	/* compute error using 4th and 5th order 
	   dormand prince coefficients
	*/        
	for (i=0; i< n_eqs; i++) Error[i] = h_used*(
				 (b1-b1_prime)*k1[i] 
				+ (b2-b2_prime)*k2[i] 
			        + (b3-b3_prime)*k3[i] 
				+ (b4-b4_prime)*k4[i] 
	                        + (b5-b5_prime)*k5[i] 
			        + (b6-b6_prime)*k6[i]);

            
        // Scaling and recalculating max computation error
        for (i=0;i<2*N;i++) {
            scale[i] = (abs(ynew[i]) + abs(ynew[i+2*N]) + epsilon);
        }
        
	/* Ecurrent cannot be zero when scaling h_used 
	   (see end of function)
	 */    
        Ecurrent = abs(Error[0]/scale[0]);             
        for (i=1; i< 2*N; i++) {
            if(abs(Error[i]/scale[i]) > Ecurrent) {
                Ecurrent = abs(Error[i]/scale[i]);
            }
        }
        } // end while loop
        
        /* If our error is very small increase 
            step size for next cycle
        */
        if (0.25>(abs(Ecurrent/Emax))) {
            h = h_used*pow(abs(Emax/Ecurrent),1./5);
        }
        else {
            h = h_used;
        }
        
        // update values for our next step
        for (i=0;i<n_eqs;i++) {
            yold[i] = ynew[i];
        }
        
        /* prevent a very large nstep from wrapping 
            around to a negative value
        */
        if ((int) ((tmax-tinit)/h) < 0) {
            nstep = 1073741824; // 2^30
        }
        else {
            /* for a new timestep we must adjust 
                the number of loop iterations
             */
            nstep = (int) ((tmax-tinit)/h);
        }

        // a change in h and nstep requires a change in istep
        istep = (int) (tn-tinit)/h;
                
        // Three-body output
        output_file << h_used << setw(15) << tn << setw(15) << ynew[0] 
		    << setw(15) << ynew[1] << setw(15) << ynew[2] 
		    << setw(15) << ynew[3] << setw(15) << ynew[4] 
		    << setw(15) << ynew[5] << endl;

        cout << "h: " << h_used << endl;
        cout << "t : " << tn << endl;
        cout << "Xe : " << ynew[0] << endl;
        cout << "Xm : " << ynew[1] << endl;
        cout << "Xm2 : " << ynew[2] << endl;
        cout << "Ye : " << ynew[3] << endl;
        cout << "Ym : " << ynew[4] << endl;
        cout << "Ym2 : " << ynew[5] << endl;

    } // end for loop over nsteps

    delete[] k1; // free memory
    return;
    
} // end adaptive_rk

int main() {
    
    // create/open output files for plotting 
    ofstream ODE;
    
    ODE.open("ODE.dat",ofstream::out);
    if (ODE.fail() || ODE.bad()) {
        cout << "Cannot open file" << endl;
        return -1;
    }

    // 3-BODY Initialization
    N = 3;
    n_eqs = 2*2*N; // n-body problem in 2D has 2*2*N equations
    // initialization
    double yold_3body[] = {0.0,0.0,-4.97e8,0.0,3.84e8,0.0,-12.593,
			   1019.0,965.0,0.0,0.0,820.0};
    tmax = 60*60*24*600.; // 600 days
    tinit = 0.0;
    h = 100.0;
    Emax = 0.1;
    adaptive_rk(yold_3body,ynew,ynew_prime,h,Emax,tinit,
		tmax,n_eqs,three_body,ODE);
    
    ODE.close();
    return(EXIT_SUCCESS);
}
