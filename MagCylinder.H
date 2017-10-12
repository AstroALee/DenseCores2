
#ifndef __MAGCYLINDER_H__
#define __MAGCYLINDER_H__

// NR routines
#include "NRRoutines/nr3.h"
#include "NRRoutines/stepper.h"
#include "NRRoutines/odeint.h"
#include "NRRoutines/stepperdopr5.h"
//#include "stepperross.h"

#include "DCores_Global.H"
#include "TheState.H"


void calcMagCylinder(Params& simP,TheState& curState,int idx);


// Integration works in units of cs = G = 1.
// We are going to find the set of (rho_c,r) such that rho(r)=1 and lambda(code) = lambda(input)

// Derivatives for OD solver
struct derivs {
public:
    derivs(double n, double b)
    {
        ncst = n;
        beta0 = b;
    }

    double ncst;
    double beta0;

    // Derivatives
    void operator() (const double t, VecDoub_I &y, VecDoub_O &dydx) {

        // Fixed parameters?
        //double ncst= nCyl;  //0.5;
        double alp= sqrt(8.0*PI/beta0); // alp^2 = 8*pi*rho^(1-2n)_intfy / beta_infty

        // Derived quantities
        double a2 = alp*alp;
        double gm = a2*ncst/4.0/PI;

        // Use r instead of t (for sanity)
        double r = t;
        double rho = y[0]; // we'll use this often
        double psi = y[1];
        double dpsi = y[2];

        // Second derivative place holder
        double ddpsi;

        // First derivative of r*A
        double drA;

        if(r==0) // Likely derived using L'Hopital
        {
            ddpsi = 2.0*PI*rho;
        }
        else
        {
            ddpsi = 4.0*PI*rho - dpsi/r;
        }


        drA = r*(alp*pow(rho,ncst));

        double drho = -rho*dpsi/(1.0+gm*pow(rho,2.*ncst-1));


        dydx[0] = drho; //y[1]; // d/dr(rho) = rho'
        //dydx[1] = ddrho; // d/dr(rho') = rho''

        dydx[1] = dpsi; // d/dr(psi) = psi'
        dydx[2] = ddpsi; // d/dr(psi') = psi''

        dydx[3] = drA;   // d/dr(r*A) = (r*A)'

        // Mass/length integral of cylinder
        dydx[4] = 2.0*PI*r*rho; // d/dr(lambda) = 2*pi*r*rho

    }

};


// Linearly interpolate on the output to find the particular r value where it occurs
inline double LIntY(VecDoub Xarray, MatDoub Yarray, double point, int size, int elem)
{
    double value = 0;
    int found = 0;
    // Find the two points surrounding 'point' in Yarray[elem][:]
    for(int i=0;i<size-1;i++)
    {
        if( (Yarray[elem][i] <= point && Yarray[elem][i+1] >= point) || (Yarray[elem][i] >= point && Yarray[elem][i+1] <= point)  )
        {
            found = i;

            double y1 = Yarray[elem][i];
            double y2 = Yarray[elem][i+1];
            double x1 = Xarray[i];
            double x2 = Xarray[i+1];
            double m = (y2-y1)/(x2-x1);

            value = x1 + (point-y1)/m;
        }
    }

    if(!found) cout << "Interpolated value was not found!" << endl;
    return value;
};


// Linearly Interpolate an array onto a point. Returns value.
inline double LInt(VecDoub Xarray, MatDoub Yarray, double point, int size,int elem)
{
    double value = 0;
    int found = 0;

    // find the two points surrounding 'point' in the X array
    // Assumes X array is arranged in ascending value
    int i;
    for(i=0;i<size-1;i++) if( Xarray[i+1] >= point)
    {
        double RightY = Yarray[elem][i+1];
        double LeftY  = Yarray[elem][i];
        double RightX = Xarray[i+1];
        double LeftX  = Xarray[i];
        // hi

        double m = (RightY-LeftY)/(RightX-LeftX);

        value = LeftY + m*(point-LeftX);
        found = 1;
        break;
    }

    if(!found) cout << "Interpolated value was not found!" << endl;

    return value;
};


#endif
