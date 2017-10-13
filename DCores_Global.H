/*
	Global header for the Dense Core code, v2.0
	Includes header inclusions and namespace declarations that are used throughout
	code.
*/

#ifndef __DCORES_GLOBAL_H__
#define __DCORES_GLOBAL_H__

#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <string.h>
#include <math.h>
#include <time.h>

// No'bdy gonna make me type std:: over and over again
using namespace std;

// Constants (wrap in parentheses to avoid subtle arithmetic issues)
#define PI (4.0*atan(1))      // Delicious with whipped cream
#define Gcst (6.67408e-8)     // Gravitational constant, cgs
#define Cbackground (18821.9544) // background sound speed (T=10 K, mu=2.33), cgs

// Unit conversions
#define Pc2Cm (3.0856775e18)
#define Sol2G (1.989e33)


// Enuemrate Declaration (Calculating Rho is easy, so memory is not being allocated for it)
#define NStates 4
enum states { Vpot, Apot, Q, dQdP };


// Definition of Parameter class
class Params
{
public:

	// Simulation parameters
	int M, N; // grid points in each dimension
	double rL, zL; // dimensions of box
	double dR, dZ; // grid size in each dimension

	int convergeLoopMax;
	double convergeLoopTol;

	// Cylinder parameters
	double betaCyl;		// plasma beta at r=0
	double nCyl;		// power-law exponent
	double lambda;		// mass per unit length
	double rTozPhys;
	double rCyl;		// radius of cylidner
	double Rbdy;		// Density at radius of cylinder
	double CylCsts[2];  // Constants needed for analytic Vpot outside cylinder

	// waistband parameters
	double wRatio;
	
	// Point mass parameters
	double mExcess;		// mass excess beyond cylinder
	int pointLoopMax;
	int pointLoopTotNum;
	double pointLoopTol;


	// Output
	string outfname;

	// Units
	double Pc2Code,Sol2Code;


	double *RhoTop;   // density at top, won't change

protected:

private:

};

// Index Generators
// Row major on interior cells : i = Ridx + Zidx * M
inline int Ridx(int i, int M)
{
    return i%(M);
}

inline int Zidx(int i, int M)
{
    return (i-Ridx(i,M))/(M);
}

inline int Sidx(int r, int z, int M )
{
    return ( r + ( M * z ) );
};

// Calculates the position of the cell in one of the directions.
// Feed in the cell size and the cell index. Cell points are at edges (i=0 corresponds to r=0 or z=0)
inline double cPos(int i, double dx) {return ((double) i)*dx;}

// Plural function, because grammar matters
inline string Plural(int i,string s)
{
    if(i==1)
    {
        string blank = "";
        return blank;
    }
    else
    {
        return s;
    }
};



#endif  // __DCORES_GLOBAL_H__
