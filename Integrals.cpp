
#include "Integrals.H"

// For some reason this has to be included here and not in the header file ... not sure why
#include <Eigen/Dense>
using namespace Eigen;

void CylinderMassCheck(Params simP, TheState curState)
{
    double curMass =  Trap2D(simP,curState);
    //double curMass =  MCInt(simP,curState);
    //double curMass =  MCIntSpline(simP,curState);
    //double curMass =  DirectIntSpline(simP,curState);
    double cylMass = 2.0*simP.zL*simP.lambda;
    cout << "Current mass in the box is " << curMass/simP.Sol2Code << " Sol" << endl;
    cout << "Cylinder mass in the box is " << cylMass/simP.Sol2Code << " Sol" << endl;
    cout << "   Difference is " << curMass/simP.Sol2Code - cylMass/simP.Sol2Code << " Sol" << endl;
    cout << "   Error (if mExcess=0) is " << 100*(curMass-cylMass)/cylMass << " %" << endl;
};


// Evaluates r*rho at (rPoint,zPoint)
double CubicSpline(double rPoint, double zPoint, Params simP, TheState curState)
{
    double val = 0;

    int M = ( (int) floor(rPoint) );
    int N = ( (int) floor(zPoint) );
    rPoint -= ((double) M);
    zPoint -= ((double) N); // these are now a number between 0 and 1

    MatrixXd coef1(4,4),coef2(4,4),funcM(4,4);
    VectorXd  rvec(4), zvec(4);

    rvec << 1.0 , pow(rPoint,1) , pow(rPoint,2) , pow(rPoint,3);
    zvec << 1.0 , pow(zPoint,1) , pow(zPoint,2) , pow(zPoint,3);

    coef1 << 1,0,0,0,
             0,0,1,0,
             -3,3,-2,-1,
             2,-2,1,1;
    coef2 << 1,0,-3,2,
             0,0,3,-2,
             0,1,-2,1,
             0,0,-1,1;
    funcM << getF(M,N,simP,curState,0) , getF(M,N+1,simP,curState,0) , getF(M,N,simP,curState,2) , getF(M,N+1,simP,curState,2),
             getF(M+1,N,simP,curState,0) , getF(M+1,N+1,simP,curState,0) , getF(M+1,N,simP,curState,2) , getF(M+1,N+1,simP,curState,2),
             getF(M,N,simP,curState,1) , getF(M,N+1,simP,curState,1) , getF(M,N,simP,curState,3) , getF(M,N+1,simP,curState,3),
             getF(M+1,N,simP,curState,1) , getF(M+1,N+1,simP,curState,1) , getF(M+1,N,simP,curState,3) , getF(M+1,N+1,simP,curState,3) ;

    val = rvec.transpose() * coef1 * funcM * coef2 * zvec ;




    return val;
};

double getF(int M,int N,Params simP,TheState curState, int which)
{
    double val = 0;
    double DeltaR = simP.dR; double DeltaZ = simP.dZ;

    double vL = 0; double vR = 0;
    double vT = 0; double vB = 0;

    int outside = 0;
    if( cPos(M,DeltaR) >= curState.VContour[N] ) outside = 1;

    switch(which)
    {
        case 0: // gets value r*rho
            val = cPos(M,DeltaR) * curState.State[Q][M][N] * exp( -curState.State[Vpot][M][N] );
            if(outside) val = 0;
            break;
        case 1: // gets d/dr ( r*rho )
            if(M==0)
            {
                val = curState.State[Q][M][N] * exp( -curState.State[Vpot][M][N] ); // rho(r=0)
            }
            else if(M==simP.M-1) // outside filament
            {
                val = 0;
            }
            else
            {
                vR = 0; vL = 0;
                if( cPos(M+1,DeltaR) < curState.VContour[N] ) vR = cPos(M+1,DeltaR) * curState.State[Q][M+1][N] * exp( -curState.State[Vpot][M+1][N] );
                if( cPos(M-1,DeltaR) < curState.VContour[N] ) vL = cPos(M-1,DeltaR) * curState.State[Q][M-1][N] * exp( -curState.State[Vpot][M-1][N] );

                val = (vR-vL)/2.0/DeltaR;
            }
            break;
        case 2: // gets d/dz (r*rho)
            if(N==0 or N==simP.N-1)
            {
                val = 0;
            }
            else
            {
                vT=0;vB=0;
                if( cPos(M,DeltaR) < curState.VContour[N+1] ) vT = cPos(M,DeltaR) * curState.State[Q][M][N+1] * exp( -curState.State[Vpot][M][N+1] );
                if( cPos(M,DeltaR) < curState.VContour[N-1] ) vB = cPos(M,DeltaR) * curState.State[Q][M][N-1] * exp( -curState.State[Vpot][M][N-1] );

                val = (vT-vB)/2.0/DeltaZ;
            }
            break;
        case 3: // gets d/dzdr (r*rho)
            if(N==0 or N==simP.N-1 or M==simP.M-1)
            {
                val = 0; // d/dz = 0 so d/dr(d/dz) = 0
            }
            else if(M==0)
            {
                vT = curState.State[Q][M][N+1] * exp( -curState.State[Vpot][M][N+1] ); // rho(r=0) values
                vB = curState.State[Q][M][N-1] * exp( -curState.State[Vpot][M][N-1] );

                val = (vT-vB)/2.0/DeltaZ;
            }
            else  // d/drdz(r rho) = drho/dz + r * d^2 rho/dzdr
            {
                // drho/dz
                vT=0;vB=0;
                if( cPos(M,DeltaR) < curState.VContour[N+1] ) vT = curState.State[Q][M][N+1] * exp( -curState.State[Vpot][M][N+1] ); // rho values
                if( cPos(M,DeltaR) < curState.VContour[N-1] ) vB = curState.State[Q][M][N-1] * exp( -curState.State[Vpot][M][N-1] );
                double drhodz = (vT-vB)/2.0/DeltaZ;

                // d^2 rho/dzdr
                vR=0;vL=0;
                if( cPos(M+1,DeltaR) < curState.VContour[N] ) vR = curState.State[Q][M+1][N] * exp( -curState.State[Vpot][M+1][N] );
                if( cPos(M-1,DeltaR) < curState.VContour[N] ) vL = curState.State[Q][M-1][N] * exp( -curState.State[Vpot][M-1][N] );
                double ddrhodzdr = (vT + vR - vL - vB)/4.0/DeltaR/DeltaZ;

                val = drhodz + ( cPos(M,DeltaR) * ddrhodzdr );

            }
            break;
        default:
            cout << "There's a problem...." << endl;
    } // end of switch


    return val;
}

double DirectIntSpline(Params simP, TheState curState)
{
    double val = 0;

    // loops over each cell, calculates the cubic spline polynomial, and directly
    // integrates it over the entire cell.


    MatrixXd coef1(4,4),coef2(4,4),funcM(4,4);
    VectorXd  rvec(4), zvec(4);

    // these don't change, only have to create once
    coef1 << 1,0,0,0,
             0,0,1,0,
             -3,3,-2,-1,
             2,-2,1,1;
    coef2 << 1,0,-3,2,
             0,0,3,-2,
             0,1,-2,1,
             0,0,-1,1;

    for(int M=0;M<simP.M-1;M++) for(int N=0;N<simP.N-1;N++)
    {
        double rPoint = cPos(M,simP.dR);
        double zPoint = cPos(N,simP.dZ);

        // integrated vectors
        rvec << (1./1.)*(pow(rPoint+simP.dR,1)-pow(rPoint,1)) ,
                (1./2.)*(pow(rPoint+simP.dR,2)-pow(rPoint,2)) ,
                (1./3.)*(pow(rPoint+simP.dR,3)-pow(rPoint,3)) ,
                (1./4.)*(pow(rPoint+simP.dR,4)-pow(rPoint,4));

        zvec << (1./1.)*(pow(zPoint+simP.dZ,1)-pow(zPoint,1)) ,
                (1./2.)*(pow(zPoint+simP.dZ,2)-pow(zPoint,2)) ,
                (1./3.)*(pow(zPoint+simP.dZ,3)-pow(zPoint,3)) ,
                (1./4.)*(pow(zPoint+simP.dZ,4)-pow(zPoint,4));

        funcM << getF(M,N,simP,curState,0) , getF(M,N+1,simP,curState,0) , getF(M,N,simP,curState,2) , getF(M,N+1,simP,curState,2),
                 getF(M+1,N,simP,curState,0) , getF(M+1,N+1,simP,curState,0) , getF(M+1,N,simP,curState,2) , getF(M+1,N+1,simP,curState,2),
                 getF(M,N,simP,curState,1) , getF(M,N+1,simP,curState,1) , getF(M,N,simP,curState,3) , getF(M,N+1,simP,curState,3),
                 getF(M+1,N,simP,curState,1) , getF(M+1,N+1,simP,curState,1) , getF(M+1,N,simP,curState,3) , getF(M+1,N+1,simP,curState,3) ;

        val = val + rvec.transpose() * (coef1 * funcM * coef2) * zvec ;
    }

    return val;
}


double MCIntSpline(Params simP, TheState curState)
{
    double sum = 0;

    int Nmax = 100*simP.M*simP.N;

    for(int n=0;n<Nmax;n++)
    {
        double randR = ((double) simP.M-1.0) * ( (double) rand() )/( (double) RAND_MAX );
        double randZ = ((double) simP.N-1.0) * ( (double) rand() )/( (double) RAND_MAX );

        double locval = CubicSpline(randR, randZ, simP, curState);

        // evaluates interpolated value of r*rho and adds to sum
        sum = sum + locval ;
    }
    sum = sum / ( (double) Nmax ); // divides by number of evaluations to get average value of rho*r

    sum = 2.0 * 2.0*PI * simP.zL * simP.rL * sum ; // 2 for 2zL, 2pi for axisym, rL*zL for Monte Carlo
    return sum;

}


double MCInt(Params simP, TheState curState)
{
    double sum = 0;

    int Nmax = 1000*simP.M*simP.N;
    int Nalter = 0;

    for(int n=0;n<Nmax;n++)
    {
        double randR = ((double) simP.M-1.0) * ( (double) rand() )/( (double) RAND_MAX );
        double randZ = ((double) simP.N-1.0) * ( (double) rand() )/( (double) RAND_MAX );
        int Mleft = ( (int) floor(randR) );
        int Nleft = ( (int) floor(randZ) );

        // Bilinear interpoltion of four corners Qrz
        double Qv[2][2];
        Qv[0][0] = cPos(Mleft,simP.dR);
        Qv[0][1] = cPos(Mleft,simP.dR);
        Qv[1][0] = cPos(Mleft+1,simP.dR);
        Qv[1][1] = cPos(Mleft+1,simP.dR); // the Q's are now just Q = r
        if(Qv[0][0] >= curState.VContour[Nleft]) Qv[0][0] = 0;
        else Qv[0][0] = Qv[0][0] * curState.State[Q][Mleft][Nleft] * exp(-curState.State[Vpot][Mleft][Nleft]);
        if(Qv[0][1] >= curState.VContour[Nleft+1]) Qv[0][1] = 0;
        else Qv[0][1] = Qv[0][1] * curState.State[Q][Mleft][Nleft+1] * exp(-curState.State[Vpot][Mleft][Nleft+1]);
        if(Qv[1][0] >= curState.VContour[Nleft]) Qv[1][0] = 0;
        else Qv[1][0] = Qv[1][0] * curState.State[Q][Mleft+1][Nleft] * exp(-curState.State[Vpot][Mleft+1][Nleft]);
        if(Qv[1][1] >= curState.VContour[Nleft+1]) Qv[1][1] = 0;
        else Qv[1][1] = Qv[1][1] * curState.State[Q][Mleft+1][Nleft+1] * exp(-curState.State[Vpot][Mleft+1][Nleft+1]); // Q = r*rho

        //if(Qv[0][0]==0 or Qv[0][1]==0 or Qv[1][0]==0 or Qv[1][1]==0 ) { Nalter += 1; continue; }

        double Rvec[2], Zvec[2];
        Rvec[0] = ((double) Mleft) + 1 - randR;
        Rvec[1] = randR - ((double) Mleft) ;
        Zvec[0] = ((double) Nleft) + 1 - randZ;
        Zvec[1] = randZ - ((double) Nleft) ;

        double Qr0 = Qv[0][0]*Zvec[0] + Qv[0][1]*Zvec[1];
        double Qr1 = Qv[1][0]*Zvec[0] + Qv[1][1]*Zvec[1];

        // evaluates interpolated value of r*rho and adds to sum
        sum = sum + ( Rvec[0]*Qr0 + Rvec[1]*Qr1 ) ;

    }
    sum = sum / ( (double) Nmax - (double) Nalter); // divides by number of evaluations to get average value of rho*r

    sum = 2.0 * 2.0*PI * simP.zL * simP.rL * sum ; // 2 for 2zL, 2pi for axisym, rL*zL for Monte Carlo
    return sum;
};

double Trap2D(Params simP, TheState curState, double s0) // Simpson
{
    double sum = 0;

    for(int i = 0; i < simP.M ; i++) for(int j = 0; j< simP.N ; j++)
    {
        double rz = curState.VContour[j];
        double sFunc = 1.0 + s0*(cPos(j,simP.dZ) - simP.zL)*(cPos(i,simP.dR) - rz);
        double curV = sFunc*curState.State[Vpot][i][j];
        double curRho = curState.State[Q][i][j] * exp(-curV);
        if(cPos(i,simP.dR)>curState.VContour[j]) { curRho = 0.0; continue; }

        double curInt = curRho * cPos(i,simP.dR);

        double curWeight = 1.0;
        if( (i==0 or i==simP.M-1) and (j==0 or j==simP.N-1) )
            {
                curWeight = 1.0;
            }
        else if( (i==0) or (i==simP.M-1)  ) // won't be called if previous line is called
            {
                if(j%2==0) curWeight = 2.0;
                else curWeight = 4.0;
            }
        else if( (j==0) or (j==simP.N-1) ) // won't be called if previous line is called
            {
                if(i%2==0) curWeight = 2.0;
                else curWeight = 4.0;
            }
        else if( i%2==0 )
            {
                if(j%2==0) curWeight = 4.0;
                else curWeight = 8.0;
            }
        else
            {
                if(j%2==0) curWeight = 8.0;
                else curWeight = 16.0;
            }

        sum = sum + ( curWeight * curInt );
    }

    sum = 2.0 * 2.0*PI * sum ; // 2 for z integral, 2pi for theta integral

    sum = (1.0/9.0) * simP.dR * simP.dZ * sum;

    return sum;
};


double Trap2D(Params simP, TheState curState) // Simpson
{
    double sum = 0;

    for(int i = 0; i < simP.M ; i++) for(int j = 0; j< simP.N ; j++)
    {
        double curRho = curState.State[Q][i][j] * exp(-curState.State[Vpot][i][j]);
        if(cPos(i,simP.dR)>curState.VContour[j]) { curRho = 0.0; continue; }

        double curInt = curRho * cPos(i,simP.dR);

        double curWeight = 1.0;
        if( (i==0 or i==simP.M-1) and (j==0 or j==simP.N-1) )
            {
                curWeight = 1.0;
            }
        else if( (i==0) or (i==simP.M-1)  ) // won't be called if previous line is called
            {
                if(j%2==0) curWeight = 2.0;
                else curWeight = 4.0;
            }
        else if( (j==0) or (j==simP.N-1) ) // won't be called if previous line is called
            {
                if(i%2==0) curWeight = 2.0;
                else curWeight = 4.0;
            }
        else if( i%2==0 )
            {
                if(j%2==0) curWeight = 4.0;
                else curWeight = 8.0;
            }
        else
            {
                if(j%2==0) curWeight = 8.0;
                else curWeight = 16.0;
            }

        sum = sum + ( curWeight * curInt );
    }

    sum = 2.0 * 2.0*PI * sum ; // 2 for z integral, 2pi for theta integral

    sum = (1.0/9.0) * simP.dR * simP.dZ * sum;

    return sum;
}


double Trap2D_wrong(Params simP, TheState curState)
{
    double sum = 0;
    // Given the density discontinuity... perhaps a Monte Carlo integration would be better
    // Sum = 0.25 * DeltaR * DeltaZ * Sum( S_mn * I_mn )
    //       [ 1  2 ----  2  1 ]
    //       [ 2  4 ----  4  2 ]
    // Smn = [ |  |  \\   4  2 ]
    //       [ 2  4 ----  4  2 ]
    //       [ 1  2 ----  2  1 ]

    // int(rho dV) = int( rho * r dr dtheta dz) = 2pi * int(rho * r dr dz; r=0,rL; z=-zL,zL)
    //             = 4*pi * int(rho*r dr dz; r=0,rL; z=0,zL)
    // I_mn = rho*r
    for(int i = 0; i < simP.M-1 ; i++) for(int j = 0; j< simP.N-1 ; j++) // one less cell being used
    {
        double Rho[2][2];
        double rhoweight = 0.25;
        Rho[0][0] = curState.State[Q][i][j] * exp(-curState.State[Vpot][i][j]);
        Rho[0][1] = curState.State[Q][i][j+1] * exp(-curState.State[Vpot][i][j+1]);
        Rho[1][0] = curState.State[Q][i+1][j] * exp(-curState.State[Vpot][i+1][j]);
        Rho[1][1] = curState.State[Q][i+1][j+1] * exp(-curState.State[Vpot][i+1][j+1]);
        for(int ii=0;ii<2;ii++) for(int jj=0;jj<2;jj++)
            if( cPos(i+ii,simP.dR) > curState.VContour[j+jj] ) {Rho[ii][jj] = 0.0; rhoweight += 0.25;}

        // Way that we average these four points
        double curRho = rhoweight * (Rho[0][0] + Rho[0][1] + Rho[1][0] + Rho[1][1]) ;

        double curInt = curRho * (cPos(i,simP.dR) + 0.5*simP.dR);

        double curWeight = 1.0;
        if( (i==0 or i==simP.M-2) and (j==0 or j==simP.N-2) )
            {curWeight = 1.0;} // cout << "i,j = " << i << "," << j << endl;}
        else if( (i==0) or (i==simP.M-2) or (j==0) or (j==simP.N-2) ) // won't be called if previous line is called
            {curWeight = 2.0;} //cout << "edge value: i,j = " << i << "," << j << endl;}
        else
            {curWeight = 4.0;}

        sum = sum + curWeight * curInt;
    }

    sum = 2.0 * 2.0*PI * sum ;

    sum = 0.25 * simP.dR * simP.dZ * sum;



    return sum;
}

double Trap2D_old(Params simP, TheState curState)
{
    double sum = 0;

    // Given the density discontinuity... perhaps a Monte Carlo integration would be better
    // Sum = 0.25 * DeltaR * DeltaZ * Sum( S_mn * I_mn )
    //       [ 1  2 ----  2  1 ]
    //       [ 2  4 ----  4  2 ]
    // Smn = [ |  |  \\   4  2 ]
    //       [ 2  4 ----  4  2 ]
    //       [ 1  2 ----  2  1 ]

    // int(rho dV) = int( rho * r dr dtheta dz) = 2pi * int(rho * r dr dz; r=0,rL; z=-zL,zL)
    //             = 4*pi * int(rho*r dr dz; r=0,rL; z=0,zL)
    // I_mn = rho*r
    for(int i = 0; i < simP.M ; i++) for(int j = 0; j< simP.N ; j++)
    {
        double curRho = curState.State[Q][i][j] * exp(-curState.State[Vpot][i][j]);
        if(cPos(i,simP.dR) > curState.VContour[j]) curRho = 0;

        double curInt = curRho * cPos(i,simP.dR);

        double curWeight = 1.0;
        if( (i==0 or i==simP.M-1) and (j==0 or j==simP.N-1) )
            {curWeight = 1.0;} // cout << "i,j = " << i << "," << j << endl;}
        else if( (i==0) or (i==simP.M-1) or (j==0) or (j==simP.N-1) ) // won't be called if previous line is called
            {curWeight = 2.0;} //cout << "edge value: i,j = " << i << "," << j << endl;}
        else
            {curWeight = 4.0;}

        sum = sum + ( curWeight * curInt );
    }

    sum = 2.0 * 2.0*PI * sum ;

    sum = 0.25 * simP.dR * simP.dZ * sum;

    return sum;
};



double Simpson2D(Params simP, TheState curState)
{
    cout << "Should not be using this..." << endl;

    double val = 0;



    return val;
};
