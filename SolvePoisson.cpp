

#include "SolvePoisson.H"


// For some reason this has to be included here and not in the header file ... not sure why
#include <Eigen/Dense>
using namespace Eigen;

void createPoissonVFixedMatrix(Params simP, TheState curState, MatrixXd& PoMatrix, VectorXd& Source);
void createPoissonCstMatrix(Params simP, TheState curState, MatrixXd& PoMatrix, VectorXd& Source);
void createPoissonMatrixFixedDen(Params simP, TheState curState, MatrixXd& PoMatrix, VectorXd& Source);
void createPoissonVFixedTopBottomRightMatrix(Params simP, TheState curState, MatrixXd& PoMatrix, VectorXd& Source, VectorXd Soln_prev, int which);

void createPoissonVFixedTopRightMatrixS0(Params simP, TheState curState, MatrixXd& PoMatrix, VectorXd& Source, double s0);



void SolvePoissonS0(Params simP,TheState curState,TheState& newState, double s0)
{
    int M = curState.M;
    int N = curState.N;

    // First allocate the matrix and vectors for the finite difference scheme, source, and solution
    MatrixXd PoMatrix = MatrixXd::Zero((M)*(N),(M)*(N));
    VectorXd Source = VectorXd::Zero((M)*(N));
    VectorXd Soln = VectorXd::Zero((M)*(N));

    cout << "Creating Poisson matrix." << endl;
    createPoissonVFixedTopRightMatrixS0( simP,  curState, PoMatrix,  Source, s0);

    // Now call LinAlgebra package to invert the matrix and obtain the new potential values
    //Soln = PoMatrix.householderQr().solve(Source);
    //Soln = PoMatrix.colPivHouseholderQr().solve(Source);
    cout << "Inverting Poisson matrix." << endl;
    Soln = PoMatrix.colPivHouseholderQr().solve(Source);

    // Test to make sure solution is actually a solution (Uses L2 norm)
    double L2error = (PoMatrix*Soln - Source).norm() / Source.norm();
    cout << "Gravity Matrix L2 error = " << L2error << endl;

    // Copy the solution to the newState vector
    for(int i=0;i<M;i++) for(int j=0;j<N;j++) newState.State[Vpot][i][j] = Soln[ Sidx(i,j,M) ];

    // Vpot has been updated, re-normalize so V = 0 in the upper-left corner
    // With a fixed V=0 filament boundary, this renormalization is worked into the scheme. No need to renormalize
    // double Vcorner = newState.State[Vpot][0][N-1];
    // for(int i=0;i<M;i++) for(int j=0;j<N;j++) newState.State[Vpot][i][j] = newState.State[Vpot][i][j] - Vcorner;

};



void SolvePoisson(Params simP,TheState curState,TheState& newState)
{
    int M = curState.M;
    int N = curState.N;

    // First allocate the matrix and vectors for the finite difference scheme, source, and solution
    MatrixXd PoMatrix = MatrixXd::Zero((M)*(N),(M)*(N));
    VectorXd Source = VectorXd::Zero((M)*(N));
    VectorXd Soln = VectorXd::Zero((M)*(N));

    cout << "Creating Poisson matrix." << endl;
    //createPoissonCstMatrix(simP,curState,PoMatrix,Source);  // Solves Poisson everywhere
    //createPoissonMatrixFixedDen(simP,curState,PoMatrix,Source);  // Solves Poisson everywhere, fixed central density
    //createPoissonVFixedMatrix(simP,curState,PoMatrix,Source); // Doesn't really solve Poisson outside of the filament
    createPoissonVFixedTopBottomRightMatrix(simP,curState,PoMatrix,Source,Soln,0);

    // Now call LinAlgebra package to invert the matrix and obtain the new potential values
    //Soln = PoMatrix.householderQr().solve(Source);
    //Soln = PoMatrix.colPivHouseholderQr().solve(Source);
    cout << "Inverting Poisson matrix." << endl;
    Soln = PoMatrix.colPivHouseholderQr().solve(Source);

    // Test to make sure solution is actually a solution (Uses L2 norm)
    double L2error = (PoMatrix*Soln - Source).norm() / Source.norm();
    cout << "Gravity Matrix L2 error = " << L2error << endl;

    // Copy the solution to the newState vector
    for(int i=0;i<M;i++) for(int j=0;j<N;j++) newState.State[Vpot][i][j] = Soln[ Sidx(i,j,M) ];

    // Vpot has been updated, re-normalize so V = 0 in the upper-left corner
    // With a fixed V=0 filament boundary, this renormalization is worked into the scheme. No need to renormalize
    //double Vcorner = newState.State[Vpot][0][N-1];
    //for(int i=0;i<M;i++) for(int j=0;j<N;j++) newState.State[Vpot][i][j] = newState.State[Vpot][i][j] - Vcorner;

};


void createPoissonVFixedTopRightMatrixS0(Params simP, TheState curState, MatrixXd& PoMatrix, VectorXd& Source, double s0)
{
    int s=0;

    double DeltaR = simP.dR; double DeltaZ = simP.dZ;
    int M = simP.M; int N = simP.N;

    double Vedge = 0; // Potential value at the filament boundary

    // First, the source vector :  4*pi*dr^2*rho = 4*pi*dr^2*Q*exp(-V)
    // here we alter the given potential by a functional form s = 1 + s0*(z-zL)*(r-rz)
    for(s=0;s<M*N;s++)
    {
        if( cPos(Ridx(s,M),DeltaR) <= curState.VContour[Zidx(s,M)] )  // In the filament
        {
            if( cPos(Ridx(s,M)+1,DeltaR) > curState.VContour[Zidx(s,M)] ) // next to the filament
            {
                Source(s) = 0.0; // Matrix will be a linear interpolation
            }
            else
            {
                double curV = curState.State[Vpot][Ridx(s,M)][Zidx(s,M)];
                double sfun = 1.0 + s0*(cPos(Zidx(s,M),DeltaZ) - simP.zL)*(cPos(Ridx(s,M),DeltaR) - curState.VContour[Zidx(s,M)]);
                Source(s) = 4.0*PI*pow(DeltaR,2)*curState.State[Q][Ridx(s,M)][Zidx(s,M)]*exp(-sfun*curV);
            }
        }
        else   // Outside the filament, rho = 0 (we will also fix V out here)
            Source(s) = 0.0;
    }


    // Second, loop through and fill up the finite difference matrix
    for(s=0;s<M*N;s++)
    {
        double wi,f,Mleft,Mright,Mup,Mdown,Mcenter;

        if( cPos(Ridx(s,M),DeltaR) == 0 ) // left edge
        {
            // Useful numbers
            wi = 0;
            f  = DeltaR/DeltaZ;

            // r=0 version, derived using L'Hopital on the original DiffEq and employing dV/dr=0
            Mleft = 2.0;
            Mright = 2.0;
            Mup = pow(f,2);
            Mdown = pow(f,2);
            Mcenter = -2.0*(2.0+pow(f,2));
        }
        else
        {
            // Useful numbers
            wi = DeltaR/2.0/cPos(Ridx(s,M),DeltaR);
            f  = DeltaR/DeltaZ;

            // Original values of matrix entries (may change on boundaries)
            Mleft = 1.0 - wi;
            Mright = 1.0 + wi;
            Mup = pow(f,2);
            Mdown = pow(f,2);
            Mcenter = -2.0*(1.0+pow(f,2));
        }

        // Now we will change some of these values to employ the boundary conditions

        // The top boundary employs a V = Vcyl + alpha correction
        /*
        if( Zidx(s,M) == N-1 )
        {
            Source(s) = curState.State[Vpot][Ridx(s,M)][N-1];
            Mcenter = 1; Mleft = 0; Mright = 0; Mdown = 0; Mup = 0;
        }
        */

        // Top employs Neumann
        if( Zidx(s,M) == N-1 ) { Mdown = Mdown + Mup; Mup = 0;} // there is no up


        // The bottom boundary employs a dV/dz = 0 condition so at z=0, V(z=-1) = V(z=1)
        if( Zidx(s,M) == 0 ) { Mup = Mup + Mdown;  Mdown = 0;} // there is no down


        // The left boundary employs a dV/dr = 0 condition, so V(r=-1) = V(r=1)
        if( Ridx(s,M) == 0 ) { Mright = Mright + Mleft ; Mleft = 0;}  // there is no left

        // If we are well beyond the filament, we are not changing V
        if( cPos(Ridx(s,M),DeltaR) >= curState.VContour[Zidx(s,M)] )
        {
            Mleft = 0.0;
            Mright = 0.0;
            Mup = 0.0;
            Mdown = 0.0;
            Mcenter = 1.0;

            Source(s) = Vedge; //curState.State[Vpot][Ridx(s,M)][Zidx(s,M)];

        }
        else if( cPos(Ridx(s,M)+1,DeltaR) > curState.VContour[Zidx(s,M)] ) // If we are less than a cell from the boundary
        {
            // Instead we will interpolate between the boundary value and the adjacent cell


            double x = DeltaR / ( DeltaR + (curState.VContour[Zidx(s,M)] - cPos(Ridx(s,M),DeltaR)) );
            Mleft = 1-x;
            Mright = 0;
            Mcenter = -1;
            Mup = 0;
            Mdown = 0;

            Source(s) = -x*Vedge; // which is just zero since Vedge = 0


        }

        // Now we fill the entries of the matrix with the values
        PoMatrix(s,s) = Mcenter;
        if(Ridx(s,M)!=0) PoMatrix(s,s-1) = Mleft;
        if(Ridx(s,M)!=M-1) PoMatrix(s,s+1) = Mright;
        if(Zidx(s,M)!=N-1) PoMatrix(s,s+M) = Mup;
        if(Zidx(s,M)!=0) PoMatrix(s,s-M) = Mdown;

    }






};



void createPoissonVFixedTopBottomRightMatrix(Params simP, TheState curState, MatrixXd& PoMatrix, VectorXd& Source, VectorXd Soln_prev, int which)
{
    int s=0;

    double DeltaR = simP.dR; double DeltaZ = simP.dZ;
    int M = simP.M; int N = simP.N;

    double Vedge = 0; // Potential value at the filament boundary

    // clean out the vectors
    for(s=0;s<M*N;s++)
    {
        for(int ss=0;ss<M*N;ss++) PoMatrix(s,ss) = 0.0;

        Source(s) = 0.0;
    }

    // First, the source vector :  4*pi*dr^2*rho = 4*pi*dr^2*Q*exp(-V)
    for(s=0;s<M*N;s++)
    {
        if( cPos(Ridx(s,M),DeltaR) <= curState.VContour[Zidx(s,M)] )  // In the filament
        {
            if( cPos(Ridx(s,M)+1,DeltaR) > curState.VContour[Zidx(s,M)] ) // next to the filament
            {
                Source(s) = 0.0; // Matrix will be a linear interpolation
            }
            else
            {
                if(which==0) Source(s) = 4.0*PI*pow(DeltaR,2)*curState.State[Q][Ridx(s,M)][Zidx(s,M)]*exp(-curState.State[Vpot][Ridx(s,M)][Zidx(s,M)]);
                else if(which==1) Source(s) = 4.0*PI*pow(DeltaR,2)*curState.State[Q][Ridx(s,M)][Zidx(s,M)]*exp(-Soln_prev(s));
                else Source(s) = 4.0*PI*pow(DeltaR,2)*curState.State[Q][Ridx(s,M)][Zidx(s,M)]*exp( -(0.8*curState.State[Vpot][Ridx(s,M)][Zidx(s,M)] + 0.2*Soln_prev(s)) );
            }
        }
        else   // Outside the filament, rho = 0 (we will also fix V out here)
            Source(s) = 0.0;
    }


    // Second, loop through and fill up the finite difference matrix
    for(s=0;s<M*N;s++)
    {
        double wi,f,Mleft,Mright,Mup,Mdown,Mcenter;

        if( cPos(Ridx(s,M),DeltaR) == 0 ) // left edge
        {
            // Useful numbers
            wi = 0;
            f  = DeltaR/DeltaZ;

            // r=0 version, derived using L'Hopital on the original DiffEq and employing dV/dr=0
            Mleft = 2.0;
            Mright = 2.0;
            Mup = pow(f,2);
            Mdown = pow(f,2);
            Mcenter = -2.0*(2.0+pow(f,2));
        }
        else
        {
            // Useful numbers
            wi = DeltaR/2.0/cPos(Ridx(s,M),DeltaR);
            f  = DeltaR/DeltaZ;

            // Original values of matrix entries (may change on boundaries)
            Mleft = 1.0 - wi;
            Mright = 1.0 + wi;
            Mup = pow(f,2);
            Mdown = pow(f,2);
            Mcenter = -2.0*(1.0+pow(f,2));
        }

        // Now we will change some of these values to employ the boundary conditions

        // The top boundary employs V = cylinder condition at z=N-1
        // if( Zidx(s,M) == N-1 ) {Source(s) = simP.VTop[Ridx(s,M)]; Mcenter = 1; Mleft = 0; Mright = 0; Mdown = 0; Mup = 0;}
        if( Zidx(s,M) == N-1 ) {Source(s) = curState.State[Vpot][Ridx(s,M)][N-1]; Mcenter = 1; Mleft = 0; Mright = 0; Mdown = 0; Mup = 0;}


        // The bottom boundary employs V = cylinder + correction
        if( Zidx(s,M) == 0 )
        {
            Mcenter = 1; Mleft = 0; Mright = 0; Mdown = 0; Mup = 0;

            // If alteration is V(r) = aR^2 + bR + c,
            // requiring V(0) = Vknob -->  c = Vknob
            // requiring dVdr(0) = 0 --> b = 0
            // requiring V(rBdy) = 0 --> a rBdy^2 + Vknob = 0 --> a = -Vknob / rBdy^2
            // where rBdy is the radius of the filament bounary at z = 0
            double Vknob = -1.0*simP.Vknob ;
            double rBdybttm = curState.VContour[0] ;
            double aCst = -Vknob/rBdybttm/rBdybttm;

            // curState.State[Vpot][Ridx(s,M)][N-1] should always remain the cylinder solution
            Source(s) = curState.State[Vpot][Ridx(s,M)][N-1] + aCst*pow(cPos(Ridx(s,M),DeltaR),2) + Vknob;

        }

        // The left boundary employs a dV/dr = 0 condition, so V(r=-1) = V(r=1)
        if( Ridx(s,M) == 0 ) { Mright = Mright + Mleft ; Mleft = 0;}  // there is no left

        // If we are well beyond the filament, we are not changing V
        if( cPos(Ridx(s,M),DeltaR) >= curState.VContour[Zidx(s,M)] )
        {
            Mleft = 0.0;
            Mright = 0.0;
            Mup = 0.0;
            Mdown = 0.0;
            Mcenter = 1.0;

            /*
            if(which==0) Source(s) = curState.State[Vpot][Ridx(s,M)][Zidx(s,M)];
            else Source(s) = Soln_prev(s);
            */

            Source(s) = Vedge;

        }
        else if( cPos(Ridx(s,M)+1,DeltaR) > curState.VContour[Zidx(s,M)] and Zidx(s,M)<N-1) // If we are less than a cell from the boundary
        {
            // Instead we will interpolate between the boundary value and the adjacent cell


            double x  = DeltaR / ( DeltaR + (curState.VContour[Zidx(s,M)] - cPos(Ridx(s,M),DeltaR)) );
            Mleft = 1-x;
            Mright = 0;
            Mcenter = -1;
            Mup = 0;
            Mdown = 0;

            Source(s) = -x*Vedge; // which is just zero since Vedge = 0


        }

        // Now we fill the entries of the matrix with the values
        //if(which==0)
        if(1)
        {
            PoMatrix(s,s) = Mcenter;
            if(Ridx(s,M)!=0) PoMatrix(s,s-1) = Mleft;
            if(Ridx(s,M)!=M-1) PoMatrix(s,s+1) = Mright;
            if(Zidx(s,M)!=N-1) PoMatrix(s,s+M) = Mup;
            if(Zidx(s,M)!=0) PoMatrix(s,s-M) = Mdown;
        }


    }



};


void createPoissonVFixedMatrix(Params simP, TheState curState, MatrixXd& PoMatrix, VectorXd& Source)
{
    int s=0;

    double DeltaR = simP.dR; double DeltaZ = simP.dZ;
    int M = simP.M; int N = simP.N;

    double Vedge = 0; // Potential value at the filament boundary

    // First, the source vector :  4*pi*dr^2*rho = 4*pi*dr^2*Q*exp(-V)
    for(s=0;s<M*N;s++)
    {
        if( cPos(Ridx(s,M),DeltaR) <= curState.VContour[Zidx(s,M)] )  // In the filament
        {
            if( cPos(Ridx(s,M)+1,DeltaR) > curState.VContour[Zidx(s,M)] ) // next to the filament
            {
                Source(s) = 0.0; // Matrix will be a linear interpolation
            }
            else
            {
                Source(s) = 4.0*PI*pow(DeltaR,2)*curState.State[Q][Ridx(s,M)][Zidx(s,M)]*exp(-curState.State[Vpot][Ridx(s,M)][Zidx(s,M)]);
            }
        }
        else   // Outside the filament, rho = 0 (we will also fix V out here)
            Source(s) = 0.0;
    }


    // Second, loop through and fill up the finite difference matrix
    for(s=0;s<M*N;s++)
    {
        double wi,f,Mleft,Mright,Mup,Mdown,Mcenter;

        if( cPos(Ridx(s,M),DeltaR) == 0 ) // left edge
        {
            // Useful numbers
            wi = 0;
            f  = DeltaR/DeltaZ;

            // r=0 version, derived using L'Hopital on the original DiffEq and employing dV/dr=0
            Mleft = 2.0;
            Mright = 2.0;
            Mup = pow(f,2);
            Mdown = pow(f,2);
            Mcenter = -2.0*(2.0+pow(f,2));
        }
        else
        {
            // Useful numbers
            wi = DeltaR/2.0/cPos(Ridx(s,M),DeltaR);
            f  = DeltaR/DeltaZ;

            // Original values of matrix entries (may change on boundaries)
            Mleft = 1.0 - wi;
            Mright = 1.0 + wi;
            Mup = pow(f,2);
            Mdown = pow(f,2);
            Mcenter = -2.0*(1.0+pow(f,2));
        }

        // Now we will change some of these values to employ the boundary conditions

        // The top boundary employs a dV/dz = 0 condition so at z=N-1, V(z=N) = V(r=N-2)
        if( Zidx(s,M) == N-1 ) { Mdown = Mdown + Mup; Mup = 0;} // there is no up


        // The bottom boundary employs a dV/dz = 0 condition so at z=0, V(z=-1) = V(z=1)
        if( Zidx(s,M) == 0 ) { Mup = Mup + Mdown;  Mdown = 0;} // there is no down


        // The left boundary employs a dV/dr = 0 condition, so V(r=-1) = V(r=1)
        if( Ridx(s,M) == 0 ) { Mright = Mright + Mleft ; Mleft = 0;}  // there is no left

        // If we are well beyond the filament, we are not changing V
        if( cPos(Ridx(s,M),DeltaR) >= curState.VContour[Zidx(s,M)] )
        {
            Mleft = 0.0;
            Mright = 0.0;
            Mup = 0.0;
            Mdown = 0.0;
            Mcenter = 1.0;

            Source(s) = curState.State[Vpot][Ridx(s,M)][Zidx(s,M)];

        }
        else if( cPos(Ridx(s,M)+1,DeltaR) > curState.VContour[Zidx(s,M)] ) // If we are less than a cell from the boundary
        {
            // Instead we will interpolate between the boundary value and the adjacent cell


            double x = DeltaR / ( DeltaR + (curState.VContour[Zidx(s,M)] - cPos(Ridx(s,M),DeltaR)) );
            Mleft = 1-x;
            Mright = 0;
            Mcenter = -1;
            Mup = 0;
            Mdown = 0;

            Source(s) = -x*Vedge; // which is just zero since Vedge = 0


        }

        // The right boundary employs a dV/dr = f(r,z) condition, so at r=M-1: V(M) = V(M-2) + 2*dr*f(r,z)
        // Assuming this is always outside the filament boundary
        #if(0) // the 'and 0' means we never actually use this
        if( Ridx(s,M) == M-1)
        {
            // Uses Poisson equation and a ghost zone
            //Mleft = Mleft + Mright;
            //Mright = 0;
            //Source(s) = Source(s) - (1.0+wi)*2.0*DeltaR*VRight[Zidx(s)];

            // Employs a 2nd order approximation to the first derivative at M-1
            Mleft = -4.0;
            Mright = 0; // doesn't exist
            Mup = 0;
            Mdown = 0;
            Mcenter = 1;
            PoMatrix(s,s-2) = 1; // Need two neighbor points to the left
            Source(s) = 6.0*DeltaR*VRight[Zidx(s,M)];

        }
        #endif

        // Now we fill the entries of the matrix with the values
        PoMatrix(s,s) = Mcenter;
        if(Ridx(s,M)!=0) PoMatrix(s,s-1) = Mleft;
        if(Ridx(s,M)!=M-1) PoMatrix(s,s+1) = Mright;
        if(Zidx(s,M)!=N-1) PoMatrix(s,s+M) = Mup;
        if(Zidx(s,M)!=0) PoMatrix(s,s-M) = Mdown;

    }





};

void createPoissonCstMatrix(Params simP, TheState curState, MatrixXd& PoMatrix, VectorXd& Source)
{
    int s=0;

    double DeltaR = simP.dR; double DeltaZ = simP.dZ;
    int M = simP.M; int N = simP.N;

    // What is the potential along the contour boundary. We will use the top row.
    // s here is the radial index, after this it'll be an index that ranges from 0 the M*N
    while(true){s++; if( cPos(s,DeltaR) >= curState.VContour[N-1] ) break; }
    double x = ( cPos(s,DeltaR)-curState.VContour[N-1] )/DeltaR;
    double Vedge = x*curState.State[Vpot][s-1][N-1] + (1-x)*curState.State[Vpot][s][N-1];

    // 0th, employ the array that captures the boundary conditions on the right side of the box
    double VRight[N];

    // V cylinder outside has the form  V = C1*ln(r/rEdge) + C2
    // C1 and C2 chosen so V is smooth at rEdge (V(rE) and DV(rE) match)
    // C1 = rE*DV(rE)     C2 = V(rE)
    // so dV/dR = C1 * rEdge/r
    for(s=0;s<N;s++)
    {
        double rEdge = curState.VContour[N-1];
        VRight[s] = simP.CylCsts[0] * rEdge/cPos(M-1,DeltaR);
    }

    for(s=0;s<simP.pointLoopTotNum;s++)
    {
        double zP = ((double) s)*2.0*simP.zL;
        double Mcode = simP.mExcess;

        // add GM/r^2 contribution to each cell
        for(int j=0;j<N;j++)
        {
            // Distance locations of point particles
            double X1 = sqrt( pow(cPos(M-1,DeltaR),2) + pow(cPos(j,DeltaZ)-zP,2) );
            double X2 = sqrt( pow(cPos(M-1,DeltaR),2) + pow(cPos(j,DeltaZ)+zP,2) );

            VRight[j] = VRight[j] + Mcode/pow(X1,2);
            if(zP > 0) VRight[j] = VRight[j] + Mcode/pow(X2,2);

        }

    }

    // First, the source vector :  4*pi*dr^2*rho = 4*pi*dr^2*Q*exp(-V)
    for(s=0;s<M*N;s++)
    {
        if( cPos(Ridx(s,M),DeltaR) <= curState.VContour[Zidx(s,M)] )  // In the filament
        {
            Source(s) = 4.0*PI*pow(DeltaR,2)*curState.State[Q][Ridx(s,M)][Zidx(s,M)]*exp(-curState.State[Vpot][Ridx(s,M)][Zidx(s,M)]);
        }
        else   // Outside the filament, rho = 0
            Source(s) = 0.0;
    }

    // Second, loop through and fill up the finite difference matrix
    for(s=0;s<M*N;s++)
    {
        double wi,f,Mleft,Mright,Mup,Mdown,Mcenter;

        if( cPos(Ridx(s,M),DeltaR) == 0 )
        {
            // Useful numbers
            wi = 0;
            f  = DeltaR/DeltaZ;

            // r=0 version, derived using L'Hopital on the original DiffEq and employing dV/dr=0
            Mleft = 2.0;
            Mright = 2.0;
            Mup = pow(f,2);
            Mdown = pow(f,2);
            Mcenter = -2.0*(2.0+pow(f,2));
        }
        else
        {
            // Useful numbers
            wi = DeltaR/2.0/cPos(Ridx(s,M),DeltaR);
            f  = DeltaR/DeltaZ;

            // Original values of matrix entries (may change on boundaries)
            Mleft = 1.0 - wi;
            Mright = 1.0 + wi;
            Mup = pow(f,2);
            Mdown = pow(f,2);
            Mcenter = -2.0*(1.0+pow(f,2));
        }

        // Now we will change some of these values to employ the boundary conditions

        // The top boundary employs a dV/dz = 0 condition so at z=N-1, V(z=N) = V(r=N-2)
        if( Zidx(s,M) == N-1 ) { Mdown = Mdown + Mup; Mup = 0;} // there is no up


        // The bottom boundary employs a dV/dz = 0 condition so at z=0, V(z=-1) = V(z=1)
        if( Zidx(s,M) == 0 ) { Mup = Mup + Mdown;  Mdown = 0;} // there is no down


        // The left boundary employs a dV/dr = 0 condition, so V(r=-1) = V(r=1)
        if( Ridx(s,M) == 0 ) { Mright = Mright + Mleft ; Mleft = 0;}  // there is no left

        // If we are well beyond the filament, we are not changing V
        if( cPos(Ridx(s,M),DeltaR) >= curState.VContour[Zidx(s,M)] )
        {
            Mleft = 0.0;
            Mright = 0.0;
            Mup = 0.0;
            Mdown = 0.0;
            Mcenter = 1.0;

            Source(s) = curState.State[Vpot][Ridx(s,M)][Zidx(s,M)];

        }
        else if( cPos(Ridx(s+1,M),DeltaR) > curState.VContour[Zidx(s,M)] ) // If we are less than a cell from the boundary
        {
            // Instead we will interpolate between the boundary value and the adjacent cell


            double x = DeltaR / ( DeltaR + (curState.VContour[Zidx(s,M)] - cPos(Ridx(s,M),DeltaR)));
            Mleft = 1-x;
            Mright = 0;
            Mcenter = -1;
            Mup = 0;
            Mdown = 0;

            Source(s) = -x*Vedge;


        }

        // The right boundary employs a dV/dr = f(r,z) condition, so at r=M-1: V(M) = V(M-2) + 2*dr*f(r,z)
        // Assuming this is always outside the filament boundary
        if( Ridx(s,M) == M-1 and 0)
        {
            // Uses Poisson equation and a ghost zone
            //Mleft = Mleft + Mright;
            //Mright = 0;
            //Source(s) = Source(s) - (1.0+wi)*2.0*DeltaR*VRight[Zidx(s)];

            // Employs a 2nd order approximation to the first derivative at M-1
            Mleft = -4.0;
            Mright = 0; // doesn't exist
            Mup = 0;
            Mdown = 0;
            Mcenter = 1;
            PoMatrix(s,s-2) = 1; // Need two neighbor points to the left
            Source(s) = 6.0*DeltaR*VRight[Zidx(s,M)];

        }

        // Now we fill the entries of the matrix with the values
        PoMatrix(s,s) = Mcenter;
        if(Ridx(s,M)!=0) PoMatrix(s,s-1) = Mleft;
        if(Ridx(s,M)!=M-1) PoMatrix(s,s+1) = Mright;
        if(Zidx(s,M)!=N-1) PoMatrix(s,s+M) = Mup;
        if(Zidx(s,M)!=0) PoMatrix(s,s-M) = Mdown;

    }




};


void createPoissonMatrixFixedDen(Params simP, TheState curState, MatrixXd& PoMatrix, VectorXd& Source)
{
    int s=0;

    double DeltaR = simP.dR; double DeltaZ = simP.dZ;
    int M = simP.M; int N = simP.N;

    // What is the potential along the contour boundary. We will use the top row.
    // s here is the radial index, after this it'll be an index that ranges from 0 the M*N
    while(true){s++; if( cPos(s,DeltaR) >= curState.VContour[N-1] ) break; }
    double x = ( cPos(s,DeltaR)-curState.VContour[N-1] )/DeltaR;
    double Vedge = x*curState.State[Vpot][s-1][N-1] + (1-x)*curState.State[Vpot][s][N-1];
    double Vupleft = curState.State[Vpot][M-1][N-1];

    // 0th, employ the array that captures the boundary conditions on the right side of the box
    double VRight[N];

    // V cylinder outside has the form  V = C1*ln(r/rEdge) + C2
    // C1 and C2 chosen so V is smooth at rEdge (V(rE) and DV(rE) match)
    // C1 = rE*DV(rE)     C2 = V(rE)
    // so dV/dR = C1 * rEdge/r
    for(s=0;s<N;s++)
    {
        double rEdge = curState.VContour[N-1];
        VRight[s] = simP.CylCsts[0] * rEdge/cPos(M-1,DeltaR);
    }

    for(s=0;s<simP.pointLoopTotNum;s++)
    {
        double zP = ((double) s)*2.0*simP.zL;
        double Mcode = simP.mExcess;

        // add GM/r^2 contribution to each cell
        for(int j=0;j<N;j++)
        {
            // Distance locations of point particles
            double X1 = sqrt( pow(cPos(M-1,DeltaR),2) + pow(cPos(j,DeltaZ)-zP,2) );
            double X2 = sqrt( pow(cPos(M-1,DeltaR),2) + pow(cPos(j,DeltaZ)+zP,2) );

            VRight[j] = VRight[j] + Mcode/pow(X1,2);
            if(zP > 0) VRight[j] = VRight[j] + Mcode/pow(X2,2);

        }
    }

    for(s=0;s<N;s++) VRight[s] = VRight[s] + Vupleft; // the upper left isn't where V = 0 anymore..

    // First, the source vector :  4*pi*dr^2*rho = 4*pi*dr^2*Q*exp(-V)
    for(s=0;s<M*N;s++)
    {
        if( cPos(Ridx(s,M),DeltaR) <= curState.VContour[Zidx(s,M)] )  // In the filament
        {
            Source(s) = 4.0*PI*pow(DeltaR,2)*curState.State[Q][Ridx(s,M)][Zidx(s,M)]*exp(-curState.State[Vpot][Ridx(s,M)][Zidx(s,M)]);
        }
        else   // Outside the filament, rho = 0
            Source(s) = 0.0;
    }

    // Second, loop through and fill up the finite difference matrix
    for(s=0;s<M*N;s++)
    {
        double wi,f,Mleft,Mright,Mup,Mdown,Mcenter;

        if( cPos(Ridx(s,M),DeltaR) == 0 )
        {
            // Useful numbers
            wi = 0;
            f  = DeltaR/DeltaZ;

            // r=0 version, derived using L'Hopital on the original DiffEq and employing dV/dr = 0
            Mleft = 2.0;
            Mright = 2.0;
            Mup = pow(f,2);
            Mdown = pow(f,2);
            Mcenter = -2.0*(2.0+pow(f,2));
        }
        else
        {
            // Useful numbers
            wi = DeltaR/2.0/cPos(Ridx(s,M),DeltaR);
            f  = DeltaR/DeltaZ;

            // Original values of matrix entries (may change on boundaries)
            Mleft = 1.0 - wi;
            Mright = 1.0 + wi;
            Mup = pow(f,2);
            Mdown = pow(f,2);
            Mcenter = -2.0*(1.0+pow(f,2));
        }

        // Now we will change some of these values to employ the boundary conditions

        // The top boundary employs a dV/dz = 0 condition so at z=N-1, V(z=N) = V(r=N-2)
        if( Zidx(s,M) == N-1 ) { Mdown = Mdown + Mup; Mup = 0;} // there is no up


        // The bottom boundary employs a dV/dz = 0 condition so at z=0, V(z=-1) = V(z=1)
        if( Zidx(s,M) == 0 ) { Mup = Mup + Mdown;  Mdown = 0;} // there is no down


        // The left boundary employs a dV/dr = 0 condition, so V(r=-1) = V(r=1)
        if( Ridx(s,M) == 0 ) { Mright = Mright + Mleft ; Mleft = 0;}  // there is no left

        // If we are well beyond the filament, we are not changing V
        if( cPos(Ridx(s,M),DeltaR) >= curState.VContour[Zidx(s,M)] )
        {
            Mleft = 0.0;
            Mright = 0.0;
            Mup = 0.0;
            Mdown = 0.0;
            Mcenter = 1.0;

            Source(s) = curState.State[Vpot][Ridx(s,M)][Zidx(s,M)];

        }
        else if( cPos(Ridx(s+1,M),DeltaR) > curState.VContour[Zidx(s,M)] ) // If we are less than a cell from the boundary
        {
            // Instead we will interpolate between the boundary value and the adjacent cell


            double x = DeltaR / ( DeltaR + (curState.VContour[Zidx(s,M)] - cPos(Ridx(s,M),DeltaR)));
            Mleft = 1-x;
            Mright = 0;
            Mcenter = -1;
            Mup = 0;
            Mdown = 0;

            Source(s) = -x*Vedge;


        }

        // Reset the lower-left
        if( Ridx(s,M)==0 and Zidx(s,M)==0 )
        {
            Mleft = 0.0;
            Mright = 0;
            Mup = 0;
            Mdown = 0;
            Mcenter = 1;

            Source(s) = curState.State[Vpot][0][N-1] - log(simP.denratio) ;
            cout << "Lower-left source value : " << Source(s) << endl;

        }

        // The right boundary employs a dV/dr = f(r,z) condition, so at r=M-1: V(M) = V(M-2) + 2*dr*f(r,z)
        // Assuming this is always outside the filament boundary
        if( Ridx(s,M) == M-1 and 0)
        {
            // Uses Poisson equation and a ghost zone
            //Mleft = Mleft + Mright;
            //Mright = 0;
            //Source(s) = Source(s) - (1.0+wi)*2.0*DeltaR*VRight[Zidx(s)];

            // Employs a 2nd order approximation to the first derivative at M-1
            Mleft = -4.0;
            Mright = 0; // doesn't exist
            Mup = 0;
            Mdown = 0;
            Mcenter = 1;
            PoMatrix(s,s-2) = 1; // Need two neighbor points to the left
            Source(s) = 6.0*DeltaR*VRight[Zidx(s,M)];

        }

        // Now we fill the entries of the matrix with the values
        PoMatrix(s,s) = Mcenter;
        if(Ridx(s,M)!=0) PoMatrix(s,s-1) = Mleft;
        if(Ridx(s,M)!=M-1) PoMatrix(s,s+1) = Mright;
        if(Zidx(s,M)!=N-1) PoMatrix(s,s+M) = Mup;
        if(Zidx(s,M)!=0) PoMatrix(s,s-M) = Mdown;

    }




};
