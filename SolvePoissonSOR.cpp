

#include "SolvePoissonSOR.H"

// For some reason this has to be included here and not in the header file ... not sure why
#include <Eigen/Dense>
using namespace Eigen;

void createPoissonVFixedMatrixSOR(Params simP, TheState curState, MatrixXd& PoMatrix, VectorXd& Source, VectorXd Soln_prev, int which);
void createPoissonVFixedTopRightMatrixSOR(Params simP, TheState curState, MatrixXd& PoMatrix, VectorXd& Source, VectorXd Soln_prev, int which);
void createPoissonVFixedTopBottomRightMatrixSOR(Params simP, TheState curState, MatrixXd& PoMatrix, VectorXd& Source, VectorXd Soln_prev, int which);

void SolvePoissonSOR(Params simP,TheState curState,TheState& newState)
{
    int M = curState.M;
    int N = curState.N;

    // First allocate the matrix and vectors for the finite difference scheme, source, and solution
    MatrixXd PoMatrix = MatrixXd::Zero((M)*(N),(M)*(N));
    VectorXd Source = VectorXd::Zero((M)*(N));
    VectorXd Soln_new  = VectorXd::Zero((M)*(N));
    VectorXd Soln_prev = VectorXd::Zero((M)*(N));


    cout << "Creating Poisson matrix for SOR." << endl;
    //createPoissonVFixedMatrixSOR(simP,curState,PoMatrix,Source,Soln_prev,0); // Doesn't really solve Poisson outside of the filament
    //createPoissonVFixedTopRightMatrixSOR(simP,curState,PoMatrix,Source,Soln_prev,0); // V=fixed top and at filament
    createPoissonVFixedTopBottomRightMatrixSOR(simP,curState,PoMatrix,Source,Soln_prev,0);

    for(int i=0;i<M*N;i++) Soln_prev(i) = curState.State[Vpot][Ridx(i,M)][Zidx(i,M)];

    double SORomega = simP.SORomega;
    // This is the relaxation loop

    cout << "Applying SOR for Poisson" << endl;
    for(int kk=0;kk<simP.SORnum;kk++)
    {
        for(int k=0;k<5000;k++)
        {
            for(int i=0;i<M*N;i++)
            {
                double val = 0;

                double a_diag = PoMatrix(i,i);
                double b_val  = Source(i);

                if(a_diag==0) cout << "ISSUE: Diagonal entry of Poisson SOR solve is zero!" << endl;

                for(int j=0;j<i;j++) val = val + PoMatrix(i,j)*Soln_new(j);
                for(int j=i+1;j<M*N;j++) val = val + PoMatrix(i,j)*Soln_prev(j);
                //for(int j=0;j<M*N;j++) {if(j==i) continue; else val = val + PoMatrix(i,j)*Soln_prev(j);}
                val = (b_val - val);
                val = SORomega*(1.0/a_diag)*val; // this alone with SORomega = 1 is Gauss-Seidel
                val = val + (1.0-SORomega)*Soln_prev(i); // interpolates with old solution (SOR)

                Soln_new(i) = val;
            }

            for(int i=0;i<M*N;i++) Soln_prev(i) = Soln_new(i);

            // Test to make sure solution is actually a solution (Uses L2 norm)
            double L2error = (PoMatrix*Soln_new - Source).norm() / Source.norm();
            if(k%100==0) cout << "SOR: Gravity Matrix L2 error = " << L2error << endl;

            if(L2error <= simP.SORtol) {cout << "SOR: Breaking with " << k << " iterations." << endl; break; }

            // else we need to update source vector and try again
            // 0 = uses curState to load source vector, at this point it hasn't been changed
            // 1 = uses Soln_prev to load source vector, which has recently been overwritten as Soln_new
            // 2 = uses a blend
            //createPoissonVFixedTopRightMatrixSOR(simP,curState,PoMatrix,Source,Soln_prev,0); // V=fixed top and at filament
            createPoissonVFixedTopBottomRightMatrixSOR(simP,curState,PoMatrix,Source,Soln_prev,0);

        }


        for(int i=0;i<M*N;i++) curState.State[Vpot][Ridx(i,M)][Zidx(i,M)] = Soln_new(i);

        // Comment out to redraw Q and dQdPhi
        //updateQ(simP, curState, 0);
        //updateDQDPHI(simP, curState);

        // Updates Source vector
        //createPoissonVFixedMatrixSOR(simP,curState,PoMatrix,Source,Soln_prev,0); // use NEW curState and Q
        // 0  = uses curState to load source, but was just written as the converged value.
        // 2 = us
        //createPoissonVFixedTopRightMatrixSOR(simP,curState,PoMatrix,Source,Soln_prev,0);
        createPoissonVFixedTopBottomRightMatrixSOR(simP,curState,PoMatrix,Source,Soln_prev,0);

        double L2error = (PoMatrix*Soln_new - Source).norm() / Source.norm();
        cout << "SOR: Gravity Matrix L2 outer loop error (" << kk << ") = " << L2error << endl;
        if(L2error <= 1.0e-3) {cout << "SOR: Breaking outer loop with " << kk << " iterations." << endl; break; }
    }


    for(int i=0;i<M*N;i++) newState.State[Vpot][Ridx(i,M)][Zidx(i,M)] = Soln_new(i);

    return;
};



void createPoissonVFixedTopBottomRightMatrixSOR(Params simP, TheState curState, MatrixXd& PoMatrix, VectorXd& Source, VectorXd Soln_prev, int which)
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



void createPoissonVFixedTopRightMatrixSOR(Params simP, TheState curState, MatrixXd& PoMatrix, VectorXd& Source, VectorXd Soln_prev, int which)
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



void createPoissonVFixedMatrixSOR(Params simP, TheState curState, MatrixXd& PoMatrix, VectorXd& Source, VectorXd Soln_prev, int which)
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
                if(which==0) Source(s) = 4.0*PI*pow(DeltaR,2)*curState.State[Q][Ridx(s,M)][Zidx(s,M)]*exp(-curState.State[Vpot][Ridx(s,M)][Zidx(s,M)]);
                else  Source(s) = 4.0*PI*pow(DeltaR,2)*curState.State[Q][Ridx(s,M)][Zidx(s,M)]*exp(-Soln_prev(s));
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

            if(which==0) Source(s) = curState.State[Vpot][Ridx(s,M)][Zidx(s,M)];
            else Source(s) = Soln_prev(s);

            Source(s) = Vedge;

        }
        else if( cPos(Ridx(s,M)+1,DeltaR) > curState.VContour[Zidx(s,M)] ) // If we are less than a cell from the boundary
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
        if(which==0)
        {
            PoMatrix(s,s) = Mcenter;
            if(Ridx(s,M)!=0) PoMatrix(s,s-1) = Mleft;
            if(Ridx(s,M)!=M-1) PoMatrix(s,s+1) = Mright;
            if(Zidx(s,M)!=N-1) PoMatrix(s,s+M) = Mup;
            if(Zidx(s,M)!=0) PoMatrix(s,s-M) = Mdown;
        }


    }





};
