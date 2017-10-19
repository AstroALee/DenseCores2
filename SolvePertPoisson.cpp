
#include "SolvePertPoisson.H"

// For some reason this has to be included here and not in the header file ... not sure why
#include <Eigen/Dense>
using namespace Eigen;

void createPertPoissonMatrix(Params simP, TheState curState, TheState prevState, MatrixXd& PpmMatrix, VectorXd& Source);

void SolvePertPoisson(Params simP,TheState curState, TheState prevState, TheState& newState)
{
    int M = curState.M;
    int N = curState.N;

    // First allocate the matrix and vectors for the finite difference scheme, source, and solution
    MatrixXd PoMatrix = MatrixXd::Zero((M)*(N),(M)*(N));
    VectorXd Source = VectorXd::Zero((M)*(N));
    VectorXd Soln = VectorXd::Zero((M)*(N));

    cout << "Creating Perturbed Poisson matrix." << endl;
    createPertPoissonMatrix(simP,curState,prevState,PoMatrix,Source);

    // Now call LinAlgebra package to invert the matrix and obtain the new potential values
    //Soln = PoMatrix.householderQr().solve(Source);
    //Soln = PoMatrix.colPivHouseholderQr().solve(Source);
    cout << "Inverting Perturbed Poisson matrix." << endl;
    Soln = PoMatrix.colPivHouseholderQr().solve(Source);

    // Test to make sure solution is actually a solution (Uses L2 norm)
    double L2error = (PoMatrix*Soln - Source).norm() / Source.norm();
    cout << "Pert Gravity Matrix L2 error = " << L2error << endl;

    // Copy the solution to the newState vector
    for(int i=0;i<M;i++) for(int j=0;j<N;j++) newState.State[Vpot][i][j] = Soln[ Sidx(i,j,M) ];

    // Vpot has been updated, re-normalize so V = 0 in the upper-left corner
    //double Vcorner = newState.State[Vpot][0][N-1];
    //for(int i=0;i<M;i++) for(int j=0;j<N;j++) newState.State[Vpot][i][j] = newState.State[Vpot][i][j] - Vcorner;

};


void createPertPoissonMatrix(Params simP, TheState curState, TheState prevState, MatrixXd& PpmMatrix, VectorXd& Source)
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
                    // get delta Q
                    double dQ = curState.State[Q][Ridx(s,M)][Zidx(s,M)] - prevState.State[Q][Ridx(s,M)][Zidx(s,M)];

                    // source term in perturbation
                    Source(s) = 4.0*PI*pow(DeltaR,2)*dQ*exp(-curState.State[Vpot][Ridx(s,M)][Zidx(s,M)]);
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
                Mcenter = -2.0*(2.0+pow(f,2)) + 4*PI*pow(DeltaR,2)*curState.State[Q][Ridx(s,M)][Zidx(s,M)]*exp(-curState.State[Vpot][Ridx(s,M)][Zidx(s,M)]);
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
                Mcenter = -2.0*(1.0+pow(f,2)) + 4.0*PI*pow(DeltaR,2)*curState.State[Q][Ridx(s,M)][Zidx(s,M)]*exp(-curState.State[Vpot][Ridx(s,M)][Zidx(s,M)]);
                    // includes an additional term from the perturbation expansion 4*pi*G*rho_0*dr^2*DeltaV
            }

            // Now we will change some of these values to employ the boundary conditions

            // The top boundary employs a dV/dz = 0 condition so at z=N-1, V(z=N) = V(r=N-2)
            if( Zidx(s,M) == N-1 ) { Mdown = Mdown + Mup; Mup = 0;} // there is no up


            // The bottom boundary employs a dV/dz = 0 condition so at z=0, V(z=-1) = V(z=1)
            if( Zidx(s,M) == 0 ) { Mup = Mup + Mdown;  Mdown = 0;} // there is no down

            // The left boundary employs a dV/dr = 0 condition, so V(r=-1) = V(r=1)
            if( Ridx(s,M) == 0 ) { Mright = Mright + Mleft ; Mleft = 0;}  // there is no left

            // If we are well beyond the filament, we are not changing V (this admittedly is kinda redundant in setting the Source(s) value)
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

            // Now we fill the entries of the matrix with the values
            PpmMatrix(s,s) = Mcenter;
            if(Ridx(s,M)!=0) PpmMatrix(s,s-1) = Mleft;
            if(Ridx(s,M)!=M-1) PpmMatrix(s,s+1) = Mright;
            if(Zidx(s,M)!=N-1) PpmMatrix(s,s+M) = Mup;
            if(Zidx(s,M)!=0) PpmMatrix(s,s-M) = Mdown;

        }

};
