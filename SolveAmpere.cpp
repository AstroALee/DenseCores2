#include "SolveAmpere.H"

// For some reason this has to be included here and not in the header file ... not sure why
#include <Eigen/Dense>
using namespace Eigen;


void createAmpereMatrix(Params simP, TheState curState, MatrixXd& AmMatrix, VectorXd& Source);
void createAmpereMatrixFixedTop(Params simP, TheState curState, MatrixXd& PoMatrix, VectorXd& Source, VectorXd Soln_prev, int which);
void createAmpereMatrixFixedTopS0(Params simP, TheState curState, MatrixXd& PoMatrix, VectorXd& Source, double s0);

void SolveAmpereS0(Params simP,TheState curState,TheState& newState, double s0)
{
        int M = curState.M;
        int N = curState.N;

        // First allocate the matrix and vectors for the finite difference scheme, source, and solution
        MatrixXd AmMatrix = MatrixXd::Zero((M)*(N),(M)*(N));
        VectorXd Source = VectorXd::Zero((M)*(N));
        VectorXd Soln = VectorXd::Zero((M)*(N));

        cout << "Creating Ampere matrix." << endl;
        createAmpereMatrixFixedTopS0(simP,curState,AmMatrix,Source,s0);

        // Now call LinAlgebra package to invert the matrix and obtain the new potential values
        //Soln = PoMatrix.householderQr().solve(Source);
        //Soln = PoMatrix.colPivHouseholderQr().solve(Source);
        cout << "Inverting Ampere matrix." << endl;
        Soln = AmMatrix.colPivHouseholderQr().solve(Source);

        // Test to make sure solution is actually a solution (Uses L2 norm)
        double L2error = (AmMatrix*Soln - Source).norm() / Source.norm();
        cout << "Ampere Matrix L2 error = " << L2error << endl;

        // Copy the solution to the newState vector
        for(int i=0;i<M;i++) for(int j=0;j<N;j++) newState.State[Apot][i][j] = Soln[ Sidx(i,j,M) ];

}


void SolveAmpere(Params simP,TheState curState,TheState& newState)
{

        int M = curState.M;
        int N = curState.N;

        // First allocate the matrix and vectors for the finite difference scheme, source, and solution
        MatrixXd AmMatrix = MatrixXd::Zero((M)*(N),(M)*(N));
        VectorXd Source = VectorXd::Zero((M)*(N));
        VectorXd Soln = VectorXd::Zero((M)*(N));

        cout << "Creating Ampere matrix." << endl;
        //createAmpereMatrix(simP,curState,AmMatrix,Source);
        createAmpereMatrixFixedTop(simP,curState,AmMatrix,Source,Soln,0);

        // Now call LinAlgebra package to invert the matrix and obtain the new potential values
        //Soln = PoMatrix.householderQr().solve(Source);
        //Soln = PoMatrix.colPivHouseholderQr().solve(Source);
        cout << "Inverting Ampere matrix." << endl;
        Soln = AmMatrix.colPivHouseholderQr().solve(Source);

        // Test to make sure solution is actually a solution (Uses L2 norm)
        double L2error = (AmMatrix*Soln - Source).norm() / Source.norm();
        cout << "Ampere Matrix L2 error = " << L2error << endl;

        // Copy the solution to the newState vector
        for(int i=0;i<M;i++) for(int j=0;j<N;j++) newState.State[Apot][i][j] = Soln[ Sidx(i,j,M) ];

}

void createAmpereMatrixFixedTopS0(Params simP, TheState curState, MatrixXd& AmMatrix, VectorXd& Source, double s0)
{
    int s;

    int M = simP.M;
    int N = simP.N;
    double DeltaR = simP.dR;
    double DeltaZ = simP.dZ;

    double Binf = sqrt(8.0*PI/simP.betaCyl)*pow(simP.Rbdy,simP.nCyl);
    double alpha = sqrt(8.0*PI/simP.betaCyl); //Binf; // alpha^2 = 8*pi/beta_inf = 8*pi*(B^2/8*pi)/(rho_*c^2)_inft -> B^2
    //double Ccst = PI*simP.RhoTop[0]/2.0/(1.0+1.0/simP.betaCyl); // RhoTop[0] should be 1


    // Source vector
    int BdyRight = 2;
    for(s=0;s<M*N;s++)
    {
        if( Ridx(s,M) == 0 )
        {
            // We will manually set A = 0 on the left hand side
            Source(s) = 0.0;
        }
        else if( Ridx(s,M) == M-1 )
        {
            if(BdyRight==0) // setting A
            {
                // We will manually set A = B_inf * r / 2
                Source(s) = Binf*cPos(M-1,DeltaR)/2.0;

                // Manually set to cylinder + empty space (the cylinder solution assumes nCyl = 0.5)
                double Redge = curState.VContour[N-1];
                //Source(s) = (0.5*alpha*sqrt(simP.RhoTop[0])/Ccst)*log(Ccst*Redge*Redge+1.0)/simP.rL + 0.5*Binf*(simP.rL-Redge*Redge/simP.rL);
            }
            else if(BdyRight==1) // setting dA/dr
            {

            }
            else // setting cros(A)
            {
                // We will impose a boundary condition (1/r)*d(rA)/dr = B_infty at this boundary
                Source(s) = 2.0*DeltaR*Binf;
            }
        }
        else
        {
            double dQdPhi = curState.State[dQdP][Ridx(s,M)][Zidx(s,M)];

            // Alters density
            double curV = curState.State[Vpot][Ridx(s,M)][Zidx(s,M)];
            double sFun = 1.0 + s0*(cPos(Zidx(s,M),DeltaZ)-curState.VContour[Zidx(s,M)])*(cPos(Zidx(s,M),DeltaZ)-simP.zL);

            if( cPos(Ridx(s,M),DeltaR) <= curState.VContour[Zidx(s,M)]) Source(s) = -4.0*PI*pow(DeltaR,2)*cPos(Ridx(s,M),DeltaR)*dQdPhi*exp(-sFun*curV);
            else Source(s) = 0;

        }
    }


    // Finite Difference Matrix
    for(s=0;s<M*N;s++)
    {
        double wi,f,Mleft,Mright,Mup,Mdown,Mcenter;

        if( Ridx(s,M) == 0 )
        {
            // Left side will normalize A = 0, so all the entries here are 0
            Mleft = 0;
            Mright = 0;
            Mup = 0;
            Mdown = 0;
            Mcenter = 1;
        }
        else if(Ridx(s,M)==M-1)
        {
            if(BdyRight==0)
            {
                // Imposing value of A
                Mleft = 0;
                Mright = 0;
                Mup = 0;
                Mdown = 0;
                Mcenter = 1.0;
            }
            else if(BdyRight==1)
            {

            }
            else if(BdyRight==2)
            {
                // Right side implemented a Robin boundary condition
                Mleft = -4.0;
                Mright = 0; // doesn't exist
                Mup = 0;
                Mdown = 0;
                Mcenter = ( 2.0*DeltaR/cPos(M-1,DeltaR) + 3.0 );
                AmMatrix(s,s-2) = 1.0;
            }
        }
        else
        {
            // Base values
            wi = DeltaR/2.0/cPos(Ridx(s,M),DeltaR);
            f = DeltaR/DeltaZ;

            Mright = 1.0+wi;
            Mleft = 1.0-wi;
            Mup = pow(f,2);
            Mdown = pow(f,2);
            Mcenter = -2.0-2.0*pow(f,2)-4.0*pow(wi,2);
        }

        // Now we will adjust to incorproate z-boundary conditions

        // at bottom, dA/dz = 0
        if( Zidx(s,M) == 0 )
        {
            Mup = Mup + Mdown;
            Mdown = 0;
        }

        // at top A = cylinder?
        if( Zidx(s,M) == N-1)
        {
            if(Ridx(s,M)<M-1) // not sure what these nested if statements are all about
            {

                // Dirchlet inside filament, Neumann outside

                if( cPos(Ridx(s,M),DeltaR) <= curState.VContour[Zidx(s,M)] ) 
                {
                  // enforce cylinder
                  Mdown = 0; Mup = 0; Mleft = 0; Mright = 0; Mcenter = 1;
                  Source(s) = simP.ATop[Ridx(s,M)];
                }
                else
                {
                  // classic Neumann
                  Mdown = Mup + Mdown;
                  Mup = 0;
                }

                //if(Ridx(s,M)==M-1) cout << "Value of ATop in upper-right = " << simP.ATop[M-1] << endl;
            }
            else
            {
                Mdown = Mup + Mdown;
                Mup = 0;
            }
        }


        // Now we fill the entries of the matrix with the values
        AmMatrix(s,s) = Mcenter;
        if(Ridx(s,M)!=0) AmMatrix(s,s-1) = Mleft;
        if(Ridx(s,M)!=M-1) AmMatrix(s,s+1) = Mright;
        if(Zidx(s,M)!=N-1) AmMatrix(s,s+M) = Mup;
        if(Zidx(s,M)!=0) AmMatrix(s,s-M) = Mdown;
    }


};


void createAmpereMatrixFixedTop(Params simP, TheState curState, MatrixXd& AmMatrix, VectorXd& Source, VectorXd Soln_prev, int which)
{
    int s;

    int M = simP.M;
    int N = simP.N;
    double DeltaR = simP.dR;
    double DeltaZ = simP.dZ;

    double Binf = sqrt(8.0*PI/simP.betaCyl)*pow(simP.Rbdy,simP.nCyl);
    double alpha = sqrt(8.0*PI/simP.betaCyl); //Binf; // alpha^2 = 8*pi/beta_inf = 8*pi*(B^2/8*pi)/(rho_*c^2)_inft -> B^2
    //double Ccst = PI*simP.RhoTop[0]/2.0/(1.0+1.0/simP.betaCyl); // RhoTop[0] should be 1


    // Source vector
    int BdyRight = 2;
    for(s=0;s<M*N;s++)
    {
        if( Ridx(s,M) == 0 )
        {
            // We will manually set A = 0 on the left hand side
            Source(s) = 0.0;
        }
        else if( Ridx(s,M) == M-1 )
        {
            if(BdyRight==0) // setting A
            {
                // We will manually set A = B_inf * r / 2
                Source(s) = Binf*cPos(M-1,DeltaR)/2.0;

                // Manually set to cylinder + empty space (the cylinder solution assumes nCyl = 0.5)
                double Redge = curState.VContour[N-1];
                //Source(s) = (0.5*alpha*sqrt(simP.RhoTop[0])/Ccst)*log(Ccst*Redge*Redge+1.0)/simP.rL + 0.5*Binf*(simP.rL-Redge*Redge/simP.rL);
            }
            else if(BdyRight==1) // setting dA/dr
            {

            }
            else // setting cros(A)
            {
                // We will impose a boundary condition (1/r)*d(rA)/dr = B_infty at this boundary
                Source(s) = 2.0*DeltaR*Binf;
            }
        }
        else
        {
            double dQdPhi = curState.State[dQdP][Ridx(s,M)][Zidx(s,M)];

            if( cPos(Ridx(s,M),DeltaR) <= curState.VContour[Zidx(s,M)]) Source(s) = -4.0*PI*pow(DeltaR,2)*cPos(Ridx(s,M),DeltaR)*dQdPhi*exp(-curState.State[Vpot][Ridx(s,M)][Zidx(s,M)]);
            else Source(s) = 0;

        }
    }


    // Finite Difference Matrix
    for(s=0;s<M*N;s++)
    {
        double wi,f,Mleft,Mright,Mup,Mdown,Mcenter;

        if( Ridx(s,M) == 0 )
        {
            // Left side will normalize A = 0, so all the entries here are 0
            Mleft = 0;
            Mright = 0;
            Mup = 0;
            Mdown = 0;
            Mcenter = 1;
        }
        else if(Ridx(s,M)==M-1)
        {
            if(BdyRight==0)
            {
                // Imposing value of A
                Mleft = 0;
                Mright = 0;
                Mup = 0;
                Mdown = 0;
                Mcenter = 1.0;
            }
            else if(BdyRight==1)
            {

            }
            else if(BdyRight==2)
            {
                // Right side implemented a Robin boundary condition
                Mleft = -4.0;
                Mright = 0; // doesn't exist
                Mup = 0;
                Mdown = 0;
                Mcenter = ( 2.0*DeltaR/cPos(M-1,DeltaR) + 3.0 );
                AmMatrix(s,s-2) = 1.0;
            }
        }
        else
        {
            // Base values
            wi = DeltaR/2.0/cPos(Ridx(s,M),DeltaR);
            f = DeltaR/DeltaZ;

            Mright = 1.0+wi;
            Mleft = 1.0-wi;
            Mup = pow(f,2);
            Mdown = pow(f,2);
            Mcenter = -2.0-2.0*pow(f,2)-4.0*pow(wi,2);
        }

        // Now we will adjust to incorproate z-boundary conditions

        // at bottom, dA/dz = 0
        if( Zidx(s,M) == 0 )
        {
            Mup = Mup + Mdown;
            Mdown = 0;
        }

        // at top A = cylinder?
        if( Zidx(s,M) == N-1)
        {
            if(Ridx(s,M)<M-1)
            {
                // enforce cylinder
                Mdown = 0; Mup = 0; Mleft = 0; Mright = 0; Mcenter = 1;
                Source(s) = simP.ATop[Ridx(s,M)];

                // classic Neumann
                //Mdown = Mup + Mdown;
                //Mup = 0;

                //if(Ridx(s,M)==M-1) cout << "Value of ATop in upper-right = " << simP.ATop[M-1] << endl;
            }
            else
            {
                Mdown = Mup + Mdown;
                Mup = 0;
            }
        }


        // Now we fill the entries of the matrix with the values
        AmMatrix(s,s) = Mcenter;
        if(Ridx(s,M)!=0) AmMatrix(s,s-1) = Mleft;
        if(Ridx(s,M)!=M-1) AmMatrix(s,s+1) = Mright;
        if(Zidx(s,M)!=N-1) AmMatrix(s,s+M) = Mup;
        if(Zidx(s,M)!=0) AmMatrix(s,s-M) = Mdown;
    }


};


void createAmpereMatrix(Params simP, TheState curState, MatrixXd& AmMatrix, VectorXd& Source)
{
    int s;

    int M = simP.M;
    int N = simP.N;
    double DeltaR = simP.dR;
    double DeltaZ = simP.dZ;

    double Binf = sqrt(8.0*PI/simP.betaCyl)*pow(simP.Rbdy,simP.nCyl);
    double alpha = sqrt(8.0*PI/simP.betaCyl); //Binf; // alpha^2 = 8*pi/beta_inf = 8*pi*(B^2/8*pi)/(rho_*c^2)_inft -> B^2
    //double Ccst = PI*simP.RhoTop[0]/2.0/(1.0+1.0/simP.betaCyl); // RhoTop[0] should be 1


    // Source vector
    int BdyRight = 2;
    for(s=0;s<M*N;s++)
    {
        if( Ridx(s,M) == 0 )
        {
            // We will manually set A = 0 on the left hand side
            Source(s) = 0.0;
        }
        else if( Ridx(s,M) == M-1 )
        {
            if(BdyRight==0) // setting A
            {
                // We will manually set A = B_inf * r / 2
                Source(s) = Binf*cPos(M-1,DeltaR)/2.0;

                // Manually set to cylinder + empty space (the cylinder solution assumes nCyl = 0.5)
                double Redge = curState.VContour[N-1];
                //Source(s) = (0.5*alpha*sqrt(simP.RhoTop[0])/Ccst)*log(Ccst*Redge*Redge+1.0)/simP.rL + 0.5*Binf*(simP.rL-Redge*Redge/simP.rL);
            }
            else if(BdyRight==1) // setting dA/dr
            {

            }
            else // setting cros(A)
            {
                // We will impose a boundary condition (1/r)*d(rA)/dr = B_infty at this boundary
                Source(s) = 2.0*DeltaR*Binf;
            }
        }
        else
        {
            double dQdPhi = curState.State[dQdP][Ridx(s,M)][Zidx(s,M)];

            if( cPos(Ridx(s,M),DeltaR) <= curState.VContour[Zidx(s,M)]) Source(s) = -4.0*PI*pow(DeltaR,2)*cPos(Ridx(s,M),DeltaR)*dQdPhi*exp(-curState.State[Vpot][Ridx(s,M)][Zidx(s,M)]);
            else Source(s) = 0;

        }
    }


    // Finite Difference Matrix
    for(s=0;s<M*N;s++)
    {
        double wi,f,Mleft,Mright,Mup,Mdown,Mcenter;

        if( Ridx(s,M) == 0 )
        {
            // Left side will normalize A = 0, so all the entries here are 0
            Mleft = 0;
            Mright = 0;
            Mup = 0;
            Mdown = 0;
            Mcenter = 1;
        }
        else if(Ridx(s,M)==M-1)
        {
            if(BdyRight==0)
            {
                // Imposing value of A
                Mleft = 0;
                Mright = 0;
                Mup = 0;
                Mdown = 0;
                Mcenter = 1.0;
            }
            else if(BdyRight==1)
            {

            }
            else if(BdyRight==2)
            {
                // Right side implemented a Robin boundary condition
                Mleft = -4.0;
                Mright = 0; // doesn't exist
                Mup = 0;
                Mdown = 0;
                Mcenter = ( 2.0*DeltaR/cPos(M-1,DeltaR) + 3.0 );
                AmMatrix(s,s-2) = 1.0;
            }
        }
        else
        {
            // Base values
            wi = DeltaR/2.0/cPos(Ridx(s,M),DeltaR);
            f = DeltaR/DeltaZ;

            Mright = 1.0+wi;
            Mleft = 1.0-wi;
            Mup = pow(f,2);
            Mdown = pow(f,2);
            Mcenter = -2.0-2.0*pow(f,2)-4.0*pow(wi,2);
        }

        // Now we will adjust to incorproate z-boundary conditions

        // at top and bottom, dA/dz = 0
        if( Zidx(s,M) == N-1)
        {
            Mdown = Mdown + Mup;
            Mup = 0;
        }

        if( Zidx(s,M) == 0 )
        {
            Mup = Mup + Mdown;
            Mdown = 0;
        }

        // Now we fill the entries of the matrix with the values
        AmMatrix(s,s) = Mcenter;
        if(Ridx(s,M)!=0) AmMatrix(s,s-1) = Mleft;
        if(Ridx(s,M)!=M-1) AmMatrix(s,s+1) = Mright;
        if(Zidx(s,M)!=N-1) AmMatrix(s,s+M) = Mup;
        if(Zidx(s,M)!=0) AmMatrix(s,s-M) = Mdown;
    }

}
