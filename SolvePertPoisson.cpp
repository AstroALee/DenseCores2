
#include "SolvePertPoisson.H"

// For some reason this has to be included here and not in the header file ... not sure why
#include <Eigen/Dense>
using namespace Eigen;

void createPertPoissonMatrix(Params simP, TheState curState, MatrixXd& PpmMatrix, VectorXd& Source);

void SolvePertPoisson(Params simP,TheState curState,TheState& newState)
{
    int M = curState.M;
    int N = curState.N;

    // First allocate the matrix and vectors for the finite difference scheme, source, and solution
    MatrixXd PoMatrix = MatrixXd::Zero((M)*(N),(M)*(N));
    VectorXd Source = VectorXd::Zero((M)*(N));
    VectorXd Soln = VectorXd::Zero((M)*(N));

    cout << "Creating Perturbed Poisson matrix." << endl;
    createPertPoissonMatrix(simP,curState,PoMatrix,Source);

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
    double Vcorner = newState.State[Vpot][0][N-1];
    for(int i=0;i<M;i++) for(int j=0;j<N;j++) newState.State[Vpot][i][j] = newState.State[Vpot][i][j] - Vcorner;





};


void createPertPoissonMatrix(Params simP, TheState curState, MatrixXd& PpmMatrix, VectorXd& Source)
{

};
