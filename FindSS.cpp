

#include "FindSS.H"


void FindTheSteadyState(Params simP, TheState& curState)
{
    // Need memory for the previous state and updated state
    TheState preState(simP.M,simP.N);
    TheState newState(simP.M,simP.N);
    preState = curState;

    for(int i = 0; i<simP.convergeLoopMax; i++)
    {
        cout << "Performing single full update " << i+1 << " (max allowed: " << simP.convergeLoopMax << ")" << endl;
        preState = curState;
        SingleFullUpdate(simP,curState,newState);
        if(convergeTest(simP, preState, curState)) break;
        if(i==simP.convergeLoopMax-1) cout << "Convergence never reached!!" << endl;
    }

};

void FindTheSteadyStateWithPert(Params simP, TheState& curState)
{
    // Need memory for the previous state and updated state
    TheState preState(simP.M,simP.N);
    TheState newState(simP.M,simP.N);
    preState = curState;

    // If this is the first time we're doing this, we need two regular updates
    // before going into perturbation scheme.
    SingleFullUpdate(simP,curState,newState);
    preState = curState;
    SingleFullUpdate(simP,curState,newState);

    // At this point, curState and prevState are different, and both are solutions from
    // regular Poisson/Ampere solves. Now we can enter the perturbation scheme to try
    // and settle toward the equilibrium solution.
    for(int i = 0; i<simP.convergeLoopMax; i++)
    {
        cout << "Performing single perturbation update " << i+1 << " (max allowed: " << simP.convergeLoopMax << ")" << endl;
        PertFullUpdate(simP,curState,preState,newState);

        if(convergeTest(simP, preState, curState)) break;
        if(i==simP.convergeLoopMax-1) cout << "Perturbation: Convergence never reached!!" << endl;
    }



};

void PertFullUpdate(Params simP,TheState& curState, TheState prevState, TheState& newState)
{
    // First we solve the perturbative Poisson equation, which might require several iterations itself
    for(int j=0; j<simP.convergeLoopMax; j++)
    {
        SolvePertPoisson(simP, curState, prevState, newState);

        // newState now has a new set of values for Delta V. Needs to be smaller than a certain threshold
        int good = DeltaVcheck(simP,newState);
        if(good) { Blend(simP,curState,newState); break;}

        // The update was not good. Blend the curState and the prevState
        cout << "Perturbation Poisson solve " << j << " was not good. Blending back and trying again." << endl;
        curState.average(0.5,prevState);
        updateQ(simP, curState, 0);
        updateDQDPHI(simP, curState);
    }

    // If here, pert Poisson has converged. newState[Vpot] has new Vpot values (not just DeltaV)
    SolveAmpere(simP,curState,newState); // newState now has new Ampere values
    curState += newState; // relaxes
    updateQ(simP, curState, 0);
    updateDQDPHI(simP, curState);

    // The full state has been updated

};

void SingleFullUpdate(Params simP,TheState curState,TheState& newState)
{
    cout << "Starting single state update." << endl;
    SolvePoisson(simP,curState,newState);
    SolveAmpere(simP,curState,newState);
    curState += newState; // relaxes
    //setVbdy(simP, curState); // comment this line out and the boundary location never changes.
    updateQ(simP, curState, 0);
    updateDQDPHI(simP, curState);
};

void Blend(Params simP,TheState curState,TheState& newState)
{
    // DeltaV is good, so the newState values of DeltaV can be added with curState to get new V values
    for(int i=0;i<simP.M;i++) for(int j=0;j<simP.N;j++) newState.State[Vpot][i][j] += curState.State[Vpot][i][j];

}

int DeltaVcheck(Params simP,TheState newState)
{
    int good = 1; // assumes we're good to go.

    double maxval = 0;
    for(int i=0;i<simP.M;i++) for(int j=0;j<simP.N;j++)
    {
        double locval = fabs(newState.State[Vpot][i][j]);

        if( locval > maxval ) maxval = locval;
    }

    cout << "Maximum DeltaV value is " << maxval << endl;
    if(maxval > 0.1) good = 0; // changes too large...

    return(good);
};

int convergeTest(Params simP, TheState prevState, TheState curState)
{
    // Test to see if the solution has stopped changing.
    int done = 0;

    // find maximum error in dV/dR for interior cells
    double DVmax = -1; int iloc,jloc;
    for(int i=2;i<simP.M-1;i++) for(int j=1;j<simP.N-1;j++)
    {
        // find local dV/dr
        double plocDVDR = (prevState.State[Vpot][i+1][j]-prevState.State[Vpot][i-1][j]); // /2.0/simP.dR;
        double locDVDR = (curState.State[Vpot][i+1][j]-curState.State[Vpot][i-1][j]); // /2.0/simP.dR;

        double locerr = 0;
        if(plocDVDR!=0 and cPos(i,simP.dR)< curState.VContour[j]) locerr = fabs(locDVDR-plocDVDR)/fabs(plocDVDR); // don't care about outside filament
        if(locerr>DVmax)
        {
            DVmax = locerr;
            iloc = i; jloc = j;
        }

    }

    cout << "Max error for dV/dR located at " << iloc << "," << jloc << " : " << DVmax << endl;


    // Now for A
    iloc = 0; jloc = 0; double DVmax2 = -1;
    for(int i=1;i<simP.M-1;i++) for(int j=1;j<simP.N-1;j++)
    {
        // find local dV/dr
        double plocDADR = (prevState.State[Apot][i+1][j]-prevState.State[Apot][i-1][j]); // /2.0/simP.dR;
        double locDADR = (curState.State[Apot][i+1][j]-curState.State[Apot][i-1][j]); // /2.0/simP.dR;

        double locerr = 0;
        if(plocDADR!=0) locerr = fabs(locDADR-plocDADR)/fabs(plocDADR);
        if(locerr>DVmax2 )
        {
            DVmax2 = locerr;
            iloc = i; jloc = j;
        }

    }

    cout << "Max error for dA/dR located at " << iloc << "," << jloc << " : " << DVmax2 << endl;


    // Now for rho
    iloc = 0; jloc = 0; double DVmax3 = -1;
    for(int i=1;i<simP.M-1;i++) for(int j=1;j<simP.N-1;j++)
    {
        // find local rho
        double plocRHO = prevState.State[Q][i][j] * exp(-prevState.State[Vpot][i][j]);
        double locRHO = curState.State[Q][i][j] * exp(-curState.State[Vpot][i][j]);

        double locerr = 0;
        if(plocRHO!=0 and cPos(i,simP.dR)< curState.VContour[j]) locerr = fabs(locRHO-plocRHO)/fabs(plocRHO);
        if(locerr>DVmax3 )
        {
            DVmax3 = locerr;
            iloc = i; jloc = j;
        }

    }

    cout << "Max error for Rho located at " << iloc << "," << jloc << " : " << DVmax3 << endl;



    cout << "Central density here is " << curState.State[Q][0][0] * exp(-curState.State[Vpot][0][0]) << endl;

    //if(DVmax <= simP.convergeLoopTol and DVmax2 <= simP.convergeLoopTol and DVmax3 <= simP.convergeLoopTol) done = 1; // V, A, and rho need to be fixed
    if(DVmax3 <= simP.convergeLoopTol) done = 1; // only rho needs to have settled down
    return done;
};


void setVbdy(Params simP, TheState& curState)
{
    // Given V everywhre, find the V contour that goes through the lambda=simP.lambda filament boundary at the top

    // radius of filament
    double Vrad = simP.rCyl;

    // What is V at this location?
    double DeltaR = simP.dR;
    int i=0; while(true){ i++; if(cPos(i,DeltaR)>=Vrad) break; }
    double Vr = curState.State[Vpot][i][simP.N-1];
    double Vl = curState.State[Vpot][i-1][simP.N-1];
    double m = (Vr-Vl)/DeltaR;
    double desV = Vl + m*(Vrad-cPos(i-1,DeltaR));
    cout << "Desired V value is " << desV << " at " << Vrad << endl;

    // Top row is done
    curState.VContour[simP.N-1] = Vrad;

    // For all other rows, find radius where V = desV
    for(int j=simP.N-2;j>=0;j--)
    {
        i=0;
        while(true) { i++; if(curState.State[Vpot][i][j]>=desV) break; }
        Vr = curState.State[Vpot][i][j];
        Vl = curState.State[Vpot][i-1][j];
        m = DeltaR/(Vr-Vl);
        curState.VContour[j] = cPos(i-1,DeltaR) + m*(desV-Vl); // overwrites VContour
        //cout << "V radius at " << j << " is " << curState.VContour[j] << endl;
    }

    cout << "Done with V boundary " << endl;

};
