

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


int convergeTest(Params simP, TheState prevState, TheState curState)
{
    // Test to see if the solution has stopped changing.
    int done = 0;

    // find maximum error in dV/dR for interior cells
    double DVmax = -1; int iloc,jloc;
    for(int i=1;i<simP.M-1;i++) for(int j=1;j<simP.N-1;j++)
    {
        // find local dV/dr
        double plocDVDR = (prevState.State[Vpot][i+1][j]-prevState.State[Vpot][i-1][j]); // /2.0/simP.dR;
        double locDVDR = (curState.State[Vpot][i+1][j]-curState.State[Vpot][i-1][j]); // /2.0/simP.dR;

        double locerr = 0;
        if(plocDVDR!=0) locerr = fabs(locDVDR-plocDVDR)/fabs(plocDVDR);
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
        if(plocRHO!=0) locerr = fabs(locRHO-plocRHO)/fabs(plocRHO);
        if(locerr>DVmax3 )
        {
            DVmax3 = locerr;
            iloc = i; jloc = j;
        }

    }

    cout << "Max error for Rho located at " << iloc << "," << jloc << " : " << DVmax3 << endl;



    cout << "Central density here is " << curState.State[Q][0][0] * exp(-curState.State[Vpot][0][0]) << endl;

    if(DVmax <= simP.convergeLoopTol and DVmax2 <= simP.convergeLoopTol and DVmax3 <= simP.convergeLoopTol) done = 1;
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
