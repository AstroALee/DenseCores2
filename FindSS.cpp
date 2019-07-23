

#include "FindSS.H"

void FindSSGlobalRescaling(Params& simP, TheState& curState)
{
    // Need memory for the previous state and updated state
    TheState preState(simP.M,simP.N);
    TheState newState(simP.M,simP.N);
    preState = curState;

    double s0 = Finds0(0.1, simP, curState);
    cout << "First guess for s0 = " << s0 << endl;
    if(isnan(s0)) s0 = 0.1;

    double cylMass = 2.0*simP.zL*simP.lambda;
    double desMass = cylMass + simP.mExcess;
    double curMass =  Trap2D(simP,curState,s0);

    SolvePoissonS0(simP,curState,newState,s0);
    SolveAmpereS0(simP,curState,newState,s0);
    curState += newState; // relaxes
    updateQ(simP, curState, 0);
    updateDQDPHI(simP, curState);
    redrawVbdyDeltaPhi(simP, curState); // redraws filament boundary based on delta phi only

    curMass = Trap2D(simP,curState,s0);
    cout << "Current mass in the box is " << curMass/simP.Sol2Code << " Sol" << endl;
    //cout << "Cylinder mass in the box is " << cylMass/simP.Sol2Code << " Sol" << endl;
    cout << "Desired mass in the box is " << desMass/simP.Sol2Code << " Sol" << endl << endl;


    for(int jj=0;jj<simP.convergeLoopMax;jj++)
    {
        s0 = Finds0(s0, simP, curState); //uses old value as first guess
        cout << "New guess for s0 = " << s0 << endl;
        if(jj<5) s0 = 10; //40.0 - 35.0*((double) jj)/5.0;
        if(isnan(s0) or s0<-1.0 or s0>70.0) s0 = 1.5;

        SolvePoissonS0(simP,curState,newState,s0);
        SolveAmpereS0(simP,curState,newState,s0);
        curState += newState; // relaxes
        updateQ(simP, curState, 0);
        updateDQDPHI(simP, curState);
        redrawVbdyDeltaPhi(simP, curState);
        curMass = Trap2D(simP,curState,s0);
        cout << "Current mass in the box is " << curMass/simP.Sol2Code << " Sol" << endl;
        //cout << "Cylinder mass in the box is " << cylMass/simP.Sol2Code << " Sol" << endl;
        cout << "Desired mass in the box is " << desMass/simP.Sol2Code << " Sol" << endl << endl;
    }

    cout << "done solving Poisson" << endl;
    s0 = Finds0(s0, simP, curState); //uses old value as first guess
    curMass = Trap2D(simP,curState,s0);
    cout << "Current mass in the box is " << curMass/simP.Sol2Code << " Sol" << endl;
    //cout << "Cylinder mass in the box is " << cylMass/simP.Sol2Code << " Sol" << endl;
    cout << "Desired mass in the box is " << desMass/simP.Sol2Code << " Sol" << endl << endl;

    PrintState(simP,curState,"_",to_string(1));
    //PrintState(simP,curState,"_",to_string(i+1));

};

double Finds0(double firsts0, Params simP, TheState curState)
{
    double s0 = firsts0; // initial guess for s0

    double cylMass = 2.0*simP.zL*simP.lambda;
    double desMass = cylMass + simP.mExcess;
    double curMass =  Trap2D(simP,curState,s0);

    // Find the initial value of s_0 that gives the right mass
    double prevs0, preverr;
    for(int i=0;i<20;i++)
    {
        double err = fabs(curMass - desMass)/desMass;
        //cout << "s0 = " << s0 << "  curMass = " << curMass/simP.Sol2Code << "  err = " << err << endl;

        if(err<0.001)
        {
            cout << "Converged: s0 = " << s0 << "  curMass = " << curMass/simP.Sol2Code << "  Mass err = " << err << endl;
            double Nos0mass = Trap2D(simP,curState,0.0);
            double MerrNos0 = fabs(Nos0mass - desMass)/desMass;
            cout << "   Mass err without s0 = " << MerrNos0 << endl;
            break;
        }


        if(i==0)
        {
            prevs0 = s0;
            preverr = err;

            s0 = 1.1*s0;
        }
        else
        {
            double derr = (err-preverr)/(s0-prevs0);
            prevs0 = s0;
            preverr = err;

            s0 = prevs0 - err/derr;
        }
        curMass = Trap2D(simP,curState,s0);
    }


    return s0;
}

void FindSSVfixedTopBot(Params& simP, TheState& curState)
{
    // Need memory for the previous state and updated state
    TheState preState(simP.M,simP.N);
    TheState newState(simP.M,simP.N);
    preState = curState;

    // We will employ the magnetized cylinder for the top
    // and an augmented cylinder for the bottom
    //updateQ(simP, curState, 0);
    //updateDQDPHI(simP, curState);
    //PrintState(simP,curState,"_","0"); // initial state

    for(int i = 0; i<simP.convergeLoopMax; i++)
    {
        //SolvePoissonSOR(simP,curState,newState);
        //SolveAmpereSOR(simP,curState,newState);
        //updateQ(simP, curState, 0);
        //updateDQDPHI(simP, curState);

        SolvePoisson(simP,curState,newState);
        SolveAmpere(simP,curState,newState);
        curState += newState; // relaxes
        updateQ(simP, curState, 0);
        updateDQDPHI(simP, curState);

        PrintState(simP,curState,"_",to_string(i+1));

    }



}

void FindSSVQloopAlpha(Params& simP, TheState& curState)
{
    // Need memory for the previous state and updated state
    TheState preState(simP.M,simP.N);
    TheState newState(simP.M,simP.N);
    preState = curState;

    // reset top V value and boundary location
    //setVbdyAlpha(simP, curState);
    //resetV(simP,curState);
    updateQ(simP, curState, 0);
    updateDQDPHI(simP, curState);

    PrintState(simP,curState,"_","0");

    for(int i = 0; i<simP.convergeLoopMax; i++)
    {
        SolveAmpereSOR(simP,curState,newState);
        curState.relaxA(newState);
        updateQ(simP, curState, 0);
        updateDQDPHI(simP, curState);

        PrintState(simP,curState,"_",to_string(i+1));


    }

    for(int i = simP.convergeLoopMax; i<2*simP.convergeLoopMax; i++)
    {
        ConvergeOnV(simP,curState,newState); // Solve V and Q for a given A, adjusting alpha until the mass is right. (also updates Q)
        updateQ(simP, curState, 0);
        updateDQDPHI(simP, curState);
        PrintState(simP,curState,"_",to_string(simP.convergeLoopMax+1));
    }

    /*
    for(int i = 0; i<simP.convergeLoopMax; i++)
    {
        cout << "Starting single full alpha update " << i+1 << " (max allowed: " << simP.convergeLoopMax << ")" << endl;

        preState = curState;
        ConvergeOnV(simP,curState,newState); // Solve V and Q for a given A, adjusting alpha until the mass is right.

        //for(int j = 0; j<simP.M;j++) for(int k =0; k<simP.N; k++) curState.State[Vpot][j][k] = Veval(j,simP.N-1,simP);

        setVbdyAlpha(simP, curState);
        PrintState(simP,curState,"_",to_string(i+1));
    }
    */

}


void ConvergeOnV(Params& simP,TheState curState,TheState& newState)
{
    double alpha_prev;
    double mEx_prev;
    for(int i = 0 ; i < 1; i++)
    {
        cout << "Converging on V iteration " << i+1 << " of 1" << endl;
        SolvePoissonSOR(simP,curState,newState);
        curState.relaxV(newState);
        //curState += newState; // relaxes
        updateQ(simP, curState, 0);
        updateDQDPHI(simP, curState);


        // Checks mass in box
        double curMass =  Trap2D(simP,curState);
        double cylMass = 2.0*simP.zL*simP.lambda;
        double curMexcess = curMass - cylMass; // shouldn't ever be negative ...
        cout << "Current mass excess with alpha " << simP.alpha_fac << " is " << curMexcess << endl;

        if(i==0) // prepare the next state manually
        {
            alpha_prev = simP.alpha_fac;
            mEx_prev = curMexcess;
            // need a new value of alpha_fac
            if(curMexcess - simP.mExcess < 0) simP.alpha_fac = 1.1*simP.alpha_fac;
            else  simP.alpha_fac = 0.9*simP.alpha_fac;

            //override
            simP.alpha_fac = 2.0;


        }
        else // NR a new value of alpha
        {
            // function to find root:  cur_ mExcess - desired_mExcess
            // new alpha = old alpha - function/function'
            double f = curMexcess - simP.mExcess;
            double fp = (curMexcess - mEx_prev)/(simP.alpha_fac - alpha_prev);
            cout << "alpha_prev = " << alpha_prev << "  mEx_prev = " << mEx_prev << endl;
            cout << "f = " << f << "  , " << "fp = " << fp << endl;

            alpha_prev = simP.alpha_fac;
            mEx_prev = curMexcess;

            simP.alpha_fac = simP.alpha_fac - f/fp;

            if(simP.alpha_fac<0) simP.alpha_fac = 0.0;

            //override
            simP.alpha_fac = 2.0;
            cout << "New alpha value = " << simP.alpha_fac << endl;

        }

        // reset top V value and boundary location
        // setVbdyAlpha(simP, curState);
        resetV(simP,curState);

        updateQ(simP, curState, 0);
        updateDQDPHI(simP, curState);

        if( fabs(curMexcess - simP.mExcess)/simP.mExcess < simP.convergeLoopTol ) break; // we win!
    }

}



void FindTheSteadyStateWithAlpha(Params& simP, TheState& curState)
{
    // Need memory for the previous state and updated state
    TheState preState(simP.M,simP.N);
    TheState newState(simP.M,simP.N);
    preState = curState;

    // reset top V value
    setVbdyAlpha(simP, curState);
    resetV(simP,curState);
    updateQ(simP, curState, 0);
    updateDQDPHI(simP, curState);

    PrintState(simP,curState,"_","0");

    double alpha_prev;
    double mEx_prev;
    for(int i = 0; i<simP.convergeLoopMax; i++)
    {
        cout << "Starting single full alpha update " << i+1 << " (max allowed: " << simP.convergeLoopMax << ")" << endl;

        preState = curState;
        SingleSORUpdate(simP,curState,newState); // updates V and A and Q and dQdPhi , relaxes from newstate to curState
        PrintState(simP,curState,"_",to_string(i+1));

        // Checks mass in box
        double curMass =  Trap2D(simP,curState);
        double cylMass = 2.0*simP.zL*simP.lambda;
        double curMexcess = curMass - cylMass; // shouldn't ever be negative ...

        cout << "Desired mass excess = " << simP.mExcess << endl;
        cout << "Masses in box (" << i+1 << "): curMass = " << curMass << " , " << cylMass << " , " << curMexcess << endl;

        if( fabs(curMexcess - simP.mExcess)/simP.mExcess < simP.convergeLoopTol ) break; // we win!

        if(0) // make 0 to fix alpha throughout simulation
        {
            if(i==0) // prepare the next state manually
            {
                alpha_prev = simP.alpha_fac;
                mEx_prev = curMexcess;
                // need a new value of alpha_fac
                if(curMexcess - simP.mExcess < 0) simP.alpha_fac = 1.1*simP.alpha_fac;
                else  simP.alpha_fac = 0.9*simP.alpha_fac;


            }
            else // NR a new value of alpha
            {
                // function to find root:  cur_ mExcess - desired_mExcess
                // new alpha = old alpha - function/function'
                double f = curMexcess - simP.mExcess;
                double fp = (curMexcess - mEx_prev)/(simP.alpha_fac - alpha_prev);
                cout << "alpha_prev = " << alpha_prev << "  mEx_prev = " << mEx_prev << endl;
                cout << "f = " << f << "  , " << "fp = " << fp << endl;

                alpha_prev = simP.alpha_fac;
                mEx_prev = curMexcess;


                simP.alpha_fac = simP.alpha_fac - f/fp;

                if(simP.alpha_fac<0) simP.alpha_fac = 0;

            }
        }
        cout << "New alpha value = " << simP.alpha_fac << endl;

        // reset top V value
        //setVbdyAlpha(simP, curState); // comment out and boundary doesn't change
        //resetV(simP,curState);
        updateQ(simP, curState, 0);
        updateDQDPHI(simP, curState);

    } // broken out of a loop, converged on mExcess


};

void FindTheSteadyStateWithSOR(Params simP, TheState& curState)
{
    // Need memory for the previous state and updated state
    TheState preState(simP.M,simP.N);
    TheState newState(simP.M,simP.N);
    preState = curState;

    for(int i = 0; i<simP.convergeLoopMax; i++)
    {
        cout << "Performing single full update " << i+1 << " (max allowed: " << simP.convergeLoopMax << ")" << endl;
        preState = curState;
        SingleSORUpdate(simP,curState,newState);

        PrintState(simP,curState,"_",to_string(i+1));

        if(convergeTest(simP, preState, curState)) break;
        if(i==simP.convergeLoopMax-1) cout << "Convergence never reached!!" << endl;
    }


};

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

void SingleSORUpdate(Params simP,TheState curState,TheState& newState)
{
    cout << "Starting single state update with SOR." << endl;

    // Solves Poisson using SOR
    SolvePoissonSOR(simP,curState,newState);
    //SolvePoisson(simP,curState,newState);

    // Solves Ampere using SOR
    SolveAmpereSOR(simP,curState,newState);
    //SolveAmpere(simP,curState,newState);

    curState += newState; // relaxes
    //setVbdy(simP, curState); // comment this line out and the boundary location never changes (unless called in SORAlpha routine)
    updateQ(simP, curState, 0);
    updateDQDPHI(simP, curState);

};


void redrawVbdyDeltaPhi(Params simP, TheState& curState)
{
    // Define the filament boundary as an analytic curve that starts at the
    // anchor point where lambda = simP.lambda at the top

    // radius of filament
    double rFil = simP.rCyl;

    // phi at this location
    int i=0;
    while(true){ i++; if(cPos(i,simP.dR)>=rFil) break; }
    double PhiR = cPos(i,simP.dR)*curState.State[Apot][i][simP.N-1];
    double PhiL = cPos(i-1,simP.dR)*curState.State[Apot][i-1][simP.N-1];
    double m = (PhiR-PhiL)/simP.dR;
    double phiAnchor = PhiL + m*(rFil - cPos(i-1,simP.dR));

    // Top row is done
    curState.VContour[simP.N-1] = rFil;

    // For all other rows, find radius where phi = function(z; delPhi)
    for(int j=simP.N-2;j>=0;j--)
    {
        double phiValue = phiAnchor + simP.delPhi*(1.0 - cPos(j,simP.dZ)/simP.zL);

        // Find the radius where phi = phiValue
        i=0;
        while(true) { i++; if(cPos(i,simP.dR)*curState.State[Apot][i][j] >= phiValue) break; }
        //cout << "Here! " << i << " " << j << endl;
        double PhiR = cPos(i,simP.dR)*curState.State[Apot][i][j];
        double PhiL = cPos(i-1,simP.dR)*curState.State[Apot][i-1][j];
        m = simP.dR/(PhiR-PhiL);
        curState.VContour[j] = cPos(i-1,simP.dR) + m*(phiValue-PhiL); // overwrites VContour
        //cout << "V radius at " << j << " is " << curState.VContour[j] << endl;
    }



    // What is V at this location?
    i=0; while(true){ i++; if(cPos(i,simP.dR)>=rFil) break; }
    double Vr = curState.State[Vpot][i][simP.N-1];
    double Vl = curState.State[Vpot][i-1][simP.N-1];
    m = (Vr-Vl)/simP.dR;
    double desV = Vl + m*(rFil-cPos(i-1,simP.dR));
    cout << "Desired V value is " << desV << endl;

    for(i=0;i<simP.M;i++)
      for(int j=0;j<simP.N;j++)
        if( curState.State[Vpot][i][j] >= desV ) curState.State[Vpot][i][j] = desV;


    cout << "Done with Phi-Boundary V boundary " << endl;





};

void SingleFullUpdate(Params simP,TheState curState,TheState& newState)
{
    cout << "Starting single state update." << endl;
    SolvePoisson(simP,curState,newState);
    SolveAmpere(simP,curState,newState);
    curState += newState; // relaxes
    setVbdy(simP, curState); // comment this line out and the boundary location never changes.
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
    for(int i=2;i<simP.M-1;i++) for(int j=1;j<simP.N-1;j++)
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

void resetV(Params simP, TheState& curState)
{
    for(int i=0;i<simP.M;i++)
        if(cPos(i,simP.dR) < curState.VContour[simP.N-1])
            curState.State[Vpot][i][simP.N-1] =  Veval(i,simP.N-1,simP);
        else
            curState.State[Vpot][i][simP.N-1] = 0.0;
};

void setVbdyAlpha(Params simP, TheState& curState)
{
    // Given V everywhre, find the V contour that goes through the lambda=simP.lambda
    // filament boundary at the top

    // radius of filament
    double Vrad = simP.rCyl;

    // value of V on filament
    double Vedge = 0.0;

    // temporary array
    double temp[simP.M][simP.N];

    // uses current (likely just updated) value of alpha_fac
    for(int i=0;i<simP.M;i++) for(int j=0;j<simP.N;j++) temp[i][j] = Veval(i,j,simP);

    // Now find the contour on this temp array

    // Top row is done
    curState.VContour[simP.N-1] = Vrad;

    // double check the value of V at desired location
    double DeltaR = simP.dR;
    int i=0; while(true){ i++; if(cPos(i,DeltaR)>=Vrad) break; }
    double Vr = temp[i][simP.N-1];
    double Vl = temp[i-1][simP.N-1];
    double m = (Vr-Vl)/DeltaR;
    double desV = Vl + m*(Vrad-cPos(i-1,DeltaR));
    cout << "CHECK: Desired V value is " << desV << " at " << Vrad << endl;

    //find radius where V = Vedge
    for(int j=simP.N-2;j>=0;j--)
    {
        int i=0;
        while(true) { i++; if(temp[i][j]>=Vedge) break; }
        double Vr = temp[i][j];
        double Vl = temp[i-1][j];
        double m = simP.dR/(Vr-Vl);
        curState.VContour[j] = cPos(i-1,simP.dR) + m*(Vedge-Vl); // overwrites VContour
        //cout << "V radius at " << j << " is " << curState.VContour[j] << endl;
    }

    cout << "Done with V boundary " << endl;

    for(int i=0;i<simP.M;i++) for(int j=0; j<simP.M;j++) if(cPos(i,simP.dR) > curState.VContour[j]) curState.State[Vpot][i][j] = Vedge;

};
