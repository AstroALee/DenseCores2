
#include "SetIC.H"

int SetICnoPoints(Params& simP, TheState& curState)
{
  // Find the cylinder radius that has the desired lambda value
  calcMagCylinder(simP,curState,0); // 0 = only sets VContour location (non-dim)

  // Contour location at top will not change, save as parameter
  simP.rCyl = curState.VContour[0];

  // The contour location where lambda = simP.lambda helps set the non-dimensional units
  unitConvert(simP);

  // Sets initial conditions of cylinder for state (sets V and A)
  calcMagCylinder(simP,curState,1);

  // Get the first filament boundary
  getVbdyDeltaPhi(simP, curState);

  // Determine Q
  updateQ(simP, curState, 0);

  // Calculate dQ/dPhi
  updateDQDPHI(simP, curState);

  // Looks like we made it...
  cout << "Done setting up initial conditions!" << endl;

  // all went well
  return 0;

};




int SetIC(Params& simP, TheState& curState)
{

    // Find the cylinder radius that has the desired lambda value
    calcMagCylinder(simP,curState,0); // 0 = only sets VContour location (non-dim)

    // Contour location at top will not change, save as parameter
    simP.rCyl = curState.VContour[0];

    // The contour location where lambda = simP.lambda helps set the non-dimensional units
    unitConvert(simP);

    // Sets initial conditions of cylinder for state (sets V and A)
    calcMagCylinder(simP,curState,1);

    // Determine Q
    //updateQ(simP, curState, 0);
    //simP.alpha_fac = 0.25;

    // Adds in points (changes only V)
    double firstsmoothlength = 100.0; // should be decently large so the filament boundary is accurate (smooth length = zL/firstsmoothlength)
    initVpoints(simP, curState, 1.0, firstsmoothlength);

    // Get the first filament boundary
    getVbdy(simP, curState);

    //initVpoints(simP, curState, -1.0, firstsmoothlength); // not needed if stretchV is used.

    // Manually set filament boundary
    //setTheVbdy(simP, curState);

    // Stretches the cylinder solution at the top from r=0 to the local filament boundary
    // (overwrites V everywhere)
    stretchV(simP,curState);

    // not used anymore
    // if(simP.wRatio>1) stretchV(simP,curState);
    // if(simP.mExcess > 0) smoothV(simP,curState,firstsmoothlength);

    // Determine Q
    updateQ(simP, curState, 0);

    // Calculate dQ/dPhi
    updateDQDPHI(simP, curState);

    // Looks like we made it...
    cout << "Done setting up initial conditions!" << endl;

    // all went well
    return 0;
};

void unitConvert(Params& simP)
{
    // We have passed in rEdge, which is the non-dimensional radius of the cylinder
    // where lambda = simP.lambda. We want to associate this with a particular physical
    // radius, which we do using the value of rTozPhys from our ini file

    double rEdge = simP.rCyl;

    // Physical radius of the cylinder (cgs)
    double Rphy = (simP.rTozPhys)*(simP.zL)*Pc2Cm;

    // Rho unit (Rcyl^2/Rphy^2) * Cs^2/G
    double RhoTopCenter = pow(rEdge/Rphy,2)*pow(Cbackground,2)/Gcst;

    // Unit conversion, multiply by this to convert parsecs to code units
    double cgsCodeLength = Cbackground/sqrt(Gcst*RhoTopCenter);  // cm/code length
    double pcCodeLength = cgsCodeLength/Pc2Cm; // pc/code length
    simP.Pc2Code = 1/pcCodeLength; // code/pc length

    // Unit conversion, multiply by this to convert solar masses to code units
    double cgsCodeMass = pow(Cbackground,3)/sqrt(RhoTopCenter*pow(Gcst,3)); // g/code
    double solCodeMass = cgsCodeMass/Sol2G; // sol/code
    simP.Sol2Code = 1/solCodeMass;

    cout << "Unit conversions:" << endl;
    cout << "    Parsec to Code = " << simP.Pc2Code << endl;
    cout << "    Solar mass to Code = " << simP.Sol2Code << endl;
    cout << "    Cs^2/G in cgs = " << Cbackground*Cbackground/Gcst  << endl; 

    // Non-dimensionalize the physical lengths
    simP.zL = simP.zL*simP.Pc2Code;
    simP.rL = simP.rL*simP.Pc2Code;

    simP.dR = simP.dR*simP.Pc2Code;
    simP.dZ = simP.dZ*simP.Pc2Code;

    // Non-dimensionalize mass
    simP.mExcess = simP.mExcess * simP.Sol2Code;

};

void initVpoints(Params& simP, TheState& curState, double undo, double smooth)
{
    int M = simP.M;
    int N = simP.N;
    double DeltaR = simP.dR;
    double DeltaZ = simP.dZ;
    double zL = simP.zL;

    double PointLoopTol = simP.pointLoopTol;
    int PointLoopMax = simP.pointLoopMax;

    // Array for the point mass chain
    double PointV[M][N];
    double PointVold[M][N];

    int i,j;
    for(i=0;i<M;i++) for(j=0;j<N;j++) PointV[i][j] = 0.0;
    for(i=0;i<M;i++) for(j=0;j<N;j++) PointVold[i][j] = 0.0;

    // Loops, adding point masses to the chain. Breaks from additional points do not change potential
    int PPidx = 0; //1;// skips one particular because GM/r blows up for the idx=0 particle in the corner.
    double zP=0;//2.0*zL;

    // Mass in code units
    double Mcode = undo*simP.mExcess; //mPoints;
    //cout << "MCode = " << Mcode << endl;

    while(true)
    {

        PPidx++;

        // Adds potential
        for(i=0;i<M;i++) for(j=0;j<N;j++)
        {
            // Distance locations of point particles
            double X1 = sqrt( pow(cPos(i,DeltaR),2) + pow(cPos(j,DeltaZ)-zP,2) );
            double X2 = sqrt( pow(cPos(i,DeltaR),2) + pow(cPos(j,DeltaZ)+zP,2) );

            if(zP==0)
            {
                // (GM/X blows up, use a softening parameter)
                double Xsoft = zL/smooth; //1.0*sqrt(pow(DeltaR,2)+pow(DeltaZ,2));
                X1 = X1 + Xsoft;
            }

            PointV[i][j] = PointV[i][j] + -Mcode/X1;
            if(zP > 0) PointV[i][j] = PointV[i][j] + -Mcode/X2;


        }


        // Calculates max error
        if(PPidx>2)
        {
            double err=0;

            for(i=0;i<M;i++) for(j=0;j<N;j++)
            {
                double locerr = fabs(PointV[i][j]-PointVold[i][j])/fabs(PointV[i][j]);
                if(locerr > err) err = locerr;
            }

            if(PPidx%100==0 || err <= PointLoopTol) cout << "Loop " << PPidx << " : Number of Points = " << 1+(PPidx-1)*2 << " (Max error = " << err << ")" << endl;

            // Is the max error small enough?
            if(err <= PointLoopTol) { PointLoopMax = PPidx; break; }
        }

        // Prepares for another iteration
        for(i=0;i<M;i++) for(j=0;j<N;j++) PointVold[i][j] = PointV[i][j];

        // New z location of point masses
        zP = zP + 2.0*zL;


        if(PPidx==PointLoopMax)
        {
            // Did not converge
            //WaterlooHeader("PrepInit:InitPoints");
            cout << "Did not converge to a solution! (PointLoopMax = " << PointLoopMax << ")" << endl;
            exit(1);
        }

    } // end of while loop

    // If we're here, we have converged.
    simP.pointLoopTotNum = PointLoopMax; // has been over-written to be the number of iterations used
    cout << "Number of loop iterations required: " << PointLoopMax << endl;

    // Normalize so that the upper-left corner has Vpot=0
    double Vnorm = PointV[0][N-1];

    // actually normalize so V at edge of cylinder is V=0
    i = 0;
    while(true){ i++; if(cPos(i,simP.dR) >= simP.rCyl) break;}
    double x = (cPos(i,simP.dR) - simP.rCyl)/simP.dR ;
    Vnorm = x*PointV[i-1][N-1] + (1-x)*PointV[i][N-1];
    cout << "Initial Vpoints normalization is " << Vnorm << endl;

    for(i=0;i<M;i++) for(j=0;j<N;j++) PointV[i][j] = PointV[i][j] - Vnorm;

    // Add solution to the potential
    for(i=0;i<M;i++) for(j=0;j<N;j++) curState.State[Vpot][i][j] = curState.State[Vpot][i][j] + PointV[i][j];

};


double radFromZed(double rT, double w, double zNorm)
{
    // zNorm should be 0 <= z/zL <= 1
    if(zNorm < 0) zNorm = 0;
    if(zNorm > 1) zNorm = 1;

    // limits of x (mapping from finite limits instead of +-infty)
    double xEdge = 3.0;

    double xVal = -xEdge + (2.0*xEdge)*zNorm ;

    double rad = rT + rT*(w-1.0)*(1.0-tanh(xVal))/2.0;

    return(rad);
};


void setTheVbdy(Params simP, TheState& curState)
{
    for(int j=0;j<simP.N;j++)
    {
        double zNorm = cPos(j,simP.dZ)/simP.zL;
        curState.VContour[j] = radFromZed(simP.rCyl, simP.wRatio, zNorm);
    }
};

void getVbdyDeltaPhi(Params simP, TheState& curState)
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

void getVbdy(Params simP, TheState& curState)
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

void stretchV(Params simP, TheState& curState)
{
    // Stretch the top gravitational potential out so the filament boundary (manually set) is
    // a potential contour.
    double Vedge = 0.0;

    double VOrig[simP.M];
    for(int i=0;i<simP.M;i++) VOrig[i] = Veval( i, simP.N-1, simP);

    //for(int i=0;i<simP.M;i++) cout << curState.State[Vpot][i][simP.N-1] << " , " << VOrig[i] << " , " << (curState.State[Vpot][i][simP.N-1]-VOrig[i])/curState.State[Vpot][i][simP.N-1] << endl;




    for(int j=0;j<simP.N;j++)
    {
        //double VOrig[simP.M];
        //for(int i=0;i<simP.M;i++) VOrig[i] = curState.State[Vpot][i][j];

        // Ratio of cylinder radius to filament boundary
        double rRat = curState.VContour[j] / simP.rCyl;

        for(int i=0;i<simP.M;i++)
        {
            if( cPos(i,simP.dR) >= curState.VContour[j] ) break; // outside the filament

            double curP = cPos(i,simP.dR);
            // Find the two stretched cells the current cell resides
            int kloc = 1;
            for(int k=1;k<simP.M;k++)
            {
                kloc = k;

                if( (curP >= rRat*cPos(k-1,simP.dR)) and (curP <= rRat*cPos(k,simP.dR)) )
                {
                    //cout << "Breaking i,j = " << i << "," << j << " at k = " << kloc << endl;
                    //cout << "     " << curP << " " << rRat*cPos(k-1,simP.dR) << " " <<  (rRat*cPos(k,simP.dR)) << endl;
                    break;
                }

            }

            double m = ( VOrig[kloc] - VOrig[kloc-1] ) / ( rRat*simP.dR );
            curState.State[Vpot][i][j] = VOrig[kloc-1] + m*( cPos(i,simP.dR) - rRat*cPos(kloc-1,simP.dR) );
            //cout << "     " << curState.State[Vpot][i][j] << " " << curState.State[Vpot][kloc-1][j] << " " <<  curState.State[Vpot][kloc][j] << endl;


        }

    }

    /*
    double Vdepth = VOrig[0];
    double lVball = cPos(simP.M/2,simP.dR);
    for(int i=0;i<simP.M;i++)
        for(int j=0;j<simP.N;j++)
        {
            double curPos = sqrt( pow(cPos(i,simP.dR),2) + pow(cPos(j,simP.dZ),2) );
            if(curPos < lVball)
            {
                curState.State[Vpot][i][j] = curState.State[Vpot][i][j] + 0.25*Vdepth*(1.0-curPos/lVball);
            }
        }
    */


    // At this point, V is stretched out in the radial direction
    // If we are using the Vknob formulation, we can also smooth the initial conditions toward the value of the knob
    // If alteration is V(r) = aR^2 + bR + c,
    // requiring V(0) = Vknob -->  c = Vknob
    // requiring dVdr(0) = 0 --> b = 0
    // requiring V(rBdy) = 0 --> a rBdy^2 + Vknob = 0 --> a = -Vknob / rBdy^2
    // where rBdy is the radius of the filament bounary at z = 0
    double Vknob = -1.0*simP.Vknob ;
    double rBdybttm = curState.VContour[0] ;
    double aCst = -Vknob/rBdybttm/rBdybttm;

    for(int i=0;i<simP.M;i++)
        for(int j=0;j<simP.N;j++)
        {
            double Vmax = aCst * pow(cPos(i,simP.dR),2) + Vknob;
            if(Vmax>0) continue;

            // No change at top, Vmax change at bottom
            // DeltaV = a + b*Z
            // DeltaV = 0 at top --> a + b*zL = 0 --> b = -a/zL
            // DeltaV = Vmax at bottom --> a = Vmax --> b = -Vmax/zL
            curState.State[Vpot][i][j] = curState.State[Vpot][i][j] + Vmax - (Vmax/simP.zL)*cPos(j,simP.dZ);

        }



    // we don't solve V outside the filament boundary, so let's just set it equal to 0 out there
    for(int i=0;i<simP.M;i++)
        for(int j=0;j<simP.N;j++)
            if(cPos(i,simP.dR)>=curState.VContour[j]) curState.State[Vpot][i][j] = Vedge;


};

void smoothV(Params simP,TheState& curState,double smooth)
{


    double mEx = simP.mExcess;

    double cylMass = 2.0*simP.zL*simP.lambda;
    double curMass = Trap2D(simP,curState);

    // undo the previous point potential
    initVpoints(simP, curState,-1.0,smooth);

    // potential is now just the cylinder
    // return here to have the potential be that of the cylinder
    // return;

    for(int i=0;i<simP.M;i++) for(int j=0;j<simP.N;j++)
    {
        double r0 = curState.VContour[0]/1.25;

        double r = sqrt( pow(cPos(i,simP.dR),2) + pow(cPos(j,simP.dZ),2) );
        if(r>=r0)
        {
            // do nothing
        }
        else
        {
            double Afac = 1.0*3.0*mEx/PI/r0;  // 1.25 factor added by hand, makes answers generally better.
            double rhoball = Afac*(1.0-r/r0);
            curState.State[Vpot][i][j] = -log( (simP.RhoTop[i] + rhoball)/curState.State[Q][i][j] );
        }
    }

    return;





    int NRcount = 0;
    int NRmax   = 100;
    double errTol = 0.001;

    double err = (curMass-mEx-cylMass)/(mEx+cylMass);

    if(fabs(err) <= errTol) return;

    // undo the previous point potential
    initVpoints(simP, curState,-1.0,smooth);
    updateQ(simP, curState, 0);

    double smooth_old = smooth;
    smooth = 0.8*smooth;

    // Else go into NR routine
    while(NRcount<NRmax)
    {
        // Add in new potential with new smoothing
        initVpoints(simP, curState,1.0,smooth);
        updateQ(simP, curState, 0);

        cylMass = 2.0*simP.zL*simP.lambda;
        curMass = Trap2D(simP,curState);

        double err_old = err;
        err = (curMass-mEx-cylMass)/(mEx+cylMass);
        cout << "NR routine " << NRcount << ": Error was " << err << endl;

        if(fabs(err) <= errTol) return;

        // undo the previous point potential
        initVpoints(simP, curState,-1.0,smooth);
        updateQ(simP, curState, 0);

        double dedsmooth = (err-err_old)/(smooth-smooth_old);
        smooth_old = smooth ;
        smooth = smooth_old - err/dedsmooth;
        cout << "NR routine " << NRcount << ": New smoothing length is " << smooth << endl;
        if(smooth <= 0 or isnan(smooth)) {cout << "Resetting smooth length." << endl; smooth = 1.0;}

        NRcount++;
    }




}
