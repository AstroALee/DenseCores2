
#include "MagCylinder.H"



void calcMagCylinder(Params& simP,TheState& curState,int dooutput)
{
    //Set for particular problem (independent variable = t)
    const int Nvar = 5;
    const double t_start = 0.0; // t here is the radius
    const double t_end   = max(3.0,1.5*simP.rL); // this should be large enough so that lambda = desiredLambda somewhere
                                                 // also should be larger than rL

    //ODE tolerances and such
    const double atolode = 1.0e-8;
    const double rtol = atolode;
    const double h1 = 0.001;
    const double hmin = 1.0e-6;
    derivs derv(simP.nCyl,simP.betaCyl);  // Only need one declaration of this

    int NRout = 10*simP.M;  // need a good number for accurate interpolations

    // Initialize intitial conditions
    VecDoub ystart(Nvar);
    // Rho start
    ystart[0] = 1.0;
    // Psi and Psi' start
    ystart[1] = 0;
    ystart[2] = 0;
    // A start
    ystart[3] = 0;
    // Mass/length start
    ystart[4] = 0;

    // Output
    Output out(NRout);

    // Integrate
    Odeint< StepperDopr5<derivs> > ode(ystart,t_start,t_end,atolode,rtol,h1,hmin,out,derv);
    ode.integrate();

    double rEdge = LIntY(out.xsave, out.ysave, simP.lambda, out.count, 4); // 4 = lambda, returns radius
    cout << "Value of rEdge is " << rEdge << endl;

    // Contour location (non-dim) of just the cylinder
    for(int i=0;i<simP.N;i++) curState.VContour[i] = rEdge;

    if(dooutput)
    {
        // The cylinder solution is good to go. Set rho (will not change) , non-dim
        for(int i=0 ; i<simP.M; i++) simP.RhoTop[i] = LInt(out.xsave,out.ysave,cPos(i,simP.dR),out.count,0); // 0 = Rho

        // Outside the filament, we use the analytic form for the potential, need V and dV/dR at this location
        double VEdge = LInt(out.xsave,out.ysave,rEdge,out.count,1); // Vpot
        double VEdgeD = LInt(out.xsave,out.ysave,rEdge,out.count,2); // d(Vpot)/dR
        double RTopEdge = LInt(out.xsave,out.ysave,rEdge,out.count,0); // Rho

        // Save this value
        simP.Rbdy = RTopEdge;

        // V outside has the form  V = C1*ln(r/rEdge) + C2
        // C1 and C2 chosen so V is smooth at rEdge (V(rE) and DV(rE) match)
        // C1 = rE*DV(rE)     C2 = V(rE)
        double C1 = VEdgeD*rEdge;
        double C2 = VEdge;
        simP.CylCsts[0] = C1;
        simP.CylCsts[1] = C2;

        int i,j;
        for(i=0;i<simP.M;i++) for(j=0;j<simP.N;j++)
        {
            double rPos = cPos(i,simP.dR);
            double zPos = cPos(j,simP.dZ);

            if(rPos <= rEdge)
            {
                // V potential
                curState.State[Vpot][i][j] = LInt(out.xsave,out.ysave,rPos,out.count,1); // 1 = Vpot

                // We solved for r*A, have to divide by r (if i=0, A = 0, as is A*r)
                if(i>0) curState.State[Apot][i][j] = LInt(out.xsave,out.ysave,rPos,out.count,3)/rPos; // 3 = Apot*r

            }
            else
            {
                // Analytic form for V
                curState.State[Vpot][i][j] = C1*log(rPos/rEdge) + C2 ;

                // Analytic form for A*r outside cylinder = Cst*ln() + 0.5*Binf*(r^2-Redge^2)
                double Binf = sqrt(8.0*PI/simP.betaCyl)*pow(RTopEdge,simP.nCyl);
                double Acyl = LInt(out.xsave,out.ysave,rEdge,out.count,3) + 0.5*Binf*( pow(rPos,2) - pow(rEdge,2) );
                Acyl = Acyl/rPos;
                curState.State[Apot][i][j] = Acyl;

            }

            // Debugging stuff
            // curState.State[Apot][i][j] = 0.97*curState.State[Apot][i][j];
        }

        // Want Vpot = 0 at top left corner
        // Should have that by default at this point
        double Vnorm = curState.State[Vpot][0][simP.N-1];

        Vnorm = VEdge; // Filament boundary is the V=0 contour.
        for(i=0;i<simP.M;i++) for(j=0;j<simP.N;j++) curState.State[Vpot][i][j] = curState.State[Vpot][i][j] - Vnorm;



    } // end of dooutput

};
