
#include "DCores_Main.H"

int main(int argc, char *argv[])
{
    cout << endl;

    // In case we need random numbers
    srand(time(NULL));

    // Set precision (give me all the numbers)
    cout << std::setprecision(16);

    // Start the clock
    clock_t startTime = clock();

    // Error flag
    int errTag = 0;

    // Read in INI information
    Params SimParams;
    errTag = ReadINI(SimParams,argc,argv);
    if(errTag) { Waterloo("main","error reading in info"); return 1; }

    errTag = CheckArgs(SimParams);
    if(errTag) { Waterloo("main","error checking arguments"); return 1;}

    // Derived units
    errTag = Derived(SimParams);
    if(errTag) { Waterloo("main","error calculating derived units"); return 1; }

    // Allocate state
    TheState curState(SimParams.M,SimParams.N);

    // Initial Conditions
    //errTag = SetIC(SimParams,curState);
    errTag = SetICnoPoints(SimParams,curState);
    if(errTag) { Waterloo("main","error setting initial conditions"); return 1; }

    // If restarting, checks that the setup is correct and overwrites state
    if(SimParams.restart)
    {
        errTag = LoadRestart(SimParams,curState);
        if(errTag) { Waterloo("main","restarting from file with different parameters"); return 1; }

        if(SimParams.squeeze)
        {
          SqueezeState(SimParams,curState);
        }
    }

    // Mass in box
    CylinderMassCheck(SimParams,curState);

    // Initial conditions set up
    PrintState(SimParams,curState,"_","0");

    // Find the steady state!
    //FindTheSteadyState(SimParams,curState);
    //FindTheSteadyStateWithPert(SimParams,curState);
    //FindTheSteadyStateWithSOR(SimParams, curState);
    //FindTheSteadyStateWithAlpha(SimParams, curState);
    //FindSSVQloopAlpha(SimParams, curState);
    //FindSSVfixedTopBot(SimParams, curState);
    FindSSGlobalRescaling(SimParams, curState);

    //PrintState(SimParams,curState,"_","1");

    // Mass in box
    CylinderMassCheck(SimParams,curState);




    // End the clock
    clock_t endTime = clock();
    ExecutionTime(startTime,endTime);

    // Let's go home!
    cout << endl;
    return 0;
};


void SqueezeState(Params& simP, TheState& curState)
{
  updateQ(simP, curState, 0);
  updateDQDPHI(simP, curState);
  // calc dMdPHI
  calcDMDPHI(simP, curState);

  // Now squeeze

  double fracSqueeze = 0.9;


  double PhiTemp[simP.M][simP.N];
  double VTemp[simP.M][simP.N];
  // create new values for Phi
  for(int i=0;i<simP.M;i++) for(int j=0;j<simP.N;j++)
  {
    double squeeze_alp = 1.0 - (( (1.0 - fracSqueeze)*curState.VContour[0]/curState.VContour[i] )*( 1.0 - cPos(i,simP.dZ)/simP.zL ));
    if(squeeze_alp<=0) squeeze_alp = 0.00000001;

    double rescaled_r = cPos(i,simP.dR)/squeeze_alp;

    // find indices for re-scaled r
    int lIdx = -1;
    for(int k=0;k<simP.M-1;k++)
    {
      lIdx = k ;
      if( cPos(k,simP.dR) >= rescaled_r) break; // if never found, k = M-2 (should be outside the bndy anyway)
    }

    // Find interpolation factor
    double int_factor = (cPos(lIdx,simP.dR)-rescaled_r)/simP.dR; // 0 <= <= 1

    // Use interpolated values in new location
    PhiTemp[i][j] =  ((1.0-int_factor)*cPos(lIdx,simP.dR)*curState.State[Apot][lIdx][j]) + (int_factor*cPos(lIdx+1,simP.dR)*curState.State[Apot][lIdx+1][j]);
    VTemp[i][j]   =  ((1.0-int_factor)*curState.State[Vpot][lIdx][j]) + (int_factor*curState.State[Vpot][lIdx+1][j]);
  }

  // Update contour location
  for(int i=0;i<simP.N;i++)
  {
    double squeeze_alp = 1.0 - (( (1.0 - fracSqueeze)*curState.VContour[0]/curState.VContour[i] )*( 1.0 - cPos(i,simP.dZ)/simP.zL ));
    curState.VContour[i] = squeeze_alp*curState.VContour[i];
  }

  // overwrite V and A (hence also phi)
  for(int i=0;i<simP.M;i++) for(int j=0;j<simP.N;j++)
  {
    if(cPos(i,simP.dR) < curState.VContour[j])
      curState.State[Vpot][i][j] = VTemp[i][j];
    else
      curState.State[Vpot][i][j] = 0.0;
    if(i>0)
      curState.State[Apot][i][j] = PhiTemp[i][j]/cPos(i,simP.dR);
    else
      curState.State[Apot][i][j] = 0.0;
  }







};


int LoadRestart(Params& simP, TheState& curState)
{
    int i,j,k;

    string ReadFile;
    ReadFile = simP.readinfname ;

    cout << "Reading curState from file " << ReadFile << endl;

    //OutFile = OutFileName.c_str() + num.c_str();

    ifstream myfile;
    myfile.open(ReadFile.c_str());

    string line;
    getline(myfile, line); // first line has various parameters
    stringstream linestream(line);

    string data;
    getline(linestream, data, ','); // Reads in M
    if( atoi(data.c_str()) != simP.M) {myfile.close(); cout << "ERROR with M!" << endl; return(1);}
    getline(linestream, data, ','); // Reads in N
    if( atoi(data.c_str()) != simP.N) {myfile.close(); cout << "ERROR with N!" << endl; return(1);}
    getline(linestream, data, ','); // Reads in zL
    if( (atof(data.c_str()) - simP.zL)/simP.zL > 1e-4) {myfile.close(); cout << "ERROR with zL!" << endl; return(1);}
    getline(linestream, data, ','); // Reads in rL
    if( (atof(data.c_str()) - simP.rL)/simP.rL > 1e-4) {myfile.close(); cout << "ERROR with rL!" << endl; return(1);}
    getline(linestream, data, ','); // Reads in mExcess
    getline(linestream, data, ','); // Reads in betaCyl
    if( (atof(data.c_str()) - simP.betaCyl)/simP.betaCyl > 1e-4) {myfile.close(); cout << "ERROR with betaCyl!" << endl; return(1);}
    getline(linestream, data, ','); // Reads in nCyl
    if( (atof(data.c_str()) - simP.nCyl)/simP.nCyl > 1e-4) {myfile.close(); cout << "ERROR with nCyl!" << endl; return(1);}
    getline(linestream, data, ','); // Reads in Rbdy
    getline(linestream, data, ','); // Reads in lambda
    if( (atof(data.c_str()) - simP.lambda)/simP.lambda > 1e-4) {myfile.close(); cout << "ERROR with lambda!" << endl; return(1);}
    getline(linestream, data, ','); // Reads in Pc2Code
    getline(linestream, data, ','); // Reads in Sol2Code

    /*
    // First print inputs as the first line
    myfile << simP.M << "," << simP.N << "," << simP.zL << "," << simP.rL << ","
        << simP.mExcess << "," << simP.betaCyl << "," << simP.nCyl << "," << simP.Rbdy << ","
        << simP.lambda << "," << simP.Pc2Code << "," << simP.Sol2Code << endl;
    */


    getline(myfile, line); // second line has contour boundary, which might change

    for(i=0;i<simP.M;i++) for(j=0;j<simP.N;j++) // this should match the file output order
    {
        getline(myfile, line); // grab next line
        stringstream linestream(line); // make it a string stream

        // parse
        // first four entries aren't needed
        getline(linestream, data, ',');getline(linestream, data, ',');
        getline(linestream, data, ',');getline(linestream, data, ',');

        // next two are meaningful
        getline(linestream, data, ','); // Vpot
        if( cPos(i,simP.dR) >= curState.VContour[j] ) curState.State[Vpot][i][j] = 0; // outside
        else if ( cPos(i+1,simP.dR) >= curState.VContour[j] ) // next to
        {
            // interpolate so Vbdy is still V = 0 contour.
            // assumes curState.State[Vpot][i-1][j] has already been read in.
            double m = (0 - curState.State[Vpot][i-1][j])/(simP.dR + (curState.VContour[j] - cPos(i,simP.dR)) );
            curState.State[Vpot][i][j] = curState.State[Vpot][i-1][j] + m*(simP.dR) ;
        }
        else curState.State[Vpot][i][j] = atof(data.c_str());

        getline(linestream, data, ','); // Apot
        curState.State[Apot][i][j] = atof(data.c_str());


    }




    myfile.close();

    updateQ(simP, curState, 0);
    updateDQDPHI(simP, curState);

    return(0);
};

int ReadINI(Params& SimParams, int argc, char** filename)
{
    if(argc==1) { return 1;} // no file name provided

    // Reading first command line input as file name
    INIReader reader(filename[1]);

    if (reader.ParseError() < 0) {
        std::cout << "Can't load input file!\n";
        return 1;
    }

    SimParams.M = reader.GetInteger("simulation","RadGrid",-1); // family, parameter, default if none found
    SimParams.N = reader.GetInteger("simulation","ZedGrid",-1);

    SimParams.rL = reader.GetReal("simulation","RadLength",-1.0);
    SimParams.zL = reader.GetReal("simulation","ZedLength",-1.0);

    SimParams.convergeLoopMax = reader.GetInteger("simulation","loopMax",-1);
    SimParams.convergeLoopTol = reader.GetReal("simulation","loopTol",-1.0);

    SimParams.relaxfrac = reader.GetReal("simulation","relaxfrac",-1.0);

    SimParams.denratio = reader.GetReal("simulation","dRatio",-1.0);

    // filament boundary (delPhi)
    SimParams.delPhi = reader.GetReal("boundary","delPhi",0.0);

    // Cylinder info
    SimParams.betaCyl = reader.GetReal("cylinder","beta",-1.0);
    SimParams.nCyl = reader.GetReal("cylinder","nCyl",-1.0);
    SimParams.lambda = reader.GetReal("cylinder","lambda",-1.0);
    SimParams.rTozPhys = reader.GetReal("cylinder","rTozPhys",-1.0);


    // Waistband
    SimParams.wRatio = reader.GetReal("waistband","wRatio",-1.0);

    // SOR
    SimParams.SORtol = reader.GetReal("SOR","looptol",-1.0);
    SimParams.SORomega = reader.GetReal("SOR","omega",1.0);
    SimParams.SORnum = reader.GetInteger("SOR","loopmax",1);


    // Point masses Vknob
    SimParams.Vknob = reader.GetReal("pointmasses","Vknob",-1.0);
    SimParams.mExcess = reader.GetReal("pointmasses","mExcess",-1.0);
    SimParams.pointLoopMax = reader.GetInteger("pointmasses","loopMax",-1);
    SimParams.pointLoopTol = reader.GetReal("pointmasses","loopTol",-1.0);
    SimParams.alpha_fac = reader.GetReal("pointmasses","alpha_fac",0.0);

    // outputs
    SimParams.outfname = reader.Get("output","filename","DEFAULTNAME");
    SimParams.restart = reader.GetInteger("output","restart",0);
    SimParams.squeeze = reader.GetInteger("output","squeeze",0);
    SimParams.readinfname = reader.Get("output","inputfn","DEFAULTNAME");


    std::cout << "Config loaded from " << filename[1] << ":" << endl
              << "M=" << SimParams.M << ", N=" << SimParams.N << endl
              << "rL=" << SimParams.rL << ", zL=" << SimParams.zL << "   (parsecs)" << endl;

    cout << "File name = " << SimParams.outfname << endl;


    // All went well, hopefully.
    return 0;
};

int CheckArgs(Params SimP)
{
    int bad = 0;

    if( SimP.M < 3 || SimP.N < 3)
    {
        cout << "The number of grid cells is weird. "
             << "You provided (M,N) = (" << SimP.M << "," << SimP.N << ")" << endl;
             bad = 1;
    }

    if( SimP.rL < 0 || SimP.zL < 0)
    {
        cout << "The size of the box is weird. "
             << "You provided (rL,zL) = (" << SimP.rL << "," << SimP.zL << ")" << endl;
             bad = 1;
    }

    // all went well
    return bad;

};

int Derived(Params& SimP)
{

    // Set cell sizes (units currently same as rL and zL -- parsecs)
    // The minus 1 allows the origin and right/top edges to be a grid point.
    SimP.dR = SimP.rL / ((double) SimP.M - 1.0);
    SimP.dZ = SimP.zL / ((double) SimP.N - 1.0);
    //SimP.dR = SimP.rL / ((double) SimP.M );
    //SimP.dZ = SimP.zL / ((double) SimP.N );

    // Allocate
    SimP.RhoTop   = new double[SimP.M];
    SimP.VTop     = new double[SimP.M];
    SimP.ATop     = new double[SimP.M];


    // All went well
    return 0;
};

void ExecutionTime(clock_t start, clock_t end )
{
    cout << "Execution took " << double( end - start ) / (double)CLOCKS_PER_SEC << " seconds." << endl;
};
