
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
    errTag = SetIC(SimParams,curState);
    if(errTag) { Waterloo("main","error setting initial conditions"); return 1; }

    // Mass in box
    CylinderMassCheck(SimParams,curState);

    // Initial conditions set up
    PrintState(SimParams,curState,"_","0");

    // Find the steady state!
    FindTheSteadyState(SimParams,curState);

    PrintState(SimParams,curState,"_","1");

    // Mass in box
    CylinderMassCheck(SimParams,curState);




    // End the clock
    clock_t endTime = clock();
    ExecutionTime(startTime,endTime);

    // Let's go home!
    cout << endl;
    return 0;
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


    // Cylinder info
    SimParams.betaCyl = reader.GetReal("cylinder","beta",-1.0);
    SimParams.nCyl = reader.GetReal("cylinder","nCyl",-1.0);
    SimParams.lambda = reader.GetReal("cylinder","lambda",-1.0);
    SimParams.rTozPhys = reader.GetReal("cylinder","rTozPhys",-1.0);

    // Waistband
    SimParams.wRatio = reader.GetReal("waistband","wRatio",-1.0);

    // Point masses
    SimParams.mExcess = reader.GetReal("pointmasses","mExcess",-1.0);
    SimParams.pointLoopMax = reader.GetInteger("pointmasses","loopMax",-1);
    SimParams.pointLoopTol = reader.GetReal("pointmasses","loopTol",-1.0);

    // outputs
    SimParams.outfname = reader.Get("output","filename","DEFAULTNAME");


    std::cout << "Config loaded from " << filename[1] << ":" << endl
              << "M=" << SimParams.M << ", N=" << SimParams.N << endl
              << "rL=" << SimParams.rL << ", zL=" << SimParams.zL << "   (parsecs)" << endl;

    cout << "File name = " << SimParams.outfname << endl;


    // All went well, hopefully.
    return 0;
};

int CheckArgs(Params SimP)
{
    if( SimP.M < 3 || SimP.N < 3)
    {
        cout << "The number of grid cells is weird. "
             << "You provided (M,N) = (" << SimP.M << "," << SimP.N << ")" << endl;
    }

    if( SimP.rL < 0 || SimP.zL < 0)
    {
        cout << "The size of the box is weird. "
             << "You provided (rL,zL) = (" << SimP.rL << "," << SimP.zL << ")" << endl;
    }

    // all went well
    return 0;

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


    // All went well
    return 0;
};

void ExecutionTime(clock_t start, clock_t end )
{
    cout << "Execution took " << double( end - start ) / (double)CLOCKS_PER_SEC << " seconds." << endl;
};


// Prints out the data to a file
void PrintState(Params simP, TheState curState, string divide, string num)
{
    int i,j,k;

    string OutFile;
    OutFile = simP.outfname + divide + num;

    cout << "Printing curState to file " << OutFile << endl;

    //OutFile = OutFileName.c_str() + num.c_str();

    ofstream myfile;
    myfile.open(OutFile.c_str());

    // First print inputs as the first line
    myfile << simP.M << "," << simP.N << "," << simP.zL << "," << simP.rL << ","
        << simP.mExcess << "," << simP.betaCyl << "," << simP.nCyl << "," << simP.Rbdy << ","
        << simP.lambda << "," << simP.Pc2Code << "," << simP.Sol2Code << endl;

    // Print VContour as the second line
    myfile << curState.VContour[0];
    for(j=1;j<simP.N;j++) myfile << "," << curState.VContour[j];
    myfile << endl;

    double DeltaR = simP.dR;
    double DeltaZ = simP.dZ;

    for(j=0; j< simP.M; j++) for(k=0;k<simP.N;k++)
    {
        double Rval = cPos(j,DeltaR);
        double Zval = cPos(k,DeltaZ);
        double curRho = curState.State[Q][j][k] * exp(-curState.State[Vpot][j][k]);
        if(cPos(j,DeltaR)>curState.VContour[k]) curRho = 0;

        myfile << j << "," << k << "," << Rval << "," << Zval;
        for(i = 0; i< NStates; i++) myfile << "," << curState.State[i][j][k];
        myfile << "," << curRho;
        myfile << endl;
    }
    myfile.close();
};
