
#include "PrintState.H"

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
