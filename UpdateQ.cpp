
#include "UpdateQ.H"

void updateQ(Params simP, TheState& curState, int type)
{

    cout << "Updating Q (type=" << type << ")" << endl;

    int i,j;
    double DeltaR = simP.dR;

    for(i=1;i<simP.M-1;i++) { if(cPos(i,DeltaR)>curState.VContour[simP.N-1]) break;}
    double Rr = simP.RhoTop[i];
    double Rl = simP.RhoTop[i-1];
    double m2 = (Rr-Rl)/DeltaR;
    double Rbdy = Rl + m2*(curState.VContour[simP.N-1]-cPos(i-1,DeltaR)); // updates global Rbdy
    //cout << "New Rbdy value = " << Rbdy << endl; // shouldn't ever change...


    if(type==1) // we know where the new boundary is and Rho has been updated
    {
        double rEdge = curState.VContour[simP.N-1];

        // Interpolate Rho to rEdge
        for(i=1;i<simP.M-1;i++) { if(cPos(i,DeltaR)>rEdge) break;}
        double Rr = simP.RhoTop[i];
        double Rl = simP.RhoTop[i-1];
        double m = (Rr-Rl)/DeltaR;
        Rbdy = Rl + m*(rEdge-cPos(i-1,DeltaR)); // updates global Rbdy
        cout << "New Rbdy value = " << Rbdy << endl;

        for(i=0;i<simP.M;i++) for(j=0;j<simP.N;j++)
        {
            if(cPos(i,DeltaR)>curState.VContour[j]) curState.State[Q][i][j] = Rbdy;
        }

        return; // we don't have to make maps and stuff
    }

    // Make map
    double PhiMap[simP.M];
    double QMap[simP.M];
    for(i=0;i<simP.M;i++) PhiMap[i] = cPos(i,DeltaR)*curState.State[Apot][i][simP.N-1]; // Phi at top
    for(i=0;i<simP.M;i++)   QMap[i] = simP.RhoTop[i] * exp( curState.State[Vpot][i][simP.N-1] ); // Q at top using RhoTop

    // What is Q at the filament boundary?
    i=0; while(true){ i++; if(cPos(i,DeltaR)>=curState.VContour[simP.N-1]) break; }
    double m = ( QMap[i] - QMap[i-1] )/DeltaR;
    double Qedge = QMap[i-1] + m*(curState.VContour[simP.N-1] - cPos(i-1,DeltaR));
    m = ( PhiMap[i] - PhiMap[i-1] )/DeltaR;
    double Phiedge = PhiMap[i-1] + m*(curState.VContour[simP.N-1] - cPos(i-1,DeltaR));

    // This will assume that the Phi contour through this point moves either vertically
    // down or inward. We will set Q = Qedge for all phi contours right of this contour.

    // Beyond the boundary, we use this Q value.
    for(i=0;i<simP.M;i++) if(cPos(i,DeltaR)>=curState.VContour[simP.N-1]) QMap[i] = Qedge;

    // Top row is done
    for(i=0;i<simP.M;i++) curState.State[Q][i][simP.N-1] = QMap[i];

    // Now for every other cell, use map to interpolate a value of Q
    for(i=0;i<simP.M;i++) for(j=0;j<simP.N-1;j++) // skipping the top row
    {
        // If right of the contour (V or Phi), who cares
        if(cPos(i,DeltaR) >= curState.VContour[j]){ curState.State[Q][i][j]=Qedge; continue; }
        if(cPos(i,DeltaR)*curState.State[Apot][i][j] >= Phiedge){ curState.State[Q][i][j]=Qedge; continue; }

        // did we succeed?
        int yay=0;

        // Local phi value we want to use with PhiMap and QMap
        double localPhi = cPos(i,DeltaR)*curState.State[Apot][i][j];
        int idx;

        // Find the two indices localPhi lies between (if not monotonic, could be weird)
        for(idx=1;idx<simP.M;idx++)
        {
            double Pleft = PhiMap[idx-1]; // cPos(idx-1,DeltaR)*curState[Apot][idx-1][j];
            double Pright= PhiMap[idx];   // cPos(idx,DeltaR)*curState[Apot][idx][j];

            if( (Pleft <= localPhi && localPhi <= Pright) || (Pleft >= localPhi && localPhi >= Pright) ) {yay=1; break;}
        }

        if(idx==simP.M-1 && yay==0)
        {
            cout << "Possibly didn't find a Phi value to interpolate in updateQ (i,j)=(" << i << "," << j <<")" << endl;
            // Possibly didn't find, but we know Phi analytically
            // Assumes we are beyond the boundary
        }

        // We now know which two Phi values the local Phi lies between. Intropolate a Q value
        double y1 = QMap[idx-1];
        double y2 = QMap[idx];
        double x1 = PhiMap[idx-1]; //cPos(idx-1,DeltaR)*curState[Apot][idx-1][j];
        double x2 = PhiMap[idx];   //cPos(idx,DeltaR)*curState[Apot][idx][j];
        double m = (y2-y1)/(x2-x1);
        curState.State[Q][i][j] = y1 + m*(localPhi-x1);

        if(isnan(curState.State[Q][i][j]))
        {
            cout << "Q has a nan! (i,j) = (" << i << "," << j << ")" << endl;
            exit(1);
        }
    }

};
