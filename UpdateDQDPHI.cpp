
#include "UpdateDQDPHI.H"

void updateDQDPHI(Params simP, TheState& curState)
{
    //                      (dQ/dr)*DeltaR + (dQ/dz)*DeltaZ       Q(r+1)-Q(r-1)+Q(z+1)-Q(z-1)
    // Estimates dQ/dPhi =  ------------------------------- = ----------------------------------
    //                      (dP/dr)*DeltaR + (dP/dz)*DeltaZ   Phi(r+1)-Phi(r-1)+Phi(z+1)-Phi(z-1)
    //
    // Note: The above is true only if you perform a centered-difference for each derivative
    // (else factors of 1/2 and such would not cancel out entirely)

    // Boundary conditions:
    // top and bottom have reflective
    // left has reflective
    // right will utilize a three stencil estimate for the derivative



    int M = simP.M;
    int N = simP.N;
    double DeltaR = simP.dR;
    double DeltaZ = simP.dZ;



    // find Q at filament boundary
    int i=0; while(true){ i++; if(cPos(i,DeltaR)>=curState.VContour[N-1]) break;}
    double Qright=curState.State[Q][i][N-1];
    double Qleft=curState.State[Q][i-1][N-1];
    double m = (Qright-Qleft)/DeltaR;
    double Qedge = Qleft + m*( curState.VContour[N-1] - cPos(i-1,DeltaR) ); // 8/23/17 - originally was cPos(i,DeltaR) .. seems incorrect?

    double Rbdy = simP.RhoTop[i-1] + (simP.RhoTop[i]-simP.RhoTop[i-1])*( curState.VContour[N-1] - cPos(i-1,DeltaR) )/DeltaR;

    int s;
    for(s=0;s<M*N;s++)
    {
        Qright = 0; Qleft = 0;
        double Qup=0, Qdown=0;
        double Phiright=0, Phileft=0, Phiup=0, Phidown=0;
        double dQdr=0,dQdz=0,dPdr=0,dPdz=0;




        if(Ridx(s,M)==0) // Left side
        {
            // We use one-sided derivatives to estimate these quantities.
            // We are dividing by DeltaR so that it cancels out the DeltaR at the end
            double Q0 = curState.State[Q][Ridx(s,M)][Zidx(s,M)];
            double Q1 = curState.State[Q][Ridx(s,M)+1][Zidx(s,M)];
            double Q2 = curState.State[Q][Ridx(s,M)+2][Zidx(s,M)];
            double Q3 = curState.State[Q][Ridx(s,M)+3][Zidx(s,M)];
            double P0 = cPos(Ridx(s,M),DeltaR)*curState.State[Apot][Ridx(s,M)][Zidx(s,M)]; // = 0
            double P1 = cPos(Ridx(s,M)+1,DeltaR)*curState.State[Apot][Ridx(s,M)+1][Zidx(s,M)];
            double P2 = cPos(Ridx(s,M)+2,DeltaR)*curState.State[Apot][Ridx(s,M)+2][Zidx(s,M)];
            double P3 = cPos(Ridx(s,M)+3,DeltaR)*curState.State[Apot][Ridx(s,M)+3][Zidx(s,M)];

            dQdr=(2*Q0-5*Q1+4*Q2-Q3)/DeltaR;
            dPdr=(2*P0-5*P1+4*P2-P3)/DeltaR;
        }
        else if(Ridx(s,M)==M-1) // Right side
        {
            // Instead of saying what Q(r+1) and Q(r-1) is, we'll
            // specify a value for Q(r+1) and set Q(r-1), where Q(r+1)
            // will estimate (dQ/dr)*DeltaR

            if(cPos(Ridx(s,M),DeltaR) > curState.VContour[Zidx(s,M)]) // on right side of contour
            {
                dQdr = 0;
            }
            else
            {
                double Qr = curState.State[Q][Ridx(s,M)][Zidx(s,M)] ;
                double Qrm1 = curState.State[Q][Ridx(s,M)-1][Zidx(s,M)] ;
                double Qrm2 = curState.State[Q][Ridx(s,M)-2][Zidx(s,M)] ;
                dQdr = ( 3*Qr - 4*Qrm1 + 1*Qrm2 )/(2*DeltaR);
            }

            // Phi is done similarly, but dPhi/dr = d(r*A)/dr = A + r*(dA/dr)
            // Phiright will estimate = (dPhi/dr)*DeltaR = A*DeltaR + r*DeltaR*(dA/dr)
            double Ar = curState.State[Apot][Ridx(s,M)][Zidx(s,M)] ;
            double Arm1 = curState.State[Apot][Ridx(s,M)-1][Zidx(s,M)] ;
            double Arm2 = curState.State[Apot][Ridx(s,M)-2][Zidx(s,M)] ;
            double dAdr = (3*Ar - 4*Arm1 + 1*Arm2 )/(2*DeltaR); // continuous across filament, no special treatment needed

            dPdr = Ar + cPos(Ridx(s,M),DeltaR)*dAdr;

        }
        else // r derivatives -- interior
        {
            Qleft = curState.State[Q][Ridx(s,M)-1][Zidx(s,M)];
            Qright = curState.State[Q][Ridx(s,M)+1][Zidx(s,M)];
            dQdr = (Qright-Qleft)/(2*DeltaR); // this is true if both points are left of the contour

            Phileft  = cPos(Ridx(s,M)-1,DeltaR)*curState.State[Apot][Ridx(s,M)-1][Zidx(s,M)];
            Phiright = cPos(Ridx(s,M)+1,DeltaR)*curState.State[Apot][Ridx(s,M)+1][Zidx(s,M)];
            dPdr = (Phiright-Phileft)/(2*DeltaR);  // this is true if both points are left of the contour

            if(cPos(Ridx(s,M),DeltaR) >= curState.VContour[Zidx(s,M)]) dQdr = 0; // cell is right of contour
            else if( cPos(Ridx(s,M)+1,DeltaR) >= curState.VContour[Zidx(s,M)] )  // cell is near the contour
            {
                // Center difference overlaps the contour. Instead interpolate

                // Need the value of Q at the contour: We need the value of V at the contour
                //double Vr = curState.State[Vpot][1+Ridx(s,M)][Zidx(s,M)];
                //double Vl = curState.State[Vpot][Ridx(s,M)][Zidx(s,M)];
                double x  = (curState.VContour[Zidx(s,M)] - cPos(Ridx(s,M),DeltaR))/(DeltaR);
                //double Vbdy = (1-x)*Vl + x*Vr;
                //double Qbdy = Rbdy * exp(Vbdy);

                double dr = (curState.VContour[Zidx(s,M)] - cPos(Ridx(s,M)-1,DeltaR)) ;
                double dQ = Qedge - Qleft;
                dQdr = dQ/dr;



                // Do the same for Phi
                Phileft = cPos(Ridx(s,M)-1,DeltaR)*curState.State[Apot][Ridx(s,M)-1][Zidx(s,M)];
                double PhiC = cPos(Ridx(s,M),DeltaR)*curState.State[Apot][Ridx(s,M)][Zidx(s,M)];
                double Pbdy = (1-x)*PhiC + x*Phiright;
                dPdr = (Pbdy - Phileft)/(curState.VContour[Zidx(s,M)] - cPos(Ridx(s,M)-1,DeltaR)) ;
            }
        }

        if(Zidx(s,M)==0 || Zidx(s,M)==N-1) // top or bottom
        {
            // Qup - Qdown = 0, similar for Phi.
            dQdz = 0;
            dPdz = 0;
        }
        else // interior
        {
            Qdown = curState.State[Q][Ridx(s,M)][Zidx(s,M)-1];
            Qup   = curState.State[Q][Ridx(s,M)][Zidx(s,M)+1];
            dQdz = (Qup-Qdown)/(2*DeltaZ); // Assumes both points are interior to contour

            Phidown  = cPos(Ridx(s,M),DeltaR)*curState.State[Apot][Ridx(s,M)][Zidx(s,M)-1];
            Phiup    = cPos(Ridx(s,M),DeltaR)*curState.State[Apot][Ridx(s,M)][Zidx(s,M)+1];
            dPdz = (Phiup-Phidown)/(2*DeltaZ); // Assumes both points are interior to contour

            if(cPos(Ridx(s,M),DeltaR) >= curState.VContour[Zidx(s,M)]) dQdz = 0; // point is right of contour
            else if( cPos(Ridx(s,M),DeltaR) >= curState.VContour[Zidx(s,M)+1] ) // point is next to contour
            {
                // Need z location of contour at cPos(Ridx(s,M),DeltaR) radius
                double m = DeltaZ/(curState.VContour[Zidx(s,M)+1]-curState.VContour[Zidx(s,M)]);
                double zCt = cPos(Zidx(s,M),DeltaZ) + m*(cPos(Ridx(s,M),DeltaR)-curState.VContour[Zidx(s,M)]);

                // Now need to get Q at that location; to do so needs V at that location
                double Vu = curState.State[Vpot][Ridx(s,M)][Zidx(s,M)+1];
                double Vd = curState.State[Vpot][Ridx(s,M)][Zidx(s,M)];
                double x = (zCt-cPos(Zidx(s,M),DeltaZ))/DeltaZ;
                double Vbdy = (1-x)*Vd + x*Vu;
                double Qbdy = Rbdy * exp(Vbdy);

                double dQ = Qbdy - Qdown;
                double dz = (zCt - cPos(Zidx(s,M)-1,DeltaZ));// + DeltaZ;
                dQdz = dQ/dz;

                // Do the same for Phi
                Phidown  = cPos(Ridx(s,M),DeltaR)*curState.State[Apot][Ridx(s,M)][Zidx(s,M)];
                double Pbdy = (1-x)*Phidown + x*Phiup;
                dPdz = (Pbdy - Phidown)/(zCt - cPos(Zidx(s,M),DeltaZ));

            }


        }

        double dQ = dQdr*DeltaR + dQdz*DeltaZ; //Qup - Qdown + Qright - Qleft;
        double dP = dPdr*DeltaR + dPdz*DeltaZ; //Phiup - Phidown + Phiright - Phileft;

        double dQdPhi = 0;
        if(dP!=0) dQdPhi = dQ/dP;

        if(isnan(dQdPhi))
        {
            cout << "dQdPhi has a nan! (i,j)=(" << Ridx(s,M) << "," << Zidx(s,M) << ")" << endl;
            cout << "    dQ = " << dQ << " , dP = " << dP << endl;
            cout << "    dQdr = " << dQdr << " , dPdr = " << dPdr << endl;
            cout << "    dQdz = " << dQdz << " , dPdz = " << dPdz << endl;
            cout << "    cPos = " << cPos(Ridx(s,M),DeltaR) << " , VCont+1 = " << curState.VContour[Zidx(s,M)+1] << endl;
            cout << "    cPos = " << cPos(Ridx(s,M),DeltaR) << " , VCont = " << curState.VContour[Zidx(s,M)] << endl;
            if(cPos(Ridx(s,M),DeltaR) >= curState.VContour[Zidx(s,M)+1]) cout << "Near a contour!" << endl;
            exit(1);
        }

        //if(isnan(dQdPhi)) dQdPhi = 0;

        //if(Zidx(s,M)==N-1) cout << dQdPhi << endl;

        if(cPos(Ridx(s,M),DeltaZ) <= curState.VContour[Zidx(s,M)]) curState.State[dQdP][Ridx(s,M)][Zidx(s,M)] = dQdPhi;
        //if(dQdPhi > 0.3) cout << "Large derivative! at " << Ridx(s,M) << " and " << Zidx(s,M) << "! (" << dQdPhi << ")" << endl;
        else curState.State[dQdP][Ridx(s,M)][Zidx(s,M)] = 0; // dQdPhi should have been set to zero anyway, but just to be safe.

    }





};
