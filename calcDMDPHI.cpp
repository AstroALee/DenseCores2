
#include "calcDMDPHI.H"

// given an array for dm/dphi, calculate Q
void QfromDMDPHI(Params simP, TheState& curState)
{

  // using top values, grab phi array
  double PhiGrid[simP.M];
  for(int i=0;i<simP.M;i++) PhiGrid[i] = cPos(i,simP.dR)*curState.State[Apot][i][simP.N-1];
  double QGrid[simP.M];
  for(int i=0;i<simP.M;i++) QGrid[i] = 0.0;

  // for each phi value, we need to evaluate the integral int(r exp(-V)*dr/dphi dz,0,zL)
  for(int curPhiIDX=0;curPhiIDX<simP.M;curPhiIDX++)
  {

    double sum = 0 ;

    // if outside boundary, who cares
    if( (cPos(curPhiIDX,simP.dR)>=curState.VContour[simP.N-1]) )
    {
      QGrid[curPhiIDX] = 0;
      continue;
    }

    // else evaluate the z-integral.
    for(int j=0;j<simP.N;j++)
    {
      // find interpolation of r and V at given phi value
      double curR, curV;
      int curRidx = 0;

      while(true)
      {
        curRidx++;
        if( cPos(curRidx,simP.dR)*curState.State[Apot][curRidx][j] >= PhiGrid[curPhiIDX] )
        {
          break;
        }
      }
      // 0 if at left cell , 1 at right cell
      double intfac = (cPos(curRidx,simP.dR)*curState.State[Apot][curRidx][j] - PhiGrid[curPhiIDX]) / (PhiGrid[curPhiIDX+1] - PhiGrid[curPhiIDX]);

      curR = (1.0-intfac)*cPos(curRidx,simP.dR) + (intfac)*cPos(curRidx+1,simP.dR);
      curV = (1.0-intfac)*curState.State[Vpot][curRidx][j] + (intfac)*curState.State[Vpot][curRidx+1][j];


      // to get dr/dphi, do center difference about the current phi value
      double rR,rL,pR,pL;
      pR = PhiGrid[curPhiIDX+1] ;
      pL = PhiGrid[curPhiIDX-1] ;

      curRidx = 0;
      while(true)
      {
        curRidx++;
        if( cPos(curRidx,simP.dR)*curState.State[Apot][curRidx][j] >= pR )
        {
          break;
        }
      }
      intfac = (cPos(curRidx,simP.dR)*curState.State[Apot][curRidx][j] - pR) / (PhiGrid[curPhiIDX+2] - pR);

      rR = (1.0-intfac)*cPos(curRidx,simP.dR) + (intfac)*cPos(curRidx+1,simP.dR);


      curRidx = 0;
      while(true)
      {
        curRidx++;
        if( cPos(curRidx,simP.dR)*curState.State[Apot][curRidx][j] >= pL )
        {
          break;
        }
      }
      intfac = (cPos(curRidx,simP.dR)*curState.State[Apot][curRidx][j] - pL) / (PhiGrid[curPhiIDX] - pL);

      rL = (1.0-intfac)*cPos(curRidx,simP.dR) + (intfac)*cPos(curRidx+1,simP.dR);

      double drdphi = (rR-rL)/(pR-pL);

      sum = sum + curR*exp(-curV)*drdphi*simP.dZ;

    }

    QGrid[curPhiIDX] = (1.0/4.0/3.14159)*curState.dMdPhi[curPhiIDX]/sum;




  }

  int Vidx = 0;
  while(true)
    if(cPos(Vidx,simP.dR) >= curState.VContour[simP.N-1]) break;
    else Vidx++;

  double slope = (QGrid[Vidx-1]-QGrid[Vidx-2])/simP.dR;
  double Qedge = QGrid[Vidx-1] + slope*( curState.VContour[simP.N-1] - cPos(Vidx-1,simP.dR) );


  // With grid created, now map to all the Phi values on the grid
  for(int i=0;i<simP.M;i++) for(int j=0;j<simP.N;j++)
  {
    if( cPos(i,simP.dR) >= curState.VContour[j] )
      curState.State[Q][i][j] = Qedge;
    else
    {
      double curPhi = cPos(i,simP.dR)*curState.State[Apot][i][j];

      // Find indx
      int Pidx = 0;
      while(true)
        if(curPhi >= PhiGrid[Pidx]) break;
        else Pidx++;

      if(Pidx==0)
      {
        curState.State[Q][i][j] = QGrid[Pidx];
      }
      else
      {
        double x = (curPhi - PhiGrid[Pidx-1])/(PhiGrid[Pidx] - PhiGrid[Pidx-1]); // = 1 if at Pidx, = 0 if at Pidx-1
        curState.State[Q][i][j] = x*QGrid[Pidx] + (1-x)*QGrid[Pidx-1];
      }
    }
  }


}

// using the current state's V,A,Q, perform the dmdphi integral and store in curState class
void calcDMDPHI(Params simP, TheState& curState)
{

  double PhiGrid[simP.M];
  double QGrid[simP.M];

  // using top values of phi, create grid of interest
  for(int i=0;i<simP.M;i++) PhiGrid[i] = cPos(i,simP.dR)*curState.State[Apot][i][simP.N-1];
  for(int i=0;i<simP.M;i++) QGrid[i] = curState.State[Q][i][simP.N-1];

  for(int curPhiIDX=0;curPhiIDX<simP.M;curPhiIDX++)
  {

    double sum = 0 ;

    // if outside boundary, who cares
    if( (curPhiIDX==0) or (cPos(curPhiIDX,simP.dR)>=curState.VContour[simP.N-1]) )
    {
      curState.dMdPhi[curPhiIDX] = 0;
      continue;
    }

    // else evaluate the z-integral.
    for(int j=0;j<simP.N;j++)
    {
      // find interpolation of r and V at given phi value
      double curR, curV;
      int curRidx = 0;

      while(true)
      {
        curRidx++;
        if( cPos(curRidx,simP.dR)*curState.State[Apot][curRidx][j] >= PhiGrid[curPhiIDX] )
        {
          break;
        }
      }
      // 0 if at left cell , 1 at right cell
      double intfac = (cPos(curRidx,simP.dR)*curState.State[Apot][curRidx][j] - PhiGrid[curPhiIDX]) / (PhiGrid[curPhiIDX+1] - PhiGrid[curPhiIDX]);

      curR = (1.0-intfac)*cPos(curRidx,simP.dR) + (intfac)*cPos(curRidx+1,simP.dR);
      curV = (1.0-intfac)*curState.State[Vpot][curRidx][j] + (intfac)*curState.State[Vpot][curRidx+1][j];


      // to get dr/dphi, do center difference about the current phi value
      double rR,rL,pR,pL;
      pR = PhiGrid[curPhiIDX+1] ;
      pL = PhiGrid[curPhiIDX-1] ;

      curRidx = 0;
      while(true)
      {
        curRidx++;
        if( cPos(curRidx,simP.dR)*curState.State[Apot][curRidx][j] >= pR )
        {
          break;
        }
      }
      intfac = (cPos(curRidx,simP.dR)*curState.State[Apot][curRidx][j] - pR) / (PhiGrid[curPhiIDX+2] - pR);

      rR = (1.0-intfac)*cPos(curRidx,simP.dR) + (intfac)*cPos(curRidx+1,simP.dR);


      curRidx = 0;
      while(true)
      {
        curRidx++;
        if( cPos(curRidx,simP.dR)*curState.State[Apot][curRidx][j] >= pL )
        {
          break;
        }
      }
      intfac = (cPos(curRidx,simP.dR)*curState.State[Apot][curRidx][j] - pL) / (PhiGrid[curPhiIDX] - pL);

      rL = (1.0-intfac)*cPos(curRidx,simP.dR) + (intfac)*cPos(curRidx+1,simP.dR);

      double drdphi = (rR-rL)/(pR-pL);

      sum = sum + curR*exp(-curV)*drdphi*simP.dZ;

    }

    curState.dMdPhi[curPhiIDX] = 4.0*3.14159*QGrid[curPhiIDX]*sum;
  }

}
