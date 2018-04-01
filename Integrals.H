
#ifndef __INTEGRALS_H__
#define __INTEGRALS_H__

#include "DCores_Global.H"
#include "TheState.H"

void CylinderMassCheck(Params simP, TheState curState);
double Trap2D(Params simP, TheState curState, double s0);
double Trap2D(Params simP, TheState curState);
double Trap2D_wrong(Params simP, TheState curState);
double Trap2D_old(Params simP, TheState curState);
double Simpson2D(Params simP, TheState curState);
double MCInt(Params simP, TheState curState);
double MCIntSpline(Params simP, TheState curState);
double DirectIntSpline(Params simP, TheState curState);
double CubicSpline(double rPoint, double zPoint, Params simP, TheState curState);
double getF(int M,int N,Params simP,TheState curState,int which);

#endif
