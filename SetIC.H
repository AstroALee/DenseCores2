
#ifndef __SETIC_H__
#define __SETIC_H__

#include "DCores_Global.H"
#include "TheState.H"
#include "MagCylinder.H"
#include "UpdateQ.H"
#include "UpdateDQDPHI.H"
#include "Integrals.H"
#include "EvalV.H"

void getVbdyDeltaPhi(Params simP, TheState& curState);
int SetICnoPoints(Params& simP, TheState& curState);
int SetIC(Params& simP, TheState& curState);
void unitConvert(Params& simP);
void initVpoints(Params& simP, TheState& curState, double undo, double smooth);
void getVbdy(Params simP, TheState& curState);
void setTheVbdy(Params simP, TheState& curState);
void stretchV(Params simP, TheState& curState);
void smoothV(Params simP, TheState& curState,double smooth);


#endif // __SETIC_H__
