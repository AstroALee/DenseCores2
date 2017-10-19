
#ifndef __FINDSS_H__
#define __FINDSS_H__

#include "DCores_Global.H"
#include "TheState.H"
#include "UpdateQ.H"
#include "UpdateDQDPHI.H"
#include "SolvePoisson.H"
#include "SolvePertPoisson.H"
#include "SolveAmpere.H"


void FindTheSteadyState(Params simP, TheState& curState);
void FindTheSteadyStateWithPert(Params simP, TheState& curState);
void SingleFullUpdate(Params simP,TheState curState,TheState& newState);
void PertFullUpdate(Params simP,TheState& curState, TheState prevState, TheState& newState);
void setVbdy(Params simP,TheState& curState);
int convergeTest(Params simP, TheState prevState, TheState curState);
int DeltaVcheck(Params simP,TheState newState);
void Blend(Params simP,TheState curState,TheState& newState);


#endif
