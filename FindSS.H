
#ifndef __FINDSS_H__
#define __FINDSS_H__

#include "DCores_Global.H"
#include "TheState.H"
#include "UpdateQ.H"
#include "UpdateDQDPHI.H"
#include "SolvePoisson.H"
#include "SolvePoissonSOR.H"
#include "SolveAmpereSOR.H"
#include "SolvePertPoisson.H"
#include "SolveAmpere.H"
#include "PrintState.H"
#include "Integrals.H"
#include "EvalV.H"

void redrawVbdyDeltaPhi(Params simP, TheState& curState);
double Finds0(double firsts0, Params simP, TheState curState);
void FindSSGlobalRescaling(Params& simP, TheState& curState);
void FindSSVfixedTopBot(Params& simP, TheState& curState);
void FindSSVQloopAlpha(Params& simP, TheState& curState);
void FindTheSteadyStateWithAlpha(Params& simP, TheState& curState);
void FindTheSteadyStateWithSOR(Params simP, TheState& curState);
void FindTheSteadyState(Params simP, TheState& curState);
void FindTheSteadyStateWithPert(Params simP, TheState& curState);
void SingleFullUpdate(Params simP,TheState curState,TheState& newState);
void SingleSORUpdate(Params simP,TheState curState,TheState& newState);
void ConvergeOnV(Params& simP,TheState curState,TheState& newState);
void PertFullUpdate(Params simP,TheState& curState, TheState prevState, TheState& newState);
void setVbdy(Params simP,TheState& curState);
void setVbdyAlpha(Params simP,TheState& curState);
int convergeTest(Params simP, TheState prevState, TheState curState);
int DeltaVcheck(Params simP,TheState newState);
void Blend(Params simP,TheState curState,TheState& newState);
void resetV(Params simP, TheState& curState);


#endif
