#ifndef __SOLVEPOISSON_H__
#define __SOLVEPOISSON_H__


#include "DCores_Global.H"
#include "TheState.H"

void SolvePoissonS0(Params simP,TheState curState,TheState& newState, double s0);
void SolvePoisson(Params simP,TheState curState,TheState& newState);


#endif
