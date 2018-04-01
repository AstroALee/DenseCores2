#ifndef __SOLVEAMPERE_H__
#define __SOLVEAMPERE_H__

#include "DCores_Global.H"
#include "TheState.H"

void SolveAmpere(Params simP,TheState curState,TheState& newState);
void SolveAmpereS0(Params simP,TheState curState,TheState& newState, double s0);

#endif
