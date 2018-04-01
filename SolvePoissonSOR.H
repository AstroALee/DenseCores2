#ifndef __SOLVEPOISSONSOR_H__
#define __SOLVEPOISSONSOR_H__


#include "DCores_Global.H"
#include "TheState.H"
#include "UpdateQ.H"
#include "UpdateDQDPHI.H"

void SolvePoissonSOR(Params simP,TheState curState,TheState& newState);


#endif
