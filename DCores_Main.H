
#ifndef __DCORES_MAIN_H__
#define __DCORES_MAIN_H__


#include "INIReader.h"
#include "DCores_Global.H"
#include "ErrorMessages.H"
#include "TheState.H"
#include "SetIC.H"
#include "Integrals.H"
#include "FindSS.H"

// How long did all of this take?
void ExecutionTime(clock_t start, clock_t end );

// Calculate derived units and arrays
int Derived(Params& SimP);

// Check that your arguments make sense
int CheckArgs(Params SimP);

// Read from the ini file
int ReadINI(Params& SimParams, int argc, char** filename);

void PrintState(Params simP, TheState curState, string divide, string num);

#endif // __DCORES_MAIN_H__
