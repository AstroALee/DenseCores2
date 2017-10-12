
#ifndef __THESTATE_H__
#define __THESTATE_H__

#include "DCores_Global.H"

class TheState
{
public:
	TheState(int M, int N);

	int M,N; // local values
	double ***State;

	double *VContour; // radial location of filament contour

	TheState& operator=(const TheState& copyfrom);
	TheState& operator+=(const TheState& relaxfrom);

protected:
private:

};

#endif // __THESTATE_H__
