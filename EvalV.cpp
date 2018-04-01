
#include "EvalV.H"


double Veval(int i, int jj, Params simP)
{
	double val = simP.VTop[i] ; // MagCylinder.cpp appropriately re-normalizes VTop to have V=0 at filament boundry

	// assumes G=1
	double a = simP.alpha_fac;
	double G = 1.0;
	double m = simP.mExcess;

	// First find Vedge so V=0 at top-right (filament boundary at rCyl,zL)
	double z1 = simP.zL;
	double z2 = simP.zL;

	double Vedge = 0;
	double rCyl = simP.rCyl;
	//cout << "Loop total = " << simP.pointLoopTotNum << endl;
	for(int j=0; j<= simP.pointLoopTotNum; j++)
	{
		double pos1 = sqrt( rCyl*rCyl + z1*z1 );
		double pos2 = sqrt( rCyl*rCyl + z2*z2 );
		if(j==0) pos1 = pos1 + simP.zL/100.0;

		Vedge = Vedge - G*a*m/pos1;
		if(j!=0) Vedge = Vedge - G*a*m/pos2;

		z1 = z1 + 2.0*simP.zL;
		z2 = z2 - 2.0*simP.zL;
	}
	//cout << "Normalization value in EvalV is " << Vedge << endl;

	// Now eval value of V at desired i,j
	z1 = cPos(jj,simP.dZ); z2 = z1;
	double r = cPos(i,simP.dR);
	for(int j=0; j<= simP.pointLoopTotNum; j++)
	{
		double pos1 = sqrt( r*r + z1*z1 );
		double pos2 = sqrt( r*r + z2*z2 );
		if(j==0) pos1 = pos1 + simP.zL/100.0; // doesn't matter what it is

		val = val - G*a*m/pos1;
		if(j!=0) val = val - G*a*m/pos2;

		z1 = z1 + 2.0*simP.zL;
		z2 = z2 - 2.0*simP.zL;

	}
	val = val - Vedge; // renormalize so Vpoints has V=0 at filament boundary



	return val;
}
