#include "TheState.H"

// Constructor inlcudes allocation
TheState::TheState(int M, int N)
{
	cout << "Creating state" << endl;

	this->M = M;
	this->N = N;

	int i,j,k;
    this->State = new double **[NStates];
	for(i=0; i< NStates; i++)
    {
        this->State[i] = new double *[M];
        for(j=0; j< M; j++) this->State[i][j] = new double[N];
    }

    this->VContour = new double[N];

	// In case it's not initialized to zero
   for(i=0;i<NStates;i++) for(j=0;j<M;j++) for(k=0;k<N;k++) State[i][j][k] = 0.0;
   for(k=0;k<N;k++) VContour[k] = 0.0;



};

void TheState::relaxA(const TheState& relaxfrom)
{
	double relaxfrac = 0.5;

	double M = relaxfrom.M;
	double N = relaxfrom.N;

	int k,j;
	//for(j=0;j<M;j++) for(k=0;k<N;k++) this->State[Vpot][j][k] = relaxfrac*(this->State[Vpot][j][k]) + (1-relaxfrac)*(relaxfrom.State[Vpot][j][k]);
	for(j=0;j<M;j++) for(k=0;k<N;k++) this->State[Apot][j][k] = relaxfrac*(this->State[Apot][j][k]) + (1-relaxfrac)*(relaxfrom.State[Apot][j][k]);

	return;// *this;
};

void TheState::relaxV(const TheState& relaxfrom)
{
	double relaxfrac = 0.5;

	double M = relaxfrom.M;
	double N = relaxfrom.N;

	int k,j;
	for(j=0;j<M;j++) for(k=0;k<N;k++) this->State[Vpot][j][k] = relaxfrac*(this->State[Vpot][j][k]) + (1-relaxfrac)*(relaxfrom.State[Vpot][j][k]);
	//for(j=0;j<M;j++) for(k=0;k<N;k++) this->State[Apot][j][k] = relaxfrac*(this->State[Apot][j][k]) + (1-relaxfrac)*(relaxfrom.State[Apot][j][k]);

	return;// *this;
};

void TheState::average(double alpha, const TheState averagewith)
{
	int M = this->M;
	int N = this->N;

	if(alpha>1) alpha=1;
	if(alpha<0) alpha=0;

	int i,j,k; // really only needs to do so for Vpot and Apot, but whatever
	for(i=0;i<NStates;i++) for(j=0;j<M;j++) for(k=0;k<N;k++) (this)->State[i][j][k] = alpha*(averagewith.State[i][j][k]) + (1-alpha)*((this)->State[i][j][k]);

};


TheState& TheState::operator=(const TheState& copyfrom)
{
	int M = copyfrom.M;
	int N = copyfrom.N;

	this->M = M;
	this->N = N;

	int i,j,k;
	for(i=0;i<NStates;i++) for(j=0;j<M;j++) for(k=0;k<N;k++) (this)->State[i][j][k] = copyfrom.State[i][j][k];
    for(k=0;k<N;k++) this->VContour[k] = copyfrom.VContour[k];

	return *this;
};

TheState& TheState::operator+=(const TheState& relaxfrom)
{
	double relaxfrac = 0.5;

	double M = relaxfrom.M;
	double N = relaxfrom.N;

	int k,j;
	for(j=0;j<M;j++) for(k=0;k<N;k++) this->State[Vpot][j][k] = relaxfrac*(this->State[Vpot][j][k]) + (1-relaxfrac)*(relaxfrom.State[Vpot][j][k]);
	for(j=0;j<M;j++) for(k=0;k<N;k++) this->State[Apot][j][k] = relaxfrac*(this->State[Apot][j][k]) + (1-relaxfrac)*(relaxfrom.State[Apot][j][k]);

	return *this;

};
