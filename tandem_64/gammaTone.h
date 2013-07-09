#ifndef GAMMATONE_H
#define GAMMATONE_H 

#include "tool.h"
#include "common.h"

#define BW_CORRECTION       1.019      			/* ERB bandwidth correction 4th order */

struct gammaTonePara
{
	double lCF, rCF;
	int nChan;
	int sf;
};

class gammaToneFilter
{		
	double rPart[4], iPart[4];
		
	double dR, dI, ddR, ddI;

	double gain, twoPiT, z, f1, f2;

	int delay;

public:
	double cf, bw;

	gammaToneFilter(double cf, double bc, int sf);

	~gammaToneFilter(){;}

	void initFilter(void);
	
	void oneStep(double input);

	void filtering(double *input, int sigLength, double *response);
	
	//void filtering(double *input, int sigLength, double *response, double *insf);
};

class gammaToneFilterBank
{	
	int memIndi[MAX_CHANNEL];

	double lowerCF, upperCF;

	gammaToneFilter *gf[MAX_CHANNEL];
	
	double HzToERBRate(double Hz){ return( 21.4*log10(Hz*0.00437 + 1.0) ); }
	
	double ERBRateToHz(double ERBRate){ return( (pow(10, ERBRate/21.4) - 1) / 0.00437 ); }

public:
	int numberChannel, sigLength, sf;
	
	double *response[MAX_CHANNEL];
	
	gammaToneFilterBank(gammaTonePara x);

	~gammaToneFilterBank();

	void filtering(double *input, int sigLength, int chan);
};

#endif
