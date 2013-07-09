#include "gammaTone.h"

gammaToneFilter::gammaToneFilter(double centerF, double bC, int sf)
{
	cf = centerF;
	bw = 24.7 * ( cf*0.00437 + 1.0) * bC;
	
	double d = 1.5/PI/bw * double(sf);
	if (fmod(d, 1)>=0.5) delay=floor(d)+1;
	else delay=floor(d);
	
	twoPiT = 2*PI/double(sf);
	
	gain = pow(twoPiT*bw, 4) / 3.0;
	
	z = exp(-twoPiT * bw);

	f1 = cos(cf * twoPiT) * z;
	f2 = sin(cf * twoPiT) * z;
}

void gammaToneFilter::initFilter(void)
{
	for (int n=0; n<4; n++)
	{
		rPart[n] = 0;
		iPart[n] = 0;
	}
}

void gammaToneFilter::oneStep(double input)
{
	double x[4], y[4];
	
	for (int i=0; i<4; i++)
	{
		x[i] = f1*rPart[i] - f2*iPart[i];
		y[i] = f2*rPart[i] + f1*iPart[i];
	}

	rPart[0] = input * f1 + x[0];
	iPart[0] = input * f2 + y[0];
		
	rPart[1] = rPart[0] + x[1];
	iPart[1] = iPart[0] + y[1];
		
	rPart[2] = rPart[1] + x[1] + x[2];
	iPart[2] = iPart[1] + y[1] + y[2];
		
	rPart[3] = rPart[2] + x[1] + 2*x[2] + x[3];
	iPart[3] = iPart[2] + y[1] + 2*y[2] + y[3];

	dR = 3 * rPart[2] - twoPiT*bw * rPart[3] - twoPiT*cf* iPart[3];
	dI = 3 * iPart[2] - twoPiT*bw * iPart[3] + twoPiT*cf* rPart[3];
}

void gammaToneFilter::filtering(double *input, int sigLength, double *response)
{
	int n, k;

	initFilter();

	for (n=0; n<sigLength; n++)
	{
		k = n-delay;
		if (k>=0) response[k] = gain * rPart[3];

		oneStep(input[n]);
	}	

	for (n=0; n<delay; n++)
	{
		k = n-delay+sigLength;
		if (k>=0) response[k] = gain * rPart[3];

		oneStep(0);
	}
}

gammaToneFilterBank::gammaToneFilterBank(gammaTonePara x)
{
	numberChannel = x.nChan;

	lowerCF = x.lCF;
	upperCF = x.rCF;

	sf = x.sf;

	double lowerERB = HzToERBRate(x.lCF);
	double upperERB = HzToERBRate(x.rCF);

	double spaceERB = (x.nChan > 1) ? (upperERB-lowerERB)/(x.nChan-1) : 0;
	
	for (int chan=0; chan<x.nChan; chan++)
	{
		double cf = ERBRateToHz(lowerERB + chan*spaceERB);
		gf[chan] = new gammaToneFilter(cf, BW_CORRECTION, x.sf);
		
		memIndi[chan]=0;
	}
}

gammaToneFilterBank::~gammaToneFilterBank()
{
	for(int chan=0; chan<numberChannel; chan++)
	{
		if(memIndi[chan]==1)
		{
			delete [] response[chan];
			memIndi[chan]=0;
		}

		delete gf[chan];
	}
}

void gammaToneFilterBank::filtering(double *input, int sigLength, int chan)
{
	if (memIndi[chan]==1) delete [] response[chan];
	response[chan] = new double[sigLength];
	memIndi[chan]=1;

	gf[chan]->filtering(input, sigLength, response[chan]);
}
