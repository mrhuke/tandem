#ifndef FILTER_H
#define FILTER_H

#include "gammaTone.h"

class filter
{
	int fs;
	int fLength, nOrder, nFFT;

	double beta;
	
	double *rValue;
	double *iValue;

	double *inputR;
	double *inputI;

	void hammingPara(double delta, double transBw, int &fLength);

	void hammingBandPass(double wc, double wn);
	void hammingDiffOneBandPass(double wc, double wn);
	void hammingDiffTwoBandPass(double wc, double wn);

	void kaiserPara(double delta, double transBw);
	
	void kaiserLowPass(double wn);
	void diffKaiserLowPass(double wn);

public:
	filter(double ws, double wn, double ts, int fs, int type);

	~filter();

	void filtering(double *input, double *output, int sigLength);
};

#endif
