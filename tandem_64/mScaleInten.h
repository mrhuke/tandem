#ifndef MSCALEINTEN_H
#define MSCALEINTEN_H

#include "filter.h"

struct scalePara
{
	gammaTonePara gtP;

	int downRate;

	int scale;

	double timeScale1[MAX_SCALE];
	double timeScale2[MAX_SCALE];

	double lP1, lP2;
};

class mScaleInten
{
	filter *lowPass;
	filter *scalePass[MAX_SCALE];

	int memIndi;

	void newAmpScale(void);
	
	void deleteAmpScale(void);


public:
	int nSample, numberChannel;

	int fs, downRate;

	int scale;

	double *ampScale[MAX_CHANNEL*MAX_SCALE];
	
	mScaleInten(scalePara x);
	
	~mScaleInten(void);

	void smooth(gammaToneFilterBank *AudiPery, int sigLength, int echo);

	void readInten(int i, int j, int snr, scalePara, int sigLength, char *drn);
	void writeInten(int i, int j, int snr, scalePara, int sigLength, char *drn);
};

#endif
