#ifndef SEGMENT_H
#define SEGMENT_H

#include "mScaleInten.h"

#define MAX_FRONT	2000
#define THETA_DIS	8
#define THETA_CROSS	0.999
#define THETA_CROSS3 0.99
#define THETA_CROSS2 0.95

struct segScale
{
	int scale;

	double freqScale[MAX_SCALE];
	int timeScale[MAX_SCALE];
};

struct chanCandidate
{
	int *pos, *cPos, *cMark;
	int nC;
};

struct front
{
	int *pos[MAX_CHANNEL], *cPos[MAX_CHANNEL], *cMark[MAX_CHANNEL];
	int *sChan, *eChan;
	int nF;

	int *extendFront1;
	int *extendFront2;
};

class segment
{
	int memIndi;

	double *inten[MAX_CHANNEL];

	double *diffInten[MAX_CHANNEL];

	double thL[MAX_CHANNEL], thH[MAX_CHANNEL];

	double theta_cross;

	chanCandidate oS1[MAX_CHANNEL], oS2[MAX_CHANNEL], fS1[MAX_CHANNEL], fS2[MAX_CHANNEL], tmpC1, tmpC2;

	front oF, fF, oF1, oF2, oF3, fF1, fF2, fF3;

	
	void newCandidate(chanCandidate *s);

	void deleteCandidate(chanCandidate *s);


	void newFront(front *f);

	void deleteFront(front *f);

	
	void freqSmooth(double freqScale, double *input[MAX_CHANNEL]);
	
	
	void determineTheta(int chan, double theta1, double theta2);

	void findCandidate(int chan);

	void chanSegment(int chan);


	void matchingRange(chanCandidate ofPos, int n, int label, double &sPos1, double &ePos1, double &sPos2, double &ePos2);

	int matchingCross(chanCandidate *ofPos, int chan, int n, int k, int label);

	void candToFront(chanCandidate *ofPos, front *Front, int label);

	void candToFront2(chanCandidate *ofPos, front *Front, int label);

	        
	void oFMatch(front *oFront, front *oFront2, front *fFront2);
	void oFMatch2(front *oFront, front *oFront2);

	
	void frontToCandidate(front *oFront, front *fFront);
	
	void adjustPos(front *oFront, front *fFront);


	int newPos(int pos, int cPos, chanCandidate ofPos);

	void accuratePos(front *Front, chanCandidate *ofPos, chanCandidate *oPos2, chanCandidate *fPos2);

	void accuratePos2(front *Front, chanCandidate *ofPos, chanCandidate *oPos2, chanCandidate *fPos2);

	void removeOFSet(front *oFront, chanCandidate *oPos, int label);

	void removeOFSet2(front *oFront, chanCandidate *oPos, int label);


	void enlongFront(front *oFront, front *oFront2, chanCandidate *oPos2);

	void mergeFront(front *oFront, int label);

	void assignFront(front *tF, front *sF, int n);

public:

	int nSample, factor;
	int numberChannel;

	int *seg[MAX_CHANNEL];

	double *cross[MAX_CHANNEL], *crossEv[MAX_CHANNEL];

	segment(int nChan){ memIndi=0; factor=4; numberChannel = nChan; };

	~segment(void) { if(memIndi==1) deleteSegPara(); memIndi=0; };

	void newSegPara(void);

	void deleteSegPara(void);

	void ofFront(mScaleInten *scaleInten, segScale s);

	void frontToSeg(void);

	void backGroundSet(front *oFront, front *fFront);
};

#endif
