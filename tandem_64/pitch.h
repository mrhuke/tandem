#include "feature.h"
#define SIDE_CHAN	4
#define SIDE_FRAME	2
#define SHU_NUM		5
#define MHU_NUM		2
#define PHU_NUM		5
#define CANDIDATE_THRESHOLD	0.9
#define THETAC		0.985
#define MAX_CONTOUR	200

struct pitchSUnitNet
{
	float IW[6][SHU_NUM];
	float LW[SHU_NUM];
	float b1[SHU_NUM];
	float b2;
};

struct pitchMUnitNet
{
	float IW[(SIDE_CHAN*2+1)*(SIDE_FRAME*2+1)][MHU_NUM];
	float LW[MHU_NUM];
	float b1[MHU_NUM];
	float b2;
};

struct maskPitchNet
{
	float IW[64*7][PHU_NUM];
	float LW[PHU_NUM];
	float b1[PHU_NUM];
	float b2;
};

struct mask
{
	double value[MAX_CHANNEL];
};

struct pitchContour
{
	int sFrame, eFrame;
	int *value;
	int *oldValue;
	int *indicate;

	mask *mProb;
};

struct pitchPara
{
	double sProb[MAX_CHANNEL][MAX_D];
	double pProb[MAX_D];
	int hNum[MAX_CHANNEL][MAX_D];
	int hNumEv[MAX_CHANNEL][MAX_D];
};

class pitchMask:public feature
{
	int memIndi;

	pitchSUnitNet sNet[MAX_CHANNEL];

	pitchMUnitNet mNet[MAX_CHANNEL];

	maskPitchNet pNet;

	void readNet(char *fname1, char *fname2, char *fname3);

	void singleUnitProb(int frame, int chan, int delay);

	double compareTwoCandidateMAP(int frame, double *mask, int p1, int p2);


protected:
	double theta_p;

	void removeContour(int nCon);
	
	void newContour(int nCon);

	void deleteContour(int nCon);

	double multiUnitProb(int frame, int chan, int nCon);


public:
	pitchContour Pitch[MAX_CONTOUR];

	int numContour;

	pitchPara *Prob;

	pitchMask(featurePara x);

	~pitchMask();

	void newPitchPara(void);

	void deletePitchPara(void);

	void pitchProb(int frame, int delay);

	void readPitch(char *fname);

	void removePitch(void);

	int maskToPitchACF(int frame, double *mask, int d1, int d2);

	int maskToPitchML(int frame, double *mask, int d1, int d2);

	int maskToPitchML2(int frame, double *mask, int d1, int d2);

	int maskToPitchMAP(int frame, double *mask, int d1, int d2);
};
