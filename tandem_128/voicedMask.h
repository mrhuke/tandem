#include "pitch.h"

struct groupPara
{
	int stream[MAX_CHANNEL];
	int label[MAX_CHANNEL];
	int mark[MAX_CHANNEL];
};

class voicedMask : public pitchMask
{
	int memIndi;
	
	int min_pitch, max_pitch;
	double thd_cross, thd_evCross;  //thresholds for corss-channel correlations used in initial estimation

	int checkPitchCon(double p1, double p2);

	int checkMaskCon(double *mask1, double *mask2, int frame);
	int checkMaskCon2(double *mask1, double *mask2, int frame);
	
	int isOverlap(int m, int n);

	int isConnected(int m, int n);

	void mergeContour(void);

	
	void reDetermineMask(int n);

	void expandMask(void);

	
	void removeDuplicate(int *p1, int *p2);

	void switchCandidate(int *p1, mask *m1, int *p2, mask *m2);

	void findContour(int *p, mask *m);
	
	void convCont(int *p1, mask *m1, int *p2, mask *m2);

	
	void reEstimatePitch(int nCon);

	void checkContour(int nCon);

	void developeContour(int nCon);


public:
//	groupPara *grp;

	voicedMask(featurePara x):pitchMask(x) { 
        	memIndi=0; 
		min_pitch=x.gtP.sf/400; 
		max_pitch=x.gtP.sf/70;
		if (numberChannel==128){
			thd_cross=.985;
			thd_evCross=.985;
		}else{
			thd_cross=.945;
                        thd_evCross=.945;
		}
	};

	~voicedMask(void){;};

	void initPitchEst(void);

	void iterativePitchEst(void);

	void maskToPitch(int nCon);
	
//	void newGroupPara(void);

//	void deleteGroupPara(void);

	void dtmPitchMask(void);

	void readMask(char *filename);
};
