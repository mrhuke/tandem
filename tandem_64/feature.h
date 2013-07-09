#include "filter.h"

struct featurePara
{
	gammaTonePara gtP;

	double bP1, bP2, bPTs;
	
	double theta_p;
};

struct acfFeature
{
	double acf[MAX_CHANNEL][MAX_D];
	double acfEv[MAX_CHANNEL][MAX_D];
	double cross[MAX_CHANNEL];
	double crossEv[MAX_CHANNEL];
	double zc[MAX_CHANNEL];
	double zcEv[MAX_CHANNEL];
};

struct acfCompute
{
	double sum[MAX_D];
	double sumS[MAX_D];

	double inputXR[MAX_ACFLEN];
	double inputXI[MAX_ACFLEN];
	double inputXWR[MAX_ACFLEN];
	double inputXWI[MAX_ACFLEN];
};

class feature
{
	int memIndiEv[MAX_CHANNEL];

	int memIndiACF;

	int acfLen;
	int acfOrder;

	filter *bandPass;

	double *envelope[MAX_CHANNEL];

	void fftACF(const double* const input, int frame);
	
	void computeCross(int n);
	
	
protected:
	int window;
	int min_delay;
	int max_delay;
	
	void newFeature(void);
	void deleteFeature(void);


public:
	int sigLength, fs;

	int numFrame, numberChannel;
	
	acfFeature *corrLgm;

	acfCompute data;

	feature(featurePara x);

	~feature(){ deleteFeature();};

	void computeFeature(const gammaToneFilterBank* const AudiPery, int sigLength, int echo);
};
