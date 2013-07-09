#include "mScaleInten.h"

mScaleInten::mScaleInten(scalePara p)
{
	numberChannel = p.gtP.nChan;

	downRate = p.downRate;

	fs = p.gtP.sf/downRate;
	
	scale = p.scale;

	lowPass = new filter(0, (p.lP1+p.lP2) / double(p.gtP.sf), (p.lP2-p.lP1) / double(p.gtP.sf), p.gtP.sf, 0);

	for(int step=0; step<scale; step++)
	{
		double wn = (p.timeScale2[step] + p.timeScale1[step]) / double(fs);
		double ts = (p.timeScale2[step] - p.timeScale1[step]) / double(fs);

		scalePass[step] = new filter(0, wn, ts, fs, 0);
	}

	memIndi=0; 
}

mScaleInten::~mScaleInten(void)
{
	delete lowPass;
	
	for(int step=0; step<scale; step++)
		delete scalePass[step];

	if (memIndi>0) deleteAmpScale();
}

void mScaleInten::newAmpScale(void)
{
	if (memIndi>0) deleteAmpScale();

	for(int step=0; step<scale; step++)
		for(int chan=0; chan<numberChannel; chan++)
			ampScale[step*numberChannel+chan] = new double[nSample];

	memIndi=1;
}

void mScaleInten::deleteAmpScale(void)
{
	if (memIndi>0)
	{
		for(int step=0; step<scale; step++)
			for(int chan=0; chan<numberChannel; chan++)
				delete [] ampScale[step*numberChannel+chan];

		memIndi=0;
	}
}

void mScaleInten::smooth(gammaToneFilterBank *AudiPery, int sigLength, int echo)
{
	int n, step;

	double *temp, *temp2;
	temp = new double[sigLength];
	temp2 = new double[sigLength];

	nSample=sigLength/downRate;

	newAmpScale();
	
	if (echo) printf("Smoothing");
	for(int chan=0; chan<numberChannel; chan++)
	{
		if ( (echo) && ((chan%16)==0) ) printf(".");

		for(n=0; n<sigLength; n++)
			temp[n] = (AudiPery->response[chan][n]>0) ? AudiPery->response[chan][n]:0;

		lowPass->filtering(temp, temp2, sigLength);
				
		for(n=0; n<nSample; n++)
			temp[n] = 20*log10(fabs(temp2[n*downRate])+1);
	
		for(step=0; step<scale; step++)
			scalePass[step]->filtering(temp, ampScale[step*numberChannel+chan], nSample);
	}

	delete [] temp;
	delete [] temp2;

	if (echo) printf("\n");
}

void mScaleInten::readInten(int i, int j, int snr, scalePara p, int sigLength, char *drn)
{
	char filename[255];
	FILE *fp;

	nSample=sigLength/downRate;

	newAmpScale();

	for(int step=0; step<p.scale; step++)
	{
		int shift=step*numberChannel;
		sprintf(filename, "%s/n%d/%d.%d.%d.dat", drn, j, i, snr, int(p.timeScale1[step]));
		fp=fopen(filename, "rb");

		for(int n=0; n<nSample; n++)
			for(int chan=0; chan<numberChannel; chan++)
				fread(&ampScale[shift+chan][n], sizeof(double), 1, fp);

		fclose(fp);
	}
}

void mScaleInten::writeInten(int i, int j, int snr, scalePara p, int sigLength, char *drn)
{
	char filename[255];
	FILE *fp;

	nSample=sigLength/downRate;

	for(int step=0; step<p.scale; step++)
	{
		int shift=step*numberChannel;
		sprintf(filename, "%s/n%d/%d.%d.%d.dat", drn, j, i, snr, int(p.timeScale1[step]));
		fp=fopen(filename, "wb");

		for(int n=0; n<nSample; n++)
			for(int chan=0; chan<numberChannel; chan++)
				fwrite(&ampScale[shift+chan][n], sizeof(double), 1, fp);

		fclose(fp);
	}
}
