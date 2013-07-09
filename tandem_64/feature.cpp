#include "feature.h"
#include <algorithm>
#include <fstream>
#include <iterator>
#include <iostream>

void feature::computeFeature(const gammaToneFilterBank* const AudiPery, int sLen, int echo)
{
	int n, chan;

	double *temp;
	temp = new double[sLen];

	sigLength=sLen;

	deleteFeature();
	newFeature();
	
	if (echo) printf("Computing Envelope");

	for(chan=0; chan<numberChannel; chan++)
	{
		if ((echo) && ((chan%16)==0)) printf(".");

		for(n=0; n<sigLength; n++)
			temp[n] = (AudiPery->response[chan][n]>0) ? AudiPery->response[chan][n]:0;
		
		bandPass->filtering(temp, envelope[chan], sigLength);
	}
	if(echo) printf("\n");

	if (echo) printf("Auto-correlation");
	for(chan=0; chan<numberChannel; chan++)
	{
		if ( (echo) && ((chan%16)==0) ) printf(".");

		for(n=0; n<numFrame; n++)
		{
			fftACF(envelope[chan], n);
			
			for(int m=0; m<max_delay; m++)
				corrLgm[n].acfEv[chan][m]=data.inputXR[m];

			corrLgm[n].zcEv[chan]=zeroCross(data.inputXR, max_delay);

			fftACF(AudiPery->response[chan], n);
			for(int m=0; m<max_delay; m++)
				corrLgm[n].acf[chan][m]=data.inputXR[m];

			corrLgm[n].zc[chan]=zeroCross(data.inputXR, max_delay);
		}
	}
	
	for(n=0; n<numFrame; n++)
		computeCross(n);

	#ifdef DEBUG	
	FILE *out=fopen("cross","w");
	for(int m=0; m<numFrame; m++)
	{
		for(int c=0; c<numberChannel; c++)
			fprintf(out, "%f ", corrLgm[m].cross[c]);			
		fprintf(out, "\n");
	}
	fclose(out);

	out=fopen("evCross","w");
        for(int m=0; m<numFrame; m++)
        {
                for(int c=0; c<numberChannel; c++)
                        fprintf(out, "%f ", corrLgm[m].crossEv[c]);
                fprintf(out, "\n");
        }
        fclose(out);
	#endif

	for(chan=0; chan<numberChannel; chan++)
	{
		delete [] envelope[chan];
		memIndiEv[chan]=0;
	}

	if (echo) printf("\n");
}

void feature::newFeature(void)
{
	int chan;

	numFrame = sigLength*100/(fs);

	window = fs/50;

	acfLen = max_delay + window;
	acfOrder = ceil(log(double(acfLen)) / log(2.0));
	acfLen = pow(2.0, acfOrder);
	
	for(chan=0; chan<numberChannel; chan++)
	{
		envelope[chan]=new double[sigLength];	
		memIndiEv[chan]=1;
	}
	
	corrLgm = new acfFeature[numFrame];
	memIndiACF=1;
}

void feature::deleteFeature()
{
	for(int chan=0; chan<numberChannel; chan++)
	{
		if (memIndiEv[chan]>0) delete [] envelope[chan];
		memIndiEv[chan]=0;
	}

	if (memIndiACF>0)
	{
		delete [] corrLgm;
		memIndiACF=0;
	}
}

void feature::fftACF(const double* const response, int frame)
{
	int step, delay;

	for(step=0; step<acfLen; step++)
		data.inputXR[step] = data.inputXI[step] = data.inputXWR[step] = data.inputXWI[step] = 0;

	for(step=0; step<(window+max_delay); step++)
	{
		int tim = (frame+2) * window/2 - (step+1);
		if ( (tim>= 0)  && (tim<sigLength) )
		{
			data.inputXR[step] = response[tim];

			if(step<window) data.inputXWR[step] = data.inputXR[step];
		}
	}

	data.sum[0]=data.sumS[0]=0;
	for(step=0; step<window; step++)
	{
	//	data.sum[0] += data.inputXR[step]/double(window);
		data.sumS[0] += data.inputXR[step]*data.inputXR[step]/double(window);
	}

	for(step=1; step<max_delay; step++)
	{
		double f1 = data.inputXR[step+window-1] - data.inputXR[step-1];
		double f2 = data.inputXR[step+window-1] + data.inputXR[step-1];
		
	//	data.sum[step] = data.sum[step-1] + f1/double(window);
		data.sumS[step] = data.sumS[step-1] + f1*f2/double(window);

		if (data.sumS[step]<0) data.sumS[step]=0;
	}

	fft(data.inputXR, data.inputXI, acfOrder, 1);
	fft(data.inputXWR, data.inputXWI, acfOrder, 1);

	for(step=0; step<acfLen; step++)
	{
		double f1 = data.inputXR[step]*data.inputXWR[step] + data.inputXI[step]*data.inputXWI[step];
		double f2 = -data.inputXR[step]*data.inputXWI[step] + data.inputXI[step]*data.inputXWR[step];
	
		data.inputXR[step]=f1;
		data.inputXI[step]=f2;
	}

	fft(data.inputXR, data.inputXI, acfOrder, -1);

	for(delay=0; delay<max_delay; delay++)
	{
		data.inputXR[delay] /= double(acfLen*window);
		data.inputXR[delay] /= double(sqrt(data.sumS[0]*data.sumS[delay]+1e-5*data.sumS[0])) / data.sumS[0]; //+1e-100));
	}
}

void feature::computeCross(int n)
{
	for(int chan=0; chan<(numberChannel-1); chan++)
	{
		double m1=0, m2=0, sum=0, sum1=0, sum2=0;

		for(int delay=0; delay<max_delay; delay++)
		{
			m1 += corrLgm[n].acf[chan][delay];
			m2 += corrLgm[n].acf[chan+1][delay];

			sum += corrLgm[n].acf[chan][delay] * corrLgm[n].acf[chan+1][delay];
			sum1 += corrLgm[n].acf[chan][delay] * corrLgm[n].acf[chan][delay];
			sum2 += corrLgm[n].acf[chan+1][delay] * corrLgm[n].acf[chan+1][delay];
		}

		m1 /= double(max_delay);
		m2 /= double(max_delay);
		corrLgm[n].cross[chan] = (sum-m1*m2) /sqrt( (sum1-m1*m1) * (sum2-m2*m2) + 1e-100 );

		m1=0, m2=0, sum=0, sum1=0, sum2=0;

		for(int delay=0; delay<max_delay; delay++)
		{
			m1 += corrLgm[n].acfEv[chan][delay];
			m2 += corrLgm[n].acfEv[chan+1][delay];

			sum += corrLgm[n].acfEv[chan][delay] * corrLgm[n].acfEv[chan+1][delay];
			sum1 += corrLgm[n].acfEv[chan][delay] * corrLgm[n].acfEv[chan][delay];
			sum2 += corrLgm[n].acfEv[chan+1][delay] * corrLgm[n].acfEv[chan+1][delay];
		}

		m1 /= double(max_delay);
		m2 /= double(max_delay);
		corrLgm[n].crossEv[chan] = (sum-m1*m2) /sqrt( (sum1-m1*m1) * (sum2-m2*m2) + 1e-100 );
	}

	corrLgm[n].cross[numberChannel-1]=0;
	corrLgm[n].crossEv[numberChannel-1]=0;
}

feature::feature(featurePara x)
{	
	fs = x.gtP.sf;
	numberChannel = x.gtP.nChan;

	min_delay = fs/500;
        max_delay = fs*3/200;

	bandPass = new filter( (x.bP1 + x.bP2) / double(fs), (x.bP2 - x.bP1) / double(fs), double(x.bPTs*2) / double(fs), fs, 0);
	
	for(int chan=0; chan<numberChannel; chan++)
		memIndiEv[chan]=0;

	memIndiACF=0;
}
