#include "filter.h"

void filter::hammingBandPass(double wc, double wn)
{
	double step, k, tim;

	for (tim=0; tim<=fLength; tim++)
	{
		step = tim - fLength/2;
		k = 2*tim/double(fLength) - 1;
		
		rValue[int(tim)] = hamming(k)*sinc(step, wn);

		if (wc>0) rValue[int(tim)] *= 2*cos(step*wc* PI);
	}
}

void filter::hammingDiffOneBandPass(double wc, double wn)
{
	double tim, step;
	double k, c1=2/double(fLength), c2=wc*PI;

	for (tim=0; tim<=fLength; tim++)
	{
		step = tim - 1/c1;
		k = tim*c1 - 1;	
		
		rValue[int(tim)] = hammingDiffOne(k)*sinc(step, wn)*c1 + hamming(k)*sincDiffOne(step, wn);

		if (wc>0)
		{
			rValue[int(tim)] *= 2*cos(step*wc* PI);
			rValue[int(tim)] += -2*c2*sin(step*wc* PI) * hamming(k)*sinc(step, wn);
		}

		rValue[int(tim)] *= fs;
	}
}

void filter::hammingDiffTwoBandPass(double wc, double wn)
{
	double tim, step;
	double k, sum=0, c1=2/double(fLength), c2=wc*PI;

	for (tim=0; tim<=fLength; tim++)
	{
		step = tim - 1/c1;
		k = tim*c1 - 1;	
		
		rValue[int(tim)] = hammingDiffTwo(k)*sinc(step, wn)*c1*c1 + 2*hammingDiffOne(k)*c1*sincDiffOne(step, wn) + hamming(k)*sincDiffTwo(step, wn);

		if (wc>0)
		{
			rValue[int(tim)] *= 2*cos(step*wc*PI);
			rValue[int(tim)] += -4*c2*sin(step*wc*PI) * ( hammingDiffOne(k)*c1*sinc(step, wn) + hamming(k)*sincDiffOne(step, wn) );
			rValue[int(tim)] += -2*c2*c2*cos(step*wc*PI) * hamming(k)*sinc(step, wn);
		}

		rValue[int(tim)] *= fs*fs;
	}
}

filter::filter(double ws, double wn, double ts, int sampF, int type)
{
	fs = sampF;

	int n;
	double fLen;

	fLen=3.3/ts;

	if (fmod(fLen, 1)<0.5)	fLength = floor(fLen); 
	else fLength = ceil(fLen); 

	if ( (fLength%2)==1 ) fLength++;
	
	nOrder = ceil( log((double)fLength) / log(2.0) ) + 2; 
	nFFT = pow(2.0, nOrder);

	rValue = new double[nFFT]; iValue = new double[nFFT]; 
	inputR = new double[nFFT]; inputI = new double[nFFT];
	
	for(n=0; n<nFFT; n++)
	{
		rValue[n]=0; iValue[n]=0;
	}

	if (type==0) hammingBandPass(ws, wn);
	else if (type==1) hammingDiffOneBandPass(ws, wn);
	else hammingDiffTwoBandPass(ws, wn);
	
	fft(rValue, iValue, nOrder, 1);
}	

filter::~filter()
{
	delete [] rValue;
	delete [] iValue;
	delete [] inputR;
	delete [] inputI;
}

void filter::filtering(double *input, double *output, int sigLength)
{
	int nB, m, n, tim;

	nB = ceil(double(sigLength + fLength)/double(nFFT - fLength));
	
	for(m=0; m<nB; m++)
	{
		for(n=0; n<nFFT; n++)
		{
			tim = m * (nFFT - fLength) + n - fLength;

			inputR[n] = ( (tim>=0) & (tim<sigLength) ) ? input[tim]: 0;
			inputI[n] = 0;
		}
		
		fft(inputR, inputI, nOrder, 1);

		for(n=0; n<nFFT; n++)
		{
			double f1 = rValue[n]*inputR[n] - iValue[n]*inputI[n];
			double f2 = rValue[n]*inputI[n] + iValue[n]*inputR[n];

			inputR[n] = f1;
			inputI[n] = f2;
		}
		
		fft(inputR, inputI, nOrder, -1);

		for(n=fLength; n<nFFT; n++)
		{
			tim = m * (nFFT - fLength) + n - fLength*3/2;

			if ( (tim>=0) && (tim<sigLength) ) output[tim] = inputR[n] / double(nFFT);
		}
	}
}

void filter::kaiserPara(double delta, double transBw)
{
	double a, len;
	
	a= -20 * log10(delta);

	if (a <= 21) beta = 0;
	else if (a<= 50) beta = 0.5842 * pow(a-21.0, 0.4) + 0.07889 * (a-21);
	else beta = 0.1102 * (a - 8.7);

	len = (a - 7.95) / 14.36 / transBw;
	fLength = int(len);
	if ((len - fLength) < 0.5) fLength++;
	else fLength+=2;

	if (fLength%2 != 0) fLength++;
}
	
void filter::kaiserLowPass(double wn)
{
	int tim, step;
	double k, sum;

	for (tim=0; tim<=fLength; tim++)
	{
		k = 2*tim/double(fLength) - 1;
		rValue[tim] = bessi0( beta*sqrt(1- k*k)) / bessi0( beta );
	}

	sum=0;
	for (tim=0; tim<=fLength; tim++)
	{
		step = tim - fLength/2;
		if (step !=0) rValue[tim] *= sin(wn * PI * step) / PI / step;
		else rValue[tim] *= wn;

		sum += rValue[tim];
	}

	for (tim=0; tim<=fLength; tim++)
		rValue[tim] /= sum;
}

void filter::diffKaiserLowPass(double wn)
{
	double tim, x, y, f1, f2, df1, df2;
	
	double shift=fLength/2;

	double sum=0;
	for (tim=-shift; tim<=shift; tim++)
	{
		x = sqrt(1 - tim*tim/shift/shift);
		y = bessi0(beta);

		f1 = bessi0(beta*x) / y;
		f2 = sinc(tim, wn);

		if (x==0) df1=0;
		else df1 = - bessi1(beta*x) * (tim/shift/shift) / x * (beta/y); 
		df2 = sincDiffOne(tim, wn);

		rValue[int(tim+shift)] = f1 * df2 + f2 * df1;

		sum += f1*f2;
	}

	for (tim=0; tim<=fLength; tim++)
		rValue[int(tim)] /= sum;
}