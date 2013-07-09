#include "tool.h"

double sinc(double tim, double wn)
{
	if(tim == 0) return(wn);
	else return( sin(wn * PI * tim) / (PI * tim) );
}

double sincDiffOne(double tim, double wn)
{
	if(tim == 0) return(0);
	else
	{
		return( cos(wn * PI * tim) * wn / tim - sin(wn * PI *tim) / (PI * tim * tim) );
	}
}

double sincDiffTwo(double tim, double wn)
{
	if(tim == 0) return( -wn*wn*wn/PI/PI);
	else
	{
		return( -wn*wn*PI*sin(wn*PI*tim)/tim - 2*wn*cos(wn*PI*tim)/tim/tim + 2*sin(wn*PI*tim)/PI/tim/tim/tim );
	}
}

void fft(double *inputR, double *inputI, int N, double direct)
{
	long sigL, i, j, k, n, period, twoPeriod;
	double tmpR, tmpI, uR, uI, wR, wI;

	sigL = long(pow(2.0, N));

	j = 1;
	for(i=1; i<sigL; i++)
	{
		if(i < j)
		{
			tmpR = inputR[j-1];
			tmpI = inputI[j-1];

			inputR[j-1] = inputR[i-1];
			inputI[j-1] = inputI[i-1];

			inputR[i-1] = tmpR;
			inputI[i-1] = tmpI;
		}

		k = sigL/2;
		while (k < j){ j -=  k;	k /= 2;	}
		j += k;
	}

	for(n=1; n<=N; n++ )
    {  
		twoPeriod = long(pow(2.0, n));
        period = twoPeriod/2;
        uR = 1.0; 
        uI = 0.0; 
        wR = double( cos( PI/period ) ); 
        wI = double( -1.0 * sin( PI/period * direct) );

        for(j=0; j<period; j++ ) 
        {  
			for(i=j; i<sigL; i+=twoPeriod)
			{
				tmpR = inputR[i+period]*uR - inputI[i+period]*uI;
                tmpI = inputR[i+period]*uI + inputI[i+period]*uR;
				
				inputR[i+period] = inputR[i] - tmpR; 
				inputI[i+period] = inputI[i] - tmpI; 
				inputR[i] += tmpR ; 
				inputI[i] += tmpI; 
			}
			tmpR = uR*wR - uI*wI; 
			tmpI = uR*wI + uI*wR; 
			uR = tmpR; 
			uI = tmpI; 
		} 
	} 
}

int maxPos(double *input, int d1, int d2)
{
	int judge=d1;
	double mv=input[d1];
	for(int step=d1+1; step<=d2; step++)
	{	
		if (input[step]>mv)
		{
			judge=step;
			mv=input[step];
		}
	}

	return(judge);
}

int minPos(double *input, int d1, int d2)
{
	int judge=d1;
	double mv=input[d1];
	for(int step=d1+1; step<=d2; step++)
	{	
		if (input[step]<mv)
		{
			judge=step;
			mv=input[step];
		}
	}

	return(judge);
}

double bessi0(double x)
{
	double ax,ans;
	double y;

	ax = fabs(x);
	if (ax < 3.75)
	{
		y = x/3.75;
		y *= y;
		ans = 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492 + y * (0.2659732 + y * (0.360768e-1 + y*0.45813e-2)))));
	}
	else
	{
		y = 3.75 / ax;
		ans = (exp(ax) / sqrt(ax)) * (0.39894228 + y * (0.1328592e-1 + y * (0.225319e-2 + y * (-0.157565e-2 + y * (0.916281e-2 + y * (-0.2057706e-1 + y * (0.2635537e-1 + y * (-0.1647633e-1 + y * 0.392377e-2))))))));
	}
	return ans;
}

double bessi1(double x)
{
	double ax, y, ans;
	if (fabs(x)<3.75)
	{
		y = x/3.75;
		y *= y;
		ans = x * ( 0.5 + y*( 0.87890594 + y*( 0.51498869 + y*( 0.15084934 + y*( 0.02658733 + y* (0.00301532 + y*0.00032411 ) ) ) ) ) );
	}
	else
	{	
		ax = fabs(x);
		y = 3.75/ax;
		ans = (exp(ax)/sqrt(ax)) * (0.39894228 + y*( -0.03988024 + y*( -0.00362018 + y*( 0.00163801 + y*( -0.01031555 + y*( 0.02282967 + y*( -0.02895312 + y*( 0.01787654 + y*(-0.00420059) ) ) ) ) ) ) ) );

		if (x<0) ans *= -1;
	}
   
    return(ans);
}

double zeroCross(double *input, int sLen)
{
	int nCross=-1;
	double currentSign=1, sp, ep;
	double zC=sLen*2;

	if (input[0]<0) currentSign=-1;
	
	for(int n=1; n<sLen; n++)
	{
		if ((input[n]*currentSign)<0)
		{
			nCross++;
			if (nCross==0) sp=n;
			else ep=n;

			currentSign *= -1;

			if ( ((nCross%2)==0) && (nCross>0))
				zC=(ep-sp)/double(nCross);
		}
	}

	if(nCross==1) zC=ep*2.0/3.0;

	return(zC);
}

int findLocalPeak(double *input, int p1, int p2, double theta)
{
        int pos=p1;
        double value=theta;

        for(int p=p1+1; p<p2; p++)
        {
                if( (input[p]>input[p-1]) && (input[p]>input[p+1]) && (input[p]>value) )
                {
                        value=input[p];
                        pos=p;
                }
        }

        return(pos);
}

int findLocalValley(double *input, int p1, int p2, double theta)
{
        int pos=p1;
        double value=theta;

        for(int p=p1+1; p<p2; p++)
        {
                if( (input[p]<input[p-1]) && (input[p]<input[p+1]) && (input[p]<value) )
                {
                        value=input[p];
                        pos=p;
                }
        }

        return(pos);
}
