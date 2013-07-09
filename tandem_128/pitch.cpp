#include "pitch.h"

pitchMask::pitchMask(featurePara x):feature(x)
{
	theta_p = x.theta_p;
 
	char fname1[255], fname2[255], fname3[255];
	memIndi=0; 
		
	sprintf(fname1, "net//su.dat");
	sprintf(fname2, "net//mu.8.2.dat");
	sprintf(fname3, "net//pu.2.dat");
	readNet(fname1, fname2, fname3);
}

pitchMask::~pitchMask(void)
{
	if (memIndi>0) deletePitchPara();
}

void pitchMask::readNet(char *fname1, char *fname2, char *fname3)
{
	FILE *fp;

	fp = fopen(fname1, "rb");
	for(int chan=0; chan<numberChannel; chan++)
	{
		for(int k=0; k<6; k++)
			fread(sNet[chan].IW[k], sizeof(float), SHU_NUM, fp);

		fread(sNet[chan].b1, sizeof(float), SHU_NUM, fp);
		fread(sNet[chan].LW, sizeof(float), SHU_NUM, fp);
		fread(&sNet[chan].b2, sizeof(float), 1, fp);
	}
	fclose(fp);

	fp = fopen(fname2, "rb");
	for(int chan=0; chan<numberChannel; chan++)
	{	
		int m=(SIDE_CHAN*2+1)*(SIDE_FRAME*2+1);

		for(int k=0; k<m; k++)
			fread(mNet[chan].IW[k], sizeof(float), MHU_NUM, fp);
		fread(mNet[chan].b1, sizeof(float), MHU_NUM, fp);
		fread(mNet[chan].LW, sizeof(float), MHU_NUM, fp);
		fread(&mNet[chan].b2, sizeof(float), 1, fp);
	}
	fclose(fp);

	fp = fopen(fname3, "rb");
	for(int k=0; k<(4); k++)
		fread(pNet.IW[k], sizeof(float), PHU_NUM, fp);
	fread(pNet.b1, sizeof(float), PHU_NUM, fp);
	fread(pNet.LW, sizeof(float), PHU_NUM, fp);
	fread(&pNet.b2, sizeof(float), 1, fp);
	fclose(fp);
}

void pitchMask::singleUnitProb(int frame, int chan, int delay)
{
	double f[6];

	f[0] = corrLgm[frame].acf[chan][delay]/corrLgm[frame].acf[chan][0];
	f[1] = corrLgm[frame].acfEv[chan][delay]/corrLgm[frame].acfEv[chan][0];

	f[2] = double(delay+1)/(corrLgm[frame].zc[chan]+1e-10)/2;
	f[3] = double(delay+1)/(corrLgm[frame].zcEv[chan]+1e-10)/2;

	f[4] = int(f[2]); if ((f[2]-f[4])>=0.5) f[4]++;
	f[5] = int(f[3]); if ((f[3]-f[5])>=0.5) f[5]++;

	f[2] = fabs(f[2]-f[4]);
	f[3] = fabs(f[3]-f[5]);

	Prob[frame].hNum[chan][delay] = f[4];
	Prob[frame].hNumEv[chan][delay] = f[5];

	f[4] *= 0.05; if (f[4]>1) f[4]=1; if (f[4]==0) f[4]=2;
	f[5] *= 0.1; if (f[5]>1) f[5]=1; if (f[4]==0) f[4]=2;

	double s1[SHU_NUM], s2=sNet[chan].b2;

	for(int n=0; n<SHU_NUM; n++)
	{
		s1[n] = sNet[chan].b1[n];

		for(int m=0; m<6; m++)
			s1[n] += f[m]*sNet[chan].IW[m][n];

		s2 += sigmoid(s1[n]) * sNet[chan].LW[n];
	}

	Prob[frame].sProb[chan][delay] = sigmoid(s2);
}

double pitchMask::multiUnitProb(int frame, int chan, int nCon)
{
	double f[(SIDE_CHAN*2+1)*(SIDE_FRAME*2+1)];
	int ch, fm;

	for(ch=chan-SIDE_CHAN; ch<=chan+SIDE_CHAN; ch++)
	{
		int offset=(SIDE_FRAME)*(SIDE_CHAN*2+1)+ch-chan+SIDE_CHAN;

		f[offset]=0;
		if ((ch>=0) && (ch<numberChannel))		
		{
			int p=Pitch[nCon].value[frame];
			if (Pitch[nCon].mProb[frame].value[ch]<-0.99) f[offset]=0;
			else
			{
				if(Prob[frame].sProb[ch][p]<-1) singleUnitProb(frame, ch, p);
				f[offset]=Prob[frame].sProb[ch][p];
			}
		}
	}

	for(fm=frame-1; fm>=frame-SIDE_FRAME; fm--)
		for(ch=chan-SIDE_CHAN; ch<=chan+SIDE_CHAN; ch++)
		{
			int offset=(fm-frame+SIDE_FRAME)*(SIDE_CHAN*2+1)+ch-chan+SIDE_CHAN;

			f[offset]=0;
			if ((fm>=0) && (fm<numFrame) && (ch>=0) && (ch<numberChannel))		
			{
				int p=Pitch[nCon].value[fm];
				if (fm<Pitch[nCon].sFrame) f[offset]=f[offset+SIDE_CHAN*2+1];
				else
				{
					if (Pitch[nCon].mProb[fm].value[chan]<-0.99) f[offset]=0;
					else
					{
						if(Prob[fm].sProb[ch][p]<-1) singleUnitProb(fm, ch, p);
						f[offset]=Prob[fm].sProb[ch][p];
					}
				}
			}
		}

	for(fm=frame+1; fm<=frame+SIDE_FRAME; fm++)
		for(ch=chan-SIDE_CHAN; ch<=chan+SIDE_CHAN; ch++)
		{
			int offset=(fm-frame+SIDE_FRAME)*(SIDE_CHAN*2+1)+ch-chan+SIDE_CHAN;

			f[offset]=0;
			if ((fm>=0) && (fm<numFrame) && (ch>=0) && (ch<numberChannel))		
			{
				int p=Pitch[nCon].value[fm];
				if (fm>Pitch[nCon].eFrame) f[offset]=f[offset-SIDE_CHAN*2-1];
				else
				{
					if (Pitch[nCon].mProb[fm].value[chan]<-0.99) f[offset]=0;
					else
					{
						if(Prob[fm].sProb[ch][p]<-1) singleUnitProb(fm, ch, p);
						f[offset]=Prob[fm].sProb[ch][p];
					}
				}
			}
		}

	double s1[MHU_NUM], s2=mNet[chan].b2;

	for(int n=0; n<MHU_NUM; n++)
	{
		s1[n] = mNet[chan].b1[n];

		for(int m=0; m<(SIDE_CHAN*2+1)*(SIDE_FRAME*2+1); m++)
			s1[n] += f[m]*mNet[chan].IW[m][n];

		s2 += sigmoid(s1[n]) * mNet[chan].LW[n];
	}

	return(sigmoid(s2));
}

double pitchMask::compareTwoCandidateMAP(int frame, double *mask, int p1, int p2)
{
	int temp, chan;
	double count;

	if (p1>p2)
	{
		temp=p2; p2=p1; p1=temp;
	}

	double f[4];

	count=0;
	f[3]=double(p2)/double(p1);
	f[1]=int(f[3]); if ((f[3]-f[1])>0.5) f[1]++;
	f[3]-=f[1];
	f[1]=1/f[1];
	f[0]=f[2]=0;

	for(chan=0; chan<numberChannel; chan++)
	{
		if(mask[chan]>0)
		{
			count++;

			if(Prob[frame].sProb[chan][p1]<-1) singleUnitProb(frame, chan, p1);
			if(Prob[frame].sProb[chan][p2]<-1) singleUnitProb(frame, chan, p2);

			f[0] += Prob[frame].sProb[chan][p1];
			f[2] += Prob[frame].sProb[chan][p2];
		}
	}

	if(count>0)
	{
		f[0] /= count;
		f[2] /= count;
	}
		
	double s1[PHU_NUM], s2=pNet.b2;

	for(int n=0; n<PHU_NUM; n++)
	{
		s1[n] = pNet.b1[n];

		for(int m=0; m<(4); m++)
			s1[n] += f[m]*pNet.IW[m][n];

		s2 += sigmoid(s1[n]) * pNet.LW[n];
	}

	return(sigmoid(s2));
}

void pitchMask::removeContour(int nCon)
{
	int k, frame, chan;
	if((nCon+1)<numContour)
	{
		for(k=nCon+1; k<numContour; k++)
		{
			Pitch[k-1].sFrame=Pitch[k].sFrame;
			Pitch[k-1].eFrame=Pitch[k].eFrame;

			for(frame=0; frame<numFrame; frame++)
			{
				if (Pitch[k].value[frame]>0)
				{
					for(chan=0; chan<numberChannel; chan++)
						Pitch[k-1].mProb[frame].value[chan]=Pitch[k].mProb[frame].value[chan];
				}
				else if (Pitch[k-1].value[frame]>0)
				{
					for(chan=0; chan<numberChannel; chan++)
						Pitch[k-1].mProb[frame].value[chan]=0;
				}

				Pitch[k-1].value[frame]=Pitch[k].value[frame];
			}						
		}
	}


	numContour--;
}

void pitchMask::newContour(int nCon)
{
	int frame, chan;

	Pitch[nCon].value = new int[numFrame];

	Pitch[nCon].oldValue = new int[numFrame];

	Pitch[nCon].indicate = new int[numFrame+2];

	Pitch[nCon].mProb = new mask[numFrame];

	for(frame=0; frame<numFrame; frame++)
	{
		Pitch[nCon].value[frame]=0;

		for(chan=0; chan<numberChannel; chan++)
			Pitch[nCon].mProb[frame].value[chan]=0;
	}
}

void pitchMask::deleteContour(int nCon)
{
	delete [] Pitch[nCon].value;

	delete [] Pitch[nCon].oldValue;

	delete [] Pitch[nCon].indicate;

	delete [] Pitch[nCon].mProb;
}

int pitchMask::maskToPitchACF(int frame, double *mask, int d1, int d2)
{
	int chan, delay, pitch, count;
	double sumF[MAX_D];
	double maxV;

	for(delay=d1; delay<=d2; delay++)
		sumF[delay]=0;

	count=0;
	for(chan=0; chan<numberChannel; chan++)
	{
		if(mask[chan]>0)
		{
			count++;

			for(delay=d1; delay<=d2; delay++)
				sumF[delay] += corrLgm[frame].acf[chan][delay]/corrLgm[frame].acf[chan][0] + corrLgm[frame].acfEv[chan][delay]/corrLgm[frame].acfEv[chan][0];
		}
	}

	if (count==0) return(0);

	pitch=d1;
	maxV=sumF[d1];
	for(delay=d1+1; delay<=d2; delay++)
	{
		if (sumF[delay]>maxV)
		{
			maxV=sumF[delay];
			pitch=delay;
		}
	}

	return(pitch);
}

int pitchMask::maskToPitchML(int frame, double *mask, int d1, int d2)
{
	int chan, delay, pitch, count;
	double sumF[MAX_D];
	double maxV;

	for(delay=d1; delay<=d2; delay++)
		sumF[delay]=0;

	count=0;
	for(chan=0; chan<numberChannel; chan++)
	{
		if(mask[chan]>0)
		{
			count++;

			for(delay=d1; delay<=d2; delay++)
			{
				if(Prob[frame].sProb[chan][delay]<-1) singleUnitProb(frame, chan, delay);
				//sumF[delay] += log(Prob[frame].sProb[chan][delay]/2+0.5+1e-10);
				sumF[delay] += Prob[frame].sProb[chan][delay];
			}
		}
	}

	if (count==0) return(0);

	pitch=d1;
	maxV=sumF[d1];
	for(delay=d1+1; delay<=d2; delay++)
	{
		if (sumF[delay]>maxV)
		{
			maxV=sumF[delay];
			pitch=delay;
		}
	}

	return(pitch);
}

int pitchMask::maskToPitchML2(int frame, double *mask, int d1, int d2)
{
	int chan, delay, pitch, count;
	double sumF[MAX_D];
	double maxV;

	for(delay=d1; delay<=d2; delay++)
		sumF[delay]=0;

	count=0;
	for(chan=0; chan<numberChannel; chan++)
	{
		if(mask[chan]>0)
		{
			count++;

			for(delay=d1; delay<=d2; delay++)
			{
				if(Prob[frame].sProb[chan][delay]<-1) singleUnitProb(frame, chan, delay);
				//sumF[delay] += log(Prob[frame].sProb[chan][delay]/2+0.5+1e-10);
				if (Prob[frame].sProb[chan][delay]>theta_p)	sumF[delay] ++;
					//sumF[delay] += Prob[frame].sProb[chan][delay];
			}
		}
	}

	if (count==0) return(0);

	pitch=d1;
	maxV=sumF[d1];
	for(delay=d1+1; delay<=d2; delay++)
	{
		if (sumF[delay]>maxV)
		{
			maxV=sumF[delay];
			pitch=delay;
		}
	}

	return(pitch);
}

int pitchMask::maskToPitchMAP(int frame, double *mask, int d1, int d2)
{
	int chan, delay, pitch, count;
	double sumF[MAX_D];

	for(delay=d1; delay<=d2; delay++)
		sumF[delay]=0;

	count=0;
	for(chan=0; chan<numberChannel; chan++)
	{
		if(mask[chan]>0)
		{
			count++;

			for(delay=d1; delay<=d2; delay++)
			{
				if(Prob[frame].sProb[chan][delay]<-1) singleUnitProb(frame, chan, delay);
				//sumF[delay] += log(Prob[frame].sProb[chan][delay]/2+0.5+1e-10);
				sumF[delay] += Prob[frame].sProb[chan][delay];
			}
		}
	}

	if (count==0) return(0);

	double mv=sumF[d1]/double(count);
	pitch=d1;
	for(delay=d1+1; delay<d2; delay++)
	{
		sumF[delay] /= double(count);
		if (sumF[delay]>mv)
		{
			mv=sumF[delay];
			pitch=delay;
		}
	}

	for(delay=d1+1; delay<d2; delay++)
	{
		//if ( (sumF[delay]>(mv*0.25)) )
		if ( (sumF[delay]>sumF[delay-1]) && (sumF[delay]>sumF[delay+1]) && (sumF[delay]>(mv-0.1)) )
		{
			if (delay>pitch)
			{
				double f=compareTwoCandidateMAP(frame, mask, pitch, delay);
				if (f>0)
				{
					pitch=delay;
				}
			}
			else
			{
				double f=compareTwoCandidateMAP(frame, mask, delay, pitch);
				if (f<0) 
				{
					pitch=delay;
				}
			}
		}
	}

	return(pitch);
}

void pitchMask::pitchProb(int frame, int delay)
{
	int chan=0;

	for(chan=0; chan<numberChannel; chan++)
	{
		if (Prob[frame].sProb[chan][delay]<-1) singleUnitProb(frame, chan, delay);
	}
}

void pitchMask::newPitchPara(void)
{
	int frame, chan, n;

	if (memIndi>0) deletePitchPara();

	Prob = new pitchPara[numFrame];

	for(frame=0; frame<numFrame; frame++)
		for(chan=0; chan<numberChannel; chan++)
			for(n=0; n<max_delay; n++)
				Prob[frame].sProb[chan][n]=-2;

	numContour=0;

	for(n=0; n<MAX_CONTOUR; n++)
		newContour(n);

	memIndi=1;
}

void pitchMask::deletePitchPara(void)
{
	if (memIndi>0) delete [] Prob;

	for(int n=0; n<MAX_CONTOUR; n++)
		deleteContour(n);

	memIndi=0;
}
