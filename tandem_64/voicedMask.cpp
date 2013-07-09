#include "voicedMask.h"
#include <iostream>
#include <fstream>
using namespace std;

int voicedMask::checkPitchCon(double p1, double p2)
{
	if ( (fabs(p1-p2)<(p1*0.2)) && (fabs(p1-p2)<(p2*0.2)) )
		return(1);
	else return(0);
}

int voicedMask::checkMaskCon(double *mask1, double *mask2, int frame)
{
	double count1=0, count2=0, count3=0, count4=0;

	for(int chan=0; chan<numberChannel; chan++)
	{
		if (mask1[chan]>0) count1++; //=sqrt(corrLgm[frame].acf[chan][0]+1);
		if (mask2[chan]>0) count2++; //=sqrt(corrLgm[frame+1].acf[chan][0]+1);

		if ((mask1[chan]>0) && (mask2[chan]>0))
		{
			count3++; //=sqrt(corrLgm[frame].acf[chan][0]+1);
			count4++; //=sqrt(corrLgm[frame+1].acf[chan][0]+1);
		}
	}

	if ( ((count3*2)>count1) && ((count4*2)>count2) )
		return(1);
	else return(0);
}

int voicedMask::checkMaskCon2(double *mask1, double *mask2, int frame)
{
	double count1=0, count2=0, count3=0, count4=0;

	for(int chan=0; chan<numberChannel; chan++)
	{
		if (mask1[chan]>0) count1++; //=sqrt(corrLgm[frame].acf[chan][0]+1);
		if (mask2[chan]>0) count2++; //=sqrt(corrLgm[frame+1].acf[chan][0]+1);

		if ((mask1[chan]>0) && (mask2[chan]>0))
		{
			count3++; //=sqrt(corrLgm[frame].acf[chan][0]+1);
			count4++; //=sqrt(corrLgm[frame+1].acf[chan][0]+1);
		}
	}

	if ( ((count3*2)>count1) && ((count4*2)>count2) )
		return(1);
	else return(0);
}

int voicedMask::isOverlap(int m, int n)
{
	int frame;

	for(frame=Pitch[n].sFrame; frame<=Pitch[n].eFrame; frame++)
	{
		if ( (Pitch[n].value[frame]==Pitch[m].value[frame]) && (Pitch[n].value[frame]>0) )
		{
			if (checkMaskCon(Pitch[n].mProb[frame].value, Pitch[m].mProb[frame].value, frame)>0)
				return(1);
		}
	}

	return(0);
}

int voicedMask::isConnected(int m, int n)
{
	if(Pitch[n].sFrame==(Pitch[m].eFrame+1))
	{
		if (checkPitchCon(Pitch[n].value[Pitch[n].sFrame], Pitch[m].value[Pitch[m].eFrame])==0) return(0);

		else
		{
			if (checkMaskCon(Pitch[m].mProb[Pitch[m].eFrame].value, Pitch[n].mProb[Pitch[n].sFrame].value, Pitch[m].eFrame)==0) return(0);
			else return(1);
		}
	}

	if(Pitch[n].eFrame==(Pitch[m].sFrame-1))
	{
		if (checkPitchCon(Pitch[n].value[Pitch[n].eFrame], Pitch[m].value[Pitch[m].sFrame])==0) return(0);

		else
		{
			if (checkMaskCon(Pitch[n].mProb[Pitch[n].eFrame].value, Pitch[m].mProb[Pitch[m].sFrame].value, Pitch[n].eFrame)==0) return(0);
			else return(1);
		}
	}
	
	return(0);
}

void voicedMask::mergeContour(void)
{
	int n=0, m, judge, chan, frame;

    while (n<(numContour-1))
	{
	    m=n+1;
		
		while (m<numContour)
		{
			judge=0;

			if (isOverlap(m, n)>0) judge=1;
			else{
				if (isConnected(m, n)>0) judge=1;
			}
				
			if (judge>0)
			{
				for(frame=Pitch[m].sFrame; frame<=Pitch[m].eFrame; frame++)
				{
					Pitch[n].value[frame]=Pitch[m].value[frame];

					for(chan=0; chan<numberChannel; chan++)
						Pitch[n].mProb[frame].value[chan]=Pitch[m].mProb[frame].value[chan];
				}

				removeContour(m);

				break;
			}

			m++;
		}

		n++;
	}
}

void voicedMask::reDetermineMask(int nCon)
{
	int frame, m, chan;

	for(frame=Pitch[nCon].sFrame; frame<=Pitch[nCon].eFrame; frame++)
	{
		if(Pitch[nCon].value[frame]>0)
		{
			for(m=0; m<numContour; m++)
			{
				if ( (m!=nCon) & (Pitch[m].value[frame]>0) )
				{
					for(chan=0; chan<numberChannel; chan++)
					{
						if ( (Pitch[nCon].mProb[frame].value[chan]>0) && (Pitch[nCon].mProb[frame].value[chan]<Pitch[m].mProb[frame].value[chan]) )
							Pitch[nCon].mProb[frame].value[chan]=-1;
					}
				}
			}
		}
	}
}

void voicedMask::expandMask(void)
{
	int n, chan, frame;

	for(n=0; n<numContour; n++)
	{
		frame=Pitch[n].sFrame-1;
		if (frame>4 && frame<numFrame-1)
		{
			for(int m=0; m<numFrame; m++)
				Pitch[numContour].value[m]=0;
			
			Pitch[numContour].sFrame=frame;
			Pitch[numContour].eFrame=frame;

			int dis=double(Pitch[n].value[frame+1])*0.2;
						
			int d1=Pitch[n].value[frame+1]-int(dis); if (d1<min_pitch) d1=min_pitch;
			int d2=Pitch[n].value[frame+1]+int(dis); if (d2>max_pitch) d2=max_pitch;

			Pitch[numContour].value[frame]=maskToPitchMAP(frame, Pitch[n].mProb[frame+1].value, d1, d2);
	
			for(chan=0; chan<numberChannel; chan++)
			{
				Pitch[numContour].mProb[frame].value[chan]=multiUnitProb(frame, chan, numContour);
				if (Pitch[n].mProb[frame+1].value[chan]<=0) Pitch[numContour].mProb[frame].value[chan]=0;
			}

			if( (checkPitchCon(Pitch[numContour].value[frame], Pitch[n].value[frame+1])>0) && 
				(checkMaskCon(Pitch[numContour].mProb[frame].value, Pitch[n].mProb[frame+1].value, frame)>0) )
			{
				Pitch[n].value[frame]=Pitch[numContour].value[frame];
				Pitch[n].sFrame=frame;

				for(chan=0; chan<numberChannel; chan++)
					Pitch[n].mProb[frame].value[chan]=Pitch[numContour].mProb[frame].value[chan];
			}
		}
	     
    
		frame=Pitch[n].eFrame+1;
		
		if (frame<numFrame)
		{
			for(int m=0; m<numFrame; m++)
				Pitch[numContour].value[m]=0;
			
			Pitch[numContour].sFrame=frame;
			Pitch[numContour].eFrame=frame;

			int dis=double(Pitch[n].value[frame-1])*0.2;
						
			int d1=Pitch[n].value[frame-1]-int(dis); if (d1<min_pitch) d1=min_pitch;
			int d2=Pitch[n].value[frame-1]+int(dis); if (d2>max_pitch) d2=max_pitch;

			Pitch[numContour].value[frame]=maskToPitchMAP(frame, Pitch[n].mProb[frame-1].value, d1, d2);
	
			for(chan=0; chan<numberChannel; chan++)
			{
				Pitch[numContour].mProb[frame].value[chan]=multiUnitProb(frame, chan, numContour);
				if (Pitch[n].mProb[frame-1].value[chan]<=0) Pitch[numContour].mProb[frame].value[chan]=0;
			}

			if( (checkPitchCon(Pitch[numContour].value[frame], Pitch[n].value[frame-1])>0) && 
				(checkMaskCon(Pitch[n].mProb[frame-1].value, Pitch[numContour].mProb[frame].value, frame-1)>0) )
			{
				Pitch[n].value[frame]=Pitch[numContour].value[frame];
				Pitch[n].eFrame=frame;

				for(chan=0; chan<numberChannel; chan++)
					Pitch[n].mProb[frame].value[chan]=Pitch[numContour].mProb[frame].value[chan];
			}
		}
	}
}

void voicedMask::removeDuplicate(int *p1, int *p2)
{
	int frame;
	double f1, f2;

	for(frame=0; frame<numFrame; frame++)
	{
		f1=p1[frame];
		f2=p2[frame];

        if ( (f1>0) && (fabs(f1-f2)<(0.05*f1)) && (fabs(f1-f2)<(0.05*f2)) )
			p2[frame]=0;
	}
}

void voicedMask::switchCandidate(int *p1, mask *m1, int *p2, mask *m2)
{
	int frame, chan, tmp1;
	double tmp2;

	for(frame=1; frame<numFrame; frame++)
	{
		if ( ( (checkPitchCon(p1[frame], p2[frame-1])>0 ) && (checkMaskCon(m2[frame-1].value, m1[frame].value, frame-1)>0) ) || 
			 ( (checkPitchCon(p1[frame-1], p2[frame])>0 ) && (checkMaskCon(m1[frame-1].value, m2[frame].value, frame-1)>0) ) )
		{
			tmp1=p1[frame]; p1[frame]=p2[frame]; p2[frame]=tmp1;

			for(chan=0; chan<numberChannel; chan++)
			{
				tmp2=m1[frame].value[chan]; 
				m1[frame].value[chan]=m2[frame].value[chan]; 
				m2[frame].value[chan]=tmp2;
			}
		}
	}
}

void voicedMask::findContour(int *p1, mask *m1)
{
	int frame=1, chan;
	while(frame<(numFrame-1))
	{
		if (p1[frame]>0)
		{			
			if ( (checkPitchCon(p1[frame], p1[frame-1])>0) && (checkMaskCon(m1[frame-1].value, m1[frame].value, frame-1)>0) &&
				 (checkPitchCon(p1[frame], p1[frame+1])>0) && (checkMaskCon(m1[frame].value, m1[frame+1].value, frame)>0) )
			{
				//newContour(numContour);

				Pitch[numContour].value[frame-1]=p1[frame-1];
				Pitch[numContour].value[frame]=p1[frame];
				Pitch[numContour].value[frame+1]=p1[frame+1];

				Pitch[numContour].sFrame=frame-1;
				Pitch[numContour].eFrame=frame+1;

				for(chan=0; chan<numberChannel; chan++)
				{
					Pitch[numContour].mProb[frame-1].value[chan]=m1[frame-1].value[chan];
					Pitch[numContour].mProb[frame].value[chan]=m1[frame].value[chan];
					Pitch[numContour].mProb[frame+1].value[chan]=m1[frame+1].value[chan];
				}

				frame += 2;

				while( (checkPitchCon(p1[frame], p1[frame-1])>0) && (checkMaskCon(m1[frame-1].value, m1[frame].value, frame-1)>0) )
				{
					Pitch[numContour].value[frame]=p1[frame];

					for(chan=0; chan<numberChannel; chan++)
						Pitch[numContour].mProb[frame].value[chan]=m1[frame].value[chan];

					Pitch[numContour].eFrame=frame;

					frame++;
				}

				numContour++;
			}
		}

		frame++;
	}
}

void voicedMask::convCont(int *p1, mask *m1, int *p2, mask *m2)
{
	removeDuplicate(p1, p2);

	switchCandidate(p1, m1, p2, m2);

	numContour=0;

	findContour(p1, m1);
	findContour(p2, m2);
}

void voicedMask::reEstimatePitch(int nCon)
{
	int frame, chan, d1, d2, dis;

	Pitch[nCon].sFrame=numFrame;
	Pitch[nCon].eFrame=0;

	for(frame=0; frame<numFrame; frame++)
	{
		Pitch[nCon].oldValue[frame]=Pitch[nCon].value[frame];

		if(Pitch[nCon].oldValue[frame]==0) Pitch[nCon].value[frame]=maskToPitchMAP(frame, Pitch[nCon].mProb[frame].value, min_pitch, max_pitch);
		else
		{
			dis=double(Pitch[nCon].value[frame])*0.2;
						
			d1=Pitch[nCon].value[frame]-int(dis); if (d1<min_pitch) d1=min_pitch;
			d2=Pitch[nCon].value[frame]+int(dis); if (d2>max_pitch) d2=max_pitch;
			Pitch[nCon].value[frame]=maskToPitchMAP(frame, Pitch[nCon].mProb[frame].value, d1, d2);
		}

		if(Pitch[nCon].value[frame]>0)
		{
			Pitch[nCon].eFrame=frame;

			if (Pitch[nCon].sFrame>frame) Pitch[nCon].sFrame=frame;
		}
	}
	
	if(Pitch[nCon].eFrame>=Pitch[nCon].sFrame)
	{
		for(frame=Pitch[nCon].sFrame; frame<=Pitch[nCon].eFrame; frame++)
		{
			for(chan=0; chan<numberChannel; chan++)
				Pitch[nCon].mProb[frame].value[chan]=multiUnitProb(frame, chan, nCon);
		}
	}
	
	reDetermineMask(nCon);
}

void voicedMask::checkContour(int nCon)
{
	int frame;
	for(frame=0; frame<(numFrame+2); frame++)
		Pitch[nCon].indicate[frame]=0;

	for(frame=1; frame<(numFrame-1); frame++)
	{
		if ((checkPitchCon(Pitch[nCon].value[frame], Pitch[nCon].value[frame-1])>0) && 
			(checkMaskCon2(Pitch[nCon].mProb[frame-1].value, Pitch[nCon].mProb[frame].value, frame-1)>0) &&
			(checkPitchCon(Pitch[nCon].value[frame], Pitch[nCon].value[frame+1])>0)  && 
			(checkMaskCon2(Pitch[nCon].mProb[frame].value, Pitch[nCon].mProb[frame+1].value, frame)>0) )

			Pitch[nCon].indicate[frame]=Pitch[nCon].indicate[frame+1]=Pitch[nCon].indicate[frame+2]=1;
	}
}			

void voicedMask::developeContour(int nCon)
{
	int judge=1, frame, chan;
	int d1, d2;
	double dis;

	while (judge>0)
	{
		judge=0;
		for(frame=Pitch[nCon].sFrame; frame<=Pitch[nCon].eFrame; frame++)
		{
			if (Pitch[nCon].indicate[frame+1]==0)
			{
				if(Pitch[nCon].indicate[frame]>0)
				{
					dis=double(Pitch[nCon].value[frame-1])*0.2;
						
					d1=Pitch[nCon].value[frame-1]-int(dis); if (d1<min_pitch) d1=min_pitch;
					d2=Pitch[nCon].value[frame-1]+int(dis); if (d2>max_pitch) d2=max_pitch;
										  
					Pitch[nCon].value[frame]=maskToPitchMAP(frame, Pitch[nCon].mProb[frame].value, d1, d2);

					Pitch[nCon].indicate[frame]=1;
					//judge=1;
				}
				else if(Pitch[nCon].indicate[frame+2]>0)
				{
					dis=double(Pitch[nCon].value[frame+1])*0.2;
						
					d1=Pitch[nCon].value[frame+1]-int(dis); if (d1<min_pitch) d1=min_pitch;
					d2=Pitch[nCon].value[frame+1]+int(dis); if (d2>max_pitch) d2=max_pitch;
										  
					Pitch[nCon].value[frame]=maskToPitchMAP(frame, Pitch[nCon].mProb[frame].value, d1, d2);

					Pitch[nCon].indicate[frame]=1;
					//judg=1;
				}
			}
		}
	}
	
	for(frame=Pitch[nCon].sFrame; frame<=Pitch[nCon].eFrame; frame++)
	{
		for(chan=0; chan<numberChannel; chan++)
			Pitch[nCon].mProb[frame].value[chan]=multiUnitProb(frame, chan, nCon);
	}

	reDetermineMask(nCon);
}

void voicedMask::initPitchEst(void)
{
	int *p1, *p2, frame, chan;
	mask *m1, *m2;

	p1 = new int[numFrame];
	p2 = new int[numFrame];

	m1 = new mask[numFrame];
	m2 = new mask[numFrame];

	#ifdef DEBUG	
	FILE *out1=fopen("m1","w");
	FILE *out2=fopen("m2","w");
	FILE *out3=fopen("cm1","w");
	FILE *out4=fopen("cm2","w");
	FILE *f1=fopen("f1","w");
	FILE *f2=fopen("f2","w");
	FILE *f3=fopen("f3","w");
	FILE *f4=fopen("f4","w");
	#endif

	for(frame=0; frame<10; frame++)
	{
		p1[frame]=0; p2[frame]=0;
	}

	for(frame=10; frame<numFrame; frame++)
	{
		for(int delay=min_pitch; delay<=max_pitch; delay++)
			pitchProb(frame, delay);

		p1[frame]=p2[frame]=0;

		int count=0;
		for(chan=0; chan<numberChannel; chan++)
		{
			if ( (corrLgm[frame].cross[chan]>thd_cross) || (corrLgm[frame].crossEv[chan]>thd_evCross) )
			{
				m1[frame].value[chan]=1;	
				count++;
			}
			else m1[frame].value[chan]=0;
		}
		#ifdef DEBUG
		for (int c=0; c<numberChannel;c++)			
		{
			fprintf(out3,"%f ",m1[frame].value[c]);	
		}
		fprintf(out3,"\n");		

		#endif
				      
		if (count>2)
		{
			p1[frame]=maskToPitchML2(frame, m1[frame].value, min_pitch, max_pitch);

			if (p1[frame]>0)
			{
				count=0;
				for(chan=0; chan<numberChannel; chan++)
				{
					if ( (m1[frame].value[chan]>0) && (Prob[frame].sProb[chan][p1[frame]]<=theta_p))
					{
						m2[frame].value[chan]=1;
						count++;
					}
					else m2[frame].value[chan]=0;
				
					m1[frame].value[chan]=Prob[frame].sProb[chan][p1[frame]];
				}

				#ifdef DEBUG
				for (int c=0; c<numberChannel;c++)			
				{
					fprintf(out1,"%f ",m1[frame].value[c]);	
				}
				fprintf(out1,"\n");		
				// response ACF feature
				ofstream os;
//				if (frame==109) os.open("prob");
				for (int c=0; c<numberChannel;c++)
		                {
//					if (frame==109){
//						for (int delay=min_pitch; delay!=max_pitch; ++delay){
//							os<<Prob[frame].sProb[c][delay]<<" "<<flush;
//						}
//						os<<endl;
//					}
                		        fprintf(f1,"%f ",corrLgm[frame].acf[c][p1[frame]]/corrLgm[frame].acf[c][0]);
                		        fprintf(f2,"%f ",corrLgm[frame].acfEv[c][p1[frame]]/corrLgm[frame].acfEv[c][0]);
                		        fprintf(f3,"%f ",(p1[frame]+1)/(corrLgm[frame].zc[c]+1e-10)/2);
                		        fprintf(f4,"%f ",(p1[frame]+1)/(corrLgm[frame].zcEv[c]+1e-10)/2);
		                }
//				os.close();
                		fprintf(f1,"\n");
                		fprintf(f2,"\n");
                		fprintf(f3,"\n");
                		fprintf(f4,"\n");
				#endif
				#ifdef DEBUG
				for (int c=0; c<numberChannel;c++)			
				{
					fprintf(out4,"%f ",m2[frame].value[c]);	
				}
				fprintf(out4,"\n");		
				#endif

				if(count>2)
				{
					p2[frame]=maskToPitchML2(frame, m2[frame].value, min_pitch, max_pitch);

					if (p2[frame]>0)
					{
						for(chan=0; chan<numberChannel; chan++)
							m2[frame].value[chan]=Prob[frame].sProb[chan][p2[frame]];
					}
				}

				#ifdef DEBUG
				for (int c=0; c<numberChannel;c++)			
				{
					fprintf(out2,"%f ",m2[frame].value[c]);	
				}
				fprintf(out2,"\n");	
				#endif	
			}
		}		
	}
	#ifdef DEBUG
	fclose(out1);
	fclose(out2);
	fclose(out3);
	fclose(out4);
        fclose(f1);
        fclose(f2);
        fclose(f3);
        fclose(f4);
	#endif

	#ifdef DEBUG	
	out1=fopen("p1","w");
	out2=fopen("p2","w");
	for(int m=0; m<numFrame; m++)
	{
		fprintf(out1, "%d ", p1[m]);			
		fprintf(out2, "%d ", p2[m]);			
	}
	fclose(out1);
	fclose(out2);
	#endif
	
	convCont(p1, m1, p2, m2);

	delete [] p1;
	delete [] p2;
	delete [] m1;
	delete [] m2;
}

void voicedMask::iterativePitchEst(void)
{
	int n, judge=1, increase=10, count=0, frame;

	for(n=0; n<numContour; n++)
		reDetermineMask(n);

	while ( (judge>0) && ( (count<20)  || (increase>0) ) && (count<50) )
	{
		count++;
		increase=0;
		judge=0;
    
		printf("%d ", count);

		expandMask();
    
		n=0;
		while(n<numContour)
		{
        		maskToPitch(n);
			
			int count1=0, count2=0;

			for(frame=0; frame<=numFrame; frame++)
			{
				 if(Pitch[n].value[frame]>0) count1++;
				 if(Pitch[n].oldValue[frame]>0) count2++;

				 if(Pitch[n].value[frame]!=Pitch[n].oldValue[frame]) judge=1;
			}

			increase += count1-count2;

			if (count1==0) removeContour(n);
			else n++;
		}		       

		for(n=0; n<numContour; n++)
			reDetermineMask(n);
    
	    mergeContour();
	}
	
	printf("\n");
}

void voicedMask::maskToPitch(int nCon)
{
	int frame, chan, judge, oldjudge;

	reEstimatePitch(nCon);

	checkContour(nCon);

	oldjudge=Pitch[nCon].eFrame;
	judge=0;
	for(frame=Pitch[nCon].sFrame; frame<=Pitch[nCon].eFrame; frame++)
	{
		if ( (Pitch[nCon].value[frame]>0) && (Pitch[nCon].indicate[frame+1]==0) )
			judge++;
	}
	
	while ( (judge<oldjudge) && (judge>0) )
	{
		oldjudge=judge;
    
		developeContour(nCon);    
    
		checkContour(nCon);
   
		judge=0;
		for(int frame=Pitch[nCon].sFrame; frame<=Pitch[nCon].eFrame; frame++)
		{
			if ( (Pitch[nCon].value[frame]>0) && (Pitch[nCon].indicate[frame+1]==0) )
				judge++;
		}
	}

	Pitch[nCon].sFrame=numFrame; Pitch[nCon].eFrame=0;
	for(frame=0; frame<numFrame; frame++)
	{
		if(Pitch[nCon].indicate[frame+1]==0)
		{
			Pitch[nCon].value[frame]=0;

			for(chan=0; chan<numberChannel; chan++)
				Pitch[nCon].mProb[frame].value[chan]=0;
		}

		else
		{
			Pitch[nCon].eFrame=frame;
			if(Pitch[nCon].sFrame>frame) Pitch[nCon].sFrame=frame;
		}
	}
}

void voicedMask::dtmPitchMask(void)
{
	initPitchEst();

	iterativePitchEst();

	// refine pitch contours
	int n=0;
        while(n<numContour)
        {
                float count=0;
                float sum=0;
                for(int frame=Pitch[n].sFrame; frame<=Pitch[n].eFrame; frame++)
                {
			// remove pitch points near the boundary
			if ( (Pitch[n].value[frame]<min_pitch*1.1 || Pitch[n].value[frame]>max_pitch*.9 ) && Pitch[n].value[frame]>0 )
                        {
                                removeContour(n);
                                n--;
                                break;
                        }
                }
                n++;
        }
}
