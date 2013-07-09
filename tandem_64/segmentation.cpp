# include "segment.h"

void segment::newCandidate(chanCandidate *s)
{
	s->pos = new int[nSample];
	s->cPos = new int[nSample];
	s->cMark = new int[nSample];
}

void segment::deleteCandidate(chanCandidate *s)
{
	delete [] s->pos;
	delete [] s->cPos;
	delete [] s->cMark;
}

void segment::newFront(front *f)
{
	for(int chan=0; chan<numberChannel; chan++)
	{
		f->pos[chan] = new int[MAX_FRONT];
		f->cPos[chan] = new int[MAX_FRONT];
		f->cMark[chan] = new int[MAX_FRONT];
	}

	f->sChan = new int[MAX_FRONT];
	f->eChan = new int[MAX_FRONT];
	f->extendFront1 = new int[MAX_FRONT];
	f->extendFront2 = new int[MAX_FRONT];
}

void segment::deleteFront(front *f)
{
	for(int chan=0; chan<numberChannel; chan++)
	{
		delete [] f->pos[chan];
		delete [] f->cPos[chan];
		delete [] f->cMark[chan];
	}

	delete [] f->sChan;
	delete [] f->eChan;
	delete [] f->extendFront1;
	delete [] f->extendFront2;
}

void segment::newSegPara(void)
{
	int chan;

	if (memIndi==1) deleteSegPara();

	newCandidate(&tmpC1);
	newCandidate(&tmpC2);

	for(chan=0; chan<numberChannel; chan++)
	{
		seg[chan] = new int[nSample/factor+10];
		cross[chan] = new double[nSample/factor+10];
		crossEv[chan] = new double[nSample/factor+10];
		inten[chan] = new double[nSample];
		diffInten[chan] = new double[nSample];

		newCandidate(&oS1[chan]);
		newCandidate(&oS2[chan]);
		newCandidate(&fS1[chan]);
		newCandidate(&fS2[chan]);
	}

	newFront(&oF);
	newFront(&oF1);
	newFront(&oF2);
	newFront(&oF3);
	newFront(&fF);
	newFront(&fF1);
	newFront(&fF2);
	newFront(&fF3);

	memIndi=1;
}

void segment::deleteSegPara(void)
{
	if (memIndi==0) return;

	deleteCandidate(&tmpC1);
	deleteCandidate(&tmpC2);

	for(int chan=0; chan<numberChannel; chan++)
	{
		delete [] seg[chan];
		delete [] cross[chan];
		delete [] crossEv[chan];
		delete [] inten[chan];
		delete [] diffInten[chan];

		deleteCandidate(&oS1[chan]);
		deleteCandidate(&oS2[chan]);
		deleteCandidate(&fS1[chan]);
		deleteCandidate(&fS2[chan]);
	}

	deleteFront(&oF);
	deleteFront(&oF1);
	deleteFront(&oF2);
	deleteFront(&oF3);
	deleteFront(&fF);
	deleteFront(&fF1);
	deleteFront(&fF2);
	deleteFront(&fF3);

	memIndi=0;
}

void segment::freqSmooth(double freqScale, double *input[MAX_CHANNEL])
{
	double c1=1/sqrt(2*PI)/freqScale;
	double c2=0.5/freqScale/freqScale;

	double temp[MAX_CHANNEL];

	double *h;
	h=new double[(int)(8*freqScale+1)];

	for(int n=0; n<=8*freqScale; n++)
		h[n]=c1*exp(-pow(n-4*freqScale, 2)*c2);

	for(int f=0; f<nSample; f++)
	{
		for(int chan=0; chan<numberChannel; chan++)
		{
			temp[chan]=0;
			for(int n=0; n<=8*freqScale; n++)
			{
				int p=(double)(chan+4*freqScale-n);
				if( (p>=0) && (p<numberChannel) ) temp[chan] += h[n]*input[p][f];
			}
		}

		for(int chan=0; chan<numberChannel; chan++)
			input[chan][f]=temp[chan];
	}

	delete [] h;
}

void segment::determineTheta(int chan, double factor1, double factor2)
{
	int n, ch;

	double sumI=0;
	double sumI2=0;

	for(ch=0; ch<numberChannel; ch++)
	for(n=20; n<(nSample-20); n++)
	{
		double f = diffInten[ch][n];

		sumI += f;

		sumI2 += f*f;
	}

	sumI /= double(nSample-40)*double(numberChannel);
	sumI2 /= double(nSample-40)*double(numberChannel);
	sumI2 -= sumI*sumI;

	thH[chan] = sumI+factor1*sqrt(sumI2);
	thL[chan] = sumI-factor2*sqrt(sumI2);
}


void segment::findCandidate(int chan)
{
	int n;

	tmpC1.nC = tmpC2.nC = 0;

	for(n=1; n<(nSample-2); n++)
	{
		if ( (diffInten[chan][n]>diffInten[chan][n-1]) && (diffInten[chan][n]>diffInten[chan][n+1]) && (diffInten[chan][n]>0) )
		{
			tmpC1.pos[tmpC1.nC]=n;
			tmpC1.cMark[tmpC1.nC]=1;

			tmpC1.nC++;

			if (diffInten[chan][n]>thH[chan])
			{
				tmpC2.pos[tmpC2.nC]=n;
				tmpC2.cMark[tmpC2.nC]=1;

				tmpC2.nC++;
			}
		}

		else if ( (diffInten[chan][n]<diffInten[chan][n-1]) && (diffInten[chan][n]<diffInten[chan][n+1])  && (diffInten[chan][n]<0)  )
		{
			if(tmpC1.nC==0)
			{
				int p=findLocalPeak(diffInten[chan], 1, n-1, -1000);
				tmpC1.pos[0]=p;
				tmpC1.cMark[0]=1;
				tmpC1.nC=1;
			}

			tmpC1.pos[tmpC1.nC]=n;
			tmpC1.cMark[tmpC1.nC]=-1;

			tmpC1.nC++;

			if (diffInten[chan][n]<thL[chan])
			{
				if(tmpC2.nC==0)
				{
					int p=findLocalPeak(diffInten[chan], 1, n-1, -1000);
					tmpC2.pos[0]=p;
					tmpC2.cMark[0]=1;
					tmpC2.nC=1;
				}

				tmpC2.pos[tmpC2.nC]=n;
				tmpC2.cMark[tmpC2.nC]=-1;

				tmpC2.nC++;
			}
		}
	}

	if(tmpC1.cMark[tmpC1.nC-1]==1)
	{
		int d=tmpC1.pos[tmpC1.nC-1];
		int p=findLocalValley(diffInten[chan], d+1, nSample-3, 1000);

		if (p==(d+1)) p=nSample-1;

		tmpC1.pos[tmpC1.nC] = p;
		tmpC1.cMark[tmpC1.nC] = -1;
		tmpC1.nC++;
	}

	if(tmpC2.cMark[tmpC2.nC-1]==1)
	{
		int d=tmpC2.pos[tmpC2.nC-1];
		int p=findLocalValley(diffInten[chan], d+1, nSample-3, 1000);

		if (p==(d+1)) p=nSample-1;

		tmpC2.pos[tmpC2.nC] = p;
		tmpC2.cMark[tmpC2.nC] = -1;
		tmpC2.nC++;
	}
}

void segment::chanSegment(int chan)
{
	int n;

	findCandidate(chan);

	oS1[chan].nC = oS2[chan].nC = fS1[chan].nC = fS2[chan].nC = 0;

	for(n=0; n<(tmpC1.nC); n++)
	{
		if (tmpC1.cMark[n]>0)
		{
			oS1[chan].pos[oS1[chan].nC] = tmpC1.pos[n];
			oS1[chan].cPos[oS1[chan].nC] = tmpC1.pos[n+1];
			oS1[chan].cMark[oS1[chan].nC] = tmpC1.cMark[n+1];

			oS1[chan].nC++;
		}
		else
		{
			fS1[chan].pos[fS1[chan].nC] = tmpC1.pos[n];
			fS1[chan].cPos[fS1[chan].nC] = tmpC1.pos[n-1];
			fS1[chan].cMark[fS1[chan].nC] = tmpC1.cMark[n-1];

			fS1[chan].nC++;
		}
	}

	for(n=0; n<(tmpC2.nC); n++)
	{
		if (tmpC2.cMark[n]>0)
		{
			oS2[chan].pos[oS2[chan].nC] = tmpC2.pos[n];
			oS2[chan].cPos[oS2[chan].nC] = tmpC2.pos[n+1];
			oS2[chan].cMark[oS2[chan].nC] = tmpC2.cMark[n+1];

			oS2[chan].nC++;
		}
		else
		{
			fS2[chan].pos[fS2[chan].nC] = tmpC2.pos[n];
			fS2[chan].cPos[fS2[chan].nC] = tmpC2.pos[n-1];
			fS2[chan].cMark[fS2[chan].nC] = tmpC2.cMark[n-1];

			fS2[chan].nC++;
		}
	}
}

void segment::matchingRange(chanCandidate ofPos, int n, int label, double &sPos1, double &ePos1, double &sPos2, double &ePos2)
{
	if (label>0)
	{
		if (n==0) sPos1 = double(ofPos.pos[n])/2;
		else sPos1 = double(ofPos.pos[n] + ofPos.pos[n-1])/2;

		if (sPos1<(ofPos.pos[n]-THETA_DIS)) sPos1 = ofPos.pos[n]-THETA_DIS;

		sPos2 = ePos1 = double(ofPos.pos[n]+ofPos.cPos[n])/2;

		if (ePos1>(ofPos.pos[n]+THETA_DIS)) ePos1 = ofPos.pos[n]+THETA_DIS;

		if (sPos2<(ofPos.cPos[n]-THETA_DIS)) sPos2 = ofPos.cPos[n]-THETA_DIS;

		if (n==(ofPos.nC-1)) ePos2 = double(ofPos.cPos[n]+nSample)/2;
		else ePos2 = double(ofPos.cPos[n] + ofPos.cPos[n+1])/2;

		if (ePos2>(ofPos.cPos[n]+THETA_DIS)) ePos2 = ofPos.cPos[n]+THETA_DIS;
	}
	else
	{
		if (n==0) sPos2 = double(ofPos.cPos[n])/2;
		else sPos2 = double(ofPos.cPos[n] + ofPos.cPos[n-1])/2;

		if (sPos2<(ofPos.cPos[n]-THETA_DIS)) sPos2 = ofPos.cPos[n]-THETA_DIS;

		sPos1 = ePos2 = double(ofPos.pos[n]+ofPos.cPos[n])/2;

		if (ePos2>(ofPos.cPos[n]+THETA_DIS)) ePos2 = ofPos.cPos[n]+THETA_DIS;

		if (sPos1<(ofPos.pos[n]-THETA_DIS)) sPos1 = ofPos.pos[n]-THETA_DIS;

		if (n==(ofPos.nC-1)) ePos1 = double(ofPos.pos[n]+nSample)/2;
		else ePos1 = double(ofPos.pos[n] + ofPos.pos[n+1])/2;

		if (ePos1>(ofPos.pos[n]+THETA_DIS)) ePos1 = ofPos.pos[n]+THETA_DIS;
	}
}

int segment::matchingCross(chanCandidate *ofPos, int chan, int n, int k, int label)
{
    int sPos, ePos;

	if (label>0)
	{
		sPos = ofPos[chan].pos[n]; if (sPos<ofPos[chan-1].pos[k]) sPos = ofPos[chan-1].pos[k];
		ePos = ofPos[chan].cPos[n]; if (ePos>ofPos[chan-1].cPos[k]) ePos = ofPos[chan-1].cPos[k];
	}
	else
	{
		sPos = ofPos[chan].cPos[n]; if (sPos<ofPos[chan-1].cPos[k]) sPos = ofPos[chan-1].cPos[k];
		ePos = ofPos[chan].pos[n]; if (ePos>ofPos[chan-1].pos[k]) ePos = ofPos[chan-1].pos[k];
	}

	if (sPos>=ePos) return(0);

	int frame;
	double count=0, v1=0, v2=0;
	for(frame=ceil(double(sPos+1)/double(factor)-1); frame<=floor(double(ePos+1)/double(factor)-1); frame++)
	{
		v1 += cross[chan-1][frame];
		v2 += crossEv[chan-1][frame];
		count ++;
	}

	if ( ((v1/count)>THETA_CROSS2) || ((v2/count)>THETA_CROSS2) ) return(1);

	double sum12=0, sum1=0, sum2=0, m1=0, m2=0;
	for(int m=sPos; m<=ePos; m++)
	{
		double f1 = inten[chan-1][m];
		double f2 = inten[chan][m];

		m1 += f1;
		m2 += f2;
	}

	m1 /= double(ePos-sPos+1);
	m2 /= double(ePos-sPos+1);

	for(int m=sPos; m<=ePos; m++)
	{
		double f1 = inten[chan-1][m]-m1;
		double f2 = inten[chan][m]-m2;

		sum12 += f1*f2;
		sum1 += f1*f1;
		sum2 += f2*f2;
	}

	sum12 /=sqrt(sum1*sum2);

	if (sum12>THETA_CROSS) return(1);
	else return(0);
}

void segment::candToFront(chanCandidate *ofPos, front *Front, int label)
{
	int m, n, k, mdis;
	int chan, *mark, *oldMark;
	double sPos1, ePos1, sPos2, ePos2;

	mark = new int[nSample];
	oldMark = new int[nSample];

	Front->nF=0;

	for(chan=0; chan<numberChannel; chan++)
	{
		for(n=0; n<nSample; n++)
			oldMark[n]=mark[n];

		for(n=0; n<ofPos[chan].nC; n++)
		{
			m=-1; mdis=THETA_DIS;

		    if (chan>0)
			{
				matchingRange(ofPos[chan], n, label, sPos1, ePos1, sPos2, ePos2);

				for(k=0; k<ofPos[chan-1].nC; k++)
				{
					if ( (ofPos[chan-1].pos[k]>sPos1) && (ofPos[chan-1].pos[k]<ePos1) &&
						 (ofPos[chan-1].cPos[k]>sPos2) && (ofPos[chan-1].cPos[k]<ePos2) &&
						 (ofPos[chan-1].cMark[k]==ofPos[chan].cMark[n]) )
					{
						m=k; break;
					}
					else if( (ofPos[chan-1].pos[k]>sPos1) && (ofPos[chan-1].pos[k]<ePos1) &&
							 (abs(ofPos[chan-1].pos[k] - ofPos[chan].pos[n])<mdis) )
					{
						if (matchingCross(ofPos, chan, n, k, label)>0)
						{
							m = k;
							mdis = abs(ofPos[chan-1].pos[k] - ofPos[chan].pos[n]);
						}
					}
				}
			}

			if (m>=0) m=oldMark[m];
			else
			{
				m=Front->nF;

				Front->nF++;
				Front->sChan[m]=chan;
			}

			mark[n]=m;

			Front->pos[chan][m] = ofPos[chan].pos[n];
			Front->cPos[chan][m] = ofPos[chan].cPos[n];
			Front->cMark[chan][m] = ofPos[chan].cMark[n];
			Front->eChan[m] = chan;
		}
	}

	delete [] oldMark;
	delete [] mark;
}

void segment::candToFront2(chanCandidate *ofPos, front *Front, int label)
{
	int m, n, k, mdis;
	int chan, *mark, *oldMark;
	double sPos1, ePos1, sPos2, ePos2;

	mark = new int[nSample];
	oldMark = new int[nSample];

	Front->nF=0;

	for(chan=0; chan<numberChannel; chan++)
	{
		for(n=0; n<nSample; n++)
			oldMark[n]=mark[n];

		for(n=0; n<ofPos[chan].nC; n++)
		{
			m=-1; mdis=THETA_DIS;

		    if (chan>0)
			{
				matchingRange(ofPos[chan], n, label, sPos1, ePos1, sPos2, ePos2);

				for(k=0; k<ofPos[chan-1].nC; k++)
				{
					if ( (ofPos[chan-1].pos[k]>sPos1) && (ofPos[chan-1].pos[k]<ePos1) &&
						 (ofPos[chan-1].cPos[k]>sPos2) && (ofPos[chan-1].cPos[k]<ePos2) &&
						 (ofPos[chan-1].cMark[k]==ofPos[chan].cMark[n]) )
					{
						m=k; break;
					}
				}
			}

			if (m>=0) m=oldMark[m];
			else
			{
				m=Front->nF;

				Front->nF++;
				Front->sChan[m]=chan;
			}

			mark[n]=m;

			Front->pos[chan][m] = ofPos[chan].pos[n];
			Front->cPos[chan][m] = ofPos[chan].cPos[n];
			Front->cMark[chan][m] = ofPos[chan].cMark[n];
			Front->eChan[m] = chan;
		}
	}

	delete [] oldMark;
	delete [] mark;
}

void segment::oFMatch(front *oFront, front *oFront2, front *fFront2)
{
	int n, chan, sChan, eChan, mark[MAX_CHANNEL];

	for(n=0; n<oFront->nF; n++)
	{
		sChan=oFront->sChan[n];
		eChan=oFront->eChan[n];

		for(chan=sChan; chan<=eChan; chan++)
			mark[chan]=1;

		int overlap=1;
		int lChan=eChan-sChan+1;
		while ( (lChan>0) && (overlap>0))
		{
			overlap = 0;
			int index=-1;

            for(int m=0; m<oFront2->nF; m++)
			{
				int nChan=0;
				for(chan=sChan; chan<=eChan; chan++)
				{
					if ( (mark[chan]>0) && (oFront->cPos[chan][n]==oFront2->pos[chan][m]) &&
						 (chan>=oFront2->sChan[m]) && (chan<=oFront2->eChan[m]) )
						nChan++;
				}

				if (nChan>overlap)
				{
					overlap=nChan;
					index=m;
				}
			}

			for(int m=0; m<fFront2->nF; m++)
			{
				int nChan=0;
				for(chan=sChan; chan<=eChan; chan++)
				{
					if ( (mark[chan]>0) && (oFront->cPos[chan][n]==fFront2->pos[chan][m]) &&
						 (chan>=fFront2->sChan[m]) && (chan<=fFront2->eChan[m]) )
						nChan++;
				}

				if (nChan>overlap)
				{
					overlap=nChan;
					index=m+oFront2->nF;
				}
			}

			if ( (index>=0) && (index<oFront2->nF))
			{
				for(chan=sChan; chan<=eChan; chan++)
				{
					if ((mark[chan]>0) && (oFront2->sChan[index]<=chan) && (oFront2->eChan[index]>=chan) )
					{
						mark[chan]=0;
						oFront->cPos[chan][n]=oFront2->pos[chan][index];
						oFront->cMark[chan][n]=1;
					}
				}
			}
			else if (index>=oFront2->nF)
			{
				index -= oFront2->nF;

				for(chan=sChan; chan<=eChan; chan++)
				{
					if ((mark[chan]>0) && (fFront2->sChan[index]<=chan) && (fFront2->eChan[index]>=chan) )
					{
						mark[chan]=0;
						oFront->cPos[chan][n]=fFront2->pos[chan][index];
						oFront->cMark[chan][n]=-1;
					}
				}
			}

			lChan=0;
			for(chan=sChan; chan<=eChan; chan++)
			{
				if(mark[chan]>0) lChan++;
			}
		}
	}
}

void segment::oFMatch2(front *oFront, front *oFront2)
{
	int n, chan, sChan, eChan, mark[MAX_CHANNEL];

	for(n=0; n<oFront->nF; n++)
	{
		sChan=oFront->sChan[n];
		eChan=oFront->eChan[n];

		for(chan=sChan; chan<=eChan; chan++)
			mark[chan]=1;

		int overlap=1;
		int lChan=eChan-sChan+1;
		while ( (lChan>0) && (overlap>0))
		{
			overlap = 0;
			int index=-1;

            for(int m=0; m<oFront2->nF; m++)
			{
				int nChan=0;
				for(chan=sChan; chan<=eChan; chan++)
				{
					if ( (mark[chan]>0) && (oFront->pos[chan][n]==oFront2->pos[chan][m]) &&
						 (chan>=oFront2->sChan[m]) && (chan<=oFront2->eChan[m]) )
						nChan++;
				}

				if (nChan>overlap)
				{
					overlap=nChan;
					index=m;
				}
			}

			if (index>=0)
			{
				for(chan=sChan; chan<=eChan; chan++)
				{
					if ((mark[chan]>0) && (oFront2->sChan[index]<=chan) && (oFront2->eChan[index]>=chan) )
					{
						mark[chan]=0;
						oFront->pos[chan][n]=oFront2->pos[chan][index];
					}
				}
			}

			lChan=0;
			for(chan=sChan; chan<=eChan; chan++)
			{
				if(mark[chan]>0) lChan++;
			}
		}
	}
}

void segment::frontToCandidate(front *oFront, front *fFront)
{
	int *index[MAX_CHANNEL], chan, n, pos;

	for(chan=0; chan<numberChannel; chan++)
	{
		index[chan] = new int[nSample+10];

		for(pos=0; pos<nSample; pos++)
			index[chan][pos]=0;
	}

	for(n=0; n<oFront->nF; n++)
		for(chan=oFront->sChan[n]; chan<=oFront->eChan[n]; chan++)
		{
			pos=oFront->pos[chan][n];
			index[chan][pos]=1;

			pos=oFront->cPos[chan][n];
			index[chan][pos]=oFront->cMark[chan][n];
		}

	for(n=0; n<fFront->nF; n++)
		for(chan=fFront->sChan[n]; chan<=fFront->eChan[n]; chan++)
		{
			pos=fFront->pos[chan][n];
			index[chan][pos]=-1;

			pos=fFront->cPos[chan][n];
			index[chan][pos]=fFront->cMark[chan][n];
		}

	for(chan=0; chan<numberChannel; chan++)
	{
		tmpC2.nC = 1;
		tmpC2.pos[0] = 0; tmpC2.cMark[0] = 1;

		for(pos=1; pos<(nSample-1); pos++)
		{
			if ( index[chan][pos]>0)
			{
				tmpC2.pos[tmpC2.nC]=pos;
				tmpC2.cMark[tmpC2.nC]=1;

				tmpC2.nC++;
			}

			else if ( index[chan][pos]<0 )
			{
				tmpC2.pos[tmpC2.nC]=pos;
				tmpC2.cMark[tmpC2.nC]=-1;

				tmpC2.nC++;
			}
		}

		tmpC2.pos[tmpC2.nC] = nSample; tmpC2.cMark[tmpC2.nC] = -1; tmpC2.nC++;

		oS2[chan].nC = fS2[chan].nC = 0;

		for(n=2; n<tmpC2.nC-1; n++)
		{
			if (tmpC2.cMark[n]>0)
			{
				oS2[chan].pos[oS2[chan].nC] = tmpC2.pos[n];
				oS2[chan].cPos[oS2[chan].nC] = tmpC2.pos[n+1];
				oS2[chan].cMark[oS2[chan].nC] = tmpC2.cMark[n+1];

				oS2[chan].nC++;
			}
			else
			{
				fS2[chan].pos[fS2[chan].nC] = tmpC2.pos[n];
				fS2[chan].cPos[fS2[chan].nC] = tmpC2.pos[n-1];
				fS2[chan].cMark[fS2[chan].nC] = tmpC2.cMark[n-1];

				fS2[chan].nC++;
			}
		}
	}

	for(chan=0; chan<numberChannel; chan++)
		delete [] index[chan];
}

void segment::adjustPos(front *oFront, front *fFront)
{
	int n, chan, m;
	frontToCandidate(oFront, fFront);

	for(n=0; n<oFront->nF; n++)
		for(chan=oFront->sChan[n]; chan<=oFront->eChan[n]; chan++)
		{
			for(m=0; m<oS2[chan].nC; m++)
			{
				if (oFront->pos[chan][n]==oS2[chan].pos[m])
				{
					oFront->cPos[chan][n]=oS2[chan].cPos[m];
					oFront->cMark[chan][n]=oS2[chan].cMark[m];
				}
			}
		}

    for(n=0; n<fFront->nF; n++)
		for(chan=fFront->sChan[n]; chan<=fFront->eChan[n]; chan++)
		{
			for(m=0; m<fS2[chan].nC; m++)
			{
				if (fFront->pos[chan][n]==fS2[chan].pos[m])
				{
					fFront->cPos[chan][n]=fS2[chan].cPos[m];
					fFront->cMark[chan][n]=fS2[chan].cMark[m];
				}
			}
		}
}

int segment::newPos(int pos, int cPos, chanCandidate ofPos)
{
	int m, index, mDis;
	double sPos, ePos;

	sPos=pos-THETA_DIS;
	ePos=pos+THETA_DIS;

	if (pos<cPos)
	{
		if (ePos>(double(pos+cPos)/2)) ePos=double(pos+cPos)/2;
	}

	else
	{
		if (sPos>(double(pos+cPos)/2)) sPos=double(pos+cPos)/2;
	}

	index=-1;
	mDis=THETA_DIS;
	for(m=0; m<ofPos.nC; m++)
	{
		if ( (ofPos.pos[m]>sPos) && (ofPos.pos[m]<ePos) && (abs(ofPos.pos[m]-pos)<mDis) )
		{
			index=m;
			mDis=abs(ofPos.pos[m]-pos);
		}
	}

	return(index);
}

void segment::accuratePos(front *Front, chanCandidate *ofPos, chanCandidate *oPos2, chanCandidate *fPos2)
{
	int n=0, chan;

	while (n<Front->nF)
	{
		int markChan[MAX_CHANNEL], index;
		for(chan=Front->sChan[n]; chan<=Front->eChan[n]; chan++)
		{
			index=newPos(Front->pos[chan][n], Front->cPos[chan][n], ofPos[chan]);

			if (index>=0)
			{
				Front->pos[chan][n]=ofPos[chan].pos[index];
				markChan[chan]=1;
			}
			else markChan[chan]=0;

			if(Front->cMark[chan][n]>0)
			{
				index=newPos(Front->cPos[chan][n], Front->pos[chan][n], oPos2[chan]);

				if(index>=0) Front->cPos[chan][n]=oPos2[chan].pos[index];
			}
			else
			{
				index=newPos(Front->cPos[chan][n], Front->pos[chan][n], fPos2[chan]);

				if(index>=0) Front->cPos[chan][n]=fPos2[chan].pos[index];
			}
		}

		chan=Front->sChan[n];
		while ( (chan<=Front->eChan[n]) && (markChan[chan]==0) )
		{
			chan++;
			Front->sChan[n]=chan;
		}

		chan=Front->eChan[n];
		while ( (chan>=Front->sChan[n]) & (markChan[chan]==0) )
		{
			chan--;
			Front->eChan[n]=chan;
		}

		if(Front->sChan[n]>Front->eChan[n])
		{
			for(int m=n+1; m<Front->nF; m++)
			{
				Front->sChan[m-1]=Front->sChan[m];
				Front->eChan[m-1]=Front->eChan[m];

				for(chan=Front->sChan[m]; chan<=Front->eChan[m]; chan++)
				{
					Front->pos[chan][m-1]=Front->pos[chan][m];
					Front->cPos[chan][m-1]=Front->cPos[chan][m];
					Front->cMark[chan][m-1]=Front->cMark[chan][m];
				}
			}

			Front->nF--;
			n--;
		}

		else
		{
			for(chan=Front->sChan[n]; chan<=Front->eChan[n]; chan++)
			{
				if(Front->pos[chan][n]<Front->cPos[chan][n])
				{
					index=newPos(Front->pos[chan][n], Front->cPos[chan][n], oPos2[chan]);

					if (index>=0) Front->pos[chan][n]=oPos2[chan].pos[index];
				}

				else
				{
					index=newPos(Front->pos[chan][n], Front->cPos[chan][n], fPos2[chan]);

					if (index>=0) Front->pos[chan][n]=fPos2[chan].pos[index];
				}
			}
		}

		n++;
	}
}

void segment::accuratePos2(front *Front, chanCandidate *ofPos, chanCandidate *oPos2, chanCandidate *fPos2)
{
	int n=0, chan;

	while (n<Front->nF)
	{
		int markChan[MAX_CHANNEL], index;
		for(chan=Front->sChan[n]; chan<=Front->eChan[n]; chan++)
		{
			index=newPos(Front->pos[chan][n], Front->cPos[chan][n], ofPos[chan]);

			if (index>=0) Front->pos[chan][n]=ofPos[chan].pos[index];


			if(Front->cMark[chan][n]>0)
			{
				index=newPos(Front->cPos[chan][n], Front->pos[chan][n], oPos2[chan]);

				if(index>=0) Front->cPos[chan][n]=oPos2[chan].pos[index];
			}
			else
			{
				index=newPos(Front->cPos[chan][n], Front->pos[chan][n], fPos2[chan]);

				if(index>=0) Front->cPos[chan][n]=fPos2[chan].pos[index];
			}
		}

		n++;
	}
}

void segment::removeOFSet(front *oFront, chanCandidate *oPos, int label)
{
	int n, m, *mark, chan;

	mark=new int[nSample];

	for(n=0; n<oFront->nF; n++)
		for(chan=oFront->sChan[n]; chan<=oFront->eChan[n]; chan++)
		{
			for(m=0; m<oPos[chan].nC; m++)
			{
				mark[m]=1;

				if (label>0)
				{
					if ( (oPos[chan].pos[m]>=oFront->pos[chan][n]) && (oPos[chan].pos[m]<=oFront->cPos[chan][n]) ) mark[m]=0;
				}
				else
				{
					if ( (oPos[chan].pos[m]>=oFront->cPos[chan][n]) && (oPos[chan].pos[m]<=oFront->pos[chan][n]) ) mark[m]=0;
				}
			}

			int nC=0;
			for(m=0; m<oPos[chan].nC; m++)
			{
				if(mark[m]>0)
				{
					oPos[chan].pos[nC]=oPos[chan].pos[m];
					oPos[chan].cPos[nC]=oPos[chan].cPos[m];
					oPos[chan].cMark[nC]=oPos[chan].cMark[m];

					nC++;
				}
			}

			oPos[chan].nC=nC;
		}

	delete [] mark;
}

void segment::removeOFSet2(front *oFront, chanCandidate *oPos, int label)
{
	int n, m, *mark, chan;

	mark=new int[nSample];

	for(n=0; n<oFront->nF; n++)
		for(chan=oFront->sChan[n]; chan<=oFront->eChan[n]; chan++)
		{
			for(m=0; m<oPos[chan].nC; m++)
			{
				mark[m]=1;

				if (label>0)
				{
					if ( (oPos[chan].pos[m]>=oFront->pos[chan][n]) && (oPos[chan].pos[m]<oFront->cPos[chan][n]) ) mark[m]=0;
				}
				else
				{
					if ( (oPos[chan].pos[m]>oFront->cPos[chan][n]) && (oPos[chan].pos[m]<=oFront->pos[chan][n]) ) mark[m]=0;
				}
			}

			int nC=0;
			for(m=0; m<oPos[chan].nC; m++)
			{
				if(mark[m]>0)
				{
					oPos[chan].pos[nC]=oPos[chan].pos[m];
					oPos[chan].cPos[nC]=oPos[chan].cPos[m];
					oPos[chan].cMark[nC]=oPos[chan].cMark[m];

					nC++;
				}
			}

			oPos[chan].nC=nC;
		}

	delete [] mark;
}

void segment::enlongFront(front *oFront, front *oFront2, chanCandidate *oPos)
{
	int n, m, sChan, eChan, chan;
	for(n=0; n<oFront->nF; n++)
	{
		oFront->extendFront1[n]=0;
		oFront->extendFront2[n]=0;

		sChan=oFront->sChan[n];
		eChan=oFront->eChan[n];
		for(m=0; m<oFront2->nF; m++)
		{
			if ( (sChan>oFront2->sChan[m]) && (sChan<=oFront2->eChan[m]) &&
				 (oFront->pos[sChan][n]==oFront2->pos[sChan][m]) &&
				 (oFront->cPos[sChan][n]==oFront2->cPos[sChan][m]) )
			{
				int judge=1;
				chan=sChan-1;

				while ( (judge>0) && (chan>=oFront2->sChan[m]) )
				{
					judge=0;
					for(int k=0; k<oPos[chan].nC; k++)
					{
						if (oPos[chan].pos[k]==oFront2->pos[chan][m])
						{
							judge=1;

							oFront->extendFront1[n]=1;

							oFront->pos[chan][n]=oFront2->pos[chan][m];
							oFront->cPos[chan][n]=oFront2->cPos[chan][m];
							oFront->cMark[chan][n]=oFront2->cMark[chan][m];

							oFront->sChan[n]=chan;
							sChan=chan;

							break;
						}
					}

					chan--;

				}
			}

			if ( (eChan>=oFront2->sChan[m]) && (eChan<oFront2->eChan[m]) &&
				 (oFront->pos[eChan][n]==oFront2->pos[eChan][m]) &&
				 (oFront->cPos[eChan][n]==oFront2->cPos[eChan][m]) )
			{
				int judge=1;
				chan=eChan+1;

				while ( (judge>0) && (chan<=oFront2->eChan[m]) )
				{
					judge=0;
					for(int k=0; k<oPos[chan].nC; k++)
					{
						if (oPos[chan].pos[k]==oFront2->pos[chan][m])
						{
							judge=1;

							oFront->extendFront2[n]=1;

							oFront->pos[chan][n]=oFront2->pos[chan][m];
							oFront->cPos[chan][n]=oFront2->cPos[chan][m];
							oFront->cMark[chan][n]=oFront2->cMark[chan][m];

							oFront->eChan[n]=chan;
							eChan=chan;

							break;
						}
					}

					chan++;

				}
			}
		}
	}
}

void segment::mergeFront(front *oFront, int label)
{
	int m, n, chan, index;
	n=0;
	while (n<(oFront->nF-1))
	{
		index=n;
		for(m=n+1; m<oFront->nF; m++)
		{
			for(chan=oFront->sChan[n]; chan<=oFront->eChan[n]; chan++)
			{
				if( (chan<=oFront->eChan[m]) && (chan>=oFront->sChan[m]) &&
					(oFront->pos[chan][n]==oFront->pos[chan][m]) &&
					( ( oFront->extendFront1[n]>0) || (oFront->extendFront2[m]>0) ) )
				{
					index=m;
					break;
				}
			}

			if(index>n) break;
		}

		if (index>n)
		{
			for(chan=oFront->sChan[index]; chan<oFront->sChan[n]; chan++)
			{
				oFront->pos[chan][n]=oFront->pos[chan][index];
				oFront->cPos[chan][n]=oFront->cPos[chan][index];
				oFront->cMark[chan][n]=oFront->cMark[chan][index];
			}

			for(chan=oFront->eChan[n]+1; chan<=oFront->eChan[index]; chan++)
			{
				oFront->pos[chan][n]=oFront->pos[chan][index];
				oFront->cPos[chan][n]=oFront->cPos[chan][index];
				oFront->cMark[chan][n]=oFront->cMark[chan][index];
			}

			if(oFront->sChan[index]<oFront->sChan[n]) oFront->sChan[n]=oFront->sChan[index];
			if(oFront->eChan[index]>oFront->eChan[n]) oFront->eChan[n]=oFront->eChan[index];

			for(m=index+1; m<oFront->nF; m++)
			{
				oFront->sChan[m-1]=oFront->sChan[m];
				oFront->eChan[m-1]=oFront->eChan[m];

				for(chan=oFront->sChan[m]; chan<=oFront->eChan[m]; chan++)
				{
					oFront->pos[chan][m-1]=oFront->pos[chan][m];
					oFront->cPos[chan][m-1]=oFront->cPos[chan][m];
					oFront->cMark[chan][m-1]=oFront->cMark[chan][m];
				}
			}

			oFront->nF--;
			n--;
		}
		n++;
	}
}

void segment::assignFront(front *tF, front *sF, int m)
{
	tF->nF=m+sF->nF;
	for(int n=m; n<tF->nF; n++)
	{
		tF->sChan[n]=sF->sChan[n-m];
		tF->eChan[n]=sF->eChan[n-m];

		for(int chan=tF->sChan[n]; chan<=tF->eChan[n]; chan++)
		{
			tF->pos[chan][n]=sF->pos[chan][n-m];
			tF->cPos[chan][n]=sF->cPos[chan][n-m];
			tF->cMark[chan][n]=sF->cMark[chan][n-m];
		}
	}
}

void segment::ofFront(mScaleInten *scaleInten, segScale s)
{
	int n, chan, scale, shift;

	for(scale=0; scale<s.scale; scale++)
	{
		shift = s.timeScale[scale] * numberChannel;

		if (s.freqScale[scale]==6) theta_cross=THETA_CROSS;
		else theta_cross=THETA_CROSS3;

		for(chan=0; chan<numberChannel; chan++)
			for(n=0; n<(nSample); n++)
				inten[chan][n] = scaleInten->ampScale[shift+chan][n];

		freqSmooth(s.freqScale[scale], inten);

		for(chan=0; chan<numberChannel; chan++)
			for(n=0; n<(nSample-1); n++)
				diffInten[chan][n] = inten[chan][n+1]-inten[chan][n];

		determineTheta(0, 1, 1);
		for(chan=0; chan<numberChannel; chan++)
		{
			thH[chan] = thH[0];
			thL[chan] = thL[0];
		}

		for(chan=0; chan<numberChannel; chan++)
			chanSegment(chan);

		if (scale==0)
		{
			candToFront(oS1, &oF1, 1);
			candToFront(oS2, &oF, 1);
			candToFront(fS1, &fF1, -1);
			candToFront(fS2, &fF, -1);

			oFMatch(&oF, &oF1, &fF1);
			oFMatch(&fF, &oF1, &fF1);
		}

	    else
		{
			if (scale<4)
			{
				accuratePos(&oF, oS2, oS1, fS1);
				accuratePos(&fF, fS2, oS1, fS1);

				candToFront(oS1, &oF1, 1);
				candToFront(fS1, &fF1, -1);

				oFMatch(&oF, &oF1, &fF1);
				oFMatch(&fF, &oF1, &fF1);

				candToFront2(oS2, &oF2, 1);
				candToFront2(fS2, &fF2, -1);

				removeOFSet(&oF, oS2, 1);
				removeOFSet(&fF, oS2, -1);
				removeOFSet(&oF, fS2, 1);
				removeOFSet(&fF, fS2, -1);

				enlongFront(&oF, &oF2, oS2);
				mergeFront(&oF, 1);
				enlongFront(&fF, &fF2, fS2);
				mergeFront(&fF, -1);

				removeOFSet(&oF, oS2, 1);
				removeOFSet(&fF, oS2, -1);
				removeOFSet(&oF, fS2, 1);
				removeOFSet(&fF, fS2, -1);

				candToFront(oS2, &oF2, 1);
				candToFront(fS2, &fF2, -1);

				assignFront(&oF, &oF2, oF.nF);
				assignFront(&fF, &fF2, fF.nF);
			}

			else
			{
				accuratePos(&oF, oS2, oS1, fS1);
				accuratePos(&fF, fS2, oS1, fS1);

				candToFront(oS1, &oF1, 1);
				candToFront(fS1, &fF1, -1);

				oFMatch(&oF, &oF1, &fF1);
				oFMatch(&fF, &oF1, &fF1);

				candToFront2(oS2, &oF2, 1);
				candToFront2(fS2, &fF2, -1);

				removeOFSet(&oF, oS2, 1);
				removeOFSet(&fF, fS2, -1);

				enlongFront(&oF, &oF2, oS2);
				mergeFront(&oF, 1);
				enlongFront(&fF, &fF2, fS2);
				mergeFront(&fF, -1);

				removeOFSet(&oF, oS2, 1);
				removeOFSet(&fF, oS2, -1);
				removeOFSet(&oF, fS2, -1);
				removeOFSet(&fF, fS2, -1);

				candToFront(oS2, &oF2, 1);
				candToFront(fS2, &fF2, -1);

				assignFront(&oF, &oF2, oF.nF);
				assignFront(&fF, &fF2, fF.nF);
			}
		}

		n=0;

		while(n<oF.nF)
		{
			if ((oF.eChan[n]-oF.sChan[n])<5)
			{
				for(int m=n+1; m<oF.nF; m++)
				{
					oF.sChan[m-1]=oF.sChan[m];
					oF.eChan[m-1]=oF.eChan[m];

					for(chan=oF.sChan[m]; chan<=oF.eChan[m]; chan++)
					{
						oF.pos[chan][m-1]=oF.pos[chan][m];
						oF.cPos[chan][m-1]=oF.cPos[chan][m];
						oF.cMark[chan][m-1]=oF.cMark[chan][m];
					}
				}

				oF.nF--;
				n--;
			}

			n++;
		}

		n=0;
		while(n<fF.nF)
		{
			if ((fF.eChan[n]-fF.sChan[n])<5)
			{
				for(int m=n+1; m<fF.nF; m++)
				{
					fF.sChan[m-1]=fF.sChan[m];
					fF.eChan[m-1]=fF.eChan[m];

					for(chan=fF.sChan[m]; chan<=fF.eChan[m]; chan++)
					{
						fF.pos[chan][m-1]=fF.pos[chan][m];
						fF.cPos[chan][m-1]=fF.cPos[chan][m];
						fF.cMark[chan][m-1]=fF.cMark[chan][m];
					}
				}

				fF.nF--;
				n--;
			}

			n++;
		}
	}

	adjustPos(&oF, &fF);
}

void segment::backGroundSet(front *oF, front *fF)
{
	int *index[MAX_CHANNEL], chan, n, pos, k;

	for(chan=0; chan<numberChannel; chan++)
	{
		index[chan] = new int[nSample];

		for(pos=0; pos<nSample; pos++)
			index[chan][pos]=0;
	}

	for(n=0; n<oF->nF; n++)
		for(chan=oF->sChan[n]; chan<=oF->eChan[n]; chan++)
		{
			for(k=oF->pos[chan][n]; k<=oF->cPos[chan][n]; k++)
				index[chan][k]=1;
		}

	for(n=0; n<fF->nF; n++)
		for(chan=fF->sChan[n]; chan<=fF->eChan[n]; chan++)
		{
			for(k=fF->cPos[chan][n]; k<=fF->pos[chan][n]; k++)
				index[chan][k]=1;
		}

	for(chan=0; chan<numberChannel; chan++)
	{
		tmpC2.nC = 1;
		tmpC2.pos[0] = 0; tmpC2.cMark[0] = 1;

		for(pos=1; pos<(nSample-1); pos++)
		{
			if ( (index[chan][pos-1]>0) && (index[chan][pos]==0) && (index[chan][pos+1]==0) )
			{
				tmpC2.pos[tmpC2.nC]=pos;
				tmpC2.cMark[tmpC2.nC]=1;

				tmpC2.nC++;
			}

			else if ( (index[chan][pos-1]==0) && (index[chan][pos]==0) && (index[chan][pos+1]==1) )
			{
				tmpC2.pos[tmpC2.nC]=pos;
				tmpC2.cMark[tmpC2.nC]=-1;

				tmpC2.nC++;
			}
		}

		tmpC2.pos[tmpC2.nC] = nSample; tmpC2.cMark[tmpC2.nC] = -1; tmpC2.nC++;

		oS2[chan].nC = fS2[chan].nC = 0;

		for(n=2; n<tmpC2.nC-1; n++)
		{
			if (tmpC2.cMark[n]>0)
			{
				oS2[chan].pos[oS2[chan].nC] = tmpC2.pos[n];
				oS2[chan].cPos[oS2[chan].nC] = tmpC2.pos[n+1];
				oS2[chan].cMark[oS2[chan].nC] = tmpC2.cMark[n+1];

				oS2[chan].nC++;
			}
			else
			{
				fS2[chan].pos[fS2[chan].nC] = tmpC2.pos[n];
				fS2[chan].cPos[fS2[chan].nC] = tmpC2.pos[n-1];
				fS2[chan].cMark[fS2[chan].nC] = tmpC2.cMark[n-1];

				fS2[chan].nC++;
			}
		}
	}
}

void segment::frontToSeg()
{
	int n, frame, chan;
	int index, index2;
	int nSeg;

	for(frame=0; frame<(nSample/factor); frame++)
		for(chan=0; chan<numberChannel; chan++)
			seg[chan][frame]=0;

	for(n=0; n<oF.nF; n++)
		for(chan=oF.sChan[n]; chan<=oF.eChan[n]; chan++)
			for(frame=floor((oF.pos[chan][n]+1)/factor); frame<=floor((oF.cPos[chan][n]+1)/factor)-1; frame++)
				seg[chan][frame]=n+1;

	nSeg=oF.nF;
	for(n=0; n<fF.nF; n++)
	{
		index=nSeg+1;
		nSeg++;

		for(chan=fF.sChan[n]; chan<=fF.eChan[n]; chan++)
			for(frame=floor((fF.cPos[chan][n]+1)/factor); frame<=floor((fF.pos[chan][n]+1)/factor)-1; frame++)
			{
				if (seg[chan][frame]==0) seg[chan][frame]=index;

				else if(index!=seg[chan][frame])
				{
					if (index>seg[chan][frame])
					{
						index2=index;
						index=seg[chan][frame];
					}
					else index2=seg[chan][frame];

					int f, ch;

					seg[chan][frame]=index;

					for(f=0; f<(nSample/factor); f++)
						for(ch=0; ch<numberChannel; ch++)
						{
							if(seg[ch][f]==index2) seg[ch][f]=index;
							if(seg[ch][f]>index2) seg[ch][f]--;
						}

					nSeg--;
				}
			}
	}
}
