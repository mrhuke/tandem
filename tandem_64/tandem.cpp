#include "voicedMask.h"
#include "segment.h"

#include <fstream>
#include <vector>
#include <iostream>
#include <iterator>
#include <functional>
#include <algorithm>
#include <numeric>
#include <sstream>
using namespace std;

void initVoicedMask(gammaToneFilterBank *&AudiPery, voicedMask *&TGroup)
{
	featurePara x;
	gammaTonePara y;

	y.lCF = 50;	 y.rCF = 8000;
	y.nChan = 64; y.sf = 20000;
	
	x.gtP = y;
	x.bP1 = 50; x.bP2 = 450; x.bPTs = 20;
	x.theta_p = .5;
	
	AudiPery = new gammaToneFilterBank(x.gtP);
	TGroup = new voicedMask(x);
}

void initOnOff(segment *&TSeg, mScaleInten *&scaleInten, segScale &sScale)
{
	TSeg = new segment(64);

	gammaTonePara y;
        scalePara x;
	y.lCF = 50; y.rCF = 8000;
        y.nChan = 64; y.sf = 20000;
        x.gtP = y;
        x.lP1 = 30;     x.lP2 = 50;
        x.scale = 6;
        x.timeScale1[0] = 4; x.timeScale2[0] = 14;
        x.timeScale1[1] = 6; x.timeScale2[1] = 16;
        x.timeScale1[2] = 8; x.timeScale2[2] = 18;
        x.timeScale1[3] = 10; x.timeScale2[3] = 19;
        x.timeScale1[4] = 12; x.timeScale2[4] = 22;
        x.timeScale1[5] = 14; x.timeScale2[5] = 24;
        x.downRate = 50;
        scaleInten = new mScaleInten(x);
        
	sScale.scale=4;
        sScale.timeScale[0]=0; sScale.freqScale[0]=3;
        sScale.timeScale[1]=2; sScale.freqScale[1]=3;
        sScale.timeScale[2]=5; sScale.freqScale[2]=2;
        sScale.timeScale[3]=5; sScale.freqScale[3]=2;        
}	

double* readInput(char *filename, int *numSample)
{
	ifstream is(filename);
	if (!is.is_open())
		return 0;

	istream_iterator<double> is_it(is);
	istream_iterator<double> eof;
	vector<double> tmp;
	double sumE=0;
	while(is_it!=eof){
		sumE += (*is_it)*(*is_it);
		tmp.push_back(*is_it++);
	}
	sumE = sqrt(sumE/tmp.size());
	double scale = 1000/sumE; // scale to 60-dB loudness level

	cout<<"Processing "<<tmp.size()<<" samples...";
	*numSample = tmp.size();
	double *signal = new double[tmp.size()];
	for (vector<double>::size_type k=0; k!=tmp.size(); ++k)
		signal[k] = scale*tmp[k];
	return signal;

	#ifdef DEBUG
	vector<double> dummy;
	transform(signal, signal+tmp.size(), signal, back_inserter(dummy), std::multiplies<double>());
	cout<<"Average loudness level is now "<<10*log10(accumulate(dummy.begin(), dummy.end(), 0)/dummy.size())<<endl;
	#endif
}


void SaveMasks(char *outfn, voicedMask *TGroup)
{
	char filename[255];
	FILE *fp1, *fp2;

	sprintf(filename, "%s.%d.pitch.dat", outfn, TGroup->numberChannel);
	fp1=fopen(filename, "w");

	sprintf(filename, "%s.%d.simul.dat", outfn, TGroup->numberChannel);
	fp2=fopen(filename, "w");

	fprintf(fp1, "%d %d\n", TGroup->numContour, TGroup->numFrame);

	for(int n=0; n<TGroup->numContour; n++)
	{
		for(int frame=0; frame<TGroup->numFrame; frame++)
		{
			fprintf(fp1, "%d ", TGroup->Pitch[n].value[frame]);
			
			if(TGroup->Pitch[n].value[frame]>0)
			{
				for(int chan=0; chan<TGroup->numberChannel; chan++)
					fprintf(fp2, "%6.4f ", TGroup->Pitch[n].mProb[frame].value[chan]/2+.5);
				fprintf(fp2, "\n");
			}
			else{
				for(int chan=0; chan<TGroup->numberChannel; chan++)
                                        fprintf(fp2, "0 ");
				fprintf(fp2, "\n");
			}
		}
		fprintf(fp1, "\n");
	}
	fclose(fp1);
	fclose(fp2);
}

void SaveMasks(char *outfn, voicedMask *TGroup, vector<vector<vector<int> > > &allMasks)
{
        char filename[255];
        FILE *fp1, *fp2;

        sprintf(filename, "%s.%d.pitchT.dat", outfn, TGroup->numberChannel);
        fp1=fopen(filename, "w");

        sprintf(filename, "%s.%d.simulT.dat", outfn, TGroup->numberChannel);
        fp2=fopen(filename, "w");

        fprintf(fp1, "%d %d\n", TGroup->numContour, TGroup->numFrame);

        for(int n=0; n<TGroup->numContour; n++)
        {
                for(int frame=0; frame<TGroup->numFrame; frame++)
                {
                        fprintf(fp1, "%d ", TGroup->Pitch[n].value[frame]);

                        for(int chan=0; chan<TGroup->numberChannel; chan++)
				fprintf(fp2, "%d ", allMasks[n][frame][chan]>0);
                        fprintf(fp2, "\n");
                }
                fprintf(fp1, "\n");
        }
        fclose(fp1);
        fclose(fp2);
}

void SaveCross(char *fn_cross, char *fn_evCross, voicedMask *TGroup)
{
	char filename1[255], filename2[255];
        FILE *fp1, *fp2;

	sprintf(filename1, "%s.%d", fn_cross, TGroup->numberChannel);
	sprintf(filename2, "%s.%d", fn_evCross, TGroup->numberChannel);
        fp1=fopen(filename1, "w");
        fp2=fopen(filename2, "w");

        for(int frame=0; frame<TGroup->numFrame; frame++){
                for(int chan=0; chan<TGroup->numberChannel; chan++){
			fprintf(fp1, "%6.4f ", TGroup->corrLgm[frame].cross[chan]);
			fprintf(fp2, "%6.4f ", TGroup->corrLgm[frame].crossEv[chan]);
                }
                fprintf(fp1, "\n");
		fprintf(fp2, "\n");
        }

        fclose(fp1);
        fclose(fp2);
}

void SaveEng(char *fn_eng, voicedMask *TGroup)
{
	char filename[255];
        FILE *fp;

	sprintf(filename, "%s.%d", fn_eng, TGroup->numberChannel);
        fp=fopen(filename, "w");

        for(int frame=0; frame<TGroup->numFrame; frame++){
                for(int chan=0; chan<TGroup->numberChannel; chan++){
                        fprintf(fp, "%6.4f ", TGroup->corrLgm[frame].acf[chan][0]);
                }
                fprintf(fp, "\n");
        }
        fclose(fp);
}

void voicedMaskEst(char *filename, gammaToneFilterBank *AudiPery, voicedMask *TGroup)
{
	int sigLength, chan;

	double *Input=readInput(filename, &sigLength);
		
	printf("%d frames\n", sigLength/(TGroup->fs/100));

	TGroup->numFrame=sigLength/(TGroup->fs/100);
	
	TGroup->newPitchPara();

	AudiPery->sigLength = sigLength;
	for(int chan=0; chan<TGroup->numberChannel; chan++)
		AudiPery->filtering(Input, sigLength, chan);
	
	TGroup->computeFeature(AudiPery, sigLength, 1);

	TGroup->dtmPitchMask();

	delete [] Input;
}

void onOffSeg(gammaToneFilterBank *AudiPery, voicedMask *TGroup, 
	      segment *TSeg, mScaleInten *scaleInten, segScale &sScale)
{
	TSeg->nSample = scaleInten->nSample = AudiPery->sigLength/(scaleInten->downRate);

        scaleInten->smooth(AudiPery, AudiPery->sigLength, 1);

        TSeg->newSegPara();

	// copy response and envelope cross correlations
	for(int m=0; m<TGroup->numFrame; m++)
                for(int c=0; c<TGroup->numberChannel; c++)
			TSeg->cross[c][m] = TGroup->corrLgm[m].cross[c];

	for(int m=0; m<TGroup->numFrame; m++)
                for(int c=0; c<TGroup->numberChannel; c++)
			TSeg->crossEv[c][m] = TGroup->corrLgm[m].crossEv[c];
        
        TSeg->ofFront(scaleInten, sScale);

        TSeg->frontToSeg();
	
	#ifdef DEBUG	
	FILE *fp=fopen("onoff.dbg", "w");
        for(int frame=0; frame<(TSeg->nSample/TSeg->factor); frame++){
                for(int chan=0; chan<TSeg->numberChannel; chan++){
			fprintf(fp, "%d ", TSeg->seg[chan][frame]);
		}				
		fprintf(fp, "\n");
	}
        fclose(fp);
	#endif
}

void expandVoicedMask(vector<vector<int> > &mask, voicedMask *TGroup, segment *TSeg)
{
	#ifdef DEBUG
        FILE *fid=fopen("tmask.dbg", "w");
        for(int chan=0; chan<mask.size(); chan++){
                for(int frame=0; frame<mask[chan].size(); frame++){
                        fprintf(fid, "%d ", mask[chan][frame]);
                }
                fprintf(fid, "\n");
        }
        fclose(fid);
        #endif

	// get the voiced-mask active frames
	vector<int> sumM = sum(mask,2);

	// convert segments into a 2-D vector
	vector<vector<int> > seg;
	seg.resize(TGroup->numberChannel);
	for (int c=0; c<TGroup->numberChannel; ++c){
		for (int m=0; m<TGroup->numFrame; ++m){
			seg[c].push_back(TSeg->seg[c][m]);
   		}
        }

	#ifdef DEBUG
        FILE *fp=fopen("onoff.dbg", "w");
        for(int chan=0; chan<seg.size(); chan++){
        	for(int frame=0; frame<seg[chan].size(); frame++){
                        fprintf(fp, "%d ", seg[chan][frame]);
                }
                fprintf(fp, "\n");
        }
        fclose(fp);
        #endif

	int nSeg = max(seg);
	// process each segment
	for (int n=1; n<=nSeg; ++n){
		for (int c=0; c<TGroup->numberChannel; ++c){
			vector<int> t_seg=(seg[c]==n);  // active units at channel c for segment n
			if (accumulate(t_seg.begin(), t_seg.end(),0)){
				// eT: segment energy
				double eT=0;
				for (vector<bool>::size_type k=0; k!=t_seg.size(); ++k)
					if (t_seg[k]) eT += TGroup->corrLgm[k].acf[c][0];

				vector<vector<int>::size_type> indexV; // frames where segment and voiced mask are both 1
				for(vector<int>::size_type k=0; k<sumM.size(); ++k)
					if (sumM[k]>0 && t_seg[k]) indexV.push_back(k);

				if (!indexV.empty()){
					double eTV=0, eTVI=0;
					for (int k=0; k<indexV.size(); ++k){
						// eTV: segment energy in voiced frames
						eTV += TGroup->corrLgm[indexV[k]].acf[c][0];
						// eTVI: segment energy overlapping with voiced mask
						if (mask[indexV[k]][c]>0) eTVI +=TGroup->corrLgm[indexV[k]].acf[c][0];
					}
					if (eTV*2>eT && eTVI*2>eTV){
						for (int m=0; m<mask.size(); ++m)
							mask[m][c] = mask[m][c]+t_seg[m];
					}
				}
			}
		}
	}

	#ifdef DEBUG
        fid=fopen("mask.dbg", "w");
        for(int chan=0; chan<mask.size(); chan++){
                for(int frame=0; frame<mask[chan].size(); frame++){
                        fprintf(fid, "%d ", mask[chan][frame]);
                }
                fprintf(fid, "\n");
        }
        fclose(fid);
        #endif
}

int main(int argc, char *argv[])
{	
	if (argc!=3 && argc!=5 && argc!=6)
	{
		printf("Usage: tandem in out\n");
		printf("       tandem in out cross evCross \n");
		printf("       tandem in out cross evCross eng\n");
		return 1;
	}

	char *filename=argv[1];
	char *outfn=argv[2];

	// initialization
	gammaToneFilterBank *AudiPery;
	voicedMask *TGroup;
	initVoicedMask(AudiPery, TGroup);
			
	// pitch and mask estimation
	voicedMaskEst(filename, AudiPery, TGroup);
	SaveMasks(outfn, TGroup); // save simultaneous streams

	// binarize simultaneous streams
	vector<vector<vector<int> > > allMasks;
	allMasks.resize(TGroup->numContour);
	for(int n=0; n<TGroup->numContour; n++){
		allMasks[n].resize(TGroup->numFrame);
                for(int frame=0; frame<TGroup->numFrame; frame++){
			allMasks[n][frame].resize(TGroup->numberChannel);
                        if(TGroup->Pitch[n].value[frame]>0){
                                for(int chan=0; chan<TGroup->numberChannel; chan++)
					allMasks[n][frame][chan] = TGroup->Pitch[n].mProb[frame].value[chan]>0;
                        }
                }
        }

	// generate T-segments
	cout<<"Onset/offset segmentation..."<<endl;
	segment *TSeg;
	mScaleInten *scaleInten;
	segScale sScale;		
	initOnOff(TSeg, scaleInten, sScale);
	onOffSeg(AudiPery, TGroup, TSeg, scaleInten, sScale);

	// incorporate T-segments
	cout<<"Grouping T-segments..."<<endl;
        for (int i=0; i<TGroup->numContour; ++i)
		expandVoicedMask(allMasks[i], TGroup, TSeg);
	
	// save T-segments incorporated simultaneous streams
	cout<<"Saving outputs..."<<endl;
	SaveMasks(outfn, TGroup, allMasks); // save simultaneous streams

	if (argc==5 || argc==6)
		SaveCross(argv[3], argv[4], TGroup);

	if (argc==6)
		SaveEng(argv[5], TGroup);
	
	delete AudiPery;
	delete TGroup;
	delete TSeg;
        delete scaleInten;

	return 0;
}
