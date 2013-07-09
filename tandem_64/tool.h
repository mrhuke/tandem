#ifndef TOOL_H
#define TOOL_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <algorithm>
using namespace std;

#define PI (3.1415926535897932384626433832795)

inline double sigmoid(double x){ return( 2/(1+exp(-2*x)) -1 ); }

inline double hamming(double x){ return( 0.54 + 0.46*cos(x*PI) ); }

inline double hammingDiffOne(double x){ return( -0.46*PI*sin(x*PI) ); }

inline double hammingDiffTwo(double x){ return( -0.46*PI*PI*cos(x*PI) ); }

double sinc(double tim, double wn);
double sincDiffOne(double tim, double wn);
double sincDiffTwo(double tim, double wn);

double bessi0(double x);

double bessi1(double x);

int maxPos(double *input, int d1, int d2);

int minPos(double *input, int d1, int d2);

void fft(double *inputR, double *inputI, int N, double direct);

double zeroCross(double *input, int sLen);

int findLocalPeak(double *input, int p1, int p2, double theta);

int findLocalValley(double *input, int p1, int p2, double theta);

/* operator overloading */

// max element in a 2D matrix
template <typename T>
T max(vector<vector<T> > &v)
{
        vector<T> row_max;
        for (typename vector<vector<T> >::iterator it=v.begin(); it!=v.end(); ++it)
                row_max.push_back(*max_element(it->begin(), it->end()));
        return *max_element(row_max.begin(), row_max.end());
}
// matrix DOT product
template <typename T>
vector<vector<T> > operator*(const vector<vector<T> > &lhs, const vector<vector<T> > &rhs)
{
        vector<vector<T> > pd;
        for (typename vector<vector<T> >::size_type i=0; i!=lhs.size(); ++i){
                vector<T> row;
                for (typename vector<T>::size_type j=0; j!=lhs[i].size(); ++j)
                        row.push_back(lhs[i][j]*rhs[i][j]);
                pd.push_back(row);
        }
        return pd;
}
// MATLAB-like sum function
template <typename T>
vector<T> sum(vector<vector<T> > &v, unsigned int i)
{
	vector<T> out;
	if(i==1){
		out.resize(v[0].size());
		for(int i=0; i<v.size(); ++i)
			for (int j=0; j<v[i].size(); ++j)
				out[j] += v[i][j];
	}
	else if (i==2){
		out.resize(v.size());
                for(int i=0; i<v.size(); ++i)
                        for (int j=0; j<v[i].size(); ++j)
                                out[i] += v[i][j];
	}
	return out;
}
// vector dot product
template <typename T>
vector<T> operator*(const vector<T> &lhs, const vector<T> &rhs) // note it is DOT product
{
        vector<T> pd;
        for (typename vector<T>::size_type i=0; i!=lhs.size(); ++i)
                pd.push_back(lhs[i]*rhs[i]);
        return pd;
}
// vector addition
template <typename T>
vector<T> operator+(const vector<T> &lhs, const vector<T> &rhs)
{
        vector<T> pd;
        for (typename vector<T>::size_type i=0; i!=lhs.size(); ++i)
                pd.push_back(lhs[i]+rhs[i]);
        return pd;
}
// MATLAB-like == operator
template <typename T>
vector<int> operator==(vector<T> &v, T val)
{
	vector<int> ind;
        for (typename vector<T>::iterator it=v.begin(); it!=v.end(); ++it)
		ind.push_back(static_cast<int>(*it==val));
	return ind;
}

#endif
