/*Program to find roots in a given range from the Bernoulli function with a given lambda and range of position (x) values.  Results are stored as a table in a text file.

Please ensure that nr.h, nrutil_nr.h, and nrtypes_nr.h are in the same folder as this file before running.*/
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "nr.h"
using namespace std;
double lambda;
double position;
DP fx(const DP);

int main(void)
{
	ofstream output;
	string filename = "bondi_roots.txt"; //edit this to change filename that the table will go to
        const int N=100,NBMAX=20;
        const DP X1=-0.5,X2=3.75; //edit this to change initial range that's checked for roots
        int i,nb=NBMAX;
        DP xacc,root;
        Vec_DP xb1(NBMAX),xb2(NBMAX);
	lambda = .25*exp(1.5); //edit this to change lambda
	output.open(filename.c_str());
	for (double k=0.0;k<=10.0;k+=0.1){ //edit "k=0" for minimum position (x), "k<=10" for maximum position, "k+=0.1" for increment value
		position = k;
        	NR::zbrak(fx,X1,X2,N,xb1,xb2,nb);
		for (i=0;i<nb;i++) {
			xacc=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
			root=NR::rtbis(fx,xb1[i],xb2[i],xacc);
			output << position << " " << root << endl;
        	}
	}
	output.close();
        return 0;
}
DP fx(const DP v)
{
        return v*v/2.0+log(lambda)-log(position*position)-log(v)-1.0/position;
}

void NR::zbrak(DP fx(const DP), const DP x1, const DP x2, const int n,
	Vec_O_DP &xb1, Vec_O_DP &xb2, int &nroot)
{
	int i;
	DP x,fp,fc,dx;

	int nb=xb1.size();
	nroot=0;
	dx=(x2-x1)/n;
	fp=fx(x=x1);
	for (i=0;i<n;i++) {
		fc=fx(x += dx);
		if (fc*fp <= 0.0) {
			xb1[nroot]=x-dx;
			xb2[nroot++]=x;
			if(nroot == nb) return;
		}
		fp=fc;
	}
}

DP NR::rtbis(DP func(const DP), const DP x1, const DP x2, const DP xacc)
{
	const int JMAX=40;
	int j;
	DP dx,f,fmid,xmid,rtb;

	f=func(x1);
	fmid=func(x2);
	if (f*fmid >= 0.0) nrerror("Root must be bracketed for bisection in rtbis");
	rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
	for (j=0;j<JMAX;j++) {
		fmid=func(xmid=rtb+(dx *= 0.5));
		if (fmid <= 0.0) rtb=xmid;
		if (fabs(dx) < xacc || fmid == 0.0) return rtb;
	}
	nrerror("Too many bisections in rtbis");
	return 0.0;
}
