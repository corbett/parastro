#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>


#include "tipsy.h"


/*
** simpleprof <filename> <rMax> <rMin> <nBins> <rMinTarget> [x|y|z]
*/
int main(int argc,char **argv) {
    const float ghaloCenter[3] = {1.332850e-03, -2.497214e-02, 1.053231e-02};
    const float b2Center[3] = {1.75830873e-03, -2.54608821e-02, 1.25542628e-02};
    const float b1Center[3] = {1.624474e-03, -2.555602e-02, 1.255183e-02};
    const float b1newCenter[3] = {1.765654e-03, -2.546322e-02, 1.265646e-02};
    const float b0Center[3] = {-2.179324e-01, -4.479892e-01, -3.571862e-01};
    const float vl2Center[3] = {3.612025e-01, 2.106574e-01, 6.877048e-03};
    const float silver0Center[3] = {8.391391e-02, -1.194611e-02, 8.308139e-02};
    TCTX in;
    unsigned np,i;
    int j,type,iRet,nBins,nb,ib,ibo;
    double dTime,dSoft;
    struct base_particle *p;
    double rel[3];
    double r2,rMax,rMin,rMinTarget,rMinTrue,dLogRbin,diLogRbin,r2max;
    double dir2Min,r2Max,rMinOffset,dir2oMin;
    double fCenter[3];
    double *dMassBin;
    double *dDensBin;
    unsigned *nNumBin;
    double *doMassBin; /* offset mass bin */
    double *doDensBin; /* offset density bin */
    unsigned *noNumBin; /* offset number of particles bin */
    int iAxis = -1.0;
    double ct2;
    double lScaleKpc;

    if (strcasestr(argv[1],"ghalo")) {
	fprintf(stderr,"using ghalo center coordinates\n");
	for (j=0;j<3;++j) {
	    fCenter[j] = ghaloCenter[j];
	    }
	lScaleKpc = 40000.0;
	}
    else if (strcasestr(argv[1],"b2")) {
	fprintf(stderr,"using b2 center coordinates\n");
	for (j=0;j<3;++j) {
	    fCenter[j] = b2Center[j];
	    }
	lScaleKpc = 40000.0;
	}
    else if (strcasestr(argv[1],"b1_new")) {
	fprintf(stderr,"using b1_new center coordinates\n");
	for (j=0;j<3;++j) {
	    fCenter[j] = b1newCenter[j];
	    }
	lScaleKpc = 40000.0;
	}
    else if (strcasestr(argv[1],"b1")) {
	fprintf(stderr,"using b1 center coordinates\n");
	for (j=0;j<3;++j) {
	    fCenter[j] = b1Center[j];
	    }
	lScaleKpc = 40000.0;
	}
    else if (strcasestr(argv[1],"b0")) {
	fprintf(stderr,"using b0 center coordinates\n");
	for (j=0;j<3;++j) {
	    fCenter[j] = b0Center[j];
	    }
	lScaleKpc = 40000.0;
	}
    else if (strcasestr(argv[1],"vl2")) {
	fprintf(stderr,"using vl2 center coordinates\n");
	for (j=0;j<3;++j) {
	    fCenter[j] = vl2Center[j];
	    }
	lScaleKpc = 40000.0;
	}
    else if (strcasestr(argv[1],"silver0")) {
	fprintf(stderr,"using silver0 center coordinates\n");
	for (j=0;j<3;++j) {
	    fCenter[j] = silver0Center[j];
	    }
	lScaleKpc = 180000.0;
	}
    else {
	fprintf(stderr,"Could not match simulation name!\n");
	exit(1);
	}

    TipsyInitialize(&in,0,argv[1]);

    rMax = atof(argv[2])/lScaleKpc;
    rMin = atof(argv[3])/lScaleKpc;
    nBins = atoi(argv[4]);
    rMinTarget = atof(argv[5])/lScaleKpc;
    if (argc == 7) {
	if (argv[6][0] == 'x') iAxis = 0;
	else if (argv[6][0] == 'y') iAxis = 1;
	else if (argv[6][0] == 'z') iAxis = 2;
	fprintf(stderr,"doing %c-axis\n",'x'+iAxis);
	}

    dLogRbin = (log(rMax) - log(rMin))/nBins;
    nb = ceil((log(rMax) - log(rMinTarget))/dLogRbin);
    rMinTrue = exp(log(rMax) - nb*dLogRbin);
    dir2Min = 1.0/(rMinTrue*rMinTrue);
    rMinOffset = exp(log(rMax) - (nb-0.5)*dLogRbin);
    dir2oMin = 1.0/(rMinOffset*rMinOffset);
    diLogRbin = 0.5/dLogRbin;
    fprintf(stderr,"rMinTrue:%.14g number of bins:%d dLogRbin:%.14g\n",rMinTrue*lScaleKpc,nb,dLogRbin);

    /*
    ** Allocate bins.
    */
    dMassBin = malloc(nb*sizeof(double));
    assert(dMassBin != NULL);
    dDensBin = malloc(nb*sizeof(double));
    assert(dDensBin != NULL);
    doMassBin = malloc(nb*sizeof(double));
    assert(doMassBin != NULL);
    doDensBin = malloc(nb*sizeof(double));
    assert(doDensBin != NULL);
    nNumBin = malloc(nb*sizeof(unsigned));
    assert(nNumBin != NULL);
    noNumBin = malloc(nb*sizeof(unsigned));
    assert(noNumBin != NULL);
    for (ib=0;ib<nb;++ib) {
	nNumBin[ib] = 0;
	noNumBin[ib] = 0;
	dMassBin[ib] = 0.0;
	doMassBin[ib] = 0.0;
	}

    np = iTipsyNumParticles(in);
    dTime = dTipsyTime(in);
    r2Max = rMax*rMax;
    ct2 = pow(cos(15.0/180.0*M_PI),2); /* 15 degree opening angle about the axis */
    for (i=0;i<np;++i) {
	p = pTipsyRead(in,&type,&dSoft);
	r2 = 0.0;
	for (j=0;j<3;++j) {
	    rel[j] = p->pos[j] - fCenter[j];
	    r2 += rel[j]*rel[j];
	    }
        if (iAxis == -1 || rel[iAxis]*rel[iAxis] > ct2*r2) {
	    if (r2 < r2Max) {
		ib = floor(log(r2*dir2Min)*diLogRbin);
		ibo = floor(log(r2*dir2oMin)*diLogRbin);
		if (ib >= 0) {
		    ++nNumBin[ib];
		    dMassBin[ib] += p->mass;
		    /*
		    ** Add the particle to the offset bin.
		    */
		    if (ibo < nb-1) {
			doMassBin[ibo] += p->mass;
			}
		    }
		}
	    }
	}
    for (ib=0;ib<nb;++ib) {
	double ri = rMinTrue*exp(ib*dLogRbin);
	double ro = rMinTrue*exp((ib+1)*dLogRbin);
	double vol = 4.0/3.0*M_PI*(ro*ro*ro - ri*ri*ri);
	dDensBin[ib] = dMassBin[ib]/vol;
	}
    /*
    ** Now calculate stuff for the offset bins.
    */
    for (ib=0;ib<(nb-1);++ib) {
	double ri = rMinOffset*exp(ib*dLogRbin);
	double ro = rMinOffset*exp((ib+1)*dLogRbin);
	double vol = 4.0/3.0*M_PI*(ro*ro*ro - ri*ri*ri);
	doDensBin[ib] = doMassBin[ib]/vol;
	}

    /*
    ** Output the profile and logarithmic slope.
    */
    for (ib=0;ib<nb;++ib) {
	double ri = rMinOffset*exp(ib*dLogRbin);
	double ro = rMinOffset*exp((ib+1)*dLogRbin);
	double rmid = rMinTrue*exp((ib+0.5)*dLogRbin);
	double rmido = rMinOffset*exp((ib+0.5)*dLogRbin);
	printf("%.14g %u %.14g ",rmid,nNumBin[ib],dDensBin[ib]);
	if (ib < (nb-1)) {
	    double dlogrho = rmido/doDensBin[ib]*((dDensBin[ib+1]-dDensBin[ib])/(ro-ri));
	    printf("%.14g %u %.14g %.14g\n",rmido,noNumBin[ib],doDensBin[ib],dlogrho);
	    }
	else {
	    /*
	    ** For the last point just duplicate the previous data so that SM does not plot anything
	    ** spurious. This is just a duplication of the previous line's data, so on a plot there
	    ** would be 2 points ontop of each other at the outer edge of the halo.
	    */
	    double ri = rMinOffset*exp((ib-1)*dLogRbin);
	    double ro = rMinOffset*exp(ib*dLogRbin);
	    double rmido = rMinOffset*exp((ib-0.5)*dLogRbin);
	    double dlogrho = rmido/doDensBin[ib-1]*((dDensBin[ib]-dDensBin[ib-1])/(ro-ri));
	    printf("%.14g %u %.14g %.14g\n",rmido,noNumBin[ib-1],doDensBin[ib-1],dlogrho);
	    }
	}
    TipsyFinish(in);
    }

