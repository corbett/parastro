#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "nrutil.h"

#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

void gaussj(double **a, int n, double **b, int m)
{
	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll;
	double big,dum,pivinv,temp;

	indxc=ivector(1,n);
	indxr=ivector(1,n);
	ipiv=ivector(1,n);
	for (j=1;j<=n;j++) ipiv[j]=0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if (ipiv[j] != 1)
				for (k=1;k<=n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					} else if (ipiv[k] > 1) nrerror("gaussj: Singular Dmatrix-1");
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
			for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) nrerror("gaussj: Singular Dmatrix-2");
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=1;l<=n;l++) a[icol][l] *= pivinv;
		for (l=1;l<=m;l++) b[icol][l] *= pivinv;
		for (ll=1;ll<=n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l=n;l>=1;l--) {
		if (indxr[l] != indxc[l])
			for (k=1;k<=n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
	free_ivector(ipiv,1,n);
	free_ivector(indxr,1,n);
	free_ivector(indxc,1,n);
}

void covsrt(double **covar, int ma, int ia[], int mfit)
{
	int i,j,k;
	double temp;

	for (i=mfit+1;i<=ma;i++)
		for (j=1;j<=i;j++) covar[i][j]=covar[j][i]=0.0;
	k=mfit;
	for (j=ma;j>=1;j--) {
		if (ia[j]) {
			for (i=1;i<=ma;i++) SWAP(covar[i][k],covar[i][j])
			for (i=1;i<=ma;i++) SWAP(covar[k][i],covar[j][i])
			k--;
		}
	}
}

#undef SWAP

void mrqcof(double x[], double y[], double sig[], int ndata, double a[], int ia[],
	int ma, double **alpha, double beta[], double *chisq,
	void (*funcs)(double, double [], double *, double [], int))
{
	int i,j,k,l,m,mfit=0;
	double ymod,wt,sig2i,dy,*dyda;

	dyda=dvector(1,ma);
	for (j=1;j<=ma;j++)
		if (ia[j]) mfit++;
	for (j=1;j<=mfit;j++) {
		for (k=1;k<=j;k++) alpha[j][k]=0.0;
		beta[j]=0.0;
	}
	*chisq=0.0;
	for (i=1;i<=ndata;i++) {
		(*funcs)(x[i],a,&ymod,dyda,ma);
		sig2i=1.0/(sig[i]*sig[i]);
		dy=y[i]-ymod;
		for (j=0,l=1;l<=ma;l++) {
			if (ia[l]) {
				wt=dyda[l]*sig2i;
				for (j++,k=0,m=1;m<=l;m++)
					if (ia[m]) alpha[j][++k] += wt*dyda[m];
				beta[j] += dy*wt;
			}
		}
		*chisq += dy*dy*sig2i;
	}
	for (j=2;j<=mfit;j++)
		for (k=1;k<j;k++) alpha[k][j]=alpha[j][k];
	free_dvector(dyda,1,ma);
}


void mrqmin(double x[], double y[], double sig[], int ndata, double a[], int ia[],
	int ma, double **covar, double **alpha, double *chisq,
	void (*funcs)(double, double [], double *, double [], int), double *alamda)
{
	void covsrt(double **covar, int ma, int ia[], int mfit);
	void gaussj(double **a, int n, double **b, int m);
	void mrqcof(double x[], double y[], double sig[], int ndata, double a[],
		int ia[], int ma, double **alpha, double beta[], double *chisq,
		void (*funcs)(double, double [], double *, double [], int));
	int j,k,l,m;
	static int mfit;
	static double ochisq,*atry,*beta,*da,**oneda;

	if (*alamda < 0.0) {
		atry=dvector(1,ma);
		beta=dvector(1,ma);
		da=dvector(1,ma);
		for (mfit=0,j=1;j<=ma;j++)
			if (ia[j]) mfit++;
		oneda=dmatrix(1,mfit,1,1);
		*alamda=0.001;
		mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,chisq,funcs);
		ochisq=(*chisq);
		for (j=1;j<=ma;j++) atry[j]=a[j];
	}
	for (j=0,l=1;l<=ma;l++) {
		if (ia[l]) {
			for (j++,k=0,m=1;m<=ma;m++) {
				if (ia[m]) {
					k++;
					covar[j][k]=alpha[j][k];
				}
			}
			covar[j][j]=alpha[j][j]*(1.0+(*alamda));
			oneda[j][1]=beta[j];
		}
	}
	gaussj(covar,mfit,oneda,1);
	for (j=1;j<=mfit;j++) da[j]=oneda[j][1];
	if (*alamda == 0.0) {
		covsrt(covar,ma,ia,mfit);
		free_dmatrix(oneda,1,mfit,1,1);
		free_dvector(da,1,ma);
		free_dvector(beta,1,ma);
		free_dvector(atry,1,ma);
		return;
	}
	for (j=0,l=1;l<=ma;l++)
		if (ia[l]) atry[l]=a[l]+da[++j];
	mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,chisq,funcs);
	if (*chisq < ochisq) {
		*alamda *= 0.1;
		ochisq=(*chisq);
		for (j=0,l=1;l<=ma;l++) {
			if (ia[l]) {
				for (j++,k=0,m=1;m<=ma;m++) {
					if (ia[m]) {
						k++;
						alpha[j][k]=covar[j][k];
					}
				}
				beta[j]=da[j];
				a[l]=atry[l];
			}
		}
	} else {
		*alamda *= 10.0;
		*chisq=ochisq;
	}
}

/*
** The 3-parameter Einasto model.
*/
void einasto(double r,double *a,double *rho,double *drhoda,int na) {
    double ra;

    assert(na == 3);
    ra = pow(r/a[2],a[3]);
    assert(a[3] > 0);
    drhoda[1] = exp(-2.0/a[3]*(ra - 1.0));
    *rho = a[1]*drhoda[1];
    drhoda[2] = (*rho)*(2.0/a[2]*ra);
    drhoda[3] = (*rho)*(2.0/(a[3]*a[3])*(ra - 1.0) - 2.0/a[3]*ra*log(r/a[2]));
/*
    printf("einasto: r:%g rho:%g drhoda:%g %g %g\n",r,*rho,
	drhoda[1],drhoda[2],drhoda[3]); 
*/
    }

void log_einasto(double r,double *a,double *rho,double *drhoda,int na) {
    einasto(r,a,rho,drhoda,na);
    drhoda[1] /= *rho;
    drhoda[2] /= *rho;
    drhoda[3] /= *rho;
    *rho = log(*rho);
    }

/*
** The 3-parameter lambda model.
*/
void lambda(double r,double *a,double *rho,double *drhoda,int na) {
    double ra;

    assert(na == 3);
    ra = r/a[2];
    drhoda[1] = exp(-a[3]*pow(log(1+ra),2));
    *rho = a[1]*drhoda[1];
    drhoda[2] = (*rho)*2*a[3]/a[2]*log(1+ra)*ra/(1+ra);
    drhoda[3] = -(*rho)*pow(log(1+ra),2);
/*
    printf("lambda: r:%g rho:%g drhoda:%g %g %g\n",r,*rho,
	drhoda[1],drhoda[2],drhoda[3]); 
*/
    }

void log_lambda(double r,double *a,double *rho,double *drhoda,int na) {
    lambda(r,a,rho,drhoda,na);
    drhoda[1] /= *rho;
    drhoda[2] /= *rho;
    drhoda[3] /= *rho;
    *rho = log(*rho);
    }

/*
** The 2-parameter lambda model.
*/
void lambda2(double r,double *a,double *rho,double *drhoda,int na) {
    double ra;

    assert(na == 2);
    ra = r/a[2];
    drhoda[1] = exp(-0.1*pow(log(1+ra),2));
    *rho = a[1]*drhoda[1];
    drhoda[2] = (*rho)*2*0.1/a[2]*log(1+ra)*ra/(1+ra);
/*
    printf("lambda: r:%g rho:%g drhoda:%g %g %g\n",r,*rho,
	drhoda[1],drhoda[2],drhoda[3]); 
*/
    }

void log_lambda2(double r,double *a,double *rho,double *drhoda,int na) {
    lambda2(r,a,rho,drhoda,na);
    drhoda[1] /= *rho;
    drhoda[2] /= *rho;
    *rho = log(*rho);
    }

/*
** The 4-parameter Prugniel-Simien model, the 4-th parameter is p, but this can be
** espressed as a function of alpha, which is done in the 3 parameter version below.
** Parameter a[3] is alpha which is defined to be 1/n, so that it becomes more 
** comparable to the einasto profile given above. Parameter a[4] = p which is the 
** power-law bit, p being the central slope, it equals zero for Einasto and is 
** like the gamma parameter of Hernquist.
*/ 
void prugnielsim4(double r,double *a,double *rho,double *drhoda,int na) {
    double ra,b;
    assert(na == 4);
    assert(a[2] > 0);
    ra = pow(r/a[2],a[3]);
    assert(a[3] > 0);
    b = 2.0/a[3] - 1.0/3.0 + 0.009876*a[3];
    drhoda[1] = pow(r/a[2],-a[4])*exp(-b*ra);
    *rho = a[1]*drhoda[1];
    drhoda[2] = (*rho)/a[2]*(a[4] + b*a[3]*ra);
    drhoda[3] = (*rho)*ra*(2.0/(a[3]*a[3]) - 0.009876 - b*log(r/a[2]));
    drhoda[4] = -(*rho)*log(r/a[2]);
    }

void log_prugnielsim4(double r,double *a,double *rho,double *drhoda,int na) {
    prugnielsim4(r,a,rho,drhoda,na);
    drhoda[1] /= *rho;
    drhoda[2] /= *rho;
    drhoda[3] /= *rho;
    drhoda[4] /= *rho;
    *rho = log(*rho);
    }

/*
** A 3-parameter version of the Prugniel-Simien model. This is close to the 
** deprojected Sersic (R/Re)^(1/n) model for a *limited* range of R and n!
** Was the favoured fit for clusters in the Merrit et. al. 2005 papers, 
** although the method of determining the simulation rho(r) was a problem
** there.
*/
void prugnielsim3(double r,double *a,double *rho,double *drhoda,int na) {
    double ap[5];
    double dap[5];

    assert(na == 3);
    ap[1] = a[1];
    ap[2] = a[2];
    ap[3] = a[3];
    ap[4] = 1.0 - a[3]*(0.6097 - 0.05463*a[3]);
    prugnielsim4(r,ap,rho,dap,4);
    drhoda[1] = dap[1];
    drhoda[2] = dap[2];
    drhoda[3] = dap[3] + dap[4]*(-0.6097 + 2.0*0.05463*a[3]);
    }

void log_prugnielsim3(double r,double *a,double *rho,double *drhoda,int na) {
    prugnielsim3(r,a,rho,drhoda,na);
    drhoda[1] /= *rho;
    drhoda[2] /= *rho;
    drhoda[3] /= *rho;
    *rho = log(*rho);
    }

/*
** The 5-parameter (alpha,beta,gamma) Hernquist model.
*/
void hernquist(double r,double *a,double *rho,double *drhoda,int na) {
    double ra;
    assert(na == 5);
    ra = pow(r/a[2],a[3]);
    drhoda[1] = pow(2.0,(a[4]-a[5])/a[3])*pow(r/a[2],-a[5])*pow(1.0+ra,(a[5]-a[4])/a[3]);
    *rho = a[1]*drhoda[1];
    drhoda[2] = (*rho)/a[2]*(a[4] + (a[5]-a[4])/(1.0+ra));
    drhoda[3] = (*rho)*(a[5]-a[4])/a[3]*(log(2.0)/a[3] - log(1.0+ra)/a[3] + ra*log(r/a[2])/(1+ra));
    drhoda[4] = (*rho)/a[3]*(log(2.0) - log(1.0+ra));
    drhoda[5] = (*rho)/a[3]*(log(2.0) - a[3]*log(r/a[2]) + log(1.0+ra));
    }

void log_hernquist(double r,double *a,double *rho,double *drhoda,int na) {
    hernquist(r,a,rho,drhoda,na);
    drhoda[1] /= *rho;
    drhoda[2] /= *rho;
    drhoda[3] /= *rho;
    drhoda[4] /= *rho;
    drhoda[5] /= *rho;
    *rho = log(*rho);
    }

/*
** The popular 2 and 3 parameter Hernquist-like (double power-law) models.
*/
void nfw(double r,double *a,double *rho,double *drhoda,int na) {
    double ah[6];
    double dah[6];

    assert(na == 2);
    ah[1] = a[1];
    ah[2] = a[2];
    ah[3] = 1.0;
    ah[4] = 3.0;
    ah[5] = 1.0;
    hernquist(r,ah,rho,dah,5);
    drhoda[1] = dah[1];
    drhoda[2] = dah[2];
    }

void log_nfw(double r,double *a,double *rho,double *drhoda,int na) {
    nfw(r,a,rho,drhoda,na);
    drhoda[1] /= *rho;
    drhoda[2] /= *rho;
    *rho = log(*rho);
    }

void moore(double r,double *a,double *rho,double *drhoda,int na) {
    double ah[6];
    double dah[6];

    assert(na == 2);
    ah[1] = a[1];
    ah[2] = a[2];
    ah[3] = 1.5;
    ah[4] = 3.0;
    ah[5] = 1.5;
    hernquist(r,ah,rho,dah,5);
    drhoda[1] = dah[1];
    drhoda[2] = dah[2];
    }

void log_moore(double r,double *a,double *rho,double *drhoda,int na) {
    moore(r,a,rho,drhoda,na);
    drhoda[1] /= *rho;
    drhoda[2] /= *rho;
    *rho = log(*rho);
    }

void generalnfw(double r,double *a,double *rho,double *drhoda,int na) {
    double ah[6];
    double dah[6];

    assert(na == 3);
    ah[1] = a[1];
    ah[2] = a[2];
    ah[3] = 1.0;
    ah[4] = 3.0;
    ah[5] = a[3];
    hernquist(r,ah,rho,dah,5);
    drhoda[1] = dah[1];
    drhoda[2] = dah[2];
    drhoda[3] = dah[5];
    }

void log_generalnfw(double r,double *a,double *rho,double *drhoda,int na) {
    generalnfw(r,a,rho,drhoda,na);
    drhoda[1] /= *rho;
    drhoda[2] /= *rho;
    drhoda[3] /= *rho;
    *rho = log(*rho);
    }

void dehnenmcl2(double r,double *a,double *rho,double *drhoda,int na) {
    double ah[6];
    double dah[6];

    assert(na == 2);
    ah[1] = a[1];
    ah[2] = a[2];
    ah[3] = 4.0/9.0;
    ah[4] = 31.0/9.0;
    ah[5] = 7.0/9.0;
    hernquist(r,ah,rho,dah,5);
    drhoda[1] = dah[1];
    drhoda[2] = dah[2];
    }

void log_dehnenmcl2(double r,double *a,double *rho,double *drhoda,int na) {
    dehnenmcl2(r,a,rho,drhoda,na);
    drhoda[1] /= *rho;
    drhoda[2] /= *rho;
    *rho = log(*rho);
    }

void dehnenmcl3(double r,double *a,double *rho,double *drhoda,int na) {
    double ah[6];
    double dah[6];

    assert(na == 3);
    ah[1] = a[1];
    ah[2] = a[2];
    ah[3] = (4.0-2.0*a[3])/9.0;
    ah[4] = (31.0-2.0*a[3])/9.0;
    ah[5] = (7.0+10.0*a[3])/9.0;
    hernquist(r,ah,rho,dah,5);
    drhoda[1] = dah[1];
    drhoda[2] = dah[2];
    drhoda[3] = dah[3]*(-2.0/9.0) + dah[4]*(-2.0/9.0) + dah[5]*(10.0/9.0);
    }

void log_dehnenmcl3(double r,double *a,double *rho,double *drhoda,int na) {
    dehnenmcl3(r,a,rho,drhoda,na);
    drhoda[1] /= *rho;
    drhoda[2] /= *rho;
    drhoda[3] /= *rho;
    *rho = log(*rho);
    }


#define FIELDS 7
#define RHO_FIELD 3
#define R_FIELD2 4
#define RHO_FIELD2 6 

/*
** profit <filename> <rMin> <rMax>
*/
int main(int argc,char **argv) {
    FILE *fp,*fpo;
    double *r;
    double *rho;
    double *sig;
    double *lrho;
    double *lsig;
    double *a;
    int *lista;
    double **covar;
    double **alpha;
    double chisq;
    double alamda;
    int ma,mfit;    
    unsigned nBin;
    int nb,iRet,i,ib;
    double fTmp,rTmp,rMin,rMax;
    double *drhoda;
    double *resEinasto,*resNFW,*resMoore,*resGenNFW,*resDehnen2,*resDehnen3,*resPrug,*resLambda2,*resLambda3;
    int iter;
    const double rhoconv = 1.0/(1.0361e-16*4e4*4e4*4e4);

    assert(argc == 4);
    rMin = atof(argv[2]);
    assert(rMin > 0);
    rMax = atof(argv[3]);
    assert(rMax > 0);
    /*
    ** First find out how many lines between rMin and rMax are in the file.
    */
    fp = fopen(argv[1],"r");
    assert(fp != NULL);
    nb = 0;
    for (iter=0;iter<200;++iter) {
	iRet = fscanf(fp,"%lf %u",&rTmp,&nBin);
	if (iRet != 2) break;
	for (i=3;i<=FIELDS;++i) {
	    iRet = fscanf(fp,"%lf",&fTmp);
	    assert(iRet == 1);
	    }
	if (rTmp >= rMin && rTmp <= rMax) ++nb;
	}
    rewind(fp);

    printf("nb:%d\n",nb);
    /*
    ** Allocate the arrays.
    */
    r = dvector(1,nb);
    rho = dvector(1,nb);
    resNFW = dvector(1,nb);
    resMoore = dvector(1,nb);
    resDehnen2 = dvector(1,nb);
    resEinasto = dvector(1,nb);
    resGenNFW = dvector(1,nb);
    resDehnen3 = dvector(1,nb);
    resPrug = dvector(1,nb);
    resLambda2 = dvector(1,nb);
    resLambda3 = dvector(1,nb);
    lrho = dvector(1,nb);
    sig = dvector(1,nb);
    lsig = dvector(1,nb);
    ib = 1;
    for (iter=0;iter<200;++iter) {
	iRet = fscanf(fp,"%lf %u",&rTmp,&nBin);
	if (iRet != 2) break;
	r[ib] = rTmp;
	for (i=3;i<=FIELDS;++i) {
	    iRet = fscanf(fp,"%lf",&fTmp);
	    assert(iRet == 1);
	    if (i == RHO_FIELD && ib <= nb) {
		rho[ib] = fTmp;
		lrho[ib] = log(fTmp);
		sig[ib] = fTmp;
		lsig[ib] = nb-3.0;
		}
	    }
	printf("%d %g %f\n",ib,rTmp,rho[ib]);
	if (rTmp >= rMin && rTmp <= rMax) ++ib;
	}

    ma = 5;
    a = dvector(1,ma);
    drhoda = dvector(1,ma);
    lista = ivector(1,ma);
    covar = dmatrix(1,ma,1,ma);
    alpha = dmatrix(1,ma,1,ma);

    /*
    ** Initial guess at parameters for the NFW profile.
    */
    ma = 2;
    a[1] = 5.7e-4*1.0361e-16*4e7*4e7*4e7;
    a[2] = 29.0/40000.0;
    chisq = -1.0;
    printf("\n%6d %.14g %.14g %.14g\n",0,a[1],a[2],chisq); 
    lista[1] = 1;
    lista[2] = 1;
    alamda = -1.0;
    i = 1;
    for (iter=0;iter<200;++iter) {
	mrqmin(r,lrho,lsig,nb,a,lista,ma,covar,alpha,&chisq,log_nfw,&alamda);
	printf("NFW %6d %.14g %.14g %g %g\n",i,a[1]*rhoconv,a[2]*40000.0,chisq,alamda); 
	++i;
	if (alamda < 1e-6) break;
	}
    for (ib=1;ib<=nb;++ib) {
	nfw(r[ib],a,&fTmp,drhoda,ma);
	resNFW[ib] = (rho[ib] - fTmp)/fTmp;
	}
    /*
    ** Initial guess at parameters for the Moore profile.
    */
    ma = 2;
    a[1] = 5.7e-4*1.0361e-16*4e7*4e7*4e7;
    a[2] = 29.0/40000.0;
    chisq = -1.0;
    printf("\n%6d %.14g %.14g %.14g\n",0,a[1],a[2],chisq); 
    lista[1] = 1;
    lista[2] = 1;
    alamda = -1.0;
    i = 1;
    for (iter=0;iter<200;++iter) {
	mrqmin(r,lrho,lsig,nb,a,lista,ma,covar,alpha,&chisq,log_moore,&alamda);
	printf("Moore %6d %.14g %.14g %g %g\n",i,a[1]*rhoconv,a[2]*40000.0,chisq,alamda); 
	++i;
	if (alamda < 1e-7) break;
	}
    for (ib=1;ib<=nb;++ib) {
	moore(r[ib],a,&fTmp,drhoda,ma);
	resMoore[ib] = (rho[ib] - fTmp)/fTmp;
	}
    /*
    ** Initial guess at parameters for the Dehnen-McLaughlin 2-parameter profile.
    */
    ma = 2;
    a[1] = 5.7e-4*1.0361e-16*4e7*4e7*4e7;
    a[2] = 29.0/40000.0;
    chisq = -1.0;
    printf("\n%6d %.14g %.14g %.14g\n",0,a[1],a[2],chisq); 
    lista[1] = 1;
    lista[2] = 1;
    alamda = -1.0;
    i = 1;
    for (iter=0;iter<200;++iter) {
	mrqmin(r,lrho,lsig,nb,a,lista,ma,covar,alpha,&chisq,log_dehnenmcl2,&alamda);
	printf("Dehnen-McLaughlin2 %6d %.14g %.14g %g %g\n",i,a[1]*rhoconv,a[2]*40000.0,chisq,alamda); 
	++i;
	if (alamda < 1e-4) break;
	}
    for (ib=1;ib<=nb;++ib) {
	dehnenmcl2(r[ib],a,&fTmp,drhoda,ma);
	resDehnen2[ib] = (rho[ib] - fTmp)/fTmp;
	}
    /*
    ** Initial guess at parameters for the lambda 2-par profile.
    */
    ma = 2;
    a[1] = 1e10;
    a[2] = 1e-7;
    chisq = -1.0;
    printf("\n%6d %.14g %.14g %.14g\n",0,a[1],a[2],chisq); 
    lista[1] = 1;
    lista[2] = 1;
    alamda = -1.0;
    i = 1;
    for (iter=0;iter<200;++iter) {
	mrqmin(r,lrho,lsig,nb,a,lista,ma,covar,alpha,&chisq,log_lambda2,&alamda);
	printf("lambda %6d %.14g %.14g %g %g\n",i,a[1]*rhoconv,a[2]*40000.0,chisq,alamda); 
	++i;
	if (alamda < 1e-6) break;
	}
    for (ib=1;ib<=nb;++ib) {
	lambda2(r[ib],a,&fTmp,drhoda,ma);
	resLambda2[ib] = (rho[ib] - fTmp)/fTmp;
	}
    /*
    ** Initial guess at parameters for the Einasto profile.
    ** These come from Doug's fit.
    */
    ma = 3;
    a[1] = 5.7e-4*1.0361e-16*4e7*4e7*4e7;
    a[2] = 29.0/40000.0;
    a[3] = 0.15;
    chisq = -1.0;
    printf("\n%6d %.14g %.14g %.14g %.14g\n",0,a[1],a[2],a[3],chisq); 
    lista[1] = 1;
    lista[2] = 1;
    lista[3] = 1;    
    alamda = -1.0;
    i = 1;
    for (iter=0;iter<200;++iter) {
	mrqmin(r,lrho,lsig,nb,a,lista,ma,covar,alpha,&chisq,log_einasto,&alamda);
	printf("Einasto %6d %.14g %.14g %.14g %g %g\n",i,a[1]*rhoconv,a[2]*40000.0,a[3],chisq,alamda); 
	++i;
	if (alamda < 1e-6) break;
	}
    for (ib=1;ib<=nb;++ib) {
	einasto(r[ib],a,&fTmp,drhoda,ma);
	resEinasto[ib] = (rho[ib] - fTmp)/fTmp;
	}
    /*
    ** Initial guess at parameters for the gen_nfw profile.
    */
    ma = 3;
    a[1] = 5.7e-4*1.0361e-16*4e7*4e7*4e7;
    a[2] = 29.0/40000.0;
    a[3] = 1.0;
    chisq = -1.0;
    printf("\n%6d %.14g %.14g %.14g %.14g\n",0,a[1],a[2],a[3],chisq); 
    lista[1] = 1;
    lista[2] = 1;
    lista[3] = 1;    
    alamda = -1.0;
    i = 1;
    for (iter=0;iter<200;++iter) {
	mrqmin(r,lrho,lsig,nb,a,lista,ma,covar,alpha,&chisq,log_generalnfw,&alamda);
	printf("General NFW %6d %.14g %.14g %.14g %g %g\n",i,a[1]*rhoconv,a[2]*40000.0,a[3],chisq,alamda); 
	++i;
	if (alamda < 1e-8) break;
	}
    for (ib=1;ib<=nb;++ib) {
	generalnfw(r[ib],a,&fTmp,drhoda,ma);
	resGenNFW[ib] = (rho[ib] - fTmp)/fTmp;
	}
    /*
    ** Initial guess at parameters for the Dehnen-McLaughlin 3-par profile.
    */
    ma = 3;
    a[1] = 5.7e-5*1.0361e-16*4e7*4e7*4e7;
    a[2] = 29.0/40000.0;
    a[3] = -0.5;
    chisq = -1.0;
    printf("\n%6d %.14g %.14g %.14g %.14g\n",0,a[1],a[2],a[3],chisq); 
    lista[1] = 1;
    lista[2] = 1;
    lista[3] = 1;    
    alamda = -1.0;
    i = 1;
    for (iter=0;iter<200;++iter) {
	mrqmin(r,lrho,lsig,nb,a,lista,ma,covar,alpha,&chisq,log_dehnenmcl3,&alamda);
	printf("Dehnen-McLaughlin3 %6d %.14g %.14g %.14g %g %g\n",i,a[1]*rhoconv,a[2]*40000.0,a[3],chisq,alamda); 
	++i;
	if (alamda < 1e-8) break;
	}
    for (ib=1;ib<=nb;++ib) {
	dehnenmcl3(r[ib],a,&fTmp,drhoda,ma);
	resDehnen3[ib] = (rho[ib] - fTmp)/fTmp;
	}
    /*
    ** Initial guess at parameters for the Prugniel-Simien 3-par profile.
    */
    ma = 3;
    a[1] = 5e4;
    a[2] = 1e-3;
    a[3] = 0.3;
    chisq = -1.0;
    printf("\n%6d %.14g %.14g %.14g %.14g\n",0,a[1],a[2],a[3],chisq); 
    lista[1] = 1;
    lista[2] = 1;
    lista[3] = 1;    
    alamda = -1.0;
    i = 1;
    for (iter=0;iter<200;++iter) {
	mrqmin(r,lrho,lsig,nb,a,lista,ma,covar,alpha,&chisq,log_prugnielsim3,&alamda);
	printf("Prugniel-Simien3 %6d %.14g %.14g %.14g %g %g\n",i,a[1]*rhoconv,a[2]*40000.0,a[3],chisq,alamda); 
	++i;
	if (alamda < 1e-7) break;
	}
    for (ib=1;ib<=nb;++ib) {
	prugnielsim3(r[ib],a,&fTmp,drhoda,ma);
	resPrug[ib] = (rho[ib] - fTmp)/fTmp;
	}
    /*
    ** Initial guess at parameters for the lambda 3-par profile.
    */
    ma = 3;
    a[1] = 1e10;
    a[2] = 1e-7;
    a[3] = 3.0/8.0/log(10.0);
    chisq = -1.0;
    printf("\n%6d %.14g %.14g %.14g %.14g\n",0,a[1],a[2],a[3],chisq); 
    lista[1] = 1;
    lista[2] = 1;
    lista[3] = 1;    
    alamda = -1.0;
    i = 1;
    for (iter=0;iter<200;++iter) {
	mrqmin(r,lrho,lsig,nb,a,lista,ma,covar,alpha,&chisq,log_lambda,&alamda);
	printf("lambda %6d %.14g %.14g %.14g %g %g\n",i,a[1]*rhoconv,a[2]*40000.0,a[3],chisq,alamda); 
	++i;
	if (alamda < 1e-7) break;
	}
    for (ib=1;ib<=nb;++ib) {
	lambda(r[ib],a,&fTmp,drhoda,ma);
	resLambda3[ib] = (rho[ib] - fTmp)/fTmp;
	}
#if 0
    /*
    ** Initial guess at parameters for the Prugniel-Simien 4-par profile.
    */
    ma = 4;
    a[1] = 5e4;
    a[2] = 3e-3;
    a[3] = 0.4;
    a[4] = 1.0;
    chisq = -1.0;
    printf("\n%6d %.14g %.14g %.14g %.14g %.14g\n",0,a[1],a[2],a[3],a[4],chisq); 
    lista[1] = 1;
    lista[2] = 1;
    lista[3] = 0;    
    lista[4] = 1;
    alamda = -1.0;
    i = 1;
    for (iter=0;iter<200;++iter) {
	mrqmin(r,lrho,lsig,nb,a,lista,ma,covar,alpha,&chisq,log_prugnielsim4,&alamda);
	printf("Prugniel-Simien4 %6d %.14g %.14g %.14g %.14g %g %g\n",i,a[1],a[2],a[3],a[4],chisq,alamda); 
	++i;
	if (alamda < 1e-5) break;
	}
    lista[1] = 1;
    lista[2] = 1;
    lista[3] = 1;    
    lista[4] = 1;
    chisq = -1.0;
    alamda = -1.0;
    i = 1;
    for (iter=0;iter<200;++iter) {
	mrqmin(r,lrho,lsig,nb,a,lista,ma,covar,alpha,&chisq,log_prugnielsim4,&alamda);
	printf("Prugniel-Simien4 %6d %.14g %.14g %.14g %.14g %g %g\n",i,a[1],a[2],a[3],a[4],chisq,alamda); 
	++i;
	if (alamda < 1e-6) break;
	}
#endif

    fpo = fopen("profit.out","w");
    assert(fpo != NULL);
    for (ib=1;ib<=nb;++ib) {
	fprintf(fpo,"%.14g %.14g",r[ib],rho[ib]);
	fprintf(fpo," %.14g",resNFW[ib]);
	fprintf(fpo," %.14g",resMoore[ib]);
	fprintf(fpo," %.14g",resDehnen2[ib]);
	fprintf(fpo," %.14g",resEinasto[ib]);
	fprintf(fpo," %.14g",resGenNFW[ib]);
	fprintf(fpo," %.14g",resDehnen3[ib]);
	fprintf(fpo," %.14g",resPrug[ib]);
	fprintf(fpo," %.14g",resLambda3[ib]);
	fprintf(fpo," %.14g",resLambda2[ib]);
	fprintf(fpo,"\n");
	}
    fclose(fpo);
    }

