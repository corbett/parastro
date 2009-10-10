#ifndef MASTER_HINCLUDED
#define MASTER_HINCLUDED
#define MASTER_H_MODULE_ID "$Id: master.h,v 1.72 2009/10/03 22:33:23 wadsley Exp $"

#include <stdint.h>

#include "param.h"
#include "pst.h"
#include "mdl.h"
#include "parameters.h"
#ifndef HAVE_CONFIG_H
#include "floattype.h"
#endif

#define MSR_INIT_E		1
#define MSR_STEP_E		0

typedef struct msrContext {
    PRM prm;
    PST pst;
    MDL mdl;
    LCL lcl;
    FLOAT fCenter[3];
    /*
    ** Parameters.
    */
    struct parameters param;
    /*
    ** Other stuff...
    */
    int nThreads;
    uint64_t N;
    uint64_t nDark;
    uint64_t nGas;
    uint64_t nStar;
    uint64_t nMaxOrder;		/* Order number of last particle */
    uint64_t nMaxOrderGas;
    uint64_t nMaxOrderDark;
    int iCurrMaxRung;
    double dCrit;
    /*
    ** Comoving coordinate variables.
    */
    double dEcosmo;
    double dUOld;
    double dTimeOld;
    /*
    ** Redshift output points.
    */
    int nMaxOuts;
    int nOuts;
    double *pdOutTime;
    int iOut;
    /*
    ** Processor mapping for one-node-output functions.
    */
    int *pMap;

    uint8_t iRungVeryActive;    /* NOTE: The first very active particle is at iRungVeryActive + 1 */

    /*
     * Domain Decomposition Done
     */
    uint64_t *nRung;
    int iLastRungRT;
    uint64_t nActive;
    int nGroups;
    int nBins;
    int bAntiGrav;

    int bSavePending;

#ifdef PLANETS
    int nPlanets; /* currently not used */
    double dEcoll;
    double dSunMass;
#endif
    } * MSR;

void msrInitialize(MSR *,MDL,int,char **);
void msrLogParams(MSR msr, FILE *fp);
void msrprintf(MSR msr, const char *Format, ... );
int msrGetLock(MSR msr);
int msrCheckForStop(MSR msr);
void msrFinish(MSR);
double msrGenerateIC(MSR);
double msrRead(MSR msr,const char *achInFile);
void msrWrite(MSR,const char *,double, int bCheckpoint );
void msrSetSoft(MSR msr,double);
void msrDomainDecomp(MSR,int iRung,int bGreater,int bSplitVA);
void msrBuildTree(MSR msr,double dTime,int bNeedEwald);
void msrBuildTreeExcludeVeryActive(MSR msr,double dTime);
void msrCalcBound(MSR msr,BND *pbnd);
void msrDomainColor(MSR);
void msrReorder(MSR);
void msrOutArray(MSR,const char *,int);
void msrOutVector(MSR,const char *,int);
void msrSmoothSetSMF(MSR msr, SMF *smf, double dTime);
void msrSmooth(MSR,double,int,int);
void msrReSmooth(MSR,double,int,int);
void msrUpdateSoft(MSR,double);
void msrGravity(MSR msr,uint8_t uRungLo, uint8_t uRungHi, double dTime,
		double dStep,int bEwald,int *piSec,uint64_t *pnActive);
void msrCalcEandL(MSR,int,double,double *,double *,double *,double *,double *);
void msrDrift(MSR,double dTime,double dDelta,uint8_t uRungLo,uint8_t uRungHi);
void msrKick(MSR,double dTime,double dDelta,uint8_t uRungLo,uint8_t uRungHi);
double msrReadCheck(MSR,int *);
void msrWriteCheck(MSR,double,int);
int msrOutTime(MSR,double);
void msrReadOuts(MSR,double);
void msrTopStepKDK(MSR msr,
		   double dStep,	/* Current step */
		   double dTime,	/* Current time */
		   double dDelta,	/* Time step */
		   int iRung,		/* Rung level */
		   int iKickRung,	/* Gravity on all rungs from iRung
					   to iKickRung */
		   int iRungVeryActive, /* rung *below which* very active particles are */
		   int iAdjust,		/* Do an adjust? */
		   double *pdActiveSum,
		   int *piSec);
void msrStepVeryActiveKDK(MSR msr, double dStep, double dTime, double dDelta,
			  int iRung);
#ifdef HERMITE
void msrTopStepHermite(MSR msr,
		       double dStep,	/* Current step */
		       double dTime,	/* Current time */
		       double dDelta,	/* Time step */
		       int iRung,		/* Rung level */
		       int iKickRung,	/* Gravity on all rungs from iRung
					   to iKickRung */
		       int iRungVeryActive,  /* current setting for iRungVeryActive */
		       int iAdjust,		/* Do an adjust? */
		       double *pdActiveSum,
		       int *piSec);
void msrStepVeryActiveHermite(MSR msr,double dStep,double dTime,double dDelta,
			      int iRung);
void msrCopy0(MSR msr,double dTime);
void msrPredictor(MSR msr,double dTime);
void msrCorrector(MSR msr,double dTime);
void msrSunCorrector(MSR msr,double dTime);
void msrPredictorInactive(MSR msr,double dTime);
void msrAarsethStep(MSR msr);
void msrFirstDt(MSR msr);
#endif /* HERMITE */

void msrBallMax(MSR msr, int iRung, int bGreater);
/*------------------*/
/* Active Functions */
/*------------------*/
void msrActiveRung(MSR msr, int iRung, int bGreater);
void msrActiveOrder(MSR msr);

/* Replacement functions */
void msrActiveMaskRung(MSR msr, unsigned int iSetMask, int iRung, int bGreater);
/*------------------*/
/* Active Functions */
/*------------------*/

void msrVelocityRung(MSR msr,int iRung,double dDelta,double dTime,int bAll);
uint64_t msrCalcWriteStart(MSR);
void msrGetNParts(MSR msr);
void msrAddDelParticles(MSR msr);
void msrGravStep(MSR msr, double dTime);
void msrAccelStep(MSR msr,uint8_t uRungLo,uint8_t uRungHi,double dTime);
void msrDensityStep(MSR msr,uint8_t uRungLo,uint8_t uRungHi,double dTime);
int msrUpdateRung(MSR msr, uint8_t uRung);

/*
** Interface functions.
*/
int msrSteps(MSR);
void msrOutput(MSR msr, int iStep, double dTime, int bCheckpoint);
char *msrOutName(MSR);
char *msrBuildName(MSR msr,char *achFile,int iStep);
char *msrBuildIoName(MSR msr,char *achFile,int iStep);
double msrDelta(MSR);
int msrLogInterval(MSR);
int msrCheckInterval(MSR);
const char *msrCheckTypes(MSR msr);
int msrOutInterval(MSR);
const char *msrOutTypes(MSR msr);
int msrRestart(MSR);
int msrComove(MSR);
double msrSoft(MSR);
int msrDoDensity(MSR);
#ifdef USE_PNG
int msrPNGResolution(MSR msr);
#endif
int msrDoGravity(MSR msr);
void msrInitStep(MSR msr);
void msrSetRung(MSR msr, uint8_t uRungLo, uint8_t uRungHi, int uRung);
int msrMaxRung(MSR msr);

void msrSwitchTheta(MSR msr,double);
uint64_t msrMaxOrder(MSR msr);

void msrFof(MSR msr, double exp);
void msrGroupMerge(MSR msr, double exp);
void msrGroupProfiles(MSR msr, double exp);
void msrOutGroups(MSR msr,const char *,int,double dTime);
void msrDeleteGroups(MSR msr);
void msrInitRelaxation(MSR msr);
void msrRelaxation(MSR msr,double dTime,double deltaT,int iSmoothType,int bSymmetric);
/* Gas routines */
void msrInitSph(MSR,double);
void msrSph(MSR msr, double dTime, double dStep);
void msrSphStep(MSR msr,uint8_t uRungLo,uint8_t uRungHi,double dTime);
void msrCoolSetup(MSR msr, double);
void msrCooling(MSR msr,double dTime,double dStep,int bUpdateState, int bUpdateTable,int bInterateDt);
void msrStarForm( MSR, double);
/* END Gas routines */

#ifdef PLANETS
double msrReadSS(MSR msr);
void msrWriteSS(MSR msr, char *pszFileName, double dTime);
void msrGravSun(MSR msr);
static char * _msrParticleLabel(MSR msr,int iColor);
void msrDoCollision(MSR msr,double dTime,double dDelta);
#ifdef SYMBA
void msrTopStepSymba(MSR msr,
		     double dStep,	/* Current step */
		     double dTime,	/* Current time */
		     double dDelta,	/* Time step */
		     int iRung,		/* Rung level */
		     int iKickRung,	/* Gravity on all rungs from iRung
					   to iKickRung */
		     int iRungVeryActive,
		     int iAdjust,		/* Do an adjust? */
		     double *pdActiveSum,
		     int *piSec);
void msrStepVeryActiveSymba(MSR msr, double dStep, double dTime, double dDelta,
			    int iRung);
void msrDrminToRung(MSR msr,int iRung);
void msrDriftSun(MSR msr,double dTime,double dDelta);
void msrKeplerDrift(MSR msr,double dDelta);
#endif  /* SYMBA */
#endif /* PLANETS */

void msrHostname(MSR msr);
void msrMemStatus(MSR msr);


void msrSelSrcAll(MSR msr);
void msrSelDstAll(MSR msr);
void msrSelSrcGas(MSR msr);
void msrSelDstGas(MSR msr);
void msrSelSrcStar(MSR msr);
void msrSelDstStar(MSR msr);
uint64_t msrSelSrcMass(MSR msr,double dMinMass,double dMaxMass,int setIfTrue,int ClearIfFalse);
uint64_t msrSelDstMass(MSR msr,double dMinMass,double dMaxMass,int setIfTrue,int ClearIfFalse);
uint64_t msrSelSrcById(MSR msr,uint64_t idStart,uint64_t idEnd,int setIfTrue,int clearIfFalse);
uint64_t msrSelDstById(MSR msr,uint64_t idStart,uint64_t idEnd,int setIfTrue,int clearIfFalse);
uint64_t msrSelSrcPhaseDensity(MSR msr,double dMinPhaseDensity,double dMaxPhaseDensity,int setIfTrue,int clearIfFalse);
uint64_t msrSelDstPhaseDensity(MSR msr,double dMinPhaseDensity,double dMaxPhaseDensity,int setIfTrue,int clearIfFalse);
uint64_t msrSelSrcBox(MSR msr,double *dCenter, double *dSize,int setIfTrue,int clearIfFalse);
uint64_t msrSelDstBox(MSR msr,double *dCenter, double *dSize,int setIfTrue,int clearIfFalse);
uint64_t msrSelSrcSphere(MSR msr,double *r, double dRadius,int setIfTrue,int clearIfFalse);
uint64_t msrSelDstSphere(MSR msr,double *r, double dRadius,int setIfTrue,int clearIfFalse);
uint64_t msrSelSrcCylinder(MSR msr,double *dP1, double *dP2, double dRadius,
		      int setIfTrue, int clearIfFalse );
uint64_t msrSelDstCylinder(MSR msr,double *dP1, double *dP2, double dRadius,
		      int setIfTrue, int clearIfFalse );

void msrDeepestPot(MSR msr,double *r, float *fPot);
double msrTotalMass(MSR msr);
void msrProfile(
    MSR msr, const PROFILEBIN **pBins, int *pnBins, double *r,
    double dMinRadius, double dLogRadius, double dMaxRadius,
    int nPerBin, int nBins, int nAccuracy );
void msrDeleteProfile(MSR msr);

void msrCalcCOM(MSR msr,const double *dCenter, double dRadius,
		double *com, double *vcm, double *L, double *M);
void msrPeakVc(MSR msr,int N,struct inPeakVc *in);

void msrInitGrid(MSR msr,int x,int y,int z);
void msrGridProject(MSR msr,double x,double y,double z);
#ifdef MDL_FFTW
void msrMeasurePk(MSR msr,double *dCenter,double dRadius,int nGrid,float *Pk);
#endif

#endif
