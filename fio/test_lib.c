#include "fio.h"
#include <stdio.h>
#include <stdlib.h>
int main(int narg, char **args)
{
	fprintf(stdout,"testing fio fioOpen\n");
	FIO grafic = \
	  fioOpen("/Users/corbett/Local/BigData/GraficIcs/ic_files", 0.01, 0.01);
  double dExpansion,dEcosmo,dTimeOld,dUOld;
  uint64_t nTot,nGas,nDark,nStar;
  
  if(!fioGetAttr(grafic,"dTime",FIO_TYPE_DOUBLE,&dExpansion)) dExpansion=0.0;
  
  if (!fioGetAttr(grafic,"dEcosmo",FIO_TYPE_DOUBLE,&dEcosmo)) dEcosmo = 0.0;
  if (!fioGetAttr(grafic,"dTimeOld",FIO_TYPE_DOUBLE,&dTimeOld)) dTimeOld = 0.0;
  if (!fioGetAttr(grafic,"dUOld",FIO_TYPE_DOUBLE,&dUOld)) dUOld = 0.0;
  nTot = fioGetN(grafic,FIO_SPECIES_ALL);
  nGas  = fioGetN(grafic,FIO_SPECIES_SPH);
  nDark = fioGetN(grafic,FIO_SPECIES_DARK);
  nStar = fioGetN(grafic,FIO_SPECIES_STAR);
  
  
  fprintf(stdout,"dExpansion: %f\n",dExpansion);
  fprintf(stdout,"dEcosmo: %f\n",dEcosmo);
  fprintf(stdout,"dTimeOld: %f\n",dTimeOld);
  fprintf(stdout,"dUOld: %f\n",dUOld);
  fprintf(stdout,"nTot: %llu\n",nTot);
  fprintf(stdout,"nGas: %llu\n",nGas);
  fprintf(stdout,"nDark: %llu\n",nDark);
  fprintf(stdout,"nStar: %llu\n",nStar);
  
  
  
  // reserve enough data for vx, vy, and vz!!!
  // WHEEEEEEEEEEEEEEEE
  double *ic_velcx = malloc(nTot*sizeof(double));
  double *ic_velcy = malloc(nTot*sizeof(double));
  double *ic_velcz = malloc(nTot*sizeof(double));
  
  fioGetAttr(grafic,"ic_velcx",FIO_TYPE_DOUBLE,ic_velcx);
  fioGetAttr(grafic,"ic_velcy",FIO_TYPE_DOUBLE,ic_velcy);
  fioGetAttr(grafic,"ic_velcz",FIO_TYPE_DOUBLE,ic_velcz);
    
  uint64_t i;
  // particle variables
  uint64_t piOrder;
  double pdPos[3];
  double pdVel[3];
  float pfMass;
  float pfSoft;
  float pfPot;
  /*

  for(i=0; i<nStar; i++) {
    fioSeek(grafic,i,FIO_SPECIES_STAR);
    fioReadDark(grafic,&piOrder,pdPos,pdVel,&pfMass,&pfSoft,&pfPot);
    // all the variables we just read!
    // these will be stored in the array!
    fprintf(stdout,"%llu\n",piOrder);
    fprintf(stdout,"%f\t%f\t%f\n",pdPos[0],pdPos[1],pdPos[2]);
    fprintf(stdout,"%f\t%f\t%f\n",pdVel[0],pdVel[1],pdVel[2]);
    fprintf(stdout,"%f\t%f\t%f\n",pfMass,pfSoft,pfPot);
  }
  for(i=0; i<nDark; i++) {
    fioSeek(grafic,i,FIO_SPECIES_DARK);
    fioReadDark(grafic,&piOrder,pdPos,pdVel,&pfMass,&pfSoft,&pfPot);
  }
  */
  for(i=nDark; i<nDark+nGas; i++) {
    fioSeek(grafic,i,FIO_SPECIES_SPH);
    fioReadDark(grafic,&piOrder,pdPos,pdVel,&pfMass,&pfSoft,&pfPot);
    // all the variables we just read!
    // these will be stored in the array!
    fprintf(stdout,"%llu\n",piOrder);
    fprintf(stdout,"%f\t%f\t%f\n",pdPos[0],pdPos[1],pdPos[2]);
    fprintf(stdout,"%f\t%f\t%f\n",pdVel[0],pdVel[1],pdVel[2]);
    fprintf(stdout,"%f\t%f\t%f\n",pfMass,pfSoft,pfPot);
  }


  // memory management, natch
  
	return 0;
}