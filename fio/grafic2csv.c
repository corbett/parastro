#include "fio.h"
#include <stdio.h>
#include <stdlib.h>
/*
 * Script takes in a directory containing grafic files, and outputs its contents to 
 * comma separated value format: type,id,x,y,z,vx,vy,vz,mass,softening,potential,
 * one particle per line.
 *
 * Author: Christine Corbett MOran
 * License: BSD
 */
int main(int argc, char **argv)
{
	if (argc != 2){
		fprintf(stderr,"usage: %s <filename>\n", argv[0]);
		return 1;
	}
	FIO grafic = \
		fioOpen(argv[1], 0.01, 0.01);
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
    
  uint64_t i;
  // particle variables
  uint64_t piOrder;
  double pdPos[3],pdVel[3];
  float pfMass,pfSoft,pfPot,pfRho,pfTemp,pfMetals,pfTform;
	// header
	fprintf(stdout,"type,id,x,y,z,vx,vy,vz,mass,softening,potential\n");
	// read star, 
  for(i=0; i<nStar; i++) {
    fioSeek(grafic,i,FIO_SPECIES_STAR);    
    fioReadStar(grafic,
								&piOrder,pdPos,pdVel,&pfMass,&pfSoft,&pfPot,&pfMetals,&pfTform);
		fprintf(stdout,"%d,%llu,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
						FIO_SPECIES_STAR,piOrder,pdPos[0],pdPos[1],pdPos[2],pdVel[0],pdVel[1],pdVel[2],pfMass,pfSoft,pfPot);
  }
  for(i=0; i<nDark; i++) {
    fioSeek(grafic,i,FIO_SPECIES_DARK);
    fioReadDark(grafic,
      &piOrder,pdPos,pdVel,&pfMass,&pfSoft,&pfPot);
    fprintf(stdout,"%d,%llu,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
						FIO_SPECIES_DARK,piOrder,pdPos[0],pdPos[1],pdPos[2],pdVel[0],pdVel[1],pdVel[2],pfMass,pfSoft,pfPot);
      
  }
  for(i=0; i<nGas; i++) {
    fioSeek(grafic,i,FIO_SPECIES_SPH);
    fioReadSph(grafic,
							 &piOrder,pdPos,pdVel,&pfMass,&pfSoft,&pfPot,&pfRho,&pfTemp,&pfMetals);
    fprintf(stdout,"%d,%llu,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
						FIO_SPECIES_SPH,piOrder,pdPos[0],pdPos[1],pdPos[2],pdVel[0],pdVel[1],pdVel[2],pfMass,pfSoft,pfPot);
  }

	return 0;
}