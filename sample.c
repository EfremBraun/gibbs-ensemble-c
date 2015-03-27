#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "conf.h"
#include "system.h"

void Sample(int I,double *En,double *Vir, FILE *FilePtr)
{
  //C     Writes Quantities To File
  int BoxID;
  double Enp[2], Press[2], Vol[2], Rho[2];
  
  for(BoxID = 0; BoxID<2;BoxID++)
    {
      Vol[BoxID] = pow(Box[BoxID],3);
      Rho[BoxID] = (double)(Npbox[BoxID])/Vol[BoxID];
      Press[BoxID] = Rho[BoxID]/Beta + Vir[BoxID]/(3.0*Vol[BoxID]);
         
      if (Npbox[BoxID]!=0)
	{
	  Enp[BoxID] = En[BoxID]/(double)(Npbox[BoxID]);
	}
      else
	{
	  Enp[BoxID] = 0.;
	}
    }
  
  fprintf(FilePtr,"%7d%11.6lf%11.6lf%11.6lf%11.6lf%11.6lf%11.6lf%7d%7d%11.3lf%11.3lf\n",I, Enp[0],Enp[1],Press[0],Press[1],Rho[0],Rho[1],Npbox[0],Npbox[1],Vol[0],Vol[1]);
  return;
}
