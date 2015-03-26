#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "system.h"
#include "conf.h"
#include "potential.h"
#include "system_2.h"

void Toterg(double *Ener,double *Vir,int BoxID)
{
  //     ---Calculates Total Energy
  
  double Xi, Yi, Zi, Eni, Viri, Tail;
  int I, Jb;
  
  *Ener = 0.0;
  *Vir = 0.0;
  for(I = 0;I<Npart-1;I++)
  {
    if(Id[I]==BoxID) 
    {
      Xi = X[I];
      Yi = Y[I];
      Zi = Z[I];
      Jb = I + 1;
      Eneri(Xi, Yi, Zi, I, Jb, &Eni,&Viri, BoxID);
      *Ener = *Ener + Eni;
      *Vir = *Vir + Viri;
    }
  }

  // Add tail corrections
  if (TruncFlag==0)
  {
    Tail = TailC(BoxID) * (double) Npbox[BoxID] * (double) Npbox[BoxID];
    //printf("Tail is %f\n", Tail / (double) Npbox[BoxID]); // For testing purposes; can compare to value on page 37 of Frenkel&Smit
    *Ener = *Ener + Tail;
    // Note to Efrem: need to add tail corrections to the pressure too, but only important for printing
  }

  return;
}
