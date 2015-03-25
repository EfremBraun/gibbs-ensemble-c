#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "system.h"
#include "conf.h"
#include "potential.h"

void Toterg(double *Ener,double *Vir,int BoxID)
{
  //     ---Calculates Total Energy
  
  double Xi, Yi, Zi, Eni, Viri, Rc3i, Tail;
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
  Rc3i = pow(Sig/Rc[BoxID], 3.0);
  Tail = 2.0/3.0 * M_PI * Npbox[BoxID] / pow(Box[BoxID], 3.0)  * Eps4 * (pow(Rc3i, 3.0) / 3.0 - Rc3i);
  *Ener = *Ener + Tail;
  return;
}
