#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "system.h"
#include "conf.h"   

int Lattice(void)
{
  // ---Place `Npbox[0]' Particles On Lattice With Volume Box[0]^3
  // ---Place `Npbox[1]' Particles On Lattice With Volume Box[1]^3
  
  int I, J, K, Itel, N, BoxID;
  double Del;

  printf("Generate Simple Cubic Lattice\n");

  for(BoxID=0; BoxID<2;BoxID++)
  {
    Del = Box[BoxID];
    N = (int)(pow(Npbox[BoxID],(1.0/3.0)) + 1);
    Del = Del/(double)N;
    Itel = 0;
  
    for(I=0;I<N;I++)
    {
      for(J=0;J<N;J++)
      {
        for(K=0;K<N;K++)
        {
          if(Itel<Npbox[BoxID])
          {
            if(BoxID==0)
            {
              X[Itel] = (double)K*Del;
              Y[Itel] = (double)J*Del;
              Z[Itel] = (double)I*Del;
              Id[Itel] = BoxID;
              Itel++;
            }
            else if(BoxID==1)
            {
              X[Itel+Npbox[0]] = (double)K*Del;
              Y[Itel+Npbox[0]] = (double)J*Del;
              Z[Itel+Npbox[0]] = (double)I*Del;
              Id[Itel+Npbox[0]] = BoxID;
              Itel++;
            }
          }
        }
      }
    }
    
    printf("Initialization on Lattice: \n"); 
    printf("%d Particles Placed on Lattice %d\n" ,Itel, BoxID); 
  }
  
  return(0);
}
