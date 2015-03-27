#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "system.h"
#include "ran_uniform.h"
#include "conf.h"  
#include "potential.h"  

int McvolMod(double En[2], double Vir[2], int *Attempt, int *Acc, double Vmax)
{
  //  Attempts To Change The Volume and Number of Particles of Box 0, keeping the total density constant and not touching Box 1

  double Ran;
  double Enn[2], Virn[2];
  double F[2], Arg, Volo[2];
  double Voln[2], Dele, Dlnv;
  double XDel, YDel, ZDel;
  int I, Idi, IndexDel;
  int BoxID=0;
  int AddFlag;
  
  (*Attempt)++;
  
  // ---Calulate New Volume
  Volo[BoxID] = pow(Box[BoxID],3);
  Ran = RandomNumber();
  if (Ran < 0.5) AddFlag=1;
  else           AddFlag=-1;
  if (AddFlag==1) Voln[BoxID] = Volo[BoxID] + DvMod;
  else
  {
    // Check that there's a particle in the box. If not, return.
    if(Npbox[BoxID]==0)
      {
        return(0);
      }
    Voln[BoxID] = Volo[BoxID] - DvMod;
  }

  // ---Calculate new box values
  Box[BoxID] = pow(Voln[BoxID],(1.0/3.0));
  F[BoxID] = Box[BoxID]/pow(Volo[BoxID],(1.0/3.0));
  Hbox[BoxID] = Box[BoxID]/2.0;
  if (Rc[BoxID] >= Hbox[BoxID])
  {
    printf("Half of box length (%lf) of box %d has become smaller than Rc (%lf)\n", Hbox[BoxID], BoxID+1, Rc[BoxID]);
    exit(1);
  }
  
  // ---Determine New Coordinates of particles in box 0
  for(I=0;I<Npart;I++)
    {
      Idi = Id[I];
      if (Idi == BoxID)
      {
        X[I] = F[Idi]*X[I];
        Y[I] = F[Idi]*Y[I];
        Z[I] = F[Idi]*Z[I];
      }
    }

  // ---Add a Particle
  if (AddFlag==1)
  {
    X[Npart] = Box[BoxID]*RandomNumber();
    Y[Npart] = Box[BoxID]*RandomNumber();
    Z[Npart] = Box[BoxID]*RandomNumber();
    Id[Npart] = BoxID;
    Npbox[BoxID]++;
    Npart++;
    if(Npart>Npmax)
    {
      printf("Error: Number Of Particles Too Large\n");
      exit(1);
    }
  }
  // --- Or, remove a particle
  else
  {
    Idi = -1;
    while(Idi!=BoxID)
      {
        IndexDel = (int)((double)(Npart)*RandomNumber());
        Idi = Id[IndexDel];
      }
    XDel = X[IndexDel];
    YDel = Y[IndexDel];
    ZDel = Z[IndexDel];
    for(I=IndexDel;I<Npart;I++)
    {
      X[I] = X[I+1];
      Y[I] = Y[I+1];
      Z[I] = Z[I+1];
      Id[I] = Id[I+1];
    }
    Npbox[BoxID]--;
    Npart--;
  }

  //    ---Calculate New Energy
  Toterg(&(Enn[BoxID]),&(Virn[BoxID]),BoxID);
     
  // ---Acceptance:
  Dele = Enn[BoxID] - En[BoxID];
  Dlnv = log(Voln[BoxID]) / (double) Npbox[BoxID]- log(Volo[BoxID]) / ((double) (Npbox[BoxID] - AddFlag));
  Arg = exp(-Beta * Dele - Dlnv) / (double) Npbox[BoxID]; // Note to Efrem: I'm excluding the de Broglie wavelength here because I don't know what to do with it
                   
  if(RandomNumber()<Arg)
  { 
    
    //    ---Accepted
    (*Acc)++;
    En[BoxID] = Enn[BoxID];        
    Vir[BoxID] = Virn[BoxID];
  }
  else
  {
    // ---Undo Adding a Particle
    if (AddFlag==1)
    {
      Npbox[BoxID]--;
      Npart--;
    }
    // --- Or, undo removing a particle
    else
    {
      for(I=(Npart+1);I>IndexDel;I--)
      {
        X[I] = X[I-1];
        Y[I] = Y[I-1];
        Z[I] = Z[I-1];
        Id[I] = Id[I-1];
      }
      X[IndexDel] = XDel;
      Y[IndexDel] = YDel;
      Z[IndexDel] = ZDel;
      Id[IndexDel] = BoxID;
      Npbox[BoxID]++;
      Npart++;
    }

    //    ---Restore The Old Configuration
    F[BoxID] = 1/F[BoxID];
    Box[BoxID] = Box[BoxID]*F[BoxID];
    Hbox[BoxID] = 0.5*Box[BoxID];
    for(I=0;I<Npart;I++)
    {
      Idi = Id[I];
      if (Idi == BoxID)
      {
        X[I] = F[Idi]*X[I];
        Y[I] = F[Idi]*Y[I];
        Z[I] = F[Idi]*Z[I];
      }
    }
  }
  
  return(0);
}
