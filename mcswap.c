#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "system.h"
#include "ran_uniform.h"
#include "conf.h"
#include "chem.h"
#include "system_2.h"

int Mcswap(double En[2], double Vir[2],int *Attempt, int *Acc)
{
  // ---Exchange A Particle Bewteen The Two Boxes
  
  double Xn, Yn, Zn, Enn, Virn, Eno, Viro, Tail=0.0; 
  double Arg, Vola, Vold; 
  double Xo, Yo, Zo, Dele;
  double DelTail[2];
  int IndexAdd, Iadd, Idel, Jb, Idi;
  
  (*Attempt)++;
  
  // ===Select A Box At Random
  
  if(RandomNumber()<0.5)
    {
      Iadd = 0;
      Idel = 1;
   }
  else
   {
     Iadd = 1;
     Idel = 0;
   }
  
  Vola = pow(Box[Iadd],3.0);
  Vold = pow(Box[Idel],3.0);
  
  // ---Add A Particle To Box Iadd
  
  Xn = Box[Iadd]*RandomNumber();
  Yn = Box[Iadd]*RandomNumber();
  Zn = Box[Iadd]*RandomNumber();
  
  // ---Calculate Energy Of This Particle
  
  Jb = 0;
  IndexAdd = Npart;
  
  Eneri(Xn, Yn, Zn, IndexAdd, Jb, &Enn, &Virn, Iadd);
  
  // ---Calculate Contibution To The Chemical Potential:
  if (TruncFlag==0)
  {
    Tail = TailC(Iadd) * (double) (Npbox[Iadd]+1); // Note to Efrem: I don't think this tail contribution to the chemical potential is right. See page 174 of Berend's book.
  }
  else if (TruncFlag==1)
  {
    Tail = 0.0;
  }
  Arg = -Beta*(Enn+Tail);
  Chp[Iadd] = Chp[Iadd] * ((double) Ichp[Iadd] / (double) (Ichp[Iadd]+1)) + Vola*exp(Arg)/(double)(Npbox[Iadd]+1) / (double) (Ichp[Iadd]+1);
  Ichp[Iadd]++;
  
  // ---Delete Particle From Box B:  
  if(Npbox[Idel]==0)
    { 
      return(0);
    }     
  Idi = -1;
  while(Idi!=Idel)
    {
      IndexAdd = (int)((double)(Npart)*RandomNumber());
      Idi = Id[IndexAdd];
    }
  Xo = X[IndexAdd];
  Yo = Y[IndexAdd];
  Zo = Z[IndexAdd];
  Eneri(Xo, Yo, Zo, IndexAdd, Jb, &Eno, &Viro, Idel);
  
  // ---Acceptance Test:
  Dele = Enn - Eno + log(Vold*(double)((Npbox[Iadd]+1))/
			 (Vola*(double)(Npbox[Idel])))/Beta;
  if (TruncFlag==0)
  {
    DelTail[Iadd] = TailC(Iadd) * (((double) (Npbox[Iadd]+1) * (double) (Npbox[Iadd]+1)) - ((double) Npbox[Iadd] * (double) Npbox[Iadd]));
    DelTail[Idel] = TailC(Idel) * (((double) (Npbox[Idel]-1) * (double) (Npbox[Idel]-1)) - ((double) Npbox[Idel] * (double) Npbox[Idel]));
    Dele = Dele + DelTail[Iadd] + DelTail[Idel];
  }
  
  if(RandomNumber()<exp(-Beta*Dele))
    { 
      //    ---Accepted:
      (*Acc)++;
      Npbox[Iadd]++;
      X[IndexAdd] = Xn;
      Y[IndexAdd] = Yn;
      Z[IndexAdd] = Zn;
      Id[IndexAdd] = Iadd;
      En[Iadd]+= Enn;
      Vir[Iadd]+= Virn;
      Npbox[Idel]--;
      En[Idel]-= Eno;
      Vir[Idel]-= Viro;
      if (TruncFlag==0)
      {
        En[Iadd]+= DelTail[Iadd];
        En[Idel]+= DelTail[Idel];
      }
    }
  return(0);
}
