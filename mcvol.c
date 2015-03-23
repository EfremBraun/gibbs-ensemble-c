#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "system.h"
#include "ran_uniform.h"
#include "conf.h"  
#include "potential.h"  

int Mcvol(double En[2], double Vir[2], int *Attempt, int *Acc, double Vmax)
{
  //  Attempts To Change The Volume

  double Enn[2], Virn[2], Yy; 
  double F[2], Arg, Volo[2], Volt, Dlnv; 
  double Voln[2], Dele1, Dele2, Dlnv1, Dlnv2; 
  double Enold;
  int I, BoxID, Idi;
  
  (*Attempt)++;
  
  // ---Calulate New Volume By Making Random Walk In Ln V
  
   Volo[0] = pow(Box[0],3);
   Volo[1] = pow(Box[1],3);
   Volt = Volo[0] + Volo[1];
   Dlnv = log(Volo[0]/Volo[1]) + (RandomNumber()-0.5)*Vmax;
   Voln[0] = exp(Dlnv)*Volt/(1.0+exp(Dlnv)); // Note to Efrem: I don't get this calculation, but I won't be making moves in ln space, so it doesn't matter much to me.
   Voln[1] = Volt - Voln[0];

   for(BoxID=0;BoxID<2;BoxID++)
     {
       Box[BoxID] = pow(Voln[BoxID],(1.0/3.0));
       F[BoxID] = Box[BoxID]/pow(Volo[BoxID],(1.0/3.0));
       Hbox[BoxID] = Box[BoxID]/2.0;
       Rc[BoxID] = F[BoxID]*Rc[BoxID];  // Note to Efrem: I see no reason why the r_cutoff should be scaled. This just looks wrong to me. Fix this before using this code for comparing my modified Gibbs Ensemble code to.
       Rc2[BoxID] = pow(Rc[BoxID],2.0); // Note to Efrem: I see no reason why the r_cutoff should be scaled. This just looks wrong to me. Fix this before using this code for comparing my modified Gibbs Ensemble code to.
     }
   
   // ---Determine New Coordinates
   for(I=0;I<Npart;I++)
     {
       Idi = Id[I];
       X[I] = F[Idi]*X[I];
       Y[I] = F[Idi]*Y[I];
       Z[I] = F[Idi]*Z[I];
     }
   
   //    ---Calculate New Energy Using Scaling // Note to Efrem: I don't know if this works or not, but if it was designed for r_cutoff scaling it shouldn't be done. It certainly won't work for modified Gibbs Ensemble.
   for(BoxID=0;BoxID<2;BoxID++)
     {
       Enold = En[BoxID];
       Yy = pow((Volo[BoxID]/Voln[BoxID]),2);
       Enn[BoxID] = Enold*Yy*(2.0-Yy) - 
	 Vir[BoxID]*Yy*(1.0-Yy)/6.0;
       Virn[BoxID] = -12.0*Enold*Yy*(Yy-1.0) + 
	 Vir[BoxID]*Yy*(2.0*Yy-1.0);
     }
      
   // ---Acceptance:
   Dele1 = Enn[0] - En[0];
   Dele2 = Enn[1] - En[1];
   Dlnv1 = log(Voln[0]/Volo[0]);
   Dlnv2 = log(Voln[1]/Volo[1]);
   Arg = exp(-Beta*(Dele1+Dele2-(double)((Npbox[0]+1))*
		    Dlnv1/Beta-(double)((Npbox[1]+1))
		    *Dlnv2/Beta));
   if(RandomNumber()<Arg)
     { 
       
       //    ---Accepted
       (*Acc)++;
       for(BoxID=0;BoxID<2;BoxID++)
	 {
	   En[BoxID] = Enn[BoxID];	  
	   Vir[BoxID] = Virn[BoxID];
	 }
     }
   else
     {
       //    ---Restore The Old Configuration
       for(BoxID=0;BoxID<2;BoxID++)
	 {
	   F[BoxID] = 1/F[BoxID];
	   Box[BoxID] = Box[BoxID]*F[BoxID];
	   Hbox[BoxID] = 0.5*Box[BoxID];
	   Rc[BoxID] = F[BoxID]*Rc[BoxID]; // Note to Efrem: same deal here as before.
	   Rc2[BoxID] = pow(Rc[BoxID],2); // Note to Efrem: same deal here as before.
	 }
       for(I=0;I<Npart;I++)
	 {
	   Idi = Id[I];
	   X[I] = F[Idi]*X[I];
	   Y[I] = F[Idi]*Y[I];
	   Z[I] = F[Idi]*Z[I];
	 }
     }
   
   return(0);
}
