#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "system.h"
#include "potential.h"
#include "conf.h"

void Readdat(int *Melt, int *Equil,int *Prod,int *Nsamp,int *Nprint, int *Ndispl,double *Dr,int *Nvol,double *Vmax,int *Nswap,double *Succ)
{
  // Reads Input Data And Model Parameters
  //
  // ---Input Parameters: File: input.settings
  // BoxInitFlag       : 0 to Initilaize From A Lattice, 1 to Read Configuration From Disk
  // Melt          : Number Of Monte Carlo Cycles During Melting
  // Equil         : Number Of Monte Carlo Cycles During Equilibration
  // Prod          : Number Of Monte Carlo Cycles During Production
  // Nsamp         : Number Of Monte Carlo Cycles Between Two Sampling Periods
  // Nprint        : Number Of Monte Carlo Cycles Between Printing of Movie PDBs
  // Dr            : Maximum Displacement
  // Vmax          : Maximum Volume Change
  // Succ          : Optimal Percentance Of Accepted Attemps
  //                 The Program Adjusts Vmax Or Dr In Just A Way That
  //                 On Average Succ% Of The Moves Are Accepted
  // Ndispl        : Number Of Attemps To Displace A Particle Per Mc Cycle
  // Nvol          : Number Of Attemps To Change The Volume  Per Mc Cycle
  // Nswap         : Number Of Attemps To Swap Particle Between The Two Boxes Per Mc Cycle
  // Npart         : Total Number of Particles
  // Temp          : Temperature
  // TempMelt      : Temperature to run melting cycles
  // Rho           : Density
  // TruncFlag     : 0 to do truncated with tail corrections, 1 to do truncated and shifted
  // ModGibbsFlag  : 0 to do normal Gibbs Ensemble, 1 to do modified Gibbs Ensemble
  
  //  ---Input Parameters: File: input.lj.model
  // Eps    = Epsilon Lennard-Jones Potential
  // Sig    = Sigma Lennard-Jones Potential
  // Mass   = Mass Of The Particle
  // Rcc    = Cut-Off Radius Of The Potential
  // Lambda = thermal de Broglie wavelength
  
  //  ---Input Parameters: File: output.restart (Restart File
  //             To Continue A Simulation From Disk)
  // Box[0]   = Length Box 0 Old Configuration
  // Hbox[0]  = Box[0]/2
  // Box[1]   = Length Box 1 Old Configuration
  // Hbox[1]  = Box[1]/2
  // Npart    = Total Number Of Particles (overrules input.settings!!)
  // Npbox[0] = Number Of Particles In Box 0
  // Npbox[1] = Number Of Particles In Box 1
  // Dr     = Optimized Maximum Displacement Old Configurations
  // Vmax   = Optimized Maximum Volume Change Old Configurations
  // X[0],Y[0],Z[0]            : Position First Particle
  //               ,Id(i)   = 0 if Particle i In Box 0
  //               ,Id(i)   = 1 if Particle i In Box 1
  //    ....
  // X(Npart),Y(Npart),Z(Npart): Position Particle Last Particle

  int BoxInitFlag, I,BoxID;
  double Eps, Rcc;
  double Vol0, Vol1;
  double TempMelt;
  FILE* fileptr;
  char line[500];
  
  fileptr=fopen("input.settings","r");
  fgets(line,300,fileptr);
  fgets(line,300,fileptr);
  sscanf(line,"%d %d %d %d %d",Melt,Equil,Prod,Nsamp,Nprint);
  fgets(line,300,fileptr);
  fgets(line,300,fileptr);
  sscanf(line,"%lf %lf %lf",Dr,Vmax,Succ);
  fgets(line,300,fileptr);
  fgets(line,300,fileptr);
  sscanf(line,"%d %d %d",Ndispl,Nvol,Nswap);
  fgets(line,300,fileptr);
  fgets(line,300,fileptr);
  sscanf(line,"%d %lf %d %lf",&Npbox[0],&Vol0,&Npbox[1],&Vol1);
  fgets(line,300,fileptr);
  fgets(line,300,fileptr);
  sscanf(line,"%lf %lf",&Temp,&TempMelt);
  fgets(line,300,fileptr);
  fgets(line,300,fileptr);
  sscanf(line,"%d",&BoxInitFlag);
  fgets(line,300,fileptr);
  fgets(line,300,fileptr);
  sscanf(line,"%d",&TruncFlag);
  fgets(line,300,fileptr);
  fgets(line,300,fileptr);
  sscanf(line,"%d",&ModGibbsFlag);
  fclose(fileptr);

  Npart = Npbox[0] + Npbox[1];

  if(Npart>Npmax)
    {
      printf("Error: Number Of Particles Too Large\n");
      exit(0);
    }
  
  //  ---Read Model Parameters
  fileptr=fopen("input.lj.model","r");
  fgets(line,300,fileptr);
  fgets(line,300,fileptr);
  sscanf(line,"%lf %lf %lf %lf %lf",&Eps,&Sig,&Mass,&Rcc,&Lambda);
  fclose(fileptr);

  //  ---Read/Generate Configuration
  Box[0]  = pow(Vol0,(1.0/3.0));
  Hbox[0] = 0.5*Box[0];
  Box[1]  = pow(Vol1,(1.0/3.0));
  Hbox[1] = Hbox[0];
  
  if(BoxInitFlag==0)
    {
      //     ---Generate Configuration Form Lattice
      Lattice();
    }
  else
    {
      printf("Read Conf From Disk\n");
      fileptr=fopen("output.restart","r");
      fscanf(fileptr,"%lf %lf %lf %lf",&(Box[0]), &(Hbox[0]), &(Box[1]), &(Hbox[1]));
      fscanf(fileptr,"%d %d %d",&Npart, &(Npbox[0]), &(Npbox[1]));
      fscanf(fileptr,"%lf %lf",Dr,Vmax);
      for(I = 0;I< Npart;I++)
	{
	  fscanf(fileptr,"%lf %lf %lf %d",&(X[I]),&(Y[I]),&(Z[I]),&(Id[I]));
	}
      //         Rewind (11)
    }
  
  //  ---Write Input Data
  printf("Number Of Melting Cycles                   : %d\n", *Melt);
  printf("Number Of Equilibration Cycles             : %d\n", *Equil);
  printf("Number Of Production Cycles                : %d\n", *Prod);
  printf("Sample Frequency                           : %d\n", *Nsamp);
  printf("Movie Printing Frequency                   : %d\n", *Nprint);
  printf("Number Of Att. To Displ. a Part. per Cycle : %d\n", *Ndispl);
  printf("Maximum Displacement                       : %lf\n" ,*Dr);
  printf("Number Of Att. to Change Volume  per Cycle : %d\n", *Nvol);
  printf("Maximum Change Volume                      : %lf\n", *Vmax);
  printf("Number Of Att. to Exch Part.per Cycle      : %d\n", *Nswap);
  printf("Number Of Particles in Box 0: %d\n",Npbox[0]);  
  printf("Volume of Box 0: %lf\n",Vol0);  
  printf("Box 0 Length: %lf\n",Box[0]);  
  printf("Density Box 0: %lf\n",(double)(Npbox[0])/pow(Box[0],3.0));  
  printf("Number Of Particles in Box 1: %d\n",Npbox[1]);  
  printf("Volume of Box 1: %lf\n",Vol1);  
  printf("Box 1 Length: %lf\n",Box[1]);  
  printf("Density Box 1: %lf\n",(double)(Npbox[1])/pow(Box[1],3.0));  
  printf("Temperature: %lf\n",Temp);  
  printf("Melting Temperature: %lf\n",TempMelt);  
  printf("Pressure: %lf\n",0.0);  
  if (TruncFlag==0) printf("Potential is: Truncated and Tail-Corrected\n");  
  else if (TruncFlag==1) printf("Potential is: Truncated and Shifted\n");  
  printf("Model Parameters\n");
  printf("Epsilon: %lf\n",Eps);  
  printf("Sigma: %lf\n",Sig);
  printf("Mass: %lf\n",Mass);
  printf("Lambda: %lf\n",Lambda);
  
  //  ---Calculate Parameters:
  Beta  = 1.0/Temp;
  BetaMelt  = 1.0/TempMelt;
  Eps4  = 4.0*Eps;
  Eps48 = 48.*Eps;
  Sig2  = Sig*Sig;
  DvMod = (Vol0+Vol1) / (double) Npart;
  
  //  ---Calculate Cut-Off Radius Potential
  for(BoxID = 0;BoxID< 2;BoxID++)
    {
      if(Rcc<Hbox[BoxID])
	{
	  Rc[BoxID] = Rcc;
	}
      else
	{
          printf("Half of both length (%lf) of box %d is smaller than Rcc (%lf)\n", Hbox[BoxID], BoxID, Rcc);
          exit(1);
	}
      Rc2[BoxID] = Rc[BoxID]*Rc[BoxID];
    }
  
  return;
}
 
