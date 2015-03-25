#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "system.h"
#include "potential.h"
#include "conf.h"


void Readdat(int *Equil,int *Prod,int *Nsamp,int *Nprint, int *Ndispl,double *Dr,int *Nvol,double *Vmax,int *Nswap,double *Succ)
{
  // Reads Input Data And Model Parameters
  //
  // ---Input Parameters: File: input.settings
  // BoxIDeg  =  0 : Initilaize From A Lattice
  //          1 : Read Configuration From Disk
  // Equil      : Number Of Monte Carlo Cycles During Equilibration
  // Prod       : Number Of Monte Carlo Cycles During Production
  // Nsamp      : Number Of Monte Carlo Cycles Between Two Sampling Periods
  // Nprint     : Number Of Monte Carlo Cycles Between Printing of Movie PDBs
  // Dr         : Maximum Displacement
  // Vmax       : Maximum Volume Change
  // Succ       : Optimal Percentance Of Accepted Attemps
  //              The Program Adjusts Vmax Or Dr In Just A Way That
  //              On Average Succ% Of The Moves Are Accepted
  // Ndispl     : Number Of Attemps To Displace A Particle Per Mc Cycle
  // Nvol       : Number Of Attemps To Change The Volume  Per Mc Cycle
  // Nswap      : Number Of Attemps To Swap Particle Between The Two Boxes Per Mc Cycle
  // Npart      : Total Number of Particles
  // Temp       : Temperature
  // Rho        : Density
  
  //  ---Input Parameters: File: input.lj.model
  // Eps    = Epsilon Lennard-Jones Potential
  // Sig    = Sigma Lennard-Jones Potential
  // Mass   = Mass Of The Particle
  // Rcc    = Cut-Off Radius Of The Potential
  
  //  ---Input Parameters: File: output.lj.restart (Restart File
  //             To Continue A Simulation From Disk)
  // Box(1)   = Length Box 1 Old Configuration
  // Hbox(1)  = Box(1)/2
  // Box(2)   = Length Box 2 Old Configuration
  // Hbox(2)  = Box(2)/2
  // Npart    = Total Number Of Particles (overrules input.settings!!)
  // Npbox(1) = Number Of Particles In Box 1
  // Npbox(2) = Number Of Particles In Box 2
  // Dr     = Optimized Maximum Displacement Old Configurations
  // Vmax   = Optimized Maximum Volume Change Old Configurations
  // X(1),Y(1),Z(1)            : Position First Particle 1
  //               ,Id(1)   = 1 Particle In Box 1
  //               ,Id(1)   = 2 Particle In Box 2
  //    ....
  // X(Npart),Y(Npart),Z(Npart): Position Particle Last Particle

  int BoxIDeg, I,BoxID;
  double Eps, Sig, Rho, Rcc;
  FILE* fileptr;
  char line[500];
  
  fileptr=fopen("input.settings","r");
  fgets(line,300,fileptr);
  fgets(line,300,fileptr);
  sscanf(line,"%d %d %d %d %d",&BoxIDeg,Equil,Prod,Nsamp,Nprint);
  fgets(line,300,fileptr);
  fgets(line,300,fileptr);
  sscanf(line,"%lf %lf %lf",Dr,Vmax,Succ);
  fgets(line,300,fileptr);
  fgets(line,300,fileptr);
  sscanf(line,"%d %d %d",Ndispl,Nvol,Nswap);
  fgets(line,300,fileptr);
  fgets(line,300,fileptr);
  sscanf(line,"%d %lf %lf",&Npart,&Temp,&Rho);
  fclose(fileptr);
  if(Npart>Npmax)
    {
      printf("Error: Number Of Particles Too Large\n");
      exit(0);
    }
  
  //  ---Read Model Parameters
  fileptr=fopen("input.lj.model","r");
  fgets(line,300,fileptr);
  fgets(line,300,fileptr);
  sscanf(line,"%lf %lf %lf %lf",&Eps,&Sig,&Mass,&Rcc);      
  fclose(fileptr);

  //  ---Read/Generate Configuration
  
  Box[0]  = pow(((double)(Npart)/(2.*Rho)),(1.0/3.0));
  Hbox[0] = 0.5*Box[0];
  Box[1]  = Box[0];
  Hbox[1] = Hbox[0];
  
  if(BoxIDeg==0)
    {
      
      //     ---Generate Configuration Form Lattice
      Lattice();
    }
  else
    {
      printf("Read Conf From Disk\n");
      fileptr=fopen("output.lj.restart","r");
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
  
  printf("Number Of Equilibration Cycles             : %d\n", *Equil);
  printf("Number Of Production Cycles                : %d\n", *Prod);
  printf("Sample Frequency                           : %d\n", *Nsamp);
  printf("Movie Printing Frequency                   : %d\n", *Nprint);
  printf("Number Of Att. To Displ. a Part. per Cycle : %d\n", *Ndispl);
  printf("Maximum Displacement                       : %lf\n" ,*Dr);
  printf("Number Of Att. to Change Volume  per Cycle : %d\n", *Nvol);
  printf("Maximum Change Volume                      : %lf\n", *Vmax);
  printf("Number Of Att. to Exch Part.per Cycle      : %d\n", *Nswap);
  printf("Number Of Particles: %d\n",Npart);  
  printf("Temperature: %lf\n",Temp);  
  printf("Pressure: %lf\n",0.0);  
  printf("Density Box 1: %lf\n",(double)(Npbox[0])/pow(Box[0],3));  
  printf("Box 1 Length: %lf\n",Box[0]);  
  printf("Density Box 2: %lf\n",(double)(Npbox[1])/pow(Box[1],3));  
  printf("Box 2 Length: %lf\n",Box[1]);  
  printf("Model Parameters\n");
  printf("Epsilon: %lf\n",Eps);  
  printf("Sigma: %lf\n",Sig);
  printf("Mass: %lf\n",Mass);
  
  //  ---Calculate Parameters:
  
  Beta  = 1.0/Temp;
  Eps4  = 4.0*Eps;
  Eps48 = 48.*Eps;
  Sig2  = Sig*Sig;
  
  //  ---Calculate Cut-Off Radius Potential
  
  for(BoxID = 0;BoxID< 2;BoxID++)
    {
      if(Rcc<Hbox[BoxID])
	{
	  Rc[BoxID] = Rcc;
	}
      else
	{
          printf("Rcc (%lf) is greater than Hbox[BoxID] (%lf)\n",Rcc, Hbox[BoxID]);
          exit(1);
	}
      Rc2[BoxID] = Rc[BoxID]*Rc[BoxID];
    }
  
  return;
}
 
