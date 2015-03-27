#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "conf.h"
#include "system.h"


void WritePdb(FILE *FilePtrBox0, FILE *FilePtrBox1)
{
  int I;
  int Countatom;
  static int Countmodel=0;

  Countmodel++;

  // Header
  fprintf(FilePtrBox0, "MODEL %4d\n", Countmodel);
  fprintf(FilePtrBox1, "MODEL %4d\n", Countmodel);
  fprintf(FilePtrBox0, "CRYST1%9.3lf%9.3lf%9.3lf%7.2lf%7.2lf%7.2lf\n", Box[0], Box[0], Box[0], 90.00, 90.00, 90.00);
  fprintf(FilePtrBox1, "CRYST1%9.3lf%9.3lf%9.3lf%7.2lf%7.2lf%7.2lf\n", Box[1], Box[1], Box[1], 90.00, 90.00, 90.00);

  // Atom
  Countatom=0;
  for(I=0;I<Npart;I++)
  {
    Countatom++;
    if (Id[I] == 0)
    {
      fprintf(FilePtrBox0,"ATOM%7d  H              %8.3lf%8.3lf%8.3lf                          H\n",
        Countatom, X[I], Y[I], Z[I]);
    }
    if (Id[I] == 1)
    {
      fprintf(FilePtrBox1,"ATOM%7d  H              %8.3lf%8.3lf%8.3lf                          H\n",
        Countatom, X[I], Y[I], Z[I]);
    }
  }

  // End
  fprintf(FilePtrBox0,"%s\n","ENDMDL");
  fprintf(FilePtrBox1,"%s\n","ENDMDL");
}
