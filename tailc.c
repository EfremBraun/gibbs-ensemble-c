#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "system.h"
#include "conf.h"
#include "potential.h"
#include "system_2.h"

double TailC(int BoxID)
{
  //     ---Calculates tail correction. Only call this function if TruncFlag is 0.
  //     Multiply by (double) Npbox[BoxID] to get results on a per particle basis.
  //     Multiply by (double) Npbox[BoxID] ^2 to get total tail correction for whole box.
  
  double Rc3i, Tail;
  Rc3i = pow(Sig/Rc[BoxID], 3.0);
  Tail = 2.0/3.0 * M_PI / pow(Box[BoxID], 3.0)  * Eps4 * (pow(Rc3i, 3.0) / 3.0 - Rc3i);
  return Tail;
}
