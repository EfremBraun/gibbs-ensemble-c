#include "system_2.h"
extern double Box[2],Hbox[2],Temp,Beta,BetaMelt,DvMod;
extern int TruncFlag, ModGibbsFlag;
extern double AvgNpbox[2], AvgBox[2], AvgDens[2];
extern int AvgCount[2];

//extern double /Sys1/,Box[2],Hbox[2],Temp,Beta;
//C
//C     Box      : Simulation Box Length
//C     Hbox     : 0.5 * Box
//C     Temp     : Temperature
//C     Beta     : 1/Temp
