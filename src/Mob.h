#ifndef MOB
#define MOB

#include <iostream>
#include <stdlib.h>

#include <math.h>

#include "Node.h"


void computeOrbit(int & start, int & end, const SubPb & sub) ;

int ComputeB(int orbitFirst, int orbitLast, const SubPb & sub) ;

void doMOB(Branching & branch, SubPb & sub) ;

#endif
