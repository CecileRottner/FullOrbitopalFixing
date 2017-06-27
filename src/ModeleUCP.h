#ifndef MODELEUCP
#define MODELEUCP

#include <ilcplex/ilocplex.h>

#include "InstanceUCP.h"


IloModel defineModel(IloEnv env, InstanceUCP* pb, const IloBoolVarArray & x, const IloBoolVarArray & u, int uNum) ;

IloModel defineModel_y(IloEnv env, InstanceUCP* pb, const IloBoolVarArray & x, const IloBoolVarArray & u) ;
IloModel defineModel_sum(IloEnv env, InstanceUCP* pb, const IloBoolVarArray & x, const IloBoolVarArray & u, int methode) ;

#endif /* MODELEUCP_INCLUDED */
