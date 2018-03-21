#ifndef MODELEUCP
#define MODELEUCP

#include <ilcplex/ilocplex.h>

#include "InstanceUCP.h"


IloModel defineModel(IloEnv env, InstanceUCP* pb, const IloBoolVarArray & x, const IloBoolVarArray & u, int uNum, int ramp) ;

void AddRSUIneq(IloModel & model, IloEnv env, InstanceUCP* pb, const IloBoolVarArray & x, const IloBoolVarArray & u, int methode) ;


IloModel defineModel_y(IloEnv env, InstanceUCP* pb, const IloBoolVarArray & x, const IloBoolVarArray & u) ;
IloModel defineModel_numberOfOnes(IloEnv env, InstanceUCP* pb, const IloBoolVarArray & x, const IloBoolVarArray & u) ;

IloModel defineModel_sum(IloEnv env, InstanceUCP* pb, const IloBoolVarArray & x, const IloBoolVarArray & u, int methode) ;
IloModel AggregatedModel(IloEnv env, InstanceUCP* pb) ;


#endif /* MODELEUCP_INCLUDED */
