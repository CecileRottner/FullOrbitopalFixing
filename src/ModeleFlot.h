#ifndef MODELEFLOT
#define MODELEFLOT

#include <ilcplex/ilocplex.h>

#include "InstanceUCP.h"


class ModeleFlot {
public:
    IloEnv env ;
    InstanceUCP* pb ;
    int T ;
    int n ;

    ModeleFlot(IloEnv enviro, InstanceUCP* pb) ;

    int Adj(int i, int t, int k) ;
    int AdjG(int g, int t, int k) ;
    double Cost(int i, int t, int k) ;
    double CostG(int g, int t, int k) ;
    int arc(int i, int t, int k) ;
    IloModel defineModelFlot(IloBoolVarArray f) ;
    IloModel AggregatedFlowModel() ;

    ~ModeleFlot() {
    }
};
//FIN CLASSE



#endif /* MODELEFLOT_INCLUDED */