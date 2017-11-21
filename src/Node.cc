#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <math.h>
#include <ctime>
#include <list>

#include "InstanceUCP.h"
#include "Process.h"
#include "Node.h"


using namespace std ;

SousMatrices::SousMatrices(IloEnv env, int n, int T, int nbG) {
    M = IloIntArray(env, n*T) ;
    ThereIsASet = IloIntArray(env, T);
    ThereIsASetForG = IloIntArray(env, nbG*T) ;
}

myNodeData::myNodeData(IloEnv env, int T,  int nbG, Methode methode) { // initialisation à la racine
    num=0 ;
    rankOf = IloIntArray(env, T*nbG) ;
    timeOf = IloIntArray(env, T*nbG) ;
    finOrdre= IloIntArray(env, nbG);
    for (int g=0 ; g < nbG ; g++) {
        for (int t=0 ; t < T; t++) {
            rankOf[g*T + t] = t ; //pour l'instant on les laisse dans l'ordre, et on permute lorsque le branchement a lieu
            timeOf[g*T + t] = t ;
        }
        if (methode.DynamicFixing()) {
            finOrdre[g] = -1 ;
        }
        else {
            finOrdre[g] = T-1 ;
        }
    }

}

myNodeData::myNodeData(IloEnv env, myNodeData* data, int group, int time, int varX, int feasibility, int T, int nbG) { // data est la donnée du noeud père. time est le t de la variable sur laquelle on branche


    finOrdre = IloIntArray(env, nbG) ;
    rankOf = IloIntArray(env, nbG*T) ;
    timeOf = IloIntArray(env, nbG*T) ;

    for (int g = 0 ; g < nbG ; g++) {
        for (int t = 0 ; t < T ; t++) {
            rankOf[g*T+t] = data->rankOf[g*T+t] ;
            timeOf[g*T+t] = data->timeOf[g*T+t] ;
        }
        finOrdre[g] = data->finOrdre[g] ;
    }

    if (varX==0 && time>0) { // on ordonne aussi time-1

        if ( rankOf[group*T + time-1] > finOrdre[group] /*&& feasibility != 2*/) { // time n'a pas encore été ordonné définitivement

            finOrdre[group] ++ ;
            int tmp_rank = rankOf[group*T + time-1] ; // le rang de time avant le réordonnancement
            int tmp_time = timeOf[group*T + finOrdre[group]] ; // le pas de temps qui se trouve à la place finOrdre

            //time à sa place
            rankOf[group*T + time-1]=finOrdre[group] ;
            timeOf[group*T + finOrdre[group]] = time-1 ;

            // tmp à sa place
            rankOf[group*T + tmp_time] = tmp_rank ;
            timeOf[group*T + tmp_rank] = tmp_time ;
        }

    }

    if ( rankOf[group*T + time] > finOrdre[group] /*&& feasibility != 2*/) { // time n'a pas encore été ordonné définitivement

        finOrdre[group] ++ ;
        int tmp_rank = rankOf[group*T + time] ; // le rang de time avant le réordonnancement
        int tmp_time = timeOf[group*T + finOrdre[group]] ; // le pas de temps qui se trouve à la place finOrdre

        //time à sa place
        rankOf[group*T + time]=finOrdre[group] ;
        timeOf[group*T + finOrdre[group]] = time ;

        // tmp à sa place
        rankOf[group*T + tmp_time] = tmp_rank ;
        timeOf[group*T + tmp_rank] = tmp_time ;
    }
}


SubPb::SubPb(IloEnv env, InstanceUCP* inst, IloBoolVarArray xx, IloBoolVarArray uu, IloNum epsilon, Methode methode) :
    n(inst->getn()), T(inst->getT()), nbG(inst->getnbG()),
    RSU(SousMatrices(env, n, T, nbG)),
    RSD(SousMatrices(env, n, T, nbG)),
    Full(SousMatrices(env, n, T, nbG))

{

   // n= inst->getn() ;

    instance=inst ;

    nbFixs = 0 ;
    timeFix = 0 ;
    unit=-1 ;
    time=-1 ;
    varX=-1 ;
    prune=0 ;

    eps=epsilon ;
    met=methode ;

    x=xx ;
    u=uu;

    UB_LB = 0 ;

    feasible = IloArray<IloCplex::ControlCallbackI::IntegerFeasibility>(env, n*T);
    //// First, FirstG, Last, LastG

    First = IloIntArray(env,n) ;
    Last = IloIntArray(env,n) ;
    Group = IloIntArray(env,n) ;
    for (int i = 0 ; i <n ; i++) { //à initialiser avec les valeurs de First et Last dans inst lorsqu'elles existeront
        First[i] = inst->getFirst(i) ; //i est le premier du groupe de symétrie
        Last[i] =  inst->getLast(i) ; //i est le dernier de son groupe de symétrie
        Group[i] = inst->getGroup(i);
    }

    cout << "Last: " << Last << endl ;
    cout << "Group: " << Group << endl;

    FirstG = IloIntArray(env, nbG) ;
    LastG = IloIntArray(env, nbG) ;
    for (int i=0 ; i < nbG ; i++) {
        FirstG[i] = inst->getFirstG(i) ;
        LastG[i] = inst->getLastG(i) ;
    }

    for (int t = 0 ; t < T ; t++) {
        Full.ThereIsASet[t]=1 ;
        for (int i=0 ; i < n ; i++) {
            Full.M[i*T+t]=1 ;
        }
        for (int g= 0 ; g < nbG ; g++) {
            Full.ThereIsASetForG[g*T+t] = 1 ;
        }

    }

    ///// rankOf, timeOf, finOrdre

    rankOf =  IloIntArray(env, nbG*T) ;
    timeOf =  IloIntArray(env, nbG*T) ;
    finOrdre =  IloIntArray(env, nbG) ;

    if (met.DynamicFixing()==0) {
        for (int g=0 ; g < nbG ; g++) {
            for (int t=0 ; t < T; t++) {
                rankOf[g*T + t] = t ; //pour l'instant on les laisse dans l'ordre, et on permute lorsque le branchement a lieu
                timeOf[g*T + t] = t ;
            }
            finOrdre[g] = T-1 ;
        }
    }

    /*if (methode!=3) {
        for (int i=0 ; i < T; i++) {
            rankOf[i] = i ; //pour l'instant on les laisse dans l'ordre, et on permute lorsque le branchement a lieu
            timeOf[i] = i ;
        }
    }
    else { //on ordonne les t par rapport à la demande
        for (int i=0 ; i < T; i++) {
            rankOf[i] = inst->getordreT(i) ;
            timeOf[inst->getordreT(i)] = i ;
        }
    }

    if (methode==4) {
        for (int g=0 ; g < nbG ; g++) {
            finOrdre[g]=-1 ;
        }
    }
    else {
        finOrdre=T-1 ;
    }*/


    ///// values, LB, UB, Xmin, Xmax

    values = IloIntArray(env, n*T) ;
    valuesU = IloIntArray(env, n*T) ;

    for (int i=0 ; i < n*T ; i++) {
        values[i] = 8 ;
    }

    LB = IloNumArray(env, n*T) ;
    UB = IloNumArray(env, n*T) ;
    LB_u = IloNumArray(env, n*T) ;
    UB_u = IloNumArray(env, n*T) ;
    x_frac = IloNumArray(env, 0);

    Xmin = IloIntArray(env, n*T) ;
    Xmax = IloIntArray(env, n*T) ;
}


void SubPb::getVar(int varID, int & unit, int & time, int & varX) {

    int nX = n*T ;
    int var=varID ;
    varX = 1 ;

    if (varID >= nX ) { // une variable u(i,t) est fixée
        var -= nX ;
        varX=0 ;
    }

    time = var % T ;
    unit = (var - time)/ T ;
}

void SubPb::update(int varID, myNodeData* data) {
    prune=0 ;

    //mise à jour unit et time
    int nX = n*T ;
    int var= varID ;


    if (varID < nX ) { // une variable x(i,t) est fixée
        varX=1 ;
    }
    else { // une variable u(i,t) est fixée
        varX=0 ;
        varU=1 ;
        var -= nX ;

        if (var >= nX) {
            var -= nX ;
        }

    }
    time = var % T ;
    unit = (var - time)/ T ;
    group = Group[unit] ;


    //il faut juste récupérer les valeurs données par myNodeData
    //dans le cas des méthodes à ordre fixe, c'est pas la peine

    node = data->num ;

    if (met.DynamicFixing()) {
        for (int g=0 ; g < nbG ; g++) {
            finOrdre[g] = data->finOrdre[g] ;
            for (int t = 0 ; t < T ; t++) {
                rankOf[g*T + t] = data->rankOf[g*T + t] ;
                timeOf[g*T + t] = data->timeOf[g*T + t] ;
            }
        }
    }

    //Mise à jour de values en fonction de l'ordre. tout est mis dans value, même les pas de temps non ordonnés. ils ne seront pas traités malgré tout par le fixing.
    for (int i=0 ; i < n ; i++) {
        int g = Group[i] ;
        for (int r=0 ; r < T ; r++) {
            // for (int t=0 ; t < T ; t++) {
            int time = timeOf[g*T + r] ;
            if (LB[i*T + time] >= 1 -eps) {
                values[i*T + r] = 1 ;
            }
            else if (UB[i*T + time] <= eps) {
                values[i*T + r] = 0 ;
            }
            else {
                values[i*T + r] = 8 ;
            }
            if (UB[i*T + time] - LB[i*T + time] < -eps ) {

                values[i*T + r] = 8 ;

                UB_LB=1 ;
                cout << "Cas où UB < LB" << endl ;
                cout << "i, t : " << i << ", " <<  time << endl ;
                cout << "at node : " << node << ", bounds : " << LB[i*T + time] << "; " << UB[i*T + time] << endl ;
                cout << "Feasibility : " << feasible[i*T + time] << endl ;
                cout << endl ;

                prune=1 ;
            }
        }
        /*for (int t= finOrdre[g] +1 ; t < T ; t++) {
            values[i*T + t] = 8 ;
        }*/
    }

}

void SubPb::computeValuesU() {

    for (int i=0 ; i < n*T; i++) {
        if (LB_u[i] >= 1 -eps) {
            valuesU[i] = 1 ;
        }
        else if (UB_u[i] <= eps) {
            valuesU[i] = 0 ;
        }
        else {
            valuesU[i] = 8 ;
        }

        if (LB_u[i] > UB_u[i]+eps) {
            cout << "Cas ou Ubu < LBu" << endl ;

        }
    }
}
