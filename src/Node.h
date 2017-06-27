#ifndef BRANCH
#define BRANCH

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

using namespace std ;

class Branching {

public :
    IloNumVarArray varLeft ;
    IloNumArray bLeft ;
    IloCplex::BranchDirectionArray dirLeft ;

    IloNumVarArray varRight ;
    IloNumArray bRight ;
    IloCplex::BranchDirectionArray dirRight;

    Branching(IloEnv env) {
        varLeft = IloNumVarArray(env, 0) ;
        bLeft = IloNumArray(env, 0) ;
        dirLeft = IloCplex::BranchDirectionArray(env, 0);

        varRight = IloNumVarArray(env, 0) ;
        bRight = IloNumArray(env, 0) ;
        dirRight = IloCplex::BranchDirectionArray(env, 0);
    }

    ~Branching() {
        varLeft.end() ;
        varRight.end() ;

        bLeft.end() ;
        bRight.end() ;

        dirLeft.end() ;
        dirRight.end() ;
    }

};


class myNodeData : public IloCplex::BranchCallbackI::NodeData {
public:
    int num ;
    //ordre des t
    IloIntArray rankOf ; //rankOf[t] : rang du pas de temps t
    IloIntArray timeOf ; //timeOf[r] : pas de temps de rang r (inverse de rankOf)
    IloIntArray finOrdre ; // les t de rang 0 à finOrdre sont dans l'ordre définitif

    myNodeData(IloEnv env, int T,  int nbG, Methode methode) ; // initialisation à la racine
    myNodeData(IloEnv env, myNodeData* data, int group, int time, int varX, int feasible, int T, int nbG) ; // data est la donnée du noeud père. time est le t de la variable sur laquelle on branche

    ~myNodeData() {
        rankOf.end();
        timeOf.end() ;
        finOrdre.end() ;
    }
};
//FIN CLASSE



class SubPb {

public :

    /////données du problème, inchangées dans l'arbre
    int n ;
    int T;

    InstanceUCP* instance ;
    IloEnv env ;

    Methode met ;

    //Donnees symétries pour le MOB
    IloIntArray First ;
    IloIntArray Last ;

    //Données symétries pour le fixing
    IloIntArray FirstG ; // tableau des premiers éléments de chaque groupe
    IloIntArray LastG ;
    IloIntArray Group ;
    IloInt nbG ; // taille de FirstG

    //ordre des t
    IloIntArray rankOf ; //rankOf[t] : rang du pas de temps t
    IloIntArray timeOf ; //timeOf[r] : pas de temps de rang r (inverse de rankOf)
    IloIntArray finOrdre ; // les t de rang 0 à finOrdre sont dans l'ordre définitif

    // variables du problème
    IloBoolVarArray x ;
    IloBoolVarArray u ;

    IloNum eps ; // tolerance

    //mis à jour dans l'arbre
    int UB_LB ;
    IloArray<IloCplex::ControlCallbackI::IntegerFeasibility> feasible ;

    /////indicateurs branchement, différents à chaque noeud

    int node ;

    //UB et LB sont là pour récupérer les valeurs fournies par Cplex. Comme ça on les déclare pas à chaque callback.
    IloNumArray UB ;
    IloNumArray LB ;
    IloNumArray UB_u ;
    IloNumArray LB_u ;
    // à partr de UB et LB, méthode update calcule values
    IloIntArray values ; //valeurs de x (0, 1 ou libre)
    IloIntArray valuesU ; //valeurs de u
    IloNumArray x_frac ;// solution fractionnaire au noeud en cours. utile que pour MOB

    //Matrices Xmin et Xmax du fixing

    IloIntArray Xmin ;
    IloIntArray Xmax ;

    //à partir de varID, méthode update calcule les indicateurs suivants:
    int unit ;
    int group ; //groupe de symétrie de l'unité unit
    int time ;
    int varX ; //=1 si la variable fixée est une variable x, 0 si u
    int varU ;
    int varW ;

    int prune;

    //temps et nombre de fixing : se cumulent d'un noeud à l'autre
    double timeFix ; //temps nécessaire aux opérations de fixing / branching
    int nbFixs ; //nombres d'opérations de fixing

    SubPb(IloEnv env, InstanceUCP* inst, IloBoolVarArray xx, IloBoolVarArray uu, IloNum epsilon, Methode methode) ;
    ~SubPb() {

        First.end() ;
        Last.end()  ;

        //Données symétries pour le fixing
        FirstG.end()  ; // tableau des premiers éléments de chaque groupe
        LastG.end()  ;
        Group.end()  ;

        //ordre des t
        rankOf.end()  ; //rankOf[t] : rang du pas de temps t
        timeOf.end()  ; //timeOf[r] : pas de temps de rang r (inverse de rankOf)
        finOrdre.end()  ; // les t de rang 0 à finOrdre sont dans l'ordre définitif

        // variables du problème
        x.end()  ;
        u.end()  ;

        feasible.end()  ;

        /////indicateurs branchement, différents à chaque noeud

        //UB et LB sont là pour récupérer les valeurs fournies par Cplex. Comme ça on les déclare pas à chaque callback.
        UB.end()  ;
        LB.end()  ;
        UB_u.end()  ;
        LB_u.end()  ;
        // à partr de UB et LB, méthode update calcule values
        values.end()  ; //valeurs de x (0, 1 ou libre)
        valuesU.end()  ; //valeurs de u
        x_frac.end()  ;// solution fractionnaire au noeud en cours. utile que pour MOB

        Xmin.end() ;
        Xmax.end() ;

        env.end() ;
    }

    void update(int varID, myNodeData* data) ;
    void computeValuesU() ;
    void getVar(int varID, int & unit, int & time, int & varX) ;

    int setXmin() ;
    int setXmax() ;

    int doFixing_side(Branching & branch, int side) ;
    void doFixing(Branching & branch, int & pruneLeft, int & pruneRight, int nNodes) ;

    int newVarU(Branching & branch, int nodes) ;
    int newVarUW(Branching & branch, int nodes) ;
    int newFeasibleVar(Branching & branch, int nodes ) ;
    int newVarFromDemand(Branching & branch, int nodes) ;

};





#endif
