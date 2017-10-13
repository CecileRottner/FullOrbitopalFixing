#ifndef INSTANCEUCP
#define INSTANCEUCP

#include <ilcplex/ilocplex.h>
#include <fstream>

class InstanceUCP {

public:
    IloEnv env ;
    int id ;

    IloInt n, T, K, S ;
    IloBoolArray Init ;
    IloIntArray L, l, nk, nS ;
    IloNumArray D, P, Pmax, cf, c0, cp ;

    IloIntArray First ; // i est le premier element de son groupe, First[i]=1
    IloIntArray Last ;

    IloIntArray First_ ; // i est le premier element de son groupe, First[i]=1
    IloIntArray Last_ ;

    IloIntArray FirstG_ ; // tableau des premiers éléments de chaque groupe
    IloIntArray LastG_ ;
    IloBoolArray Init_ ;
    IloIntArray L_, l_, nk_, nS_ ;
    IloNumArray P_, Pmax_, cf_, c0_ , cp_;

    IloIntArray C_ ;

    IloInt SommePmax ;

    //indicateurs symétrie
    int nbGroupes ; //sans compter les unités seules
    int MaxSize ; //taille du plus gros groupe de symétrie
    double MeanSize ; //taille moyenne d'un groupe

    IloIntArray FirstG ; // tableau des premiers éléments de chaque groupe
    IloIntArray LastG ;
    IloIntArray SizeG ;
    IloIntArray Group ; //group[i] = numéro du groupe de symétrie de l'unité i
    int nbG ; //nb de groupes en comptant les unités seules
    int nbG2 ; //nb de groupes de taille >= 2

    IloIntArray repartition_tailles ;

    //t triés dans l'ordre de la demande décroissante
    IloIntArray ordreT ;


    //pour le branchement maison
    int tmax1 ;
    int tmax2 ;
    int tmin ;
    int Dmax ;
    int Dmin ;

    IloIntArray ordre_ratio ;
    int nb ; // unités de base dans l'ordre ratio
    int np ; // unités de pic dans l'ordre ratio

    IloIntArray isIncreasing ;


public:

    InstanceUCP(IloEnv envir, const char* file) ;
    ~InstanceUCP() {

        Init.end() ;

        L.end() ;
        l.end() ;
        nk.end() ;
        nS.end() ;

        D.end() ;
        P.end() ;
        Pmax.end() ;
        cf.end() ;
        c0.end() ;
        cp.end() ;

        First.end() ; // i est le premier element de son groupe, First[i]=1
        Last.end() ;

        Init_.end() ;
        L_.end() ;
        l_.end() ;
        nk_.end() ;
        nS_.end() ;

        P_.end() ;
        Pmax_.end() ;
        cf_.end() ;
        c0_.end() ;
        cp_.end() ;

        C_.end() ;


        FirstG.end() ; // tableau des premiers éléments de chaque groupe
        LastG.end() ;
        Group.end() ; //group[i] = numéro du groupe de symétrie de l'unité i

        //t triés dans l'ordre de la demande décroissante
        ordreT.end() ;

        ordre_ratio.end() ;
        isIncreasing.end() ;
    }


    void Lecture(const char* file);
    void Initialise() ;

    IloInt getn() ;
    IloInt getT() ;
    IloBool getInit(IloInt i) ;
    IloInt getL(IloInt i) ;
    IloInt getl(IloInt i) ;
    IloNum getD(IloInt i) ;
    IloNum getP(IloInt i) ;
    IloNum getPmax(IloInt i) ;
    IloNum getcf(IloInt i) ;
    IloNum getc0(IloInt i) ;
    IloNum getcp(IloInt i) ;

    IloNum getnbG();
    IloNum getFirst(IloInt i);
    IloNum getLast(IloInt i);
    IloNum getFirstG(IloInt i);
    IloNum getLastG(IloInt i);
    IloNum getSizeG(IloInt g);
    IloNum getGroup(IloInt i);
    IloNum getordreT(IloInt t) ;

    IloInt getSommePmax() ;

    IloInt partition(IloNumArray const & ordre, IloIntArray & indices, IloInt p, IloInt q) ;

    void quickSort(IloNumArray const & ordre, IloIntArray & indices, IloInt p, IloInt q) ;

};

#endif /* INSTANCEUCP_INCLUDED */
