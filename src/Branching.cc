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
#include "ModeleUCP.h"

#include "Node.h"


int SubPb::newVarFromFractionalGroup(Branching & branch, int nodes) {

    int new_time=-1 ;
    int new_unit=-1 ;

    int distance = 1 ;

    for (int g=0 ; g < nbG ; g++) {

        int first = FirstG[g];

        for (int t=0 ; t< T ; t++) {
            int val = 0 ;
            for (int i = first ; i <= LastG[g] ; i++) {
                val += x_frac[i*T+t] ;
            }

            int IntPart = floor(val) ;
            int FracPart = val - IntPart ;

            if (fabs(FracPart - 0.5) < distance) {

                int possibleUnit= first + IntPart ;

                if (feasible[possibleUnit*T+t] !=2) {
                    new_unit=possibleUnit ;
                    new_time=t ;

                    distance =fabs(FracPart - 0.5) ;
                }

                else {
                    for (int i = first ; i <= LastG[g] ; i++) {

                        if (feasible[i*T+t] != 2) {
                            new_unit=i ;
                            new_time=t ;
                            distance = fabs(FracPart - 0.5) ;
                        }
                    }
                }
            }

        }
    }

    if (new_unit<0) {
        return 0 ;
    }


    unit = new_unit ;
    time = new_time ;
    (branch.varLeft).add( x[unit*T+time] ) ;
    //bound
    (branch.bLeft).add( 1 ) ;
    //direction
    (branch.dirLeft).add( IloCplex::BranchUp) ;

    //coté droit
    //variable
    (branch.varRight).add( x[unit*T+time] );
    //bound
    (branch.bRight).add(0) ;
    //direction
    (branch.dirRight).add(IloCplex::BranchDown) ;
    return 1 ;

}


int SubPb::newVarFromSymmetryGroup(Branching & branch, int nodes) {

    int new_time=-1 ;
    int new_unit=-1 ;

    int g = Group[unit] ;
    int first = FirstG[g];

    int distance = 1 ;
    for (int t=0 ; t< T ; t++) {
        int val = 0 ;
        for (int i = first ; i <= LastG[g] ; i++) {
            val += x_frac[i*T+t] ;
        }

        int IntPart = floor(val) ;
        int FracPart = val - IntPart ;

        if (fabs(FracPart - 0.5) < distance) {

            int possibleUnit= first + IntPart ;

            if (feasible[possibleUnit*T+t] !=2) {
                new_unit=possibleUnit ;
                new_time=t ;

                distance =fabs(FracPart - 0.5) ;
            }

            else {
                for (int i = first ; i <= LastG[g] ; i++) {

                    if (feasible[i*T+t] !=2) {
                        new_unit=i ;
                        new_time=t ;
                        distance =fabs(FracPart - 0.5) ;
                    }
                }
            }
        }

    }

    if (new_unit<0) {
        return 0 ;
    }


    unit = new_unit ;
    time = new_time ;
    (branch.varLeft).add( x[unit*T+time] ) ;
    //bound
    (branch.bLeft).add( 1 ) ;
    //direction
    (branch.dirLeft).add( IloCplex::BranchUp) ;

    //coté droit
    //variable
    (branch.varRight).add( x[unit*T+time] );
    //bound
    (branch.bRight).add(0) ;
    //direction
    (branch.dirRight).add(IloCplex::BranchDown) ;
    return 1 ;

}


int SubPb::newFeasibleVar(Branching & branch, int nodes) {

    cout << "unit : " << unit << endl ;
    cout << "time-1 : " << time-1 << endl ;

    if (feasible[unit*T+time] == 2) { // la variable branchée est un u, le x associé est ImpliedFeasible

        if ((time>0) && feasible[unit*T+time-1] != 2) {
            int g = Group[unit];
            //if ( lb < ub -eps ) {
            time=time-1 ;

            cout << "new var" << endl ;
            (branch.varLeft).add( u[unit*T+time] ) ;
            //bound
            (branch.bLeft).add( 1 ) ;
            //direction
            (branch.dirLeft).add( IloCplex::BranchUp) ;

            //coté droit
            //variable
            (branch.varRight).add( u[unit*T+time] );
            //bound
            (branch.bRight).add(0) ;
            //direction
            (branch.dirRight).add(IloCplex::BranchDown) ;

            return 1;
            // }

            cout << "ici" << endl ;


        }
    }
    return 0 ;
}

int SubPb::newVarFromDemand(Branching & branch, int nodes) {

    if (nodes > 10) {
        return 0;
    }


    int new_time=-1 ;
    int new_unit=-1 ;
    int demand=0 ;
    for (int t=0 ; t< T ; t++) {
        int i = 0;
        if (instance->getD(t) > demand) {
            while (i < n && feasible[i*T+t] != 1) {
                i++ ;
            }
            if (i < n) {
                new_unit=i ;
                new_time=t ;
            }
        }
    }

    if (unit<0) {
        return 0 ;
    }


    unit = new_unit ;
    time = new_time ;
    (branch.varLeft).add( x[unit*T+time] ) ;
    //bound
    (branch.bLeft).add( 1 ) ;
    //direction
    (branch.dirLeft).add( IloCplex::BranchUp) ;

    //coté droit
    //variable
    (branch.varRight).add( x[unit*T+time] );
    //bound
    (branch.bRight).add(0) ;
    //direction
    (branch.dirRight).add(IloCplex::BranchDown) ;
    return 1 ;

}


int SubPb::newVarU(Branching & branch, int nodes) {

    if (varX == 1 && nodes < 100) {

        int l = instance->getl(unit);
        int L = instance->getL(unit) ;
        int g = Group[unit];

        for (int k = fmax(0,time - l) ; k <= time - 1 ; k++) {
            if (values[unit*T + rankOf[g*T+k]] == 1) {
                return 0 ;
            }
        }

        for (int k = time + 1 ; k < fmin(T,time + L) ; k++) {
            if (values[unit*T + rankOf[g*T +k]] == 0) {
                return 0 ;
            }
        }
        if (valuesU[unit*T+time] == 0) {
            return 0 ;
        }
        //si on est là, c'est que u(unit, time) peut être fixé à 1

        varX = 0 ;
        varU = 1 ;

        //branch est mis à jour

        ///Manière qui marche pas
        /*branch.varLeft = IloNumVarArray(env, 0) ;
        (branch.varLeft).add( u[unit*T+time] );

        branch.varRight = IloNumVarArray(env, 0) ;
        (branch.varRight).add( u[unit*T+time] ) ;*/

        ///Manière qui marche pas
        /*  branch.varLeft[0] = u[unit*T+time] ;
        branch.varRight[0] = u[unit*T+time] ;*/

        ///Manière qui marche
        (branch.varLeft).add( u[unit*T+time] ) ;
        //bound
        (branch.bLeft).add( 1 ) ;
        //direction
        (branch.dirLeft).add( IloCplex::BranchUp) ;

        //coté droit
        //variable
        (branch.varRight).add( u[unit*T+time] );
        //bound
        (branch.bRight).add(0) ;
        //direction
        (branch.dirRight).add(IloCplex::BranchDown) ;
        return 1 ;
    }

    // sinon, pas possible
    return 0 ;

}


int SubPb::newVarUW(Branching & branch, int nodes) {

    if (varX == 1 && nodes < 100) {

        int l = instance->getl(unit);
        int L = instance->getL(unit) ;

        int g = Group[unit];

        if (instance->isIncreasing[time]) {

            for (int k = fmax(0,time - l) ; k <= time - 1 ; k++) {
                if (values[unit*T + rankOf[g*T+k]] == 1) {
                    return 0 ;
                }
            }

            for (int k = time + 1 ; k < fmin(T,time + L) ; k++) {
                if (values[unit*T + rankOf[g*T +k]] == 0) {
                    return 0 ;
                }
            }

            if (valuesU[unit*T+time] == 0) {
                return 0 ;
            }

            //si on est là, c'est que u(unit, time) peut être fixé à 1

            varX=0 ;
            varU=1 ;
            varW=0 ;

            //branch est mis à jour
            (branch.varLeft).add( u[unit*T+time] ) ;
            //bound
            (branch.bLeft).add( 1 ) ;
            //direction
            (branch.dirLeft).add( IloCplex::BranchUp) ;

            //coté droit
            //variable
            (branch.varRight).add( u[unit*T+time] );
            //bound
            (branch.bRight).add(0) ;
            //direction
            (branch.dirRight).add(IloCplex::BranchDown) ;
        }

        else {

        }
    }

    // sinon, pas possible
    return 0 ;

}
