#ifndef MOB
#define MOB

#include <iostream>
#include <stdlib.h>

#include <math.h>

#include "Node.h"
#include "Mob.h"



void computeOrbit(int & start, int & end, const SubPb & sub) {

    int i = sub.unit ;

    int T = sub.T ;

    int stop = sub.Last[i] ;// si i est le dernier de son groupe de symétrie, on ne fait rien à droite
    int j = i ;

    end=i ;

    while (!stop) {

        j++ ;

        int t=0 ;
        while ((t < T) & !stop) {
            if (sub.values[i*T+t] != sub.values[j*T+t]) {
                stop= 1 ;
            }
            t++ ;
        }
        if (!stop) {
            end=j;
        }
        if (sub.Last[j]) {
            stop=1 ;
        }

    }

    stop=sub.First[i] ;
    j=i ;
    start=i ;

    while (!stop) {

        j-- ;

        int t=0 ;
        while ((t < T) & !stop) {
            if (sub.values[i*T+t] != sub.values[j*T+t]) {
                stop= 1 ;
            }
            t++ ;
        }

        if (!stop) {
            start=j;
        }

        if (sub.First[j]) {
            stop=1 ;
        }
    }

}

int ComputeB(int orbitFirst, int orbitLast, const SubPb & sub) {
    double b = 0 ;
    int T = sub.T ;
    int t = sub.time ;

    for (int j = orbitFirst ; j <= orbitLast ; j++) {
        b += sub.x_frac[j*T + t] ;
    }

    // cout <<"b before ceil : " << b << endl ;
    if (b - floor(b) < sub.eps) { //si b est entier. dans ce cas ceil(b) = b+1 :/
        return fmax(1, b) ;
    }

    b= ceil(b) ;

    return fmax(1, b);
}

void doMOB(Branching & branch, SubPb & sub) {

    // il faut vérifier varX=1

    int T= sub.T ;
    int n = sub.n ;

    int time = sub.time ;
    int unit = sub.unit ;


    //cout << "at node "<< sub.node <<", x(13,16) has bounds " << sub.LB[13*T + 16]  << ", " << sub.UB[13*T + 16] << endl ;


    //Si la variable de branchement est u(i,t), alors x(i,t) ou x(i,t-1) est libre : faux

    if (sub.varX == 0) {

        cout << "branchement sur u(" << sub.unit << ", " << sub.time << ")" << endl ;
        cout << "at node : " << sub.node << endl ;

        if ((abs(sub.values[unit*T + time] - 8 ) < sub.eps) && (sub.feasible[unit*T + time] != 2)) { //si x(i,t) est libre, on branche dessus
            //on branche dessus
            cout << "x(i,t) est libre" << endl ;
            sub.varX=1 ;

        }
        else if ( (time > 0) && (abs(sub.values[unit*T + time-1] - 8 ) < sub.eps) && (sub.feasible[unit*T + time-1] != 2) ) { //sinon si x(i,t-1) est libre, on branche dessus
            //on branche sur x(i,t-1)
            time-- ;
            sub.time-- ;
            sub.varX=1 ;
            cout << "x(i,t-1) est libre : " << sub.LB[unit*T + time-1]  << ", " << sub.UB[unit*T + time-1] << endl ;

        }

        else { //sinon on branche sur la première variable libre qu'on trouve

            //cout << "first found" << endl ;
            //cout << "on cherche une autre variable" << endl ;

            // on cherche une variable libre
            int index = 0 ;
            int stop=0 ;

            while (!stop && index < n*T) {
                if (abs(sub.values[index] - 8)  < sub.eps && (sub.feasible[index] != 2)) {
                    stop=1 ; // on a trouvé une variable libre
                }
                else {
                    index++ ;
                }
            }

            if (index < n*T) {
                time = index % T ;
                unit = (index - time)/ T ;
                sub.time=time ;
                sub.unit=unit ;
                sub.varX = 1 ;

                cout << "On choisit plutot la variable x: time, unit: " << time << ", " << unit << endl ;
            }
            else {
                cout << "pas de variable x libre" << endl ;
            }

        }
    }


    int orbiteFirst ;
    int orbiteLast;

    computeOrbit(orbiteFirst, orbiteLast, sub) ;

    int b = ComputeB(orbiteFirst, orbiteLast, sub) ;
    int sizeOrb = orbiteLast-orbiteFirst+1 ;


    /*cout << "time : " << time << endl ;
         cout << "unit : " << unit << endl ;
         cout << "x(i,t) : " << sub.values[unit*T + time] << endl ;

         cout << "lb, ub : " << sub.LB[unit*T + time] <<", " <<sub.UB[unit*T + time] << endl ;

         cout << "orbite : " << orbiteFirst << ", " << orbiteLast << endl ;
         cout << "b : " << b << endl ; */

    for (int j = 0 ; j < b ; j++ ) {

        //cout << j+orbiteFirst << " à 1 ; values : " << sub.values[(j+orbiteFirst)*T + time] << endl ;

        //variables de l'orbit
        if (sub.feasible[(j+orbiteFirst)*T+time] !=2 ) {

            (branch.varLeft).add( sub.x[(j+orbiteFirst)*T+time] ) ;
            sub.nbFixs ++ ;

            //bounds
            (branch.bLeft).add( 1 ) ;

            //directions
            (branch.dirLeft).add( IloCplex::BranchUp) ;
        }
    }


    for (int j = 0 ; j < sizeOrb-b+1 ; j++ ) {

        //cout << j+orbiteFirst+b-1 << " à 0 ;  values : " << sub.values[(j+orbiteFirst+b-1)*T + time] << endl ;


        if (sub.feasible[(j+orbiteFirst+b-1)*T+time] !=2 ) {

            //variables de l'orbit
            (branch.varRight).add( sub.x[(j+orbiteFirst+b-1)*T+time] );
            sub.nbFixs ++ ;

            //bounds
            (branch.bRight).add(0) ;

            //directions
            (branch.dirRight).add(IloCplex::BranchDown) ;
        }

    }
    sub.nbFixs -= 2 ; //on retire les variables qui auraient été fixées de toute façon


}

#endif
