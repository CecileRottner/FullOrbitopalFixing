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
#include "Fixing.h"

int SubPb::setXmin() {

    for (int g = 0 ; g < nbG ; g++) {
        int last_g = LastG[g] ;
        int first_g = FirstG[g] ;

        //colonne LastG[g] de Xmin
        for (int t=0 ; t < T ;t++) {
            if (values[last_g*T+ t] < 2) {
                Xmin[last_g*T+ t] = values[last_g*T+ t] ;
            }
            else {
                Xmin[last_g*T+ t] = 0 ;
            }
        }

        for (int i = last_g-1 ; i >=first_g ; i--) {

            //colonnes i vs i+1
            int ld = -1 ; //last discriminating row

            int stop=0 ;
            int t = 0 ;

            while (!stop && (t < T)) {

                if ((values[i*T+t] == 1) && (Xmin[(i+1)*T+t] == 0)) {
                    Xmin[(i)*T+t] = values[i*T+t] ;
                    stop=1 ;
                }

                if ((values[i*T+t] == 0) && (Xmin[(i+1)*T+t] == 1)) {
                    stop = 1 ;
                    if (ld==-1) {
                        return 1 ; // 1= on élague
                    }
                    else { //on met la dernière ligne discriminante à 1
                        Xmin[i*T+ld] = 1 ;

                        // on met à 0 tout ce qui est libre en dessous de ld et au dessus de t (ce qui est en dessous de t est traité après le while)
                        for (int s=ld+1 ; s <= t ; s++) {
                            if (values[i*T+s] == 8) {
                                Xmin[i*T+s] = 0 ;
                            }
                            else {
                                Xmin[i*T+s] = values[i*T+s] ;
                            }
                        }

                    } // fin else
                }

                if (!stop) {
                    if (values[i*T+t] == 8) {
                        if (Xmin[(i+1)*T+t] == 1) {
                            Xmin[i*T+t] = 1 ;
                        }
                        else {
                            Xmin[i*T+t] = 0 ;
                            ld = t ;
                        }
                    }
                    else {
                        Xmin[i*T+t] = values[i*T+t] ;
                    }
                }
                t++ ;
            }


            if (stop) { // si on s'est arrêté, c'est que la ligne "t" (maintenant t-1) est à différence fixe
                //on met à 0 tout ce qui est en dessous}
                for (int s=t-1 ; s < T ; s++) {
                    if (values[i*T+s] == 8) {
                        Xmin[i*T+s] = 0 ;
                    }
                    else {
                        Xmin[i*T+s] = values[i*T+s];
                    }
                }
            }
            // fin i vs i+1
        }
    }
    return 0 ;
}

int SubPb::setXmax() {

    for (int g = 0 ; g < nbG ; g++) {

        int last_g = LastG[g] ;
        int first_g = FirstG[g] ;

        //colonne first_g de Xmax
        for (int t=0 ; t <T ;t++) {
            if (values[first_g*T + t] < 2) {
                Xmax[first_g*T + t] = values[first_g*T + t] ;
            }
            else {
                Xmax[first_g*T + t] = 1 ;
            }
        }

        for (int i = first_g + 1 ; i <= last_g ; i++) {

            //colonnes i-1 vs i
            int ld = -1 ; //last discriminating row

            int stop=0 ;
            int t = 0 ;

            while (!stop && (t < T)) {

                if ((Xmax[(i-1)*T+t] == 1 ) && (values[(i)*T+t] == 0)) {
                    stop=1 ;
                }

                if ((Xmax[(i-1)*T+t] == 0 ) && (values[(i)*T+t] == 1)) {
                    stop = 1 ;
                    if (ld==-1) {
                        return 1 ; // 1= on élague
                    }
                    else { //on met la dernière ligne discriminante à 0
                        Xmax[i*T+ld] = 0 ;

                        // on met à 1 tout ce qui est libre en dessous de ld et au dessus de t (ce qui est en dessous de t est traité après le while)
                        for (int s=ld+1 ; s <= t ; s++) {
                            if (values[i*T+s] == 8) {
                                Xmax[i*T+s] = 1 ;
                            }

                            else {
                                Xmax[i*T+s] = values[i*T+s] ;
                            }
                        }

                    } // fin else
                }

                if (!stop) {
                    if (values[i*T+t] == 8) {
                        if (Xmax[(i-1)*T+t] == 0) {
                            Xmax[i*T+t] = 0 ;
                        }
                        else {
                            Xmax[i*T+t] = 1 ;
                            ld = t ;
                        }
                    }

                    else { // si valeur(i,t) est fixée alors on la met dans Xmax
                        Xmax[i*T+t] = values[i*T+t] ;
                    }
                }

                t++ ;
            }

            if (stop) { // si on s'est arrêté, c'est que la ligne "t" est à différence fixe
                //on met à 1 tout ce qui est en dessous
                for (int s=t-1 ; s < T ; s++) {
                    if (values[i*T+s] == 8) {
                        Xmax[i*T+s] = 1 ;
                    }
                    else {
                        Xmax[i*T+s] = values[i*T+s];
                    }
                }

            }
        }

    }

    return 0 ;

}

int SubPb::doFixing_side(Branching & branch, int side) { // la classe SubPb contient la variable sur laquelle le branchement va être fait, branch aussi

    int rang_t = rankOf[group*T + time];
    int rang_t_1 = -1 ;

    if (time>0) {
        rang_t_1 = rankOf[group*T + time-1] ;
    }

    if (varX) { // une variable x est fixée
        values[ unit*T + rang_t ] = side ;
    }

    else if (!varX && side) { // si une variable u est fixée à 1

        values[unit*T +  rang_t ] = 1 ;
        if ((time>0) && (rang_t_1 <= finOrdre[group]) ) {
            values[unit*T + rang_t_1 ] = 0;
        }
    }

    else if  (!varX && !side)  { // un u est fixé à 0
        if ((time>0) && (rang_t_1 <= finOrdre[group]) && (values[unit*T + rang_t_1] == 0)) { //x(i,t-1) est fixé à 0
            values[unit*T + rang_t ] = 0 ;
        }

        if ((time>0) && (rang_t_1 <= finOrdre[group]) && (values[unit*T + rang_t] == 1)) { //x(i,t) est fixé à 1
            values[unit*T + rang_t_1 ] = 1 ;
        }
    }



    int elague = setXmax() ;
    if (elague) {
        return 1 ;
    }

    elague = setXmin() ;
    if (elague) {
        return 1 ;
    }

    int fix=1 ;


    //sinon, il y a une solution

    if (!elague) {
        // on détermine de quel côté (left ou right) il faut ajouter nos variables
        int left = 0 ; // left = 1 si on ajoute à gauche, 0 sinon
        if (side) {
            if (branch.bLeft[0] == 1) {
                left = 1 ;
            }
            else {
                left= 0 ;
            }
        }

        else { // side= 0
            if (branch.bLeft[0] == 0) {
                left = 1 ;
            }
            else {
                left= 0 ;
            }
        }



        for (int i=0 ; i < n ; i++) {

            int stop=0 ;
            int t=0 ;

            int g = Group[i] ;


            while (!stop && t <= finOrdre[g] ) {

                int time_t = timeOf[g*T + t] ;
                //cout << "t : " << t << endl ;
                //cout << "min-max : " << endl ;
                if (Xmax[i*T+t] == Xmin[i*T+t]) {

                    int bound = Xmax[i*T+t] ;
                    if (values[i*T+t] == 8 && feasible[i*T+time_t] != 2) {

                        /*if (feasible[i*T+time_t] == 2) {
                            cout << "i, t: " << i << " " << time_t << " is 2-feasible" << endl ;
                            cout << "LB, UB : " << LB[i*T+time_t] << ", " << UB[i*T + time_t] << endl ;
                        }*/

                        if (left) {

                            //affichage
                            if (!fix) {
                                cout << "Fixing a " <<  node*left + (node+1)*!left << endl ;
                                cout << "Variable branchée : "<< "unit, time : " << unit << ", " << time  << "= " << side << endl ;
                                cout << " x ? : " << varX << endl ;

                                cout << "Fixing of unit i, time t : " << i << ", " << time_t << endl ;
                                cout << "Group : " << FirstG[g] << ", " << LastG[g] << endl ;
                                cout << "at bound : " << bound << endl ;
                                cout << "values[(i)*T+ t] = " << values[(i)*T+t] << endl ;

                                cout << "ordre des t : " ;
                                for (int s=0 ; s <= finOrdre[g] ; s++) {
                                    cout << timeOf[g*T + s] << " " ;
                                }
                                cout << endl ;

                                for (int s= 0 ; s <= finOrdre[g] ; s++) {
                                    for (int k = FirstG[g] ; k <= LastG[g] ; k++) {
                                        cout << Xmin[k*T+s] << " " ;
                                    }
                                    cout << "          " ;
                                    for (int k = FirstG[g] ; k <= LastG[g] ; k++) {
                                        cout << values[k*T+s] << " " ;
                                    }
                                    cout << "          " ;
                                    for (int k = FirstG[g] ; k <= LastG[g] ; k++) {
                                        cout << Xmax[k*T+s] << " " ;
                                    }
                                    cout << "          " ;
                                    cout << endl ;
                                }

                                /*cout << "u : " << endl ;
                                for (int s= 0 ; s <= finOrdre[g] ; s++) {
                                    for (int k = FirstG[g] ; k <= LastG[g] ; k++) {

                                        cout << valuesU[k*T+timeOf[g*T +s]] << " " ;
                                    }
                                    cout << endl ;
                                }

                                cout << endl ;
                                cout << endl ;*/
                            }
                            fix=1 ;


                            (branch.varLeft).add( x[i*T + time_t] ) ;
                            nbFixs ++ ;

                            (branch.bLeft).add( bound ) ;
                            if (bound == 1) {
                                (branch.dirLeft).add( IloCplex::BranchUp ) ;
                            }
                            else {
                                (branch.dirLeft).add( IloCplex::BranchDown ) ;
                            }
                        }


                        else {


                            if (!fix) {
                                cout << "Fixing a " << node*left + (node+1)*!left << endl ;
                                cout << "on the right" << endl ;
                                cout << "Variable branchée : "<< "unit, time : " << unit << ", " << time  << "= " << side << endl ;
                                cout << " x ? : " << varX << endl ;

                                cout << "Fixing of unit i, time t : " << i << ", " << time_t << endl ;
                                cout << "Group : " << FirstG[g] << ", " << LastG[g] << endl ;
                                cout << "at bound : " << bound << endl ;
                                cout << "values[(i)*T+ t] = " << values[(i)*T+t] << endl ;

                                cout << "ordre des t : " ;
                                for (int s=0 ; s <= finOrdre[g] ; s++) {
                                    cout << timeOf[g*T + s] << " " ;
                                }
                                cout << endl ;

                                for (int s= 0 ; s <= finOrdre[g] ; s++) {
                                    for (int k = FirstG[g] ; k <= LastG[g] ; k++) {
                                        cout << Xmin[k*T+s] << " " ;
                                    }
                                    cout << "          " ;
                                    for (int k = FirstG[g] ; k <= LastG[g] ; k++) {
                                        cout << values[k*T+s] << " " ;
                                    }
                                    cout << "          " ;
                                    for (int k = FirstG[g] ; k <= LastG[g] ; k++) {
                                        cout << Xmax[k*T+s] << " " ;
                                    }
                                    cout << "          " ;
                                    cout << endl ;
                                }

                               /*cout << "u : " << endl ;
                                for (int s= 0 ; s <= finOrdre[g] ; s++) {
                                    for (int k = FirstG[g] ; k <= LastG[g] ; k++) {

                                        cout << valuesU[k*T+timeOf[g*T +s]] << " " ;
                                    }
                                    cout << endl ;
                                }

                                cout << endl ;
                                cout << endl ;*/
                            }
                            fix=1 ;

                            (branch.varRight).add( x[i*T + time_t] );
                            nbFixs ++ ;
                            (branch.bRight).add( bound ) ;

                            if (bound == 1) {
                                (branch.dirRight).add( IloCplex::BranchUp ) ;
                            }
                            else {
                                (branch.dirRight).add( IloCplex::BranchDown ) ;
                            }
                        }
                    }
                }

                else { // min et max sont différents
                    if (Xmin[i*T+t] > Xmax[i*T+t] + eps) {
                        cout << "contradiction min et max pour i,t = " << i << ", " << t  << endl ;
                        return 1 ;
                    }
                    stop= 1 ;
                }
                t++ ;
            } // fin while
        }
    }

    return 0 ;
}

void SubPb::doFixing(Branching & branch, int & pruneLeft, int & pruneRight, int nNodes) {

    int pruneSide0 = doFixing_side(branch, 0) ; // fixe le côté var=0, à faire en premier car pour u=1 on modifie 2 valeurs et pour u=0 on modifie au pire 1
    int pruneSide1 = doFixing_side(branch, 1) ; // fixe le côté var=1

    int left = 0 ; // left = 1 si c'est la variable de gauche qui est à 1, 0 sinon
    if (branch.bLeft[0] == 1) {
        left = 1 ;
    }
    else {
        left= 0 ;
    }


    if (left) {
        pruneLeft = pruneSide1 ;
        pruneRight = pruneSide0 ;
    }
    else {
        pruneLeft = pruneSide0 ;
        pruneRight = pruneSide1 ;
    }
    if (pruneLeft) {
        nbFixs ++ ;
    }
    if (pruneRight) {
        nbFixs ++ ;
    }
}



