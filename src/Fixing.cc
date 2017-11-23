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


// on a vérifié au préalable qu'il y avait des groupes au temps t
int SubPb::setXmin(const SubMatrices & SubM, int tsub)  {

    for (int g = 0 ; g < nbG ; g++) {
        if (SubM.SetForG(g,tsub)) {
            int last_g = LastG[g] ;
            int first_g = FirstG[g] ;

            //redéfinition de last_g

            int i=last_g ;
            while ( !SubM.inGroup(i,tsub) && i > first_g ) {
                i-- ;
            }
            last_g=i ;

            //colonne last_g de Xmin
            for (int r=0 ; r <= finOrdre[g] ;r++) {
                if (timeOf[g*T+r] >= tsub) {
                    if (values[last_g*T+ r] < 2) { // variable non libre
                        Xmin[last_g*T+ r] = values[last_g*T+ r] ;
                    }
                    else {
                        Xmin[last_g*T+ r] = 0 ;
                    }
                }
            }

            int previous_unit = last_g ;
            for (int i = last_g-1 ; i >=first_g ; i--) {

                if (SubM.inGroup(i, tsub)) { // si i est dans le groupe considéré

                    //colonnes i vs i+1
                    int ld = -1 ; //last discriminating row

                    int stop=0 ;
                    int r = 0 ;

                    while (!stop && (r <= finOrdre[g])) {

                        if (timeOf[g*T+r] >= tsub) {
                            if ((values[i*T+r] == 1) && (Xmin[previous_unit*T+r] == 0)) {
                                Xmin[(i)*T+r] = values[i*T+r] ;
                                stop=1 ;
                            }

                            if ((values[i*T+r] == 0) && (Xmin[previous_unit*T+r] == 1)) {
                                stop = 1 ;
                                if (ld==-1) {
                                    return 1 ; // 1= on élague
                                }
                                else { //on met la dernière ligne discriminante à 1
                                    Xmin[i*T+ld] = 1 ;

                                    // on met à 0 tout ce qui est libre en dessous de ld et au dessus de t (ce qui est en dessous de t est traité après le while)
                                    for (int s=ld+1 ; s <= r ; s++) {
                                        if (timeOf[g*T+s] >= tsub) {
                                            if (values[i*T+s] == 8) {
                                                Xmin[i*T+s] = 0 ;
                                            }
                                            else {
                                                Xmin[i*T+s] = values[i*T+s] ;
                                            }
                                        }
                                    }

                                } // fin else
                            }

                            if (!stop) {
                                if (values[i*T+r] == 8) {
                                    if (Xmin[previous_unit*T+r] == 1) {
                                        Xmin[i*T+r] = 1 ;
                                    }
                                    else {
                                        Xmin[i*T+r] = 0 ;
                                        ld = r ;
                                    }
                                }
                                else {
                                    Xmin[i*T+r] = values[i*T+r] ;
                                }
                            }

                        } //fin if timeOf[g*T+r] >= tsub

                        r++ ;

                    } // fin while


                    if (stop) { // si on s'est arrêté, c'est que la ligne "t" (maintenant t-1) est à différence fixe
                        //on met à 0 tout ce qui est en dessous}
                        for (int s=r-1 ; s <= finOrdre[g] ; s++) {
                            if (timeOf[g*T+s] >= tsub) {
                                if (values[i*T+s] == 8) {
                                    Xmin[i*T+s] = 0 ;
                                }
                                else {
                                    Xmin[i*T+s] = values[i*T+s];
                                }
                            }
                        }
                    }// fin i vs i+1

                    previous_unit=i ;

                }//fin if i in group
            }
        }
    }
    return 0 ;
}

int SubPb::setXmax(const SubMatrices & SubM, int tsub)  {

    for (int g = 0 ; g < nbG ; g++) {
        if (SubM.SetForG(g,tsub)) {
            int last_g = LastG[g] ;
            int first_g = FirstG[g] ;

            //redéfinition de first_g
            int i=first_g ;
            while ( !SubM.inGroup(i,tsub) && i < last_g) {
                i++ ;
            }
            first_g=i ;


            //colonne first_g de Xmax
            for (int r=0 ; r <= finOrdre[g] ; r++) {
                if (timeOf[g*T+r] >= tsub) {
                    if (values[first_g*T + r] < 2) {
                        Xmax[first_g*T + r] = values[first_g*T + r] ;
                    }
                    else {
                        Xmax[first_g*T + r] = 1 ;
                    }
                }
            }

            int previous_unit = first_g ;
            for (int i = first_g + 1 ; i <= last_g ; i++) {

                if (SubM.inGroup(i, tsub)) {

                    //colonnes i-1 vs i
                    int ld = -1 ; //last discriminating row

                    int stop=0 ;
                    int r = 0 ;

                    while (!stop && (r <= finOrdre[g])) {

                        if (timeOf[g*T+r] >= tsub) {

                            if ((Xmax[previous_unit*T+r] == 1 ) && (values[(i)*T+r] == 0)) {
                                stop=1 ;
                            }

                            if ((Xmax[previous_unit*T+r] == 0 ) && (values[(i)*T+r] == 1)) {
                                stop = 1 ;
                                if (ld==-1) {
                                    return 1 ; // 1= on élague
                                }
                                else { //on met la dernière ligne discriminante à 0
                                    Xmax[i*T+ld] = 0 ;

                                    // on met à 1 tout ce qui est libre en dessous de ld et au dessus de t (ce qui est en dessous de t est traité après le while)
                                    for (int s=ld+1 ; s <= r ; s++) {
                                        if (timeOf[g*T+s] >= tsub) {
                                            if (values[i*T+s] == 8) {
                                                Xmax[i*T+s] = 1 ;
                                            }

                                            else {
                                                Xmax[i*T+s] = values[i*T+s] ;
                                            }
                                        }
                                    }

                                } // fin else
                            }

                            if (!stop) {
                                if (values[i*T+r] == 8) {
                                    if (Xmax[previous_unit*T+r] == 0) {
                                        Xmax[i*T+r] = 0 ;
                                    }
                                    else {
                                        Xmax[i*T+r] = 1 ;
                                        ld = r ;
                                    }
                                }

                                else { // si valeur(i,t) est fixée alors on la met dans Xmax
                                    Xmax[i*T+r] = values[i*T+r] ;
                                }
                            }
                        }

                        r++ ;
                    } // fin while

                    if (stop) { // si on s'est arrêté, c'est que la ligne "r" est à différence fixe
                        //on met à 1 tout ce qui est en dessous
                        for (int s=r-1 ; s <= finOrdre[g] ; s++) {
                            if (timeOf[g*T+s] >= tsub) {
                                if (values[i*T+s] == 8) {
                                    Xmax[i*T+s] = 1 ;
                                }
                                else {
                                    Xmax[i*T+s] = values[i*T+s];
                                }
                            }
                        }
                    }

                    previous_unit=i ;
                } //fin "if subGroup(i,tsub)"
            }

        }
    }

    return 0 ;

}

void SubPb::updateValues(int side) { // on considère la matrice de la solution partielle sur la branche où la variable branchée est égale à "side". values est mis à jour en conséquence.
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

}

int SubPb::ComputeSubMatrixFixing(Branching & branch, int left, SubMatrices SubM, int tsub) { // la classe SubPb contient la variable sur laquelle le branchement va être fait, branch aussi

    int elague = setXmin(Full,0) ;
    if (elague) {
        return 1 ;
    }

    elague = setXmax(Full,0) ;
    if (elague) {
        return 1 ;
    }

    int fix=1 ;


    //sinon, il y a une solution

    if (!elague) {


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
                                cout << "Variable branchée : "<< "unit, time : " << unit << ", " << time << endl ;
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
                                cout << "Variable branchée : "<< "unit, time : " << unit << ", " << time  <<  endl ;
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

void SubPb::doFixing(Branching & branch, int & pruneLeft, int & pruneRight) {

    int left; // est-ce que la branche considérée (côté 0 ou 1) est la branche de gauche

    ////////////////////////////// Côté branche à 0 //////////////////////////////
    if (branch.bLeft[0] == 0) {
        left = 1 ;
    }
    else {
        left= 0 ;
    }

    updateValues(0) ;
    int pruneSide0 = ComputeSubMatrixFixing(branch, left, Full, 0) ; // fixe le côté var=0, à faire en premier car pour u=1 on modifie 2 valeurs et pour u=0 on modifie au pire 1


    ////////////////////////////// Côté branche à 1 //////////////////////////////
    left=!left ;
    resetValues() ;
    updateValues(1) ;
    int pruneSide1 = ComputeSubMatrixFixing(branch, left, Full, 0) ; // fixe le côté var=1


    ///////////////////////////////// Elagage ? /////////////////////////////////

    // maintenant left = 1 si c'est la variable de gauche qui est à 1, 0 sinon

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



