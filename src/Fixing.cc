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


void SubPb::affichage(int left, int t, int i, int bound, const SubMatrices & SubM, int tsub) {

    int g = Group[i] ;
    int time_of_t = timeOf[g*T + t] ;

    if (tsub >0) {
        cout << "fixing de sous symétries" << endl ;
        cout << "units in group: " ;
        for (int i=FirstG[g] ; i < LastG[g] ; i++) {
            if (SubM.inGroup(i,tsub)) {
                cout << i << " " ;
            }
            cout << endl ;
        }

        if (values[FirstG[g]*T + rankOf[g*T + t-1]] == 1)  { // à quelle valeur est fixée la première unité de g au pas de temps t-1

        }


        cout << "Fixing a " <<  node*left + (node+1)*!left << endl ;
        cout << "Variable branchée : "<< "unit, time : " << unit << ", " << time << endl ;
        cout << " x ? : " << varX << endl ;

        cout << "Fixing of unit i, time t : " << i << ", " << time_of_t << endl ;
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
        cout << endl ;
    }
}

int SubPb::ComputeSubMatrixFixing(Branching & branch, int left, const SubMatrices & SubM, int tsub) { // la classe SubPb contient la variable sur laquelle le branchement va être fait, branch aussi

    int PrintFix=0 ;

    int elague = setXmin(SubM,tsub) ;
    if (elague) {
        return 1 ;
    }
    elague = setXmax(SubM,tsub) ;
    if (elague) {
        return 1 ;
    }

    //si on n'élague pas, on peut faire du fixing
    if (!elague) {
        for (int g = 0 ; g < nbG ; g++) {
            if (SubM.SetForG(g,tsub)) {
                int last_g = LastG[g] ;
                int first_g = FirstG[g] ;

                for (int i = first_g ; i <= last_g ; i++) {
                    if (SubM.inGroup(i, tsub)) {
                        int stop=0 ;
                        int t=0 ;

                        while (!stop && t <= finOrdre[g] ) {
                            if (timeOf[g*T+t] >= tsub) {
                                int time_of_t = timeOf[g*T + t] ;

                                if (Xmax[i*T+t] == Xmin[i*T+t]) {
                                    int bound = Xmax[i*T+t] ;
                                    if (values[i*T+t] == 8 && feasible[i*T+time_of_t] != 2) {

                                        values[i*T+t] = bound ;// Mise à jour de values

                                        if (tsub==0) {
                                            nbFixs ++ ;
                                        }
                                        else {
                                            nbSubFixs++ ;
                                        }

                                        if (PrintFix) {
                                            affichage(left, t, i, bound, SubM, tsub) ;
                                        }
                                        PrintFix=0 ;

                                        if (left) {
                                            //////// Mise à jour de branch left ////////
                                            (branch.varLeft).add( x[i*T + time_of_t] ) ;


                                            (branch.bLeft).add( bound ) ;
                                            if (bound == 1) {
                                                (branch.dirLeft).add( IloCplex::BranchUp ) ;
                                            }
                                            else {
                                                (branch.dirLeft).add( IloCplex::BranchDown ) ;
                                            }
                                        }
                                        else {

                                            //////// Mise à jour de branch right ////////
                                            (branch.varRight).add( x[i*T + time_of_t] );

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

                            }
                            t++ ;
                        } // fin while t
                    } // fin if i in group
                } // for for i
            } // fin if ya un sous ensemble d'unités symétriques dans le groupe g
        }// fin for group g
    }
    return 0 ;
}

void SubPb::doFixing(Branching & branch, int & pruneLeft, int & pruneRight, int subfixing) {

    if (subfixing) {
        updateRSU_RSD() ;
    }

    int left; // est-ce que la branche considérée (côté 0 ou 1) est la branche de gauche

    int subprune ;
    ////////////////////////////// Côté branche à 0 //////////////////////////////
    if (branch.bLeft[0] == 0) {
        left = 1 ;
    }
    else {
        left= 0 ;
    }

    updateValues(0) ;
    // Fixing de la matrice entière
    int pruneSide0 = ComputeSubMatrixFixing(branch, left, Full, 0) ; // fixe le côté var=0, à faire en premier car pour u=1 on modifie 2 valeurs et pour u=0 on modifie au pire 1

    if (pruneSide0) {
        nbFixs ++ ;
    }

    //Fixing des sous matrices
    if (subfixing) {
        int t=1 ;
        while (!pruneSide0 && t < T) {

            if (RSU.SetAtT(t)) {
                subprune = ComputeSubMatrixFixing(branch, left, RSU, t) ;
                pruneSide0 = pruneSide0 || subprune ;
                if (subprune) {
                    nbSubFixs++ ;
                }
            }
            if (!pruneSide0 && RSD.SetAtT(t)) {
                subprune  = ComputeSubMatrixFixing(branch, left, RSD, t) ;
                pruneSide0 = pruneSide0 || subprune ;
                if (subprune) {
                    nbSubFixs++ ;
                }
            }
            t++ ;
        }
    }

    ////////////////////////////// Côté branche à 1 //////////////////////////////
    left=!left ;
    resetValues() ;
    updateValues(1) ;
    // Fixing de la matrice entière
    int pruneSide1 = ComputeSubMatrixFixing(branch, left, Full, 0) ; // fixe le côté var=1

    if (pruneSide1) {
        nbFixs ++ ;
    }


    //Fixing des sous matrices
    if (subfixing) {
        int t=1 ;
        while (!pruneSide1 && t < T) {

            if (RSU.SetAtT(t)) {
                subprune = ComputeSubMatrixFixing(branch, left, RSU, t) ;
                pruneSide1 = pruneSide1 || subprune ;
                if (subprune) {
                    nbSubFixs++ ;
                }
            }
            if (!pruneSide1  && RSD.SetAtT(t)) {
                subprune = ComputeSubMatrixFixing(branch, left, RSD, t) ;
                pruneSide1 = pruneSide1 || subprune ;
                if (subprune) {
                    nbSubFixs++ ;
                }
            }
            t++ ;
        }
    }



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

    if (pruneRight) {
        nbFixs ++ ;
    }
}



