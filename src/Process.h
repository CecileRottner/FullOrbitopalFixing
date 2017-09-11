#ifndef PROCESSS
#define PROCESSS

#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <math.h>
#include <ctime>
#include <list>

using namespace std ;

class Methode {
private:

    int num ; // numéro de référence de la méthode
    // méthodes de résolution
    int doCplex ;
    int doIneqSum ;
    int doIneqVarY ;
    int doDynamicFixing ;
    int doStaticFixing ;
    int doMob ;
    int doIup ; // UP+IUP


    //paramètres
    int stickToCplexFirstBranchingDecision ; // = 1 si Cplex n'a pas le droit de changer de décision de branchement après le branch callback (pour les cas où ne fait rien sur le branchement, afin de comparer)

    int UseLazyCallback ;
    int UseBranchCallback ;
    int UseEmptyBranchCallback ;
    int UseCutCallback ;

    int doSpecialBranching ; // méthode de branchement sur basant sur les pics de demande

    int useNumU ; //u definie comme variable continue dans le modèle

public:

    int allGroups ;
    int StopNode ;

    //Constructeur (par défault c'est Cplex Default)
    Methode() {// crée Cplex default
        doCplex=1 ;
        num=-1 ;

        stickToCplexFirstBranchingDecision = 0 ;
        UseLazyCallback = 0;
        UseBranchCallback = 0 ;
        UseEmptyBranchCallback = 0 ;
        UseCutCallback  = 0;
        doIneqSum = 0;
        doIneqVarY = 0;
        doStaticFixing = 0;
        doDynamicFixing = 0;
        doSpecialBranching = 0 ;
        useNumU = 0 ;
        doMob = 0 ;

        //indicateurs à l'arrache
        StopNode=10 ;
        allGroups=0 ;

    }


    // modification des paramètres
    void setNum(int nb) { num=nb ; }
    void DontUseLazyCB() {UseLazyCallback=0 ;}
    void UseLazyCB() {UseLazyCallback=1 ; }
    void UseEmptyBranchCB() {UseEmptyBranchCallback=1;}

    void DontUseBranchCB() { UseBranchCallback=0 ; }
    void UseBranchCB() {UseBranchCallback=1 ;}

    void DontUseCutCB() { UseCutCallback=0 ; }
    void UseCutCB() {UseCutCallback=1 ;}

    void DoStickToBranchingDecisions() { stickToCplexFirstBranchingDecision=1;}
    void DontStickToBranchingDecisions() { stickToCplexFirstBranchingDecision=0;}



    //Méthodes de résolution
    //paramètres utilisés = le minimum nécessaire

    void CplexCallback(int stick, int lazy, int branch, int cut) { // Cplex callback
        num=0 ;
        stickToCplexFirstBranchingDecision = stick ;
        UseLazyCallback = lazy;
        UseBranchCallback = branch ;
        UseEmptyBranchCallback = branch && !stick ;
        UseCutCallback  = cut;
    }

    void useMOB() {
        num=1 ;
        doCplex=0 ;
        UseBranchCallback = 1 ;
        UseLazyCallback = 1 ;
        stickToCplexFirstBranchingDecision=0 ;
        doMob=1 ;
        useNumU = 1 ;
    }

    void UseStaticFixing() {
        num=2 ;
        doCplex=0 ;
        doStaticFixing=1 ;
        UseBranchCallback = 1 ;
        UseLazyCallback = 1 ;
        stickToCplexFirstBranchingDecision=0 ;
    }

    void UseDynamicFixing() {
        num=4 ;
        doCplex=0 ;
        doDynamicFixing=1 ;
        UseBranchCallback = 1 ;
        UseLazyCallback = 1 ;
        stickToCplexFirstBranchingDecision=0 ;
    }

    void UseIneqSum() { // on peut choisir ou non d'utiliser un lazy callback dans ce cas (même si fixing static en plus)
        num=-3 ;
        doCplex=0 ;
        doIneqSum=1 ;
    }

    void AddIneqSum() { // on peut choisir ou non d'utiliser un lazy callback dans ce cas (même si fixing static en plus)
        doCplex=0 ;
        doIneqSum=1 ;
    }

    void UseIUP() {
        num=10 ;
        doIup=1 ;
        UseCutCallback= 1 ;
    }

    void UseSpecialBranching() {
        UseBranchCallback = 1 ;
        stickToCplexFirstBranchingDecision=0 ; // pas forcément nécessaire en DemandBranching seul, mais peu d'intérêt sans
        doSpecialBranching = 1 ;
    }


    //accès aux paramètres
    int getNum() {return num ;}
    int CplexOnly() {return doCplex;}
    int IneqSum() {return doIneqSum;}
    int IneqVarY() {return doIneqVarY;}
    int DynamicFixing() {return doDynamicFixing;}
    int StaticFixing() {return doStaticFixing;}
    int Mob() {return doMob;}
    int Iup() {return doIup;}
    int LazyCB() {return UseLazyCallback;}
    int BranchCB() {return UseBranchCallback;}
    int CutCB() {return UseCutCallback;}
    int SpecialBranching() {return doSpecialBranching;}
    int StickToBranchingDecision() {return stickToCplexFirstBranchingDecision;}
    int NumU() {return useNumU;}
    int EmptyBranchCB() {return UseEmptyBranchCallback;}

};

class InstanceProcessed {
public :

    //donnees
    int n ;
    int T ;
    int bloc ;
    int demande ;
    int symetrie ;
    int cat ;
    int intra ;
    int id ;

    //Methode
    string localisation ;

    InstanceProcessed(int n_, int T_, int bloc_, int demande_, int symetrie_, int cat_, int intra_, int id_, string loc_)
    {
        n=n_ ;
        T=T_ ;
        bloc = bloc_ ;
        demande = demande_ ;
        symetrie = symetrie_ ;
        cat = cat_ ;
        intra = intra_ ;
        id = id_ ;
        localisation = loc_ ;


    }


    string createName() const ;

    string fileName() {

        string nom = createName() ;
        string fileI = localisation + nom;
        string fileS = fileI + ".txt" ;
        cout << fileS << endl ;
        return fileS ;
    }



};

#endif 
