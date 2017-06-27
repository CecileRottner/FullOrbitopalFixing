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
#include "Process.h"
#include "Fixing.h"
#include "Mob.h"

using namespace std ;




ILOSTLBEGIN

ILOBRANCHCALLBACK4(BCallBack, Methode &, methode, SubPb &, sub, myNodeData*, data, int, count) {


    clock_t start;

    int nBranches = getNbranches() ;

    if (nBranches > 0 && !methode.EmptyBranchCB() ) {

        start = clock();

        if (nBranches >2) {
            cout<< "nb branches > 2 " << endl ;
        }
        if (nBranches ==1) {
            cout<< "nb branches = 1 " << endl ;
        }

        //////////Récupération variable branchée par Cplex //////////////:
        IloNumVarArray vars = IloNumVarArray(getEnv(), 0) ;
        IloNumArray bounds = IloNumArray(getEnv(), 0) ;
        IloCplex::BranchDirectionArray dirs = IloCplex::BranchDirectionArray(getEnv(), 0);
        getBranch(vars, bounds, dirs, 0) ; // on récupère le premier fils, sachant qu'en pratique le deuxième fils c'est la même variable, fixée à la valeur opposée.

        int varID = vars[0].getId() ;

        int T=sub.T ;
        int nX = sub.n*sub.T ;
        int var=varID ;
        int varX = 1 ;
        int varP = 0 ;

        if (varID >= nX ) { // une variable u(i,t) est fixée
            var -= nX ;
            varX=0 ;
        }
        if (var >= nX) {
            var -= nX ;
            varP=1 ;
        }
        int time = var % sub.T ;
        int unit = (var - time)/ sub.T ;



        //////////// Mise à jour de nodeData (ordre des variables) /////////////
        if (getNnodes() != 0 && methode.DynamicFixing() ) { // si on est à la racine, pas de data dans le noeud courant, data est "vide. // si on utilise une methode à ordre fixe on peut garder le même nodeData tout le long
            data = dynamic_cast <myNodeData *> (getNodeData()) ; // sinon on récupère le nodeData du noeud courant
        }

        int father = data->num ;
        if (father==0) {
            getFeasibilities(sub.feasible, sub.x) ;
        }
        //implied feasible variables

        myNodeData* dataFils = data ;
        myNodeData* dataFils2 = data ;

        if ( methode.DynamicFixing() ) {
            dataFils = new myNodeData(getEnv(), data, sub.Group[unit], time, varX, sub.feasible[unit*T+time], sub.T, sub.nbG) ;
            dataFils2 = new myNodeData(getEnv(), data, sub.Group[unit], time, varX, sub.feasible[unit*T+time], sub.T, sub.nbG) ;

        }
        count++ ;
        dataFils->num=count ;
        count++ ;
        dataFils2->num=count ;


        //////////// Mise à jour de sub ///////////////
        getUBs(sub.UB, sub.x) ;
        getLBs(sub.LB, sub.x) ;
        sub.update(varID, dataFils);

        /////////// Méthodes ////////////
        Branching branch= Branching(getEnv()) ;
        int pruneLeft = 0;
        int pruneRight = 0;

        if ( methode.StickToBranchingDecision() ) {
            getBranch(branch.varLeft, branch.bLeft, branch.dirLeft, 0) ;
            getBranch(branch.varRight, branch.bRight, branch.dirRight, 1) ;
        }

        if ( methode.Mob() ) { // MOB
            getValues(sub.x_frac,sub.x); // MOB : on a besoin des valeurs fractionnaires de x
            doMOB(branch, sub) ;
        }


        if ( methode.StaticFixing() || methode.DynamicFixing() ) { // fixing

            if ( !methode.DemandBranching() ) {

            // Fixing avec branchement cplex (static ou dynamic)
            getBranch(branch.varLeft, branch.bLeft, branch.dirLeft, 0) ;
            getBranch(branch.varRight, branch.bRight, branch.dirRight, 1) ;
            sub.doFixing(branch, pruneLeft, pruneRight, getNnodes()) ;
            }

            else {
                int newVar = sub.newVarFromDemand(branch, getNnodes()) ;
                if (!newVar) { // la variable choisie par cplex est conservée
                    getBranch(branch.varLeft, branch.bLeft, branch.dirLeft, 0) ;
                    getBranch(branch.varRight, branch.bRight, branch.dirRight, 1) ;
                }
                sub.doFixing(branch, pruneLeft, pruneRight, getNnodes()) ;
            }
        }


        int print = 0 ; //que pour le fixing dynamic
        int printRoot=0 ;


        if (father==0 && printRoot) {
            for (int t=0 ; t < T ; t++) {
                cout << "t=" << t << ", D= " << sub.instance->getD(t) <<"    " ;
                for (int i=0 ; i <sub.n; i++) {
                    cout << "(" << sub.feasible[i*T + t] << ", " <<getFeasibility(sub.u[i*T + t]) << ") " ;
                    if (sub.Last[i]) {
                        cout << "       " ;
                    }
                }
                cout << endl ;
            }
            cout << endl;
        }

        if (!pruneLeft ) {
            if (print) {
                cout << "Noeud " << father << " pere de " << dataFils->num << endl ;
                cout << "Variable branchée : " ;
                if (sub.varX) {
                    cout << "x" ;
                }
                else {
                    cout << "u" ;
                }
                cout << "(i,t) = " <<  sub.unit << ", " << sub.time  << " = " << branch.bLeft[0] << endl ;
                cout << "varLeft : " << branch.varLeft << endl ;
                cout << "feasibility u : " << getFeasibility(sub.u[unit*T + time]) << endl ;
                cout << "feasibility x : " << getFeasibility(sub.x[unit*T + time]) << endl  ;
                cout << "group : " << sub.FirstG[sub.Group[sub.unit]] << ", " << sub.LastG[sub.Group[sub.unit]] << endl ;
                cout << endl ;
            }

            if (methode.DynamicFixing() ) {
                makeBranch(branch.varLeft, branch.bLeft, branch.dirLeft, getObjValue(), dataFils) ;
            }
            else {
                makeBranch(branch.varLeft, branch.bLeft, branch.dirLeft, getObjValue()) ;
            }
        }


        if (!pruneRight) {
            if (print) {
                cout << "Noeud " << father << " pere de " << dataFils2->num << endl ;
                cout << "Variable branchée : " ;
                if (sub.varX) {
                    cout << "x" ;
                }
                else {
                    cout << "u" ;
                }
                cout << "(i,t) = " << sub.unit << ", " << sub.time << " = " << branch.bRight[0] << endl ;
                cout << "varRight : " << branch.varRight << endl ;
                cout << "feasibility u : " << getFeasibility(sub.u[unit*T + time]) << endl ;
                cout << "feasibility x : " << getFeasibility(sub.x[unit*T + time]) << endl  ;
                cout << "group : " << sub.FirstG[sub.Group[sub.unit]] << ", " << sub.LastG[sub.Group[sub.unit]] << endl ;
                /*cout << "LB, UB : " << getLB(sub.u[unit*T + sub.time]) << ", " << getUB(sub.u[unit*T + sub.time]) << endl ;
                cout << "xt : "<< getLB(sub.x[unit*T + sub.time]) << ", " << getUB(sub.x[unit*T + sub.time]) << endl ;
                cout << "xt-1 : "<< getLB(sub.x[unit*T + sub.time-1]) << ", " << getUB(sub.x[unit*T + sub.time-1]) << endl ;*/
                cout << endl ;
            }
            if ( methode.DynamicFixing() ) {
                makeBranch(branch.varRight, branch.bRight, branch.dirRight, getObjValue(), dataFils2) ;
            }
            else {
                makeBranch(branch.varRight, branch.bRight, branch.dirRight, getObjValue()) ;
            }
        }

        if (pruneRight && pruneLeft) {
            prune() ;
        }

        vars.end() ;
        bounds.end() ;
        dirs.end() ;

        sub.timeFix += ( clock() - start ) / (double) CLOCKS_PER_SEC;
    }




}
//FIN BRANCH CALLBACK


ILOLAZYCONSTRAINTCALLBACK1(LazyCB, IloInt, fake) {}

int process(InstanceProcessed I, ofstream & fichier, double & time, Methode met, IloEnv env) {

    string nom = I.fileName() ;
    const char* file = nom.c_str() ;


    int id=I.id ;

    //IloEnv env ;
    InstanceUCP* inst = new InstanceUCP(env, file) ;

    int n = inst->n ;
    int T = inst->T ;

    IloBoolVarArray x(env,n*T);
    IloBoolVarArray u(env,n*T);

    IloModel model ;
    if (met.IneqVarY()) {
        model = defineModel_y(env,inst,x,u) ;
    }
    else if (met.IneqSum()) {
        model = defineModel_sum(env,inst, x,u, -3) ;
    }
    else {
        model = defineModel(env,inst,x,u,met.NumU()) ;
    }

    IloCplex cplex(model) ;
    IloNum eps = cplex.getParam(IloCplex::Param::Simplex::Tolerances::Feasibility) ;

    SubPb sub = SubPb(env, inst, x, u, eps, met) ;

    myNodeData* dataNode ;
    dataNode = new myNodeData(env, T, sub.nbG, met) ;
    int count=0 ;

    if (met.BranchCB()) {
        cplex.use(BCallBack(env, met, sub, dataNode, count)) ;
    }
    if (met.LazyCB()) {
        cplex.use(LazyCB(env, 0)) ;
    }

    //Paramètres
    cplex.setParam(IloCplex::Param::ClockType, 1); //1 : CPU TIME
    cplex.setParam(IloCplex::Param::Threads, 1);
    cplex.setParam(IloCplex::EpGap, 0.0000001) ;
    cplex.setParam(IloCplex::Param::TimeLimit, 3600) ;


    //Résolution et affichage de la solution
    cplex.solve();

    /*cplex.getValues(sub.x_frac, x) ;

    for (int t=0 ; t < T ; t++) {
        cout << "t=" << t << "    " ;
        for (int i= 0 ; i <sub.n ; i++) {
            cout << abs(sub.x_frac[i*T + t]) << " " ;
            if (sub.Last[i]) {
                cout << "         " ;
            }
        }
        cout << endl ;
    }
    cout << endl;*/

    double t = cplex.getCplexTime();
    double opt = cplex.getObjValue() ;

    fichier << met.getNum() <<  " & " << n << " & " << T  << " & " << I.symetrie << " & " << inst->nbG  << " & " << inst->MaxSize << " & " << inst->MeanSize  << " & " << id ;
    fichier << " & " << cplex.getObjValue()  ; //Optimal value
    fichier << " & " << cplex.getMIPRelativeGap() << " \\% " ; //approx gap
    fichier << " & " << cplex.getNnodes() ;
    fichier << " & " << sub.nbFixs ;
    fichier << " & " << sub.timeFix ;
    fichier << " & " << t - time ;
    if (sub.UB_LB) {
        fichier << "        UB_LB     "  ;
    }
    fichier <<" \\\\ " << endl ;

    time = t ;

    //Destructeurs
    delete inst ;
    delete dataNode ;
    // env.end() ;
    return opt;
}


int
main(int argc,char**argv)
{
    srand(time(NULL));

    //définition des méthodes de résolution
    Methode DefaultCplex ;

    Methode CBCplex ;
    CBCplex.CplexCallback(1,0,1,0);

    Methode IneqPures;
    IneqPures.UseIneqSum();

    Methode IneqCB ;
    IneqCB.UseIneqSum() ;
    IneqCB.UseBranchCB();
    //IneqCB.DoStickToBranchingDecisions() ;
    IneqCB.UseEmptyBranchCB() ;
    IneqCB.setNum(0) ;


    Methode StaticFix ;
    StaticFix.UseStaticFixing() ;
    StaticFix.AddIneqSum() ;
    StaticFix.DontUseLazyCB();

    Methode DynamicFix ;
    DynamicFix.UseDynamicFixing() ;


    /////////////////////// SI ARGUMENTS //////////////////////
    if (argc>1) {
        ofstream fichier("result.txt", std::ofstream::out | std::ofstream::app);

        int met= atoi(argv[1]);
        string localisation = argv[2] ;

        int n = atoi(argv[3]);
        int T = atoi(argv[4]);
        int bloc = atoi(argv[5]);
        int demande = atoi(argv[6]);
        int sym = atoi(argv[7]);

        int cat01 = atoi(argv[8]);
        int intra = atoi(argv[9]);

        int id = atoi(argv[10]);

        InstanceProcessed Instance = InstanceProcessed(n, T, bloc, demande, sym, cat01, intra, id, localisation) ;

        double time = 0 ;



        fichier << endl ;

    }


    ////////////////////////// SI PAS D'ARGUMENTS ////////////////////
    if (argc==1) {
        ofstream fichier("result.txt");

        //fichier << "Instance & n & T & Sym & nG & max & mean & OptVal & RootRelax & ApproxGap &  USCuts & IUSCuts & SepTime & Prof & Nodes & CplexCuts & CPU \\\\ " << endl;
        fichier << "Methode & n & T & OptVal & Nodes & NbFixs & TimeFix & CPU \\\\ " << endl;

        double time = 0 ;

        //Paramètres de l'instance

        int T = 96;
        int n = 30 ;
        int sym = 2 ;
        int demande = 3;
        int cat01 = 0;
        int bloc = 1;
        int intra = 0 ;
        string localisation = "data/Litt_Real/" ;
        InstanceProcessed Instance = InstanceProcessed(n, T, bloc, demande, sym, cat01, intra, 0, localisation) ;


        localisation = "data/Litt_Real/" ;
        fichier << localisation << endl ;
        Instance.localisation = localisation ;

        n=60 ;
        T=48;
        Instance.n=n;
        Instance.T=T ;
        IloEnv env ;

        for (sym= 4 ; sym >= 2; sym--) {
            Instance.symetrie = sym ;
            for (int id=9; id <=20; id++) {
                Instance.id = id ;

                env=IloEnv() ;
                process(Instance, fichier, time, DefaultCplex, env) ;
                env.end() ;

                env=IloEnv() ;
                process(Instance, fichier, time, IneqPures, env) ;
                env.end() ;

                env=IloEnv() ;
                process(Instance, fichier, time, IneqCB, env) ;
                env.end() ;

                env=IloEnv() ;
                process(Instance, fichier, time, StaticFix, env) ;
                env.end() ;

                fichier << endl ;
            }
            fichier << endl ;
            fichier << endl ;
        }
    }



    return 0 ;
}
