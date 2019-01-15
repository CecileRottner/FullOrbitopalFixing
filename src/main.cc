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
#include "ModeleFlot.h"
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

        int T = sub.T ;

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

        int unit;
        int time;
        int varX ;
        sub.getVar(varID, unit, time, varX) ;


        //////////// Mise à jour de nodeData (ordre des variables) /////////////
        if (getNnodes() != 0 && methode.DynamicFixing() ) { // si on est à la racine, pas de data dans le noeud courant, data est "vide. // si on utilise une methode à ordre fixe on peut garder le même nodeData tout le long
            data = dynamic_cast <myNodeData *> (getNodeData()) ; // sinon on récupère le nodeData du noeud courant
        }


        int father = data->num ;
        if (father==0 && !methode.CplexOnly() ) {
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

        if ( !methode.CplexOnly() ) {
            getUBs(sub.UB, sub.x) ;
            getLBs(sub.LB, sub.x) ;
            sub.update(varID, dataFils);
        }

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

            if ( !methode.SpecialBranching() ) {

                // Fixing avec branchement cplex (static ou dynamic)

                getBranch(branch.varLeft, branch.bLeft, branch.dirLeft, 0) ;
                getBranch(branch.varRight, branch.bRight, branch.dirRight, 1) ;
                sub.doFixing(branch, pruneLeft, pruneRight, methode.SubFixing()) ;
            }

            else {
                int newVar = 0 ;
                int stopNode = 10 ;
                if (getNnodes() < stopNode ) {

                    getValues(sub.x_frac,sub.x);
                    if (methode.allGroups) {
                        newVar = sub.newVarFromFractionalGroup(branch, getNnodes()) ;
                    }
                    else  {
                        newVar = sub.newVarSameTimePeriod(branch, getNnodes()) ;
                    }
                }
                if (!newVar) { // la variable choisie par cplex est conservée
                    getBranch(branch.varLeft, branch.bLeft, branch.dirLeft, 0) ;
                    getBranch(branch.varRight, branch.bRight, branch.dirRight, 1) ;
                }

                sub.doFixing(branch, pruneLeft, pruneRight, methode.SubFixing()) ;

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

    cout  << "ici" << endl ;

    // cout << "ici : " << met.getNum() << endl ;

    string nom = I.fileName() ;
    const char* file = nom.c_str() ;


    int id=I.id ;

    //IloEnv env ;
    InstanceUCP* inst = new InstanceUCP(env, file) ;

    int n = inst->n ;
    int T = inst->T ;

    IloBoolVarArray x(env,n*T);
    IloBoolVarArray u(env,n*T);

    IloBoolVarArray f(env,n*(T+2)*(T+2));

    IloModel model ;

    int ramp = met.Ramping();


    if (met.IneqVarY()) {
        if (!met.IneqSum()) {
            model = defineModel_y(env,inst,x,u) ;
        }
        else {
            model = defineModel_sum(env,inst, x,u, -5) ;
        }
    }
    else if (met.IneqSum()) {
        model = defineModel_sum(env,inst, x,u, -3) ;

    }

    else if (met.NumberOfOnes()) {
        model = defineModel_numberOfOnes(env,inst, x,u) ;
    }

    else if (met.AggregatedModel()) {
        model = AggregatedModel(env, inst) ;
    }

    else if (met.ModeleFlot()) {
        ModeleFlot flot(env, inst) ;
        model = flot.AggregatedFlowModel();
    }

    else {
        model = defineModel(env,inst,x,u,met.NumU(), ramp) ;
        if (met.RSUonly()) {
            AddRSUIneq(model, env, inst, x, u,0);
        }
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

    cout << "feasible : " << cplex.isPrimalFeasible() << endl ;
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

    ///Affichage solution optimale
    /*if (met.ModeleFlot()) {
        IloNumArray flots(env, 0) ;
        cplex.getValues(flots, f) ;

        for (int i=0; i<n; i++) {
            for (int t=0 ; t < T+2 ; t++) {
                for (int k = 0 ; k < T+2; k++) {
                    if (flots[t + (T+2)*k + (T+2)*(T+2)*i] >0) {
                        cout << "arc: " << i << ", " << t << ", " << k << endl ;
                    }
                }
            }
        }
    }

    else {
        cplex.getValues(sub.x_frac, x) ;
        for (int i=0; i<n; i++) {
            cout << "unité " << i << " : " ;
            for (int t=0 ; t < T ; t++) {
                cout << sub.x_frac[i*T +t] << " " ;
            }
        }
        cout << endl ;
    }
*/
    double t = cplex.getCplexTime();
    double opt = cplex.getObjValue() ;

    fichier << met.getNum() <<  " & " << n << " & " << T  << " & " << I.symetrie << " & " << inst->nbG  << " & " << inst->MaxSize << " & " << inst->MeanSize  << " & " << id ;
    fichier << " & " << cplex.getObjValue()  ; //Optimal value
    fichier << " & " << cplex.getMIPRelativeGap() << " \\% " ; //approx gap
    fichier << " & " << cplex.getNnodes() ;
    fichier << " & " << sub.nbFixs ;
    fichier << " & " << sub.nbSubFixs ;
    fichier << " & " << sub.timeFix ;
    fichier << " & " << t - time ;
    /*if (sub.UB_LB) {
        fichier << "        UB_LB     "  ;
    }*/
    fichier <<" \\\\ " << endl ;



    time = t ;



    //Destructeurs
   // delete inst ;
   // delete dataNode ;
    // env.end() ;

    cout << "ici" << endl ;

    return 1;
}


int
main(int argc,char**argv)
{
    srand(time(NULL));

    //définition des méthodes de résolution
    Methode DefaultCplex ;

    Methode CBCplex ;
    CBCplex.CplexCallback(1,1,1,0);

    Methode IneqPures;
    IneqPures.UseIneqSum();

    Methode RampModel;
    RampModel.UseRampConstraints();

    Methode RampIneqRSU;
    RampIneqRSU.UseRampConstraints();
    RampIneqRSU.AddIneqRSU() ;

    Methode IneqVarY;
    IneqVarY.UseIneqVarY();

    Methode Flot;
    Flot.UseModeleFlot();

    Methode IneqNumberOfOnes;
    IneqNumberOfOnes.UseNumberOfOnes();

    Methode IneqSumAndVarY;
    IneqSumAndVarY.UseIneqVarY();
    IneqSumAndVarY.UseIneqSum();
    IneqSumAndVarY.setNum(-5);

    Methode AggregModel;
    AggregModel.UseAggregatedModel();

    Methode IneqCB ;
    IneqCB.UseIneqSum() ;
    IneqCB.UseBranchCB();
    //IneqCB.DoStickToBranchingDecisions() ;
    IneqCB.UseEmptyBranchCB() ;
    IneqCB.setNum(0) ;

    Methode Mob ;
    Mob.useMOB() ;

    Methode StaticFix ;
    StaticFix.UseStaticFixing() ;
    //StaticFix.AddIneqSum() ;
    //StaticFix.DontUseLazyCB();

    Methode StaticSubFix ;
    StaticSubFix.UseStaticFixing() ;
    StaticSubFix.UseSubFixing() ;
    StaticSubFix.setNum(22);

    Methode StaticFixWithIneq ;
    StaticFixWithIneq.UseStaticFixing() ;
    StaticFixWithIneq.AddIneqSum() ;
    StaticFixWithIneq.DontUseLazyCB();

    Methode StaticFixWithBranching ;
    StaticFixWithBranching.UseStaticFixing() ;
    StaticFixWithBranching.AddIneqSum() ;
    StaticFixWithBranching.DontUseLazyCB();
    StaticFixWithBranching.UseSpecialBranching();

    Methode StaticFixWithBranching_50 ;
    StaticFixWithBranching_50.UseStaticFixing() ;
    StaticFixWithBranching_50.AddIneqSum() ;
    StaticFixWithBranching_50.DontUseLazyCB();
    StaticFixWithBranching_50.UseSpecialBranching();
   // StaticFixWithBranching_50.StopNode = 50;


    Methode StaticFixWithBranching_all ;
    StaticFixWithBranching_all.UseStaticFixing() ;
    StaticFixWithBranching_all.AddIneqSum() ;
    StaticFixWithBranching_all.DontUseLazyCB();
    StaticFixWithBranching_all.UseSpecialBranching();
    StaticFixWithBranching_all.allGroups=1 ;

    Methode DynamicFix ;
    DynamicFix.UseDynamicFixing() ;

    Methode DynamicSubFix ;
    DynamicSubFix.UseDynamicFixing() ;
    DynamicSubFix.UseSubFixing() ;
    DynamicSubFix.setNum(44);

    Methode DynamicSubFixWithRamps ;
    DynamicSubFixWithRamps.UseRampConstraints();
    DynamicSubFixWithRamps.UseDynamicFixing() ;
    DynamicSubFixWithRamps.UseSubFixing() ;
    DynamicSubFixWithRamps.setNum(55);

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
        IloEnv env ;

        if (met==1) {

            env=IloEnv() ;
            process(Instance, fichier, time, DefaultCplex, env) ;
            env.end() ;

            env=IloEnv() ;
            process(Instance, fichier, time, CBCplex, env) ;
            env.end() ;

            env=IloEnv() ;
            process(Instance, fichier, time, Mob, env) ;
            env.end() ;

            env=IloEnv() ;
            process(Instance, fichier, time, DynamicFix, env) ;
            env.end() ;
        }

        if (met==2) {
            env=IloEnv() ;
            process(Instance, fichier, time, DefaultCplex, env) ;
            env.end() ;

            env=IloEnv() ;
            process(Instance, fichier, time, IneqPures, env) ;
            env.end() ;

            env=IloEnv() ;
            process(Instance, fichier, time, StaticFixWithIneq, env) ;
            env.end() ;
        }

        if (met==3) {

            env=IloEnv() ;
            int opt= process(Instance, fichier, time, DefaultCplex , env) ;
            env.end() ;

            env=IloEnv() ;
            process(Instance, fichier, time, DynamicFix , env) ;
            env.end() ;

            env=IloEnv() ;
            int solution_fixing_static = process(Instance, fichier, time, StaticSubFix , env) ;
            env.end() ;

            env=IloEnv() ;
            int solution_fixing = process(Instance, fichier, time, DynamicSubFix , env) ;
            env.end() ;

            if (fabs(opt-solution_fixing) > 0.0000001 ) {
                fichier << "ERREUR" << endl ;
            }
            if (fabs(opt-solution_fixing_static) > 0.0000001 ) {
                fichier << "ERREUR" << endl ;
            }
        }

        if (met==4) {
            env=IloEnv() ;
            process(Instance, fichier, time, RampIneqRSU , env) ;
            env.end() ;
        }
        if (met==5) {
            env=IloEnv() ;
            process(Instance, fichier, time, DynamicSubFixWithRamps, env) ;
            env.end() ;
        }

        /*env=IloEnv() ;
        process(Instance, fichier, time, IneqPures, env) ;
        env.end() ;

        env=IloEnv() ;
        process(Instance, fichier, time, IneqCB, env) ;
        env.end() ;

        env=IloEnv() ;
        process(Instance, fichier, time, StaticFix, env) ;
        env.end() ;*/

        /* env=IloEnv() ;
        process(Instance, fichier, time, StaticFixWithBranching, env) ;
        env.end() ;*/


        /*   env=IloEnv() ;
        process(Instance, fichier, time, StaticFixWithBranching_50, env) ;
        env.end() ;


        env=IloEnv() ;
        process(Instance, fichier, time, StaticFixWithBranching_all, env) ;
        env.end() ; */


    }


    ////////////////////////// SI PAS D'ARGUMENTS ////////////////////
    if (argc==1) {
        ofstream fichier("result.txt");

        //fichier << "Instance & n & T & Sym & nG & max & mean & OptVal & RootRelax & ApproxGap &  USCuts & IUSCuts & SepTime & Prof & Nodes & CplexCuts & CPU \\\\ " << endl;
        fichier << "Methode & n & T & OptVal & Nodes & NbFixs & TimeFix & CPU \\\\ " << endl;

        double time = 0 ;

        //Paramètres de l'instance

        int T = 48;
        int n = 60 ;
        int sym = 2 ;
        int demande = 3;
        int cat01 = 0;
        int bloc = 1;
        int intra = 0 ;
        string localisation = "data/smaller/" ;
        InstanceProcessed Instance = InstanceProcessed(n, T, bloc, demande, sym, cat01, intra, 0, localisation) ;

        fichier << localisation << endl ;
        Instance.localisation = localisation ;

        n=20;
        T=24;
        Instance.n=n;
        Instance.T=T ;
        IloEnv env ;

        for (sym= 3; sym >=2 ; sym--) {
            Instance.symetrie = sym ;

            for (int id=1; id <=20; id++) {
                Instance.id = id ;




                env=IloEnv() ;
                cout <<"start ramp model" << endl ;
                process(Instance, fichier, time, RampIneqRSU , env) ;
                cout <<"end ramp model" << endl ;
                env.end() ;


                env=IloEnv() ;
                process(Instance, fichier, time, DynamicSubFixWithRamps, env) ;
                env.end() ;


                /* env=IloEnv() ;
                int solution_fixing = process(Instance, fichier, time, DynamicSubFix , env) ;
                env.end() ;*/

                /*  if (fabs(opt-solution_fixing) > 0.0000001 ) {
                    fichier << "ERREUR" << endl ;
                }*/

                /*env=IloEnv() ;
                process(Instance, fichier, time, IneqCB, env) ;
                env.end() ;
*/
                /* env=IloEnv() ;
                process(Instance, fichier, time, StaticFixWithBranching, env) ;
                env.end() ;


                env=IloEnv() ;
                process(Instance, fichier, time, StaticFixWithBranching_50, env) ;
                env.end() ;


                env=IloEnv() ;
                process(Instance, fichier, time, StaticFixWithBranching_all, env) ;
                env.end() ;*/

                fichier << endl ;
            }

            fichier << endl ;
            fichier << endl ;
        }
    }



    return 0 ;
}
