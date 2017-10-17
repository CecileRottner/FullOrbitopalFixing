[1mdiff --git a/src/InstanceUCP.cc b/src/InstanceUCP.cc[m
[1mindex 2feb4f6..9f4fa07 100644[m
[1m--- a/src/InstanceUCP.cc[m
[1m+++ b/src/InstanceUCP.cc[m
[36m@@ -268,6 +268,8 @@[m [mvoid InstanceUCP::Initialise() {[m
         D += Pmax[np] ;[m
         np++ ;[m
     }*/[m
[32m+[m
[32m+[m[32m    cout << "Pmax: "<< Pmax << endl ;[m
 }[m
 [m
 void InstanceUCP::Lecture(const char* file) {[m
[1mdiff --git a/src/Process.h b/src/Process.h[m
[1mindex 41a9ab0..182156d 100644[m
[1m--- a/src/Process.h[m
[1m+++ b/src/Process.h[m
[36m@@ -63,6 +63,7 @@[m [mpublic:[m
         doSpecialBranching = 0 ;[m
         useNumU = 0 ;[m
         doMob = 0 ;[m
[32m+[m[32m        doIup = 0 ;[m
 [m
         //indicateurs à l'arrache[m
         StopNode=10 ;[m
[36m@@ -163,6 +164,24 @@[m [mpublic:[m
         doSpecialBranching = 1 ;[m
     }[m
 [m
[32m+[m[32m    void printParam()  {[m
[32m+[m[32m        cout << "method nb: " << num << endl ;[m
[32m+[m[32m        cout << "Cplex only: " << doCplex << endl ;[m
[32m+[m[32m        cout << "Use sub symmetry ineq: " << doIneqSum << endl ;[m
[32m+[m[32m        cout << "Use variables y: " << doIneqVarY<< endl ;[m
[32m+[m[32m        cout << "Agregated model: " << doAggregatedModel << endl ;[m
[32m+[m[32m        cout << "Dynamic fixing: " << doDynamicFixing << endl ;[m
[32m+[m[32m        cout << "Static fixing: " << doStaticFixing << endl ;[m
[32m+[m[32m        cout << "MOB: " << doMob << endl ;[m
[32m+[m[32m        cout << "Use IUP: " << doIup << endl ;[m
[32m+[m[32m        cout << "Lazy Callback: " << UseLazyCallback << endl ;[m
[32m+[m[32m        cout << "Branch Callback: " << UseBranchCallback << endl ;[m
[32m+[m[32m        cout << "Use empty Branch Callback: " << UseEmptyBranchCallback << endl ;[m
[32m+[m[32m        cout << "Cut Callback: " << UseCutCallback << endl ;[m
[32m+[m[32m        cout << "Special branching: " << doSpecialBranching << endl ;[m
[32m+[m[32m        cout << "Stick to Cplex first branching decision: " << stickToCplexFirstBranchingDecision << endl ;[m
[32m+[m[32m        cout << "Use continuous u variable: " << useNumU << endl ;[m
[32m+[m[32m    }[m
 [m
     //accès aux paramètres[m
     int getNum() {return num ;}[m
[1mdiff --git a/src/main.cc b/src/main.cc[m
[1mindex 0761298..6a5c3d6 100644[m
[1m--- a/src/main.cc[m
[1m+++ b/src/main.cc[m
[36m@@ -255,6 +255,7 @@[m [mint process(InstanceProcessed I, ofstream & fichier, double & time, Methode met,[m
     }[m
     else {[m
         model = defineModel(env,inst,x,u,met.NumU()) ;[m
[32m+[m[32m        cout << "ici modèle" << endl ;[m
     }[m
 [m
     IloCplex cplex(model) ;[m
[36m@@ -413,6 +414,14 @@[m [mmain(int argc,char**argv)[m
         if (met==1) {[m
 [m
             env=IloEnv() ;[m
[32m+[m[32m            process(Instance, fichier, time, DefaultCplex, env) ;[m
[32m+[m[32m            env.end() ;[m
[32m+[m
[32m+[m[32m            env=IloEnv() ;[m
[32m+[m[32m            process(Instance, fichier, time, CBCplex, env) ;[m
[32m+[m[32m            env.end() ;[m
[32m+[m
[32m+[m[32m            env=IloEnv() ;[m
             process(Instance, fichier, time, StaticFix, env) ;[m
             env.end() ;[m
 [m
[36m@@ -501,15 +510,21 @@[m [mmain(int argc,char**argv)[m
 [m
         for (sym= 4 ; sym >= 2 ; sym--) {[m
             Instance.symetrie = sym ;[m
[31m-            for (int id=1; id <=20; id++) {[m
[32m+[m[32m            for (int id=3; id <=20; id++) {[m
                 Instance.id = id ;[m
 [m
[31m-                env=IloEnv() ;[m
[32m+[m[32m               /*env=IloEnv() ;[m
                 process(Instance, fichier, time, DefaultCplex, env) ;[m
[32m+[m[32m                env.end() ;*/[m
[32m+[m
[32m+[m[32m                DynamicFix.printParam() ;[m
[32m+[m
[32m+[m
[32m+[m[32m                 env=IloEnv() ;[m
[32m+[m[32m                process(Instance, fichier, time, DynamicFix, env) ;[m
                 env.end() ;[m
[31m-                env=IloEnv() ;[m
[31m-                process(Instance, fichier, time, AggregModel, env) ;[m
[31m-                env.end() ;[m
[32m+[m
[32m+[m
 [m
                 /*env=IloEnv() ;[m
                 process(Instance, fichier, time, IneqVarY, env) ;[m
