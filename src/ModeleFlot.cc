#include "ModeleFlot.h"
#include "InstanceUCP.h"

#include <ilcplex/ilocplex.h>

using namespace std ;

ModeleFlot::ModeleFlot(IloEnv enviro, InstanceUCP* pbm) {
    env=enviro ;
    pb=pbm ;
    n = pb->getn();
    T = pb->getT() ;
}

///// Description du graphe Gi
/// arc (t,k) correspond à :
///     - (ut, vk) si t < k
///     - (vk, ut) si t > k
/// sachant que s (source) considérée comme u0 si etat initial 'up' et comme v0 si etat initial 'down'
/// et que p (puits) est u(T+1) et v(T+1) à la fois
/// note: l'arc (s,p) est noté (0, T+1) si etat initial 'up' (arc u0 -> v(T+1)) et noté (T+1,0) si état initial 'down' (arc v0 -> u(T+1))


int ModeleFlot::Adj(int i, int t, int k) { // matrice d'adjacence du graphe : est-ce que l'arc (t,k) existe dans Gi ?
    if (t==k) {
        return 0 ;
    }

    // Etats initiaux
    if (!pb->getInit(i)) { // état initial à 0
        if (t==0) { // A(i, 0, k) = 0 pour tout k \in [0, T+1]
            return 0 ;
        }
    }
    else { // état initial à 1
        if (k==0) { // A(i,t,0) = 0 pour tout t \in [0, T+1]
            return 0 ;
        }
    }

    // Temps min de marche ; pas d'arc si k \in [t+1, t+Li-1]
    if (k >= t+1 && k <= t + pb->getL(i) - 1) {
        return 0 ;
    }

    //Temps min d'arrêt ; pas d'arc si t \in [k+1, k+li-1]
    if (t >= k+1 && t <= k + pb->getl(i) - 1 ) {
        return 0 ;
    }

    return 1 ;
}

double ModeleFlot::Cost(int i, int t, int k) { // Cout de l'arc (t,k) dans Gi (0 si l'arc n'existe pas)
    if (t==k) {return 0 ;}

    else if (t > k) { // cas d'un arc (vk, ut)
        if (k <= T) { // on démarre entre 1 et T (exlut T+1 qui correspond au puits)
            return pb->getc0(i)*Adj(i,t,k) ;
        }
    }

    else { // t < k, cas d'un arc (ut, vk)
        if (t==0) { // cas où état initial est up, on ne compte pas le prix lié au pas de temps 0
            return pb->getcf(i)*(k-t-1)*Adj(i,t,k) ;
        }
        else {
            return pb->getcf(i)*(k-t)*Adj(i,t,k) ;
        }
    }
    return 0 ;
}

int ModeleFlot::arc(int i, int t, int k) { // numéro de l'arc (t,k) de Gi
    return (t + (T+2)* k + (T+2)*(T+2)* i) ;
}



IloModel ModeleFlot::defineModelFlot() {

    cout << "start model" << endl ;
    IloModel model ;

    IloBoolVarArray f(env, n*(T+2)*(T+2)) ;
    IloNumVarArray pp(env, n*T, 0.0, 1000);




    //Vecteur x
    IloExprArray x(env) ;
    for (int i = 0 ; i < n ; i++) {
        for (int t=0 ; t < T ; t++) {
            IloExpr xit(env) ;

            // calcul de xit : unité i en marche à t+1
            for (int s = 0 ; s <= t+1 ; s++) {
                for (int z = t+2 ; z <= T+1 ; z++) {
                    xit +=  Adj(i,s,z)*f[arc(i,s,z)]   ;
                }
            }
            x.add(xit) ;
            xit.end() ;
        }
    }

    cout << "expr xit" << endl ;


    // Objective Function: Minimize Cost
    IloExpr cost(env) ;
    for (int i=0; i<n; i++) {
        for (int t=0 ; t < T ; t++) {
            cost += (pp[i*T + t]+pb->getP(i)*x[i*T + t])*(pb->getcp(i)) ;
        }
        for (int t=0 ; t < T+2 ; t++) {
            for (int k = 0 ; k < T+2; k++) {
                cost += Cost(i,t,k)*f[arc(i,t,k)] ;
            }
        }
        cout << "cout pour unité " << i << endl ;
    }
    model.add(IloMinimize(env, cost));

    cout << "fonction objectif." << endl ;


    //pour ut, t \in [1,T] : flot entrant = flot sortant
    for (int i=0; i<n; i++) {
        for (int t=1 ; t <= T ; t++) {

            IloExpr FlowUt(env) ;

            for (int s = 0 ; s <= t-1 ; s++) {
                FlowUt += Adj(i,t,s)*f[arc(i,t,s)] ;
            }
            for (int s = t+1 ; s <= T ; s++) {
                FlowUt -= Adj(i,t,s)*f[arc(i,t,s)] ;
            }

            model.add(FlowUt == 0) ;
            FlowUt.end() ;
        }
    }

    cout << "contrainte de flot, ut " << endl ;

    //pour vt, t \in [1,T] : flot entrant = flot sortant
    for (int i=0; i<n; i++) {
        for (int t=1 ; t <= T ; t++) {

            IloExpr FlowVt(env) ;

            for (int s = 0 ; s <= t-1 ; s++) {
                FlowVt += Adj(i,s,t)*f[arc(i,s,t)] ;
            }
            for (int s = t+1 ; s <= T ; s++) {
                FlowVt -= Adj(i,s,t)*f[arc(i,s,t)] ;
            }

            model.add(FlowVt == 0) ;
            FlowVt.end() ;
        }
    }
    cout << "contrainte de flot, vt " << endl ;

    //pour s: qté de flot sortant =1
    for (int i=0; i<n; i++) {
        IloExpr FlowS(env) ;

        for (int s = 1 ; s <= T+1 ; s++) {
            FlowS += Adj(i,s,0)*f[arc(i,s,0)] + Adj(i,0,s)*f[arc(i,0,s)] ;
        }

        model.add(FlowS==1) ;
        FlowS.end() ;
    }
    cout << "contrainte de flot, s " << endl ;

    //Limite de production
    for (int i=0; i<n; i++) {
        for (int t=0 ; t < T ; t++) {
            model.add(pp[i*T + t] <= (pb->getPmax(i)-pb->getP(i))*x[i*T + t]);
            model.add(pp[i*T + t] >= 0);
        }
    }


    //Demande
    for (int t=0; t < T ; t++) {
        IloExpr Prod(env) ;
        for (int i=0; i<n; i++) {
            Prod += pp[i*T + t] + pb->getP(i)*x[i*T + t];
        }
        model.add(pb->getD(t) <= Prod);
        Prod.end() ;
    }

    cout << "end model" << endl ;
    return model ;

}

