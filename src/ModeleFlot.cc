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
        else {
            if (k==0) {
                return 1 ; // arc de s à n'importe quel autre sommet ut, vt, t\in [0,T]
            }
        }
    }
    else { // état initial à 1
        if (k==0) { // A(i,t,0) = 0 pour tout t \in [0, T+1]
            return 0 ;
        }
        else {
            if (t==0) {
                return 1 ;
            }
        }
    }

    if (k==T+1 && t<T+1) {
        return 1 ;
    }

    if (t==T+1 && k<T+1) {
        return 1 ;
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

int ModeleFlot::AdjG(int g, int t, int k) {
    int i = pb->getFirstG(g) ;
    return Adj(i,t,k) ;
}

double ModeleFlot::Cost(int i, int t, int k) { // Cout de l'arc (t,k) dans Gi (0 si l'arc n'existe pas)
    if (t==k) {return 0 ;}

    else if (t > k) { // cas d'un arc (vk, ut)
        if (t <= T) { // on démarre entre 1 et T (exlut T+1 qui correspond au puits)
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

double ModeleFlot::CostG(int g, int t, int k) {
    int i = pb->getFirstG(g) ;
    return Cost(i,t,k) ;
}

int ModeleFlot::arc(int i, int t, int k) { // numéro de l'arc (t,k) de Gi
    return (t + (T+2)* k + (T+2)*(T+2)* i) ;
}



IloModel ModeleFlot::defineModelFlot(IloBoolVarArray f) {

    IloModel model(env) ;

    ///Vérifications
    /*cout << "arc s -> p, unité 1 : " << Adj(0,0,T+1) << endl ;
    cout << "arc s -> p, unité 2 : " << Adj(1, 0, T+1) << endl ;
    cout << "arc v24 p, unité 3 : " << Adj(2,25,24) << endl ;
    cout << "cost: " << Cost(2,25,24) << endl ;*/


    //IloBoolVarArray f(env, n*(T+2)*(T+2)) ;
    IloNumVarArray pp(env, n*T, 0.0, 1000);

    //Vecteur x
    IloExprArray x(env, n*T) ;
    for (int i = 0 ; i < n ; i++) {
        for (int t=0 ; t < T ; t++) {

            x[i*T+t] = IloExpr(env) ;
            // calcul de xit : unité i en marche à t+1
            for (int s = 0 ; s <= t+1 ; s++) {
                for (int z = t+2 ; z <= T+1 ; z++) {
                    x[i*T+t] +=  Adj(i,s,z)*f[arc(i,s,z)]   ;
                }
            }

        }
    }

    for (int i=0; i<n; i++) {
        for (int t=0 ; t < T+2 ; t++) {
            for (int k = 0 ; k < T+2; k++) {
                if (Adj(i,t,k) == 0) {
                    model.add(f[arc(i,t,k)] == 0) ;
                }
            }
        }
    }

    // Objective Function: Minimize Cost
    IloExpr cost(env) ;
    for (int i=0; i<n; i++) {
        for (int t=0 ; t < T ; t++) {
            cost += (pp[i*T + t]+ pb->getP(i)*x[i*T+t] )*(pb->getcp(i)) ;
        }
        for (int t=0 ; t < T+2 ; t++) {
            for (int k = 0 ; k < T+2; k++) {
                cost += Cost(i,t,k)*f[arc(i,t,k)] ;
            }
        }
    }
    model.add(IloMinimize(env, cost));

    //pour ut, t \in [1,T] : flot entrant = flot sortant
    for (int i=0; i<n; i++) {
        for (int t=1 ; t <= T ; t++) {

            IloExpr FlowUt(env) ;

            for (int s = 0 ; s <= t-1 ; s++) {
                FlowUt += Adj(i,t,s)*f[arc(i,t,s)] ;
            }
            for (int s = t+1 ; s <= T+1 ; s++) {
                FlowUt -= Adj(i,t,s)*f[arc(i,t,s)] ;
            }

            model.add(FlowUt == 0) ;
            FlowUt.end() ;
        }
    }

    //pour vt, t \in [1,T] : flot entrant = flot sortant
    for (int i=0; i<n; i++) {
        for (int t=1 ; t <= T ; t++) {

            IloExpr FlowVt(env) ;

            for (int s = 0 ; s <= t-1 ; s++) {
                FlowVt += Adj(i,s,t)*f[arc(i,s,t)] ;
            }
            for (int s = t+1 ; s <= T+1 ; s++) {
                FlowVt -= Adj(i,s,t)*f[arc(i,s,t)] ;
            }

            model.add(FlowVt == 0) ;
            FlowVt.end() ;
        }
    }

    //pour s: qté de flot sortant =1
    for (int i=0; i<n; i++) {
        IloExpr FlowS(env) ;

        for (int s = 1 ; s <= T+1 ; s++) {
            FlowS += Adj(i,s,0)*f[arc(i,s,0)] + Adj(i,0,s)*f[arc(i,0,s)] ;
        }

        model.add(FlowS==1) ;
        FlowS.end() ;
    }

    //Limite de production
    for (int i=0; i<n; i++) {
        for (int t=0 ; t < T ; t++) {
            model.add(pp[i*T + t] <= (pb->getPmax(i)-pb->getP(i))*x[i*T+t] );
            model.add(pp[i*T + t] >= 0);
        }
    }

    //Demande
    for (int t=0; t < T ; t++) {
        IloExpr Prod(env) ;
        for (int i=0; i<n; i++) {
            Prod += pp[i*T + t] + pb->getP(i)*x[i*T+t];
        }
        model.add(pb->getD(t) <= Prod);
        Prod.end() ;
    }

    return model ;

}

IloModel ModeleFlot::AggregatedFlowModel() {

    IloModel model = IloModel(env);

    int nbG = pb->getnbG() ;



    IloNumVarArray pp(env, nbG*T, 0.0, 10000);
    IloIntVarArray f(env, nbG*(T+2)*(T+2), 0, n);


    for (int t=0 ; t < T+2 ; t++) {
        for (int k=0 ; t < T+2 ; t++) {
            for (int g=0; g<nbG; g++) {
                model.add(f[arc(g,t,k)] <= pb->getSizeG(g)) ;
            }
        }
    }

    //Vecteur x
    IloExprArray x(env, nbG*T) ;
    for (int g=0; g<nbG; g++) {
        for (int t=0 ; t < T ; t++) {

            x[g*T+t] = IloExpr(env) ;
            // calcul de xit : unité i en marche à t+1
            for (int s = 0 ; s <= t+1 ; s++) {
                for (int z = t+2 ; z <= T+1 ; z++) {
                    x[g*T+t] +=  AdjG(g,s,z)*f[arc(g,s,z)]   ;
                }
            }

        }
    }

    for (int g=0; g<nbG; g++) {
        for (int t=0 ; t < T+2 ; t++) {
            for (int k = 0 ; k < T+2; k++) {
                if (AdjG(g,t,k) == 0) {
                    model.add(f[arc(g,t,k)] == 0) ;
                }
            }
        }
    }

    // Objective Function: Minimize Cost
    IloExpr cost(env) ;
    for (int g=0; g<nbG; g++) {
        int i = pb->getFirstG(g) ;

        for (int t=0 ; t < T ; t++) {
            cost += (pp[g*T + t]+ pb->getP(i)*x[g*T+t] )*(pb->getcp(i)) ;
        }
        for (int t=0 ; t < T+2 ; t++) {
            for (int k = 0 ; k < T+2; k++) {
                cost += CostG(g,t,k)*f[arc(g,t,k)] ;
            }
        }
    }
    model.add(IloMinimize(env, cost));

    //pour ut, t \in [1,T] : flot entrant = flot sortant
    for (int g=0; g<nbG; g++) {
        for (int t=1 ; t <= T ; t++) {

            IloExpr FlowUt(env) ;

            for (int s = 0 ; s <= t-1 ; s++) {
                FlowUt += AdjG(g,t,s)*f[arc(g,t,s)] ;
            }
            for (int s = t+1 ; s <= T+1 ; s++) {
                FlowUt -= AdjG(g,t,s)*f[arc(g,t,s)] ;
            }

            model.add(FlowUt == 0) ;
            FlowUt.end() ;
        }
    }

    //pour vt, t \in [1,T] : flot entrant = flot sortant
    for (int g=0; g<nbG; g++) {
        for (int t=1 ; t <= T ; t++) {

            IloExpr FlowVt(env) ;

            for (int s = 0 ; s <= t-1 ; s++) {
                FlowVt += AdjG(g,s,t)*f[arc(g,s,t)] ;
            }
            for (int s = t+1 ; s <= T+1 ; s++) {
                FlowVt -= AdjG(g,s,t)*f[arc(g,s,t)] ;
            }

            model.add(FlowVt == 0) ;
            FlowVt.end() ;
        }
    }

    //pour s: qté de flot sortant =1
    for (int g=0; g<nbG; g++) {
        IloExpr FlowS(env) ;

        for (int s = 1 ; s <= T+1 ; s++) {
            FlowS += AdjG(g,s,0)*f[arc(g,s,0)] + AdjG(g,0,s)*f[arc(g,0,s)] ;
        }

        model.add(FlowS==pb->getSizeG(g)) ;
        FlowS.end() ;
    }

    //Limite de production
    for (int g=0; g<nbG; g++) {

        int i = pb->getFirstG(g) ;

        for (int t=0 ; t < T ; t++) {
            model.add(pp[g*T + t] <= (pb->getPmax(i)-pb->getP(i))*x[g*T+t] );
            model.add(pp[g*T + t] >= 0);
        }
    }

    //Demande
    for (int t=0 ; t < T ; t++) {
        IloExpr Prod(env) ;

        for (int g=0; g<nbG; g++) {
            int i = pb->getFirstG(g) ;
            Prod += pp[g*T + t] + pb->getP(i)*x[g*T+t];
        }
        model.add(pb->getD(t) <= Prod);
        Prod.end() ;
    }

    return model ;

}



IntervalModel::IntervalModel(IloEnv enviro, InstanceUCP* pbm, Methode & m) {
    env=enviro ;
    pb=pbm ;
    met=m ;
    n = pb->getn();
    T = pb->getT() ;

    Lmin = pb->getL(0) ;
    for (int i=1 ; i < n ; i++) {
        if (pb->getL(i) < Lmin) {
            Lmin = pb->getL(i) ;
        }
    }
}

double IntervalModel::FixedCost(int g, int a, int b) {

    //    cout << "Fixed cost: intervalle " << a <<", " << b+Lmin-1 << endl ;
    //    cout << "g="<< g << endl ;


    int first = pb->getFirstG(g) ;
    //    cout << "c0: " <<pb->getc0(first) << endl;
    //    cout << "cf: " << pb->getcf(first) << endl;
    double cost = 0 ;
    if (a>0) {
        cost += pb->getc0(first) ;
    }
    cost += (b+Lmin-a)*pb->getcf(first) ;
    return cost ;
}

int IntervalModel::Pindex(int g, int a, int b, int t) {
    int index = g*(T)*(T - Lmin+1)*T + a*(T-Lmin+1)*T + b*T + t ;
    return index;
}

int IntervalModel::Yindex(int g, int a, int b) {
    int index = g*(T)*(T - Lmin+1) + a*(T-Lmin+1) + b ;
    return index;
}

int IntervalModel::inUpInterval(int a, int b, int t) { // renvoie 1 si t appartient à l'intervalle  [a, b+Lmin-1]
    if (t < a) {
        return 0 ;
    }
    if (t > b+Lmin-1) {
        return 0 ;
    }
    return 1 ;
}

int IntervalModel::inCliqueInterval(int a, int b, int t, int i) { // renvoie 1 si t appartient à l'intervalle  [a, b+Lmin-1]
    int l = pb->getl(i) ;
    if (t < a) {
        return 0 ;
    }
    if (t >= b+Lmin+l) {
        return 0 ;
    }
    return 1 ;
}

IloModel IntervalModel::defineIntervalModel(IloIntVarArray Y) {
    cout << "in the model" << endl;

    double Pmaxmax = pb->getPmax(0) ;
    for (int i=1 ; i < n ; i++) {
        if (pb->getPmax(i) > Pmaxmax) {
            Pmaxmax = pb->getPmax(i) ;
        }
    }

    int nbG = pb->getnbG() ;
    //  Y = IloIntVarArray(env, nbG*(T)*(T-Lmin+1), 0, n) ;
    IloNumVarArray p(env, nbG*(T)*(T-Lmin+1)*T, 0.0, n*Pmaxmax);



    IloModel model = IloModel(env) ;



    //Borne sup sur Y et sur P pour chaque groupe. Borne à 0 lorsque l'intervalle est trop petit pour le temps min de marche, où lorsque t est en dehors (pour p).
    for (int g = 0 ; g < nbG ; g++) {

        int first = pb->getFirstG(g) ;
        int last = pb->getLastG(g) ;
        int nb = last-first+1 ;

        int L = pb->getL(first) ;

        for (int a = 0 ; a < T ; a++) {
            for (int b=0 ; b <= fmin(a + L - Lmin - 1,T-Lmin-1) ; b++) {
                //                if (b>=a) {
                //                    cout << "intervalle à zero: [" << a << ", " << b+Lmin-1 << "], L= " << pb->getL(first) << endl;
                //                }
                model.add(Y[Yindex(g,a,b)] == 0) ;

                for (int t = 0 ; t < T ; t++) {
                    model.add(p[Pindex(g,a,b,t)] == 0) ;
                }
            }

            for (int b = fmin(a + L - Lmin, T-Lmin) ; b <= T - Lmin ; b++) {

                // Borne Y
                model.add(Y[Yindex(g,a,b)] <= nb) ;

                //Borne P
                for (int t = 0 ; t < T ; t++) {
                    int inside = inUpInterval(a,b,t) ;
                    if (!inside) {
                        model.add(p[Pindex(g,a,b,t)] == 0) ;
                    }
                    else {
                        model.add(p[Pindex(g,a,b,t)] <= (pb->getPmax(first) - pb->getP(first)) *Y[Yindex(g,a,b)]) ;
                    }
                }
            }
        }
    }

    // Objective Function: Minimize Cost
    IloExpr cost(env) ;

    for (int g=0; g<nbG; g++) {

        int i = pb->getFirstG(g) ;

        for (int a = 0 ; a < T ; a++) {
            for (int b =0 ; b < T - Lmin + 1 ; b++) {

                cost += FixedCost(g,a,b)*Y[Yindex(g,a,b)] ;

                for (int t=0 ; t < T ; t++) {
                    int inside= inUpInterval(a,b,t) ;
                    if (inside) {
                        cost += (p[Pindex(g,a,b,t)] + pb->getP(i)*Y[Yindex(g,a,b)] )*(pb->getcp(i)) ;
                    }
                }
                //cout << "production cost of g="<< g << " for [" << a << ", " << b+Lmin-1 <<"] : " << pb->getP(i)*pb->getcp(i) << endl ;
            }
        }
    }

    model.add(IloMinimize(env, cost));

    // Demand constraints
    for (int t = 0 ; t < T ; t++) {
        IloExpr Prod(env) ;

        for (int g=0; g<nbG; g++) {

            int i = pb->getFirstG(g) ;

            for (int a = 0 ; a < T ; a++) {
                for (int b =0 ; b < T - Lmin + 1 ; b++) {
                    int inside = inUpInterval(a,b,t) ;
                    if (inside) {
                        Prod += p[Pindex(g,a,b,t)] + pb->getP(i)*Y[Yindex(g,a,b)];
                    }
                }
            }
        }

        model.add(pb->getD(t) <= Prod);
        Prod.end() ;
    }


    //Ramp constraints

    if (met.Ramping()) {
        for (int g=0; g<nbG; g++) {
            int i = pb->getFirstG(g) ;
            int SU = 0 ;
            double RU = (pb->getPmax(i)-pb->getP(i))/3;
            double RD = (pb->getPmax(i)-pb->getP(i))/2 ;

            for (int a = 0 ; a < T ; a++) {
                for (int b = 0 ; b < T - Lmin + 1 ; b++) {

                    if (a > 0) {
                        model.add(p[Pindex(g,a,b,a)] <= SU*Y[Yindex(g,a,b)] ) ;
                    }
                    if (b+Lmin-1 < T-1) {
                        model.add(p[Pindex(g,a,b,b+Lmin-1)] <= SU*Y[Yindex(g,a,b)] ) ;
                    }

                    for (int t=a+1 ; t < b+Lmin ; t++) {
                        model.add( p[Pindex(g,a,b,t)] - p[Pindex(g,a,b,t-1)] <= RU*Y[Yindex(g,a,b)] ) ;
                        model.add( p[Pindex(g,a,b,t-1)] - p[Pindex(g,a,b,t)] <= RD*Y[Yindex(g,a,b)] ) ;
                    }
                }
            }
        }
    }

    // Contrainte de clique
    for (int t = 0 ; t < T ; t++) {
        for (int g=0; g<nbG; g++) {

            IloExpr SumY(env) ;
            int i = pb->getFirstG(g) ;
            int j = pb->getLastG(g) ;
            int nb = j-i+1 ;
            int L = pb->getL(i) ;

            // cout << "Clique t=" << t << ", "<< "l= " << pb->getl(i) << " " ;

            for (int a = 0 ; a < T ; a++) {
                for (int b =fmin(a+L-Lmin, T-Lmin) ; b < T - Lmin + 1 ; b++) {
                    int inClique = inCliqueInterval(a,b,t,i) ;
                    if (inClique) {
                        // cout << " [" << a << ", " << b+Lmin-1 << "]"<< endl;
                        SumY += Y[Yindex(g,a,b)] ;
                    }
                }
            }

            model.add(SumY <= nb) ;
            SumY.end() ;
        }
    }

    //    IloCplex cplex(model) ;
    //    cplex.solve();
    //    IloNumArray ysol(env,0) ;
    //    IloNumArray psol(env, 0) ;
    //    cplex.getValues(Y, ysol) ;
    //    cplex.getValues(p, psol) ;

    //    for (int g = 0 ; g < nbG ; g++) {
    //        for (int a = 0 ; a < T ; a++) {
    //            for (int b =0 ; b < T - Lmin + 1 ; b++) {
    //                if ( ysol[Yindex(g,a,b)] > 0.99999) {
    //                    cout << "Groupe " << g << ", intervalle [" << a+1 << ", " << b+Lmin << "] ; y = " << ysol[Yindex(g,a,b)] << endl ;
    //                    cout << "Power: " ;
    //                    for (int t=a ; t < b+ Lmin ; t++) {
    //                        cout << psol[Pindex(g,a,b,t)] << " ";
    //                    }
    //                    cout << endl ;
    //                    cout << "Fixed cost: " << FixedCost(g,a,b) << endl ;
    //                    cout << endl ;
    //                }
    //            }
    //        }
    //    }


    //    cout << "obj value: " << cplex.getObjValue() << endl ;


    cout << "out the model" << endl;

    return model;
}

