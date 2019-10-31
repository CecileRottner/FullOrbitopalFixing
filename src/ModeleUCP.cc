#include "ModeleUCP.h"
#include "InstanceUCP.h"

#include <ilcplex/ilocplex.h>

using namespace std ;


IloModel defineModel(IloEnv env, InstanceUCP* pb, const IloBoolVarArray & x, const IloBoolVarArray & u, int uNum, int ramp) {


    IloModel model = IloModel(env);

    int n = pb->getn();
    int T = pb->getT() ;

    IloInt t ;
    IloInt i ;
    IloInt k ;



    IloNumVarArray pp(env, n*T, 0.0, 1000);


    // Objective Function: Minimize Cost
    IloExpr cost(env) ;
    for (t=0 ; t < T ; t++) {
        for (i=0; i<n; i++) {
            cost += x[i*T + t]*pb->getcf(i) + pb->getc0(i)*u[i*T + t] + (pp[i*T + t]+pb->getP(i)*x[i*T + t])*(pb->getcp(i)) ;

        }
    }

    model.add(IloMinimize(env, cost));

    cost.end() ;


    // Conditions initiales
//    for (i=0; i<n; i++) {
//        model.add(u[i*T] >= x[i*T] - pb->getInit(i) ) ;
//    }

//    for (i=0; i<n; i++) {
//        IloExpr sum(env) ;
//        for (k= 0; k < pb->getl(i) ; k++) {
//            sum += u[i*T + k] ;
//        }
//        model.add(sum <= 1 - pb->getInit(i) ) ;
//        sum.end() ;
//    }

    // Min up constraints
    for (i=0; i<n; i++) {
        for (t=pb->getL(i) -1 ; t < T ; t++) {
            IloExpr sum(env) ;
            for (k= t - pb->getL(i) + 1; k <= t ; k++) {
                sum += u[i*T + k] ;
            }
            model.add(sum <= x[i*T + t]) ;
            sum.end() ;
        }
    }


    // Min down constraints
    for (i=0; i<n; i++) {
        for (t=pb->getl(i) ; t < T ; t++) {
            IloExpr sum(env) ;
            for (k= t - pb->getl(i) + 1; k <= t ; k++) {
                sum += u[i*T + k] ;
            }
            model.add(sum <= 1 - x[i*T + t - pb->getl(i)]) ;
            sum.end() ;
        }
    }

    //Relation entre u et x
    for (i=0; i<n; i++) {
        for (t=1 ; t < T ; t++) {
            model.add(x[i*T + t] - x[i*T + t-1] <= u[i*T + t]);
        }
    }


    //Limite de production
    for (i=0; i<n; i++) {
        for (t=0 ; t < T ; t++) {
            model.add(pp[i*T + t] <= (pb->getPmax(i)-pb->getP(i))*x[i*T + t]);
            model.add(pp[i*T + t] >= 0);
        }
    }


    //Demande
    for (t=0; t < T ; t++) {
        IloExpr Prod(env) ;
        for (i=0; i<n; i++) {
            Prod += pp[i*T + t] + pb->getP(i)*x[i*T + t];
        }
        model.add(pb->getD(t) <= Prod);
        Prod.end() ;
    }

    if (uNum) {
        cout << "conversion de u" << endl ;
        model.add(IloConversion(env, u, IloNumVar::Float) ) ;
    }

    /* model.add(IloConversion(env, x, IloNumVar::Float) ) ;
    model.add(IloConversion(env, u, IloNumVar::Float) ) ;*/

    // pp.end() ;

    if (ramp==1) {
        for (i = 0 ; i <n ; i++) {
            for (t = 1 ; t < T ; t++) {
                model.add(pp[i*T + t] - pp[i*T + t-1] <= (pb->getPmax(i)-pb->getP(i))*x[i*T + t-1]/3 );
                model.add(pp[i*T + t-1] - pp[i*T + t] <= (pb->getPmax(i)-pb->getP(i))*x[i*T + t]/2 );
            }
        }
    }

    return model ;

}



void AddRSUIneq(IloModel & model, IloEnv env, InstanceUCP* pb, const IloBoolVarArray & x, const IloBoolVarArray & u, int methode) {

    // methode=-4 : inégalité pour chaque couple (i,j)
    // par défaut : seulement couple (i,i+1)

    // ajout des inégalités Ready to Start Up
    int T = pb->getT() ;

    for (int g=0 ; g < pb->nbG ; g++) {

        int first = pb->FirstG[g] ;
        int last = pb->LastG[g] ;

        for (int i=first ; i < last ; i++) {

            int l = pb->getl(i) ;
            int L = pb->getL(i) ;

            int ub_j = i+1 ;
            if (methode == -4) {
                ub_j=last ;
            }

            for (int t = l ; t < T ; t++) {

                IloExpr rhs(env) ;
                rhs += x[i*T+t] + x[i*T+t - l];

                for (int k=t -l + 1 ; k < t ; k++) {
                    rhs+= u[i*T+k] ;
                }

                for (int j=i+1 ; j <= ub_j ; j++) {
                    model.add(u[j*T + t] <= rhs) ;
                }
            }

            //if (pb->getInit(i)==0) {
                for (int j=i+1 ; j <= ub_j ; j++) {
                    for (int t = 0 ; t < L ; t++) {
                        model.add(x[j*T+t] <= x[i*T+t]) ;
                    }
                }
           // }
        }
    }
}

IloModel defineModel_y(IloEnv env, InstanceUCP* pb, const IloBoolVarArray & x, const IloBoolVarArray & u, int ramp) {


    IloModel model = defineModel(env, pb, x, u, 0,ramp) ;


    int n = pb->getn();
    int T = pb->getT() ;


    IloBoolVarArray y(env, n*T) ;
    //IloNumVarArray y(env, n*T,0,1) ;

    for (int g=0 ; g < pb->nbG ; g++) {

        int first = pb->FirstG[g] ;
        int last = pb->LastG[g] ;

        for (int i=first ; i < last ; i++) {

            model.add(x[i*T] >= x[(i+1)*T]) ;

            model.add(y[i*T] == 1 - x[i*T] + x[(i+1)*T]) ;

            for (int t = 1 ; t < T ; t++) {
                model.add(y[i*T + t] <= y[i*T + t-1]) ;
                model.add(y[i*T + t] + x[i*T+t] - x[(i+1)*T +t] <= 1) ;

                model.add(-y[i*T + t] + y[i*T + t-1] + x[(i+1)*T +t] <= 1) ;
                model.add(y[i*T + t] - y[i*T + t-1] + x[i*T+t] >= 0) ;

                //ordre lexico
                model.add(1 - 2*y[i*T + t-1] + y[i*T + t]+ x[i*T+t]  >= x[(i+1)*T + t] ) ;

            }
        }
    }

    return model ;
}

IloModel defineModel_numberOfOnes(IloEnv env, InstanceUCP* pb, const IloBoolVarArray & x, const IloBoolVarArray & u, int ramp) {

    int T = pb->getT() ;
    IloModel model  ;
    int k ;
    model = defineModel(env, pb, x, u, 0, ramp) ;

    for (int g=0 ; g < pb->nbG ; g++) {

        int first = pb->FirstG[g] ;
        int last = pb->LastG[g] ;

        if (first < last) {
            for (int i=first ; i < last ; i++) {
                IloExpr numberOfOnes(env) ;

                for (int t = 0 ; t < T ; t++) {
                    numberOfOnes += x[i*T+t] - x[(i+1)*T +t] ;
                }
                model.add(numberOfOnes >= 0) ;

            }
        }

    }
    return model ;
}

IloModel defineModel_sum(IloEnv env, InstanceUCP* pb, const IloBoolVarArray & x, const IloBoolVarArray & u, int methode, int ramp) {

    IloModel model  ;
    int k ;

    if (methode==-5) {
        model = defineModel_y(env, pb, x, u, ramp) ;
    }

    else{
        model = defineModel(env, pb, x, u, 0,ramp) ;
    }

    int T = pb->getT() ;


    //inégalités symétries
    for (int g=0 ; g < pb->nbG ; g++) {

        int first = pb->FirstG[g] ;
        int last = pb->LastG[g] ;

        for (int i=first ; i <= last ; i++) {

            int ub_j = i+1 ;
            if (methode == -4) {
                ub_j=last ;
            }

            int l = pb->getl(i) ;
            int L = pb->getL(i) ;



            if (i < last) {

                //initialisation
                for (int j=i+1 ; j <= ub_j ; j++) {
                    for (int t = 0 ; t < fmax(l,L) ; t++) {
                        model.add(x[j*T+t] <= x[i*T+t]) ;
                    }
                }


                for (int t = l ; t < T ; t++) {

                    IloExpr rhs(env) ;
                    rhs += x[i*T+t] + x[i*T+t - l];

                    for (k=t -l + 1 ; k < t ; k++) {
                        rhs+= u[i*T+k] ;
                    }

                    if (methode==5)  {
                        if (t< T-1) {
                            rhs+= 1 - x[i*T+t+1] ;
                            for (k=t + 2 ; k < fmin(T,t+L) ; k++) {
                                rhs+= x[i*T+k-1] - x[i*T+k] + u[i*T+k] ;

                            }
                        }

                        /*for (k=t + 1 ; k < fmin(T,t+L+l) ; k++) {
                                                                rhs+= x[i*T+k-1] - x[i*T+k] + u[i*T+k] ;

                                                        }*/
                    }



                    for (int j=i+1 ; j <= ub_j ; j++) {
                        model.add(u[j*T + t] <= rhs) ;
                    }
                }


            }

            if (i>first) {

                int lb_j = i-1 ;
                if (methode == -4) {
                    lb_j=first ;
                }

                for (int t = L ; t < T ; t++) {

                    IloExpr rhs_w(env) ;
                    rhs_w += 2 - x[i*T+t] - x[i*T+t - L];

                    for (k=t -L + 1 ; k < t ; k++) {
                        rhs_w += x[i*T + k-1] - x[i*T+k] + u[i*T+k] ;
                    }

                    if (methode==5)  {
                        if (t< T-1) {
                            rhs_w +=x[i*T+t+1];
                            for (k=t + 2 ; k < fmin(t+l,T) ; k++) {
                                rhs_w += u[i*T+k] ;
                            }
                        }
                    }


                    for (int j=i-1 ; j >= lb_j ; j--) {
                        model.add(x[j*T + t-1] - x[j*T+t] + u[j*T + t] <= rhs_w) ;
                    }
                }
            }
        }
    }

    return model ;
}

IloModel AggregatedModel(IloEnv env, InstanceUCP* pb) {

    IloModel model = IloModel(env);

    int n = pb->getn();
    int T = pb->getT() ;
    int nbG = pb->getnbG() ;

    IloInt t ;
    IloInt i ;
    IloInt k ;



    IloNumVarArray pp(env, nbG*T, 0.0, 10000);
    IloIntVarArray xx(env, nbG*T, 0, n);
    IloIntVarArray uu(env, nbG*T, 0, n);

    for (t=0 ; t < T ; t++) {
        for (int g=0; g<nbG; g++) {
            model.add(xx[g*T+t] <= pb->getSizeG(g)) ;
            model.add(uu[g*T+t] <= pb->getSizeG(g)) ;
        }
    }

    // Objective Function: Minimize Cost
    IloExpr cost(env) ;
    for (t=0 ; t < T ; t++) {
        for (int g=0; g<nbG; g++) {
            int i = pb->getFirstG(g) ;
            cost += xx[g*T + t]*pb->getcf(i) + pb->getc0(i)*uu[g*T + t] + (pp[g*T + t]+pb->getP(i)*xx[g*T + t])*(pb->getcp(i)) ;

        }
    }

    model.add(IloMinimize(env, cost));

    cost.end() ;


    // Conditions initiales
    for (int g=0; g<nbG; g++) {
        int i = pb->getFirstG(g) ;
        model.add(uu[g*T] >= xx[g*T] - pb->getSizeG(g)*pb->getInit(i) ) ;
    }

    for (int g=0; g<nbG; g++) {
        int i = pb->getFirstG(g) ;
        IloExpr sum(env) ;
        for (k= 0; k < pb->getl(i) ; k++) {
            sum += uu[g*T + k] ;
        }
        model.add(sum <= pb->getSizeG(g) *(1 - pb->getInit(i) ) );
        sum.end() ;
    }

    // Min up constraints
    for (int g=0; g<nbG; g++) {
        int i = pb->getFirstG(g) ;
        for (t=pb->getL(i) -1 ; t < T ; t++) {
            IloExpr sum(env) ;
            for (k= t - pb->getL(i) + 1; k <= t ; k++) {
                sum += uu[g*T + k] ;
            }
            model.add(sum <= xx[g*T + t]) ;
            sum.end() ;
        }
    }


    // Min down constraints
    for (int g=0; g<nbG; g++) {
        int i = pb->getFirstG(g) ;
        for (t=pb->getl(i) ; t < T ; t++) {
            IloExpr sum(env) ;
            for (k= t - pb->getl(i) + 1; k <= t ; k++) {
                sum += uu[g*T + k] ;
            }
            model.add(sum <= pb->getSizeG(g) - xx[g*T + t - pb->getl(i)]) ;
            sum.end() ;
        }
    }

    //Relation entre u et x
    for (int g=0; g<nbG; g++) {
        for (t=1 ; t < T ; t++) {
            model.add(xx[g*T + t] - xx[g*T + t-1] <= uu[g*T + t]);
        }
    }


    //Limite de production
    for (int g=0; g<nbG; g++) {
        int i = pb->getFirstG(g) ;
        for (t=0 ; t < T ; t++) {
            model.add(pp[g*T + t] <= (pb->getPmax(i)-pb->getP(i))*xx[g*T + t]);
            model.add(pp[g*T + t] >= 0);
        }
    }


    //Demande
    for (t=0; t < T ; t++) {
        IloExpr Prod(env) ;
        for (int g=0; g<nbG; g++) {
            int i = pb->getFirstG(g) ;
            Prod += pp[g*T + t] + pb->getP(i)*xx[g*T + t];
        }
        model.add(pb->getD(t) <= Prod);
        Prod.end() ;
    }


    return model ;

}
