#include "PricingAlgo.h"
#include <iostream>
#include <ctime>

using namespace std;


DynProgPricingAlgo::DynProgPricingAlgo(InstanceUCP* inst, Master_Model* M, const Parameters & par, int s) : Param(par) {
    //env=IloEnv() ;
    Master=M ;
    Site = s ;

    int T = inst->getT() ;

    branchingDecisions.resize(T, 8) ;

    Bellman.resize(2*T, 0) ;
    Prec.resize(2*T, 0) ;
}

//void DynProgPricingAlgo::updateObjCoefficients(InstanceUCP* inst, const Parameters & Param, const DualCosts & Dual, bool Farkas) {

//    int T = inst->getT() ;
//    int nbS = Param.nbUnits(Site);

//    cout << "nbS = " << nbS << endl ;
//    int first = Param.firstUnit(Site) ;

//    for (int i=0 ; i<nbS ; i++) {
//        for (int t=0 ; t < T ; t++) {

//            ObjCoefX.at(i*T+t) = 0 ;
//            if (!Farkas) {
//                ObjCoefX.at(i*T + t) += BaseObjCoefX.at(i) ;
//            }
//            ObjCoefX.at(i*T + t) +=  - inst->getPmin(first+i)*Dual.Mu[t] - (inst->getPmax(first+i) - inst->getPmin(first+i))*Dual.Nu[(first+i)*T+t] ;
//        }
//    }


//    cout << "c0: " << inst->getc0(first) << endl ;

//}

bool DynProgPricingAlgo::checkTransition(int prec_time, int current_time, int prec_status, int current_status) {
    if ( branchingDecisions.at(current_time) != current_status && branchingDecisions.at(current_time) !=8 ) {
        return false;
    }
    if (current_time>0) {
        if ( branchingDecisions.at(prec_time) != prec_status && branchingDecisions.at(prec_time) != 8) {
            return false;
        }
        for (int t= prec_time+1 ; t < current_time ; t++) {
            if ( branchingDecisions.at(t) != current_status &&  branchingDecisions.at(t) != 8) {
                return false;
            }
        }
    }
    return true ;
}

bool DynProgPricingAlgo::findImprovingSolution(InstanceUCP* inst, const DualCosts & Dual, double& objvalue) {

    int T = inst->getT() ;

    // cout << "branching decisions: " ;
    // for (int t=0 ; t < T ; t++) {
    //     cout << branchingDecisions.at(t) << " " ;
    // }
    // cout << endl ;



    int i = Param.firstUnit(Site) ;

    int l = inst->getl(i);
    int L = inst->getL(i);

    //cout << "L: " << L << endl ;

    //initialisation
    for (int t=0; t < 2*T ; t++) {
        Bellman.at(t) = std::numeric_limits<double>::infinity();
    }


    if ( checkTransition( -1, 0, -1, 0) ) {
        Bellman.at(0*T+ 0) = 0 ;
    }
    if ( checkTransition( -1, 0, -1, 1) ) {
        Bellman.at(1*T+ 0) = (Dual.ObjCoefX).at(i*T + 0) ;
    }



    Prec.at(0*T+ 0) = -1 ;
    Prec.at(1*T+ 0) = -1;

    //cout << "obj coef: " << endl ;

    for (int t = 1 ; t < T ; t++) {

        //cout << (Dual.ObjCoefX).at(i*T+t) << endl;

        ///// mise à jour de V(t, up) /////
        if ( t <= L-1 ) {
            if (checkTransition(t-1,t, 1, 1)) {
                Bellman.at(1*T+t) = Bellman.at(1*T+t-1) + (Dual.ObjCoefX).at(i*T + t);
                Prec.at(1*T+t) = 1*T+t-1;
            }
        }



        else if (t==T-1) {
            double up_prec = std::numeric_limits<double>::infinity();
            if (checkTransition(t-1, t, 1, 1))  {
                up_prec= Bellman.at(1*T+t-1) + (Dual.ObjCoefX).at(i*T + t);
            }

            int from = t-1 ;


            double somme_obj = (Dual.ObjCoefX).at(i*T + t) ;


            double down_prec = std::numeric_limits<double>::infinity();
            if (checkTransition(t-1, t, 0, 1)) {
                down_prec = Bellman.at(0*T+t-1) + (Dual.ObjCoefU).at(i*T + t) + somme_obj;
            }

            double best_down_prec = down_prec ;

            for (int k=t-2 ; k >= t - L ; k--) {
                somme_obj += (Dual.ObjCoefX).at(i*T + k+1) ;
                if (checkTransition(k, t, 0, 1)) {
                    down_prec = Bellman.at(0*T+k) + (Dual.ObjCoefU).at(i*T + k+1) + somme_obj ;

                    if (down_prec < best_down_prec) {
                        best_down_prec =down_prec;
                        from=k;
                    }
                }
            }
            if (up_prec < best_down_prec) {
                Bellman.at(1*T+t) = up_prec ;
                Prec.at(1*T+t) = 1*T+t-1;
            }
            else {
                if (checkTransition(from, t, 0, 1)) {
                    Bellman.at(1*T+t) =best_down_prec ;
                    Prec.at(1*T+t) = 0*T+from;
                }
            }

        }

        else {
            double up_prec = std::numeric_limits<double>::infinity();
            if ( checkTransition(t-1, t, 1, 1) ) {
                up_prec = Bellman.at(1*T+t-1) + (Dual.ObjCoefX).at(i*T + t);
            }

            double down_prec = std::numeric_limits<double>::infinity();
            if ( checkTransition(t-L, t, 0, 1) ) {
               // cout << "transition ok from down to up, from time " << t-L << " to " << t << endl ;
                down_prec = Bellman.at(0*T+t-L) + (Dual.ObjCoefU).at(i*T + t-L+1) + (Dual.ObjCoefX).at(i*T + t);

                for (int k=t-L+1 ; k < t ; k++) {
                    down_prec += (Dual.ObjCoefX).at(i*T + k);
                }
            }

            if (up_prec < down_prec) {
                if ( checkTransition(t-1, t, 1, 1) ) {
                    Bellman.at(1*T+t) = up_prec ;
                    Prec.at(1*T+t) = 1*T+t-1;
                }
            }
            else {
                if ( checkTransition(t-L, t, 0, 1) ) {
                    Bellman.at(1*T+t) = down_prec ;
                    Prec.at(1*T+t) = 0*T+t-L;
                }
            }
        }

        ///// mise à jour de V(t, down) //////
        if ( t <= l-1 ) {
            if (checkTransition(t-1, t, 0, 0)) {
            Bellman.at(0*T+t) = Bellman.at(0*T+t-1)  ;
            Prec.at(0*T+t) = 0*T+t-1;
            }
        }


        //Si t==T-1, il y a plus de prédécesseurs (pas de min-down à satisfaire)
        else if (t==T-1) {
            double down_prec =  std::numeric_limits<double>::infinity();
            if ( checkTransition(t-1, t, 0, 0) ) {
                down_prec = Bellman.at(0*T+t-1) ;
            }

            int from = t-1 ;

            double up_prec =  std::numeric_limits<double>::infinity();
            if ( checkTransition(t-1, t, 1, 0) ) {
                up_prec = Bellman.at(1*T+t-1) ;
            }


            double best_up_prec = up_prec ;

            for (int k=t-2 ; k >= t - l ; k--) {
                if ( checkTransition(k, t, 1, 0) ) {
                    up_prec = Bellman.at(1*T+k);

                    //cout << "for k = " << k << ", up_prec: " << Bellman.at(1*T+k) << endl ;
                    if (up_prec < best_up_prec) {
                        best_up_prec =up_prec;
                        from=k;
                    }
                }
            }

//            cout << "best up: " << best_up_prec << endl ;
//            cout << "down prec: " << down_prec << endl ;
            if (best_up_prec < down_prec) {
                Bellman.at(0*T+t) = best_up_prec ;
                Prec.at(0*T+t) = 1*T+from;
            }
            else {
                if ( checkTransition(t-1, t, 0, 0) ) {
                Bellman.at(0*T+t) = down_prec ;
                Prec.at(0*T+t) = 0*T+t-1;
               }
            }

            double V_up_1 = Bellman.at(1*T+T-1) ;
            double V_down_2 = Bellman.at(0*T+T-1) ;
            //cout << "Bellman t=T-1 : " << V_up_1 << " " << V_down_2 << endl ;

        }
        else {
            double up_prec = Bellman.at(1*T+t-l)  ;
            double down_prec = Bellman.at(0*T+t-1) ;
            if (up_prec < down_prec) {
                if ( checkTransition(t-l, t, 1, 0) ) {
                    Bellman.at(0*T+t) = up_prec ;
                    Prec.at(0*T+t) = 1*T+t-l;
                }
            }
            else {
                if ( checkTransition(t-1, t, 0, 0) ) {
                Bellman.at(0*T+t) = down_prec ;
                Prec.at(0*T+t) = 0*T+t-1;
                }
            }
        }
    }

    double V_up = Bellman.at(1*T+T-1) ;
    double V_down = Bellman.at(0*T+T-1) ;

    int print=0;
    if (print) cout << "Bellman: " ;
    for (int i=0 ; i <= 1 ; i++) {
        for (int t=0 ; t  < T ; t++) {
           if (print) cout << Bellman.at((!i)*T +t) << " " ;
       }
       if (print) cout << endl;
   }

//    cout << "valeur sans sigma: " <<fmin(V_up, V_down) << endl ;

   objvalue = fmin(V_up, V_down) - Dual.Sigma[Site] ;

    if (objvalue <= std::numeric_limits<double>::max() / 2) {
        return true ;
    }
    return false ;

}

void DynProgPricingAlgo::getUpDownPlan(InstanceUCP* inst, IloNumArray UpDownPlan) {

    int T = inst->getT() ;


    for (int t=0 ; t < T ; t++) {
        UpDownPlan[t] = 1 ;
    }


    double V_up = Bellman.at(1*T+T-1) ;
    double V_down = Bellman.at(0*T+T-1) ;

//    cout << "V_up : " << V_up << endl ;
//    cout << "V_down : " << V_down << endl ;

    int current = 0;
    int current_time = T-1 ;
    int status = 0 ;

    if (V_up < V_down) {
        current = 1*T+T-1;
        UpDownPlan[T-1]=1 ;
        status = 1;
    }
    else {
        current = 0*T+T-1;
        UpDownPlan[T-1]=0 ;
        status = 0;
    }

    while (current_time > 0) {
        int prec = Prec.at(current) ;
        int new_time = prec % T ;
        //cout << "prec: " << prec << endl;
        int new_status = prec / T ; // division entière

        UpDownPlan[new_time] = new_status ;

        for (int k= new_time + 1;  k < current_time ; k++ ) {
            UpDownPlan[k] = status ;
        }
        current = prec ;
        status=new_status ;
        current_time = new_time;
    }


    if (current_time == 0) {
        UpDownPlan[0] = status;
    }

    ///CHECK

//    int first = Param.firstUnit(Site) ;
//    double  cost = 0 ;
//    for (int t = 0 ; t < T ; t++) {
//        if (UpDownPlan[t]) {
//            cost += (Dual.ObjCoefX).at(t) ;

//            if (t> 0 && !UpDownPlan[t-1] ) {
//                cost += (Dual.ObjCoefU).at(t) ;  ;
//            }
//        }
//    }

//    cout << "check cost: " << cost << endl ;

}
