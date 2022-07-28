#include "PricingAlgo.h"
#include <iostream>
#include <ctime>

using namespace std;

DualCostsCustomer::DualCostsCustomer(InstanceRCFLP* inst) {
    int I = inst->getI() ;
    int J = inst->getJ() ;

    Sigma.resize(J, 0) ;

    Omega1.resize(I*J, 0) ;
    Omega2.resize(I*J, 0) ;
}



CplexPricingAlgoTime::CplexPricingAlgoTime(InstanceUCP* inst, const Parameters & par, int t) : Param(par) {
    //env=IloEnv() ;

    time=t;

    int n = inst->getn() ;
    model = IloModel(env) ;

    x = IloBoolVarArray(env, n) ;

    obj = IloAdd(model, IloMinimize(env, 0.0));

    if (Param.powerPlanGivenByMu || !Param.doubleDecompo){
        p = IloNumVarArray(env, n, 0.0, 1000.0) ;
        //Limite de production
        for (int i=0; i<n; i++) {
            model.add(p[i] <= (inst->getPmax(i)-inst->getPmin(i))*x[i] );
            model.add(p[i] >= 0);
        }

        //Demande
        IloExpr Prod(env) ;
        for (int i=0; i<n; i++) {
            Prod += p[i] + inst->getPmin(i)*x[i];
        }
        model.add(inst->getD(time) <= Prod);
        Prod.end() ;

        cplex = IloCplex(model);
        //cplex.setParam(IloCplex::EpGap, 0.01) ;
    }
    else{
        //Demande
        IloExpr Prod(env) ;
        for (int i=0; i<n; i++) {
            Prod += inst->getPmax(i)*x[i];
        }
        model.add(inst->getD(time) <= Prod);
        Prod.end() ;

        cplex = IloCplex(model);
        cplex.setParam(IloCplex::EpGap, 0.01) ;
    }


    //Initialisation des coefficients objectifs (primaux) de x
    BaseObjCoefX.resize(n, 0) ;
    for (int i=0 ; i <n ; i++) {
        BaseObjCoefX.at(i) = inst->getcf(i) + (inst->getPmin(i))*inst->getcp(i) ;
    }
}

void CplexPricingAlgoTime::updateObjCoefficients(InstanceUCP* inst, const Parameters & Param, const DualCostsTime & Dual, bool Farkas) {

    int n = inst->getn();
    int T = inst->getT() ;

    for (int i=0 ; i<n ; i++) {

        int L= inst->getL(i);
        int l= inst->getl(i);

        if (Param.doubleDecompo) {

            double dual_coef = 0 ;
            dual_coef += Dual.Omega.at(i*T+time);
            
            if (Param.minUpDownDouble) {
                if (time>0) {
                    dual_coef += - Dual.Mu.at(i*T + time) ;
                }
                if (time< T-1) {
                    dual_coef += Dual.Mu.at(i*T + time+1) ;
                }
                if (time>=L) {
                    dual_coef += - Dual.Nu.at(i*T+ time) ;
                }

                if (time<=T-l-1) {
                    dual_coef += - Dual.Xi.at(i*T + time + l) ;
                }
            }
            
            if (!Farkas) {
                if (Param.powerPlanGivenByMu){
                    obj.setLinearCoef(p[i], inst->getcp(i)) ;
                    if (Param.PminOnLambda){
                        obj.setLinearCoef(x[i], (1 - Param.costBalancingPricer.at(i)) * BaseObjCoefX.at(i)  + dual_coef );
                    }
                    else if (Param.PmaxOnLambda){
                        obj.setLinearCoef(x[i], BaseObjCoefX.at(i) - Param.costBalancingPricer.at(i) * ( inst->getcf(i) + inst->getcp(i) * inst->getPmax(i) ) + dual_coef );
                    }
                    else{
                        obj.setLinearCoef(x[i], BaseObjCoefX.at(i) - Param.costBalancingPricer.at(i) * inst->getcf(i) + dual_coef );
                    }
                }
                else{
                    obj.setLinearCoef(x[i], (1 - Param.costBalancingPricer.at(i)) * BaseObjCoefX.at(i)  + dual_coef );
                }
            }

        }
        else {

            //// Calcul du cout réduit de x
            double dual_coef = 0 ;
            if (time>0) {
                dual_coef += - Dual.Mu.at(i*T + time) ;
            }
            if (time< T-1) {
                dual_coef += Dual.Mu.at(i*T + time+1) ;
            }
            if (time>=L) {
                dual_coef+= - Dual.Nu.at(i*T+ time) ;
            }

            if (time<=T-l-1) {
                dual_coef += - Dual.Xi.at(i*T + time + l) ;
            }

            /// Mise à jour de la fonction objectif
            if (!Farkas) {
                obj.setLinearCoef(x[i],BaseObjCoefX.at(i)  + dual_coef );
                obj.setLinearCoef(p[i], inst->getcp(i)) ;
            }
            else {
                obj.setLinearCoef(x[i], dual_coef );
                obj.setLinearCoef(p[i], 0.0) ;
            }
        }
    }

}


bool CplexPricingAlgoTime::findImprovingSolution(InstanceUCP* inst, const DualCostsTime & Dual, double& objvalue, double & temps_resolution, int exact) {
    //returns True if a solution has been found

    ofstream LogFile("LogFile.txt");
    cplex.setOut(LogFile);

    if (exact) {
        cplex.setParam(IloCplex::EpGap, Param.Epsilon) ;
    }
    else {
        cplex.setParam(IloCplex::EpGap, 0.1) ;
    }

    cplex.setParam(IloCplex::Param::Threads, 1);
    //cplex.setParam(IloCplex::Param::Parallel, -1);

    clock_t start;
    start = clock();

    if ( !cplex.solve() ) {
        env.error() << "Failed to optimize Pricer with Cplex" << endl;
        exit(1);
    }

    temps_resolution = ( clock() - start ) / (double) CLOCKS_PER_SEC;

    if (cplex.getStatus()==CPX_STAT_INFEASIBLE){
        cout<<"NO SOLUTION TO PRICER"<<endl;
        cout<<endl<<" ************************* END PRICER with CPLEX"<<endl<<endl;

        return false;
    }
    else {
        objvalue = cplex.getObjValue() - Dual.Sigma[time] ;
    }
    return true ;
}

void CplexPricingAlgoTime::getUpDownPlan(InstanceUCP* inst, const DualCostsTime & Dual, IloNumArray UpDownPlan, IloNumArray PowerPlan, double& realCost, double & totalProd, bool Farkas) {


    int n = inst->getn();
    int T = inst->getT() ;

    cplex.getValues(x, UpDownPlan) ;

    if (!Param.doubleDecompo || (Param.doubleDecompo && Param.powerPlanGivenByMu)){
        IloNumArray prod(env, n) ;
        cplex.getValues(p, prod) ;
        for (int i=0 ; i <n ; i++) {
            if (UpDownPlan[i] > 1 - Param.Epsilon) {
                PowerPlan[i] = inst->getPmin(i) + prod[i] ;
            }
            else{
                PowerPlan[i] = 0;
            }
            totalProd += PowerPlan[i];
        }
    }
}