#include "PricingAlgo.h"
#include <iostream>

using namespace std;

DualCostsFacility::DualCostsFacility(InstanceRCFLP* inst) {
    int I = inst->getI() ;
    int J = inst->getJ() ;

    Sigma.resize(J, 0) ;

    Omega1.resize(I*J, 0) ;
    Omega2.resize(I*J, 0) ;
}


CplexPricingAlgoFacility::CplexPricingAlgoFacility(InstanceRCFLP* inst, const Parameters & par, int f) : Param(par) {

    //env = IloEnv() ;
    facility = f;

    int I = inst->getI() ;

    x = IloBoolVarArray(env, I) ;
    y = IloBoolVarArray(env, 1) ;

    BaseObjCoefX.resize(I, 0) ;
    BaseObjCoefY.resize(1, 0) ;

    cpuTime=0 ;

    model = IloModel(env) ;
    obj = IloAdd(model, IloMinimize(env, 0.0));


    // Capacity constraint

    IloExpr sum(env) ;
    for (int i=0; i<I; i++) {
        sum += x[i] * inst->getd(i);
    }
    if (Param.compactCapacityConstraints){
        model.add(sum <= y[0] * inst->getc(facility)) ;
    }
    else{
        model.add(sum <= inst->getc(facility)) ;
    }
    sum.end() ;
    

    cplex = IloCplex(model);
    cplex.setParam(IloCplex::EpGap, 0.00001) ;

    cplex.setParam(IloCplex::Param::Threads, 1);
    //cplex.setParam(IloCplex::Param::Parallel, -1);

    //Initialisation des coefficients objectifs (primaux)
    if (Param.doubleDecompo){
        for (int i=0 ; i < I ; i++) {
            BaseObjCoefX.at(i) =  inst->geta(facility,i);
        }
        BaseObjCoefY.at(0) = 0;
    }
    else {
        for (int i=0 ; i < I ; i++) {
            BaseObjCoefX.at(i) =  inst->geta(facility,i);
        }
        BaseObjCoefY.at(0) = inst->getb(facility);
    }
}

void CplexPricingAlgoFacility::updateObjCoefficients(InstanceRCFLP* inst, const Parameters & Param, const DualCostsFacility & Dual, bool Farkas) {

    int I = inst->getI() ;

    for (int i=0 ; i<I ; i++) {

            obj.setLinearCoef(x[i], BaseObjCoefX.at(i));

    }
    
    obj.setLinearCoef(y[0], BaseObjCoefY.at(0)) ;
}



bool CplexPricingAlgoFacility::findImprovingSolution(InstanceRCFLP* inst, const DualCostsFacility & Dual, double& objvalue) {
    //returns True if a solution has been found

    ofstream LogFile("LogFile.txt");
    cplex.setOut(LogFile);

    if ( !cplex.solve() ) {
        env.error() << "Failed to optimize Pricer with Cplex" << endl;
        exit(1);
    }

    if (cplex.getStatus()==CPX_STAT_INFEASIBLE){
        cout<<"NO SOLUTION TO PRICER"<<endl;
        cout<<endl<<" ************************* END PRICER with CPLEX"<<endl<<endl;

        return false;
    }
    else {
        /*cout << "for facility " << facility << "; " << endl ;
       cout << "obj value without sigma: " << cplex.getObjValue() << endl;*/
        objvalue = cplex.getObjValue() - Dual.Sigma[facility] ;
    }

    return true;
}


void CplexPricingAlgoFacility::getSolution(InstanceRCFLP* inst, const DualCostsFacility & Dual, IloNumArray xPlan, IloNumArray yPlan, bool Farkas) {

    cplex.getValues(x, xPlan) ;

    if (Param.compactCapacityConstraints){
        cplex.getValues(y, yPlan) ;
    }

}