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



CplexPricingAlgoCustomer::CplexPricingAlgoCustomer(InstanceRCFLP* inst, const Parameters & par, int c) : Param(par) {
    //env=IloEnv() ;

    customer = c;

    model = IloModel(env) ;

    int J = inst->getJ() ;

    x = IloBoolVarArray(env, J) ;
    y = IloBoolVarArray(env, J) ;

    BaseObjCoefX.resize(J, 0) ;
    BaseObjCoefY.resize(J, 0) ;

    obj = IloAdd(model, IloMinimize(env, 0.0));


    // Assignment of customer to a facility
    IloExpr sum(env) ;
    for (int j=0 ; j < J ; j++) {
        sum += x[j];
    }
    model.add(sum == 1) ;
    sum.end() ;

    // Reliability constraint
    for (int j=0 ; j < J ; j++) {
        IloExpr sum(env) ;
        for (int indice=0; indice<inst->getv(); indice++){
            sum += y[inst->getV(j)[indice]] * inst->getc(inst->getV(j)[indice]);
        }
        model.add(sum >= x[j] * inst->getd(customer) * inst->getK()) ;
    }


    cplex = IloCplex(model);
    cplex.setParam(IloCplex::EpGap, 0.00001) ;

    cplex.setParam(IloCplex::Param::Threads, 1);
    //cplex.setParam(IloCplex::Param::Parallel, -1);

    //Initialisation des coefficients objectifs (primaux)
    if (Param.doubleDecompo){
        for (int j=0 ; j < J ; j++) {
            BaseObjCoefX.at(j) = 0 ;
            BaseObjCoefY.at(j) = inst->getb(j) / inst->getI() ;
        }
    }
    else {
        for (int j=0 ; j < J ; j++) {
            BaseObjCoefX.at(j) = inst->geta(j,customer) ;
            BaseObjCoefY.at(j) = inst->getb(j) / inst->getI() ;
        }
    }
}

void CplexPricingAlgoCustomer::updateObjCoefficients(InstanceRCFLP* inst, const Parameters & Param, const DualCostsCustomer & Dual, bool Farkas) {

    int I = inst->getI() ;
    int J = inst->getJ() ;

    for (int j=0 ; j<J ; j++) {

            obj.setLinearCoef(x[j], BaseObjCoefX.at(j) );
            obj.setLinearCoef(y[j], BaseObjCoefY.at(j) ) ;

    }

}


bool CplexPricingAlgoCustomer::findImprovingSolution(InstanceRCFLP* inst, const DualCostsCustomer & Dual, double& objvalue) {
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
        objvalue = cplex.getObjValue() - Dual.Sigma[customer] ;
    }
    return true ;
}

void CplexPricingAlgoCustomer::getSolution(InstanceRCFLP* inst, const DualCostsCustomer & Dual, IloNumArray xPlan, IloNumArray yPlan, bool Farkas) {

    cplex.getValues(x, xPlan) ;
    cplex.getValues(y, yPlan) ;

}