#include "PricingAlgo.h"
#include <iostream>
#include <ctime>

using namespace std;

DynProgPricingAlgoCustomer::DynProgPricingAlgoCustomer(InstanceRCFLP* inst, const Parameters & par, int s) : Param(par) {
    //env=IloEnv() ;
}

void DynProgPricingAlgoCustomer::updateObjCoefficients(InstanceRCFLP* inst, const Parameters & Param, const DualCostsCustomer & Dual, bool Farkas) {


}


bool DynProgPricingAlgoCustomer::findImprovingSolution(InstanceRCFLP* inst, const DualCostsCustomer & Dual, double& objvalue) {

    return(false) ;

}


void DynProgPricingAlgoCustomer::getSolution(InstanceRCFLP* inst, const DualCostsCustomer & Dual, IloNumArray xPlan, IloNumArray yPlan, bool Farkas){

}