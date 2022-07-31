#include "PricingAlgo.h"
#include <iostream>
#include <ctime>

using namespace std;


DynProgPricingAlgoFacility::DynProgPricingAlgoFacility(InstanceRCFLP* inst, const Parameters & par, int s) : Param(par) {
    //env=IloEnv() ;
}

void DynProgPricingAlgoFacility::updateObjCoefficients(InstanceRCFLP* inst, const Parameters & Param, const DualCostsFacility & Dual, bool Farkas) {


}


bool DynProgPricingAlgoFacility::findImprovingSolution(InstanceRCFLP* inst, const DualCostsFacility & Dual, double& objvalue) {

    return(false) ;

}


void DynProgPricingAlgoFacility::getSolution(InstanceRCFLP* inst, const DualCostsFacility & Dual, IloNumArray xPlan, IloNumArray yPlan, bool Farkas){

}