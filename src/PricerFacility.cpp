#include "Pricer.h"
#include "scip/cons_linear.h"
#include <map>
#include <vector>
#include <iostream>


#define OUTPUT_PRICER
// à décommenter pour l'affichage de debug

using namespace std;
using namespace scip;


/** Constructs the pricer object with the data needed
 *
 *  An alternative is to have a problem data class which allows to access the data.
 */
ObjPricerFacility::ObjPricerFacility(
        SCIP*                                scip,          /**< SCIP pointer */
        const char*                         pp_name,      /**< name of pricer */
        MasterFacility_Model*                        M,
        InstanceRCFLP*                        instance,
        const Parameters &                  param
        ):
    ObjPricerRCFLP(scip, pp_name, instance, param)
{
    int J = inst->getJ();
    Master=M ;
    AlgoCplex = vector<CplexPricingAlgoFacility*>(J, NULL) ;
    AlgoDynProg = vector<DynProgPricingAlgoFacility*>(J, NULL) ;

    if (!Param.DynProgFacility) {
        for (int j = 0 ; j < J ; j++) {
            AlgoCplex[j] = new CplexPricingAlgoFacility(inst, param, j) ;
        }
    }
    else {
        for (int j = 0 ; j < J ; j++) {
            AlgoDynProg.at(j) = new DynProgPricingAlgoFacility(inst, param, j) ;
        }
    }
}


/** Destructs the pricer object. */
ObjPricerFacility::~ObjPricerFacility()
{
    cout<<"Destructeur du pricer"<<endl;
}

/** initialization method of variable pricer (called after problem was transformed)
 *
 *  Because SCIP transformes the original problem in preprocessing, we need to get the references to
 *  the variables and constraints in the transformed problem from the references in the original
 *  problem.
 */
SCIP_DECL_PRICERINIT(ObjPricerFacility::scip_init)
{
    //cout<<"**************PRICER INIT************ "<<endl;

    int J = inst->getJ() ;

    //convexity constraints
    for (int j = 0 ; j < J ; j++) {
        SCIPgetTransformedCons(scip, Master->convexity_cstr[j], &(Master->convexity_cstr[j]));
    }

    // TODO : other constraints

    cout<<"**************FIN PRICER INIT************ "<<endl ;

    return SCIP_OKAY;
}


/** Pricing of additional variables if LP is feasible.
 *
 *  - get the values of the dual variables you need
 *  - construct the reduced units costs from dual costs
 *  - find the cheapest production plan for each site
 *  - if this plan has negative reduced cost, add it to the LP
 *
 *  possible return values for *result:
 *  - SCIP_SUCCESS    : at least one improving variable was found, or it is ensured that no such variable exists
 *  - SCIP_DIDNOTRUN  : the pricing process was aborted by the pricer, there is no guarantee that the current LP solution is optimal
 */
//SCIP_DECL_PRICERREDCOST(ObjPricerFacility::scip_redcost)
//{
//    SCIPdebugMsg(scip, "call scip_redcost ...\n");
//    /* set result pointer, see above */
//    *result = SCIP_SUCCESS;
//    /* call pricing routine */
//    pricingRCFLP(scip,0);
//    return SCIP_OKAY;
//}

SCIP_RETCODE ObjPricerFacility::scip_redcost(SCIP* scip, SCIP_PRICER* pricer, SCIP_Real* lowerbound, SCIP_Bool* stopearly, SCIP_RESULT* result)
{

    SCIPdebugMsg(scip, "call scip_redcost ...\n");

    if( Param.PriceAndBranch && SCIPgetDepth(scip) != 0 )
    {
        *result = SCIP_SUCCESS;
        return SCIP_OKAY;
    }

    /* set result pointer, see above */
    *result = SCIP_SUCCESS;

    /* call pricing routine */
    pricingRCFLP(scip,0);

    return SCIP_OKAY;

}

SCIP_RETCODE ObjPricerFacility::scip_farkas( SCIP* scip, SCIP_PRICER* pricer, SCIP_RESULT* result ){

    SCIPdebugMsg(scip, "call scip_farkas ...\n");

    if( Param.PriceAndBranch && SCIPgetDepth(scip) != 0 )
    {
        *result = SCIP_SUCCESS;
        return SCIP_OKAY;
    }

    /* set result pointer, see above */
    *result = SCIP_SUCCESS;

    /* call pricing routine */
    pricingRCFLP(scip,1);

    return SCIP_OKAY;

}



void ObjPricerFacility::updateDualCosts(SCIP* scip, DualCostsFacility & dual_cost, bool Farkas) {
    ///// RECUPERATION DES COUTS DUAUX

    int print = 0 ;
    int J = inst->getJ() ;

    //couts duaux contrainte convexité
    for (int j = 0 ; j < J ; j++) {
        if (!Farkas) {
            dual_cost.Sigma[j] = SCIPgetDualsolLinear(scip, Master->convexity_cstr[j]);
        }
        else{
            dual_cost.Sigma[j] = SCIPgetDualfarkasLinear(scip, Master->convexity_cstr[j]);
        }
        if (print) cout << "sigma: " << dual_cost.Sigma[j] <<endl;
    }

    // TODO : autres contraintes

}

void ObjPricerFacility::pricingRCFLP( SCIP*              scip  , bool Farkas             /**< SCIP data structure */)
{
#ifdef OUTPUT_PRICER
    cout<<"**************PRICER************ "<< endl ;
    // SCIPprintBestSol(scip, NULL, FALSE);
#endif

    int print = 1 ;
    iteration++;

//    /// PMR courant et sa solution
   // SCIPwriteTransProblem(scip, NULL, NULL, FALSE);

//   // cout << "solution du PMR:" << endl ;
//    SCIPprintSol(scip, NULL, NULL, FALSE);

//    //cout << "solution réalisable:" << endl ;
//    SCIPprintBestSol(scip, NULL, FALSE);

    int I = inst->getI() ;
    int J = inst->getJ() ;

    facilityVarsToAdd.clear() ;

    // Cout duaux
    DualCostsFacility dual_cost = DualCostsFacility(inst) ;
    updateDualCosts(scip, dual_cost, Farkas);
    
    IloNumArray xPlan  ;
    IloNumArray yPlan  ;
    double redcost = 0;
    double objvalue = 0;
    bool solutionFound ;

    // Recherche par facility

    for (int j = 0 ; j < J ; j++) {

        if (print) cout << "facility "<< j << endl;

        if (!Param.DynProgFacility) {
            (AlgoCplex[j])->updateObjCoefficients(inst, Param, dual_cost, Farkas) ;
            xPlan = IloNumArray((AlgoCplex[j])->env, I) ;
            if (Param.compactCapacityConstraints){
                yPlan = IloNumArray((AlgoCplex[j])->env, 1) ;
            }
            solutionFound = (AlgoCplex[j])->findImprovingSolution(inst, dual_cost, objvalue) ;
        }

        else { // résolution par programmation dynamique
            (AlgoDynProg[j])->updateObjCoefficients(inst, Param, dual_cost, Farkas) ;
            xPlan = IloNumArray((AlgoDynProg[j])->env, I) ;
            if (Param.compactCapacityConstraints){
                yPlan = IloNumArray((AlgoDynProg[j])->env, 1) ;
            }
            solutionFound = (AlgoDynProg[j])->findImprovingSolution(inst, dual_cost, objvalue) ;
        }
        
        if (!solutionFound) {
            // Pricer detected an infeasibility : we should immediately stop pricing
            infeasibilityDetected = true ;
            break ;
        }

        // TODO : currentLowerBound += objvalue + dual_cost.Sigma[s] ;

        if (objvalue < - Param.Epsilon) {

            if (!Param.DynProgFacility) {
                (AlgoCplex[j])->getSolution(inst, dual_cost, xPlan, yPlan, Farkas);
            }
            else{
                (AlgoDynProg[j])->getSolution(inst, dual_cost, xPlan, yPlan, Farkas);
            }

            MasterFacility_Variable* lambda = new MasterFacility_Variable(j, xPlan, yPlan);
            if (print) cout << "Plan found for facility " << j << " with reduced cost = " << objvalue << " "  << endl ;

            totalDualCost += objvalue;
            facilityVarsToAdd.push_back(lambda) ;
        }
    }


    // If we detected an infeasibility in one of the pricers, we add no variables to the Master
    // This will force termination of column generation, the master problem will remain infeasible
    // And this will cause the node to be pruned

    if (!infeasibilityDetected){

        MasterFacility_Variable* lambdaFacility ;
        
        while (!facilityVarsToAdd.empty()){

            lambdaFacility = facilityVarsToAdd.front() ;

            //// CREATION D'UNE NOUVELLE VARIABLE DANS LE MASTER
            Master->initMasterFacilityVariable(scip, lambdaFacility) ;

            /* add new variable to the list of variables to price into LP (score: leave 1 here) */
            SCIP_RETCODE ajout = SCIPaddPricedVar(scip, lambdaFacility->ptr, 1.0);
            cout << "ajout var par facility: " << ajout << endl;

            ///// ADD COEFFICIENTS TO CONVEXITY and TIME/SITE EQUALITY CONSTRAINTS
            Master->addCoefsToConstraints(scip, lambdaFacility) ;

            facilityColumns++;

            facilityVarsToAdd.pop_front() ;
        }
    }

    infeasibilityDetected = false ;

#ifdef OUTPUT_PRICER
    SCIPwriteTransProblem(scip, "RCFLP.lp", "lp", FALSE);
    cout<<"************END PRICER******************"<<endl;
#endif

}


void ObjPricerFacility::addVarBound(SCIP_ConsData* consdata) {

}

void ObjPricerFacility::removeVarBound(SCIP_ConsData* consdata) {


}


