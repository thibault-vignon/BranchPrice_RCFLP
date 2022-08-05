#include "Pricer.h"
#include "Master.h"
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
ObjPricerCustomer::ObjPricerCustomer(
        SCIP*                                scip,          /**< SCIP pointer */
        const char*                         pp_name,      /**< name of pricer */
        MasterCustomer_Model*                        M,
        InstanceRCFLP*                        instance,
        const Parameters &                  param
        ):
    ObjPricerRCFLP(scip, pp_name, instance, param)
{
    Master=M;

    int I = inst->getI() ;
    AlgoCplex = vector<CplexPricingAlgoCustomer*>(I, NULL) ;
    AlgoDynProg = vector<DynProgPricingAlgoCustomer*>(I, NULL) ;

    if (!Param.DynProgCustomer) {
        for (int i=0 ; i < I ; i++) {
            AlgoCplex[i] = new CplexPricingAlgoCustomer(inst, param, i) ;
        }
    }
    else {
        for (int i=0 ; i < I ; i++) {
            AlgoDynProg.at(i) = new DynProgPricingAlgoCustomer(inst, param, i) ;
        }
    }

    cout << "ici fin constructeur" << endl ;
}


/** Destructs the pricer object. */
ObjPricerCustomer::~ObjPricerCustomer()
{
    cout<<"Destructeur du pricer"<<endl;
}

/** initialization method of variable pricer (called after problem was transformed)
 *
 *  Because SCIP transformes the original problem in preprocessing, we need to get the references to
 *  the variables and constraints in the transformed problem from the references in the original
 *  problem.
 */
SCIP_DECL_PRICERINIT(ObjPricerCustomer::scip_init)
{

    cout<<"**************PRICER INIT************ "<<endl;
    int I = inst->getI() ;
    int J = inst->getJ() ;

    //convexity constraint
    for (int i=0 ; i < I ; i++) {
        SCIPgetTransformedCons( scip, Master->convexity_cstr.at(i), &(Master->convexity_cstr.at(i)) );
    }

    //capacity constraints
    for (int j = 0 ; j < J ; j++){
        SCIPgetTransformedCons( scip, Master->capacity_cstr.at(j), &(Master->capacity_cstr.at(j)) );
    }

    //equality constraints
    for (int j = 0 ; j < J ; j++){
        for (int i=0 ; i < I-1 ; i++) {
            SCIPgetTransformedCons( scip, Master->equality_cstr.at(j*(I-1) + i), &(Master->equality_cstr.at(j*(I-1) + i)) );
        }
    }

    cout<<"**************FIN PRICER INIT************ "<<endl;

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
//SCIP_DECL_PRICERREDCOST(ObjPricerRCFLP::scip_redcost)
//{
//    SCIPdebugMsg(scip, "call scip_redcost ...\n");
//    /* set result pointer, see above */
//    *result = SCIP_SUCCESS;
//    /* call pricing routine */
//    pricingRCFLP(scip,0);
//    return SCIP_OKAY;
//}

SCIP_RETCODE ObjPricerCustomer::scip_redcost(SCIP* scip, SCIP_PRICER* pricer, SCIP_Real* lowerbound, SCIP_Bool* stopearly, SCIP_RESULT* result)
{

    cout << "call scip_redcost ..." << endl ;

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

SCIP_RETCODE ObjPricerCustomer::scip_farkas( SCIP* scip, SCIP_PRICER* pricer, SCIP_RESULT* result ){

    cout << "call scip_farkas ..." << endl ;

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



void ObjPricerCustomer::updateDualCosts(SCIP* scip, DualCostsCustomer & dual_cost, bool Farkas) {
    ///// RECUPERATION DES COUTS DUAUX

    int print = 0 ;
    int I = inst->getI() ;
    int J = inst->getJ() ;

    //cout << "solution duale :" << endl ;

    //couts duaux "convexity constraint"
    for (int i = 0 ; i < I ; i++) {
        if (!Farkas) {
            dual_cost.Sigma.at(i) = SCIPgetDualsolLinear(scip, Master->convexity_cstr.at(i));
        }
        else{
            dual_cost.Sigma.at(i) = SCIPgetDualfarkasLinear(scip, Master->convexity_cstr.at(i));
        }
        if (print)
            cout << "sigma(" << i <<") = " << dual_cost.Sigma[i] <<endl;
    }

    //couts duaux "capacity constraint"
    for (int j = 0 ; j < J ; j++){
        if (!Farkas) {
            dual_cost.Gamma.at(j) = SCIPgetDualsolLinear(scip, Master->capacity_cstr.at(j));
        }
        else{
            dual_cost.Gamma.at(j) = SCIPgetDualfarkasLinear(scip, Master->capacity_cstr.at(j));
        }
        if (print)
            cout << "gamma(" << j <<") = " << dual_cost.Gamma.at(j) <<endl;
    }

    //couts duaux "equality constraint"
    for (int j = 0 ; j < J ; j++){
        for (int i = 0 ; i < I-1 ; i++) {
            if (!Farkas) {
                dual_cost.Omega2.at(j*I + i+1) = SCIPgetDualsolLinear(scip, Master->equality_cstr.at(j*(I-1) + i));
            }
            else{
                dual_cost.Omega2.at(j*I + i+1) = SCIPgetDualfarkasLinear(scip, Master->equality_cstr.at(j*(I-1) + i));
            }
            if (print)
                cout << "omega2(" << j << "," << i+1 << ") = " << dual_cost.Omega2.at(j*I + i+1) <<endl;
        }
    }
}

void ObjPricerCustomer::pricingRCFLP( SCIP*              scip  , bool Farkas             /**< SCIP data structure */)
{
#ifdef OUTPUT_PRICER
   //cout<<"**************PRICER************ " << endl ;
    // SCIPprintBestSol(scip, NULL, FALSE);
#endif

    int print = 1 ;

    cout<<"**************PRICER************ " << endl ;

    
    iteration++;

    cout << "solution du PMR:" << endl ;
    SCIPprintSol(scip, NULL, NULL, FALSE);

    cout << "solution réalisable:" << endl ;
    SCIPprintBestSol(scip, NULL, FALSE);

    SCIPwriteLP(scip, "debug.lp");

    int I = inst->getI() ;
    int J = inst->getJ() ;

    customerVarsToAdd.clear() ;

    // Cout duaux
    DualCostsCustomer dual_cost = DualCostsCustomer(inst) ;
    updateDualCosts(scip, dual_cost, Farkas);

    IloNumArray xPlan  ;
    IloNumArray yPlan  ;
    double redcost = 0;
    double objvalue = 0;
    bool solutionFound ;

    // Recherche par client

    for (int i=0 ; i < I ; i++) {

        if (print) cout << "client "<< i << endl;

        if (!Param.DynProgCustomer) {
            (AlgoCplex[i])->updateObjCoefficients(inst, Param, dual_cost, Farkas) ;
            xPlan = IloNumArray((AlgoCplex[i])->env, J) ;
            yPlan = IloNumArray((AlgoCplex[i])->env, J) ;
            solutionFound = (AlgoCplex[i])->findImprovingSolution(inst, dual_cost, objvalue) ;
        }

        else { // résolution par programmation dynamique
            (AlgoDynProg[i])->updateObjCoefficients(inst, Param, dual_cost, Farkas) ;
            xPlan = IloNumArray((AlgoDynProg[i])->env, J) ;
            yPlan = IloNumArray((AlgoCplex[i])->env, J) ;
            solutionFound = (AlgoDynProg[i])->findImprovingSolution(inst, dual_cost, objvalue) ;
            }

        if (!solutionFound) {
            // Pricer detected an infeasibility : we should immediately stop pricing
            infeasibilityDetected = true ;
            break ;
        }

        // TODO : currentLowerBound += objvalue + dual_cost_time.Sigma[t] ;

        if (objvalue < - Param.Epsilon) {

            if (!Param.DynProgFacility) {
                (AlgoCplex[i])->getSolution(inst, dual_cost, xPlan, yPlan, Farkas);
            }
            else{
                (AlgoDynProg[i])->getSolution(inst, dual_cost, xPlan, yPlan, Farkas);
            }

            MasterCustomer_Variable* lambda = new MasterCustomer_Variable(i, xPlan, yPlan);
            if (print) cout << "Plan found for client " << i << " with reduced cost = " << objvalue << " "  << endl ;

            totalDualCost += objvalue;
            customerVarsToAdd.push_back(lambda) ;
        }
    }

    // If we detected an infeasibility in one of the pricers, we add no variables to the Master
    // This will force termination of column generation, the master problem will remain infeasible
    // And this will cause the node to be pruned

    cout << "infeasibility detected: " << infeasibilityDetected << endl;
    if (!infeasibilityDetected){

        MasterCustomer_Variable* lambdaCustomer ;

        while (!customerVarsToAdd.empty()){

            lambdaCustomer = customerVarsToAdd.front() ;

            //// CREATION D'UNE NOUVELLE VARIABLE
            Master->initMasterCustomerVariable(scip, lambdaCustomer) ;

            /* add new variable to the list of variables to price into LP (score: leave 1 here) */
            SCIP_RETCODE ajout = SCIPaddPricedVar(scip, lambdaCustomer->ptr, 1.0);
            cout << "ajout var par temps: " << ajout << endl;

            ///// ADD COEFFICIENTS TO DEMAND, POWER LIMITS and CONVEXITY CONSTRAINTS
            Master->addCoefsToConstraints(scip, lambdaCustomer) ;

            customerColumns++;

            customerVarsToAdd.pop_front() ;

            redcost = SCIPgetVarRedcost(scip, lambdaCustomer->ptr) ;

            if (redcost >= 0){
                cout << "BUG : " << endl ;
                cout << "redcost : " << redcost << endl;
            }
        }
    }

    infeasibilityDetected = false ;

#ifdef OUTPUT_PRICER
    //SCIPwriteTransProblem(scip, "RCFLP.lp", "lp", FALSE);
    cout<<"************END PRICER******************"<<endl;
#endif

}

void ObjPricerCustomer::addVarBound(SCIP_ConsData* consdata) {

        //A implémenter
}

void ObjPricerCustomer::removeVarBound(SCIP_ConsData* consdata) {

 
        //A implémenter
}
