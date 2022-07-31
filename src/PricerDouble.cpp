#include "Pricer.h"
#include "scip/cons_linear.h"
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>


#define OUTPUT_PRICER
// à décommenter pour l'affichage de debug

using namespace std;
using namespace scip;




/** Constructs the pricer object with the data needed
 *
 *  An alternative is to have a problem data class which allows to access the data.
 */
ObjPricerDouble::ObjPricerDouble(
        SCIP*                                scip,          /**< SCIP pointer */
        const char*                         pp_name,      /**< name of pricer */
        MasterDouble_Model*                        M,
        InstanceRCFLP*                        instance,
        const Parameters &                  param
        ):
    ObjPricerRCFLP(scip, pp_name, instance, param)
{
    Master=M ;

    int J = inst->getJ();
    int I = inst->getI();

    AlgoCplex_facility = vector<CplexPricingAlgoFacility*>(J, NULL) ;
    AlgoDynProg_facility = vector<DynProgPricingAlgoFacility*>(J, NULL) ;

    if (!Param.DynProgFacility) {
        for (int j = 0 ; j < J ; j++) {
            AlgoCplex_facility.at(j) = new CplexPricingAlgoFacility(inst, param, j) ;
        }
    }
    else {
        for (int j = 0 ; j < J ; j++) {
            AlgoDynProg_facility.at(j) = new DynProgPricingAlgoFacility(inst, param, j) ;
        }
    }


    AlgoCplex_customer = vector<CplexPricingAlgoCustomer*>(I, NULL) ;
    AlgoDynProg_customer = vector<DynProgPricingAlgoCustomer*>(I, NULL) ;

    if (!Param.DynProgCustomer) {
        for (int i=0 ; i < I ; i++) {
            AlgoCplex_customer[i] = new CplexPricingAlgoCustomer(inst, param, i) ;
        }
    }
    else {
        for (int i=0 ; i < I ; i++) {
            AlgoDynProg_customer.at(i) = new DynProgPricingAlgoCustomer(inst, param, i) ;
        }
    }
}


/** Destructs the pricer object. */
ObjPricerDouble::~ObjPricerDouble()
{
    cout<<"Destructeur du pricer"<<endl;
}

/** initialization method of variable pricer (called after problem was transformed)
 *
 *  Because SCIP transformes the original problem in preprocessing, we need to get the references to
 *  the variables and constraints in the transformed problem from the references in the original
 *  problem.
 */
SCIP_DECL_PRICERINIT(ObjPricerDouble::scip_init)
{
    cout<<"**************PRICER INIT************ "<<endl;

    int J = inst->getJ();
    int I = inst->getI();

    //cout<<"facility convexity"<<endl;
    //facility convexity constraints
    for (int j = 0 ; j < J ; j++) {
        SCIPgetTransformedCons(scip, Master->conv_lambda_facility.at(j), &(Master->conv_lambda_facility.at(j)));
    }

    //cout<<"customer convexity"<<endl;
    //customer convexity constraints
    for (int i=0 ; i < I ; i++) {
        SCIPgetTransformedCons(scip, Master->conv_lambda_customer.at(i), &(Master->conv_lambda_customer.at(i)));
    }

    //cout<<"equality"<<endl;
    // equality customer / facility on x
    for (int j = 0 ; j < J ; j++) {
        for (int i=0 ; i < I ; i++) {
            SCIPgetTransformedCons(scip, Master->eq_customer_facility_x.at(j*I + i), &(Master->eq_customer_facility_x.at(j*I + i)));
        }
    }

    if (Param.compactCapacityConstraints){
        // equality customer / facility on y
        for (int j = 0 ; j < J ; j++) {
            for (int i=0 ; i < I ; i++) {
                SCIPgetTransformedCons(scip, Master->eq_customer_facility_y.at(j*I + i), &(Master->eq_customer_facility_y.at(j*I + i)));
            }
        }
    }
    else{
        // equality between facilities on y
        for (int j = 0 ; j < J ; j++) {
            for (int i=1 ; i < I ; i++) {
                SCIPgetTransformedCons(scip, Master->eq_customer_facility_y.at(j*I + i), &(Master->eq_customer_facility_y.at(j*I + i)));
            }
        }
    }

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
//SCIP_DECL_PRICERREDCOST(ObjPricerSite::scip_redcost)
//{
//    SCIPdebugMsg(scip, "call scip_redcost ...\n");
//    /* set result pointer, see above */
//    *result = SCIP_SUCCESS;
//    /* call pricing routine */
//    pricingRCFLP(scip,0);
//    return SCIP_OKAY;
//}

SCIP_RETCODE ObjPricerDouble::scip_redcost(SCIP* scip, SCIP_PRICER* pricer, SCIP_Real* lowerbound, SCIP_Bool* stopearly, SCIP_RESULT* result)
{

    SCIPdebugMsg(scip, "call scip_redcost ...\n");

    if( Param.PriceAndBranch && SCIPgetDepth(scip) != 0 )
    {
        *result = SCIP_SUCCESS;
        return SCIP_OKAY;
    }

    /* set result pointer, see above */
    *result = SCIP_SUCCESS;

    //check if convergence between upper and lower bounds has been reached
    if (( (SCIPgetSolOrigObj(scip,NULL) - currentLowerBound)/(0.0000000001 + SCIPgetSolOrigObj(scip,NULL)) > 0.0001) || !Param.useLowerBound){
        /* call pricing routine */
        pricingRCFLP(scip,0);
    }

    return SCIP_OKAY;

}

SCIP_RETCODE ObjPricerDouble::scip_farkas( SCIP* scip, SCIP_PRICER* pricer, SCIP_RESULT* result ){

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



void ObjPricerDouble::updateDualCosts_facility(SCIP* scip, DualCostsFacility & dual_cost, bool Farkas) {
    ///// RECUPERATION DES COUTS DUAUX


    int print = 0 ;
    int J = inst->getJ();
    int I = inst->getI();

    //cout << "solution duale :" << endl ;
    //couts duaux contrainte égalité en x
    for (int j = 0 ; j < J ; j++) {
        for (int i=1 ; i < I ; i++) {
            if (!Farkas) {
                dual_cost.Omega1[j*I + i] = SCIPgetDualsolLinear(scip, Master->eq_customer_facility_x.at(j*I + i) );
            }
            else{
                dual_cost.Omega1[j*I + i] = SCIPgetDualfarkasLinear(scip, Master->eq_customer_facility_x.at(j*I + i) );
            }

            if (print)
                cout << "omega(" << j <<"," << i <<") = " << dual_cost.Omega1[j*I + i] <<endl;
        }
    }

    //cout << "solution duale :" << endl ;
    //couts duaux contrainte égalité en y
    for (int j = 0 ; j < J ; j++) {
        for (int i=1 ; i < I ; i++) {
            if (!Farkas) {
                dual_cost.Omega2[j*I + i] = SCIPgetDualsolLinear(scip, Master->eq_customer_facility_y.at(j*I + i) );
            }
            else{
                dual_cost.Omega2[j*I + i] = SCIPgetDualfarkasLinear(scip, Master->eq_customer_facility_y.at(j*I + i) );
            }

            if (print)
                cout << "omega(" << j <<"," << i <<") = " << dual_cost.Omega2[j*I + i] <<endl;
        }
    }

    //couts duaux contrainte convexité facility
    for (int j = 0 ; j < J ; j++) {
        if (!Farkas) {
            dual_cost.Sigma[j] = SCIPgetDualsolLinear(scip, Master->conv_lambda_facility[j]);
        }
        else{
            dual_cost.Sigma[j] = SCIPgetDualfarkasLinear(scip, Master->conv_lambda_facility[j]);
        }
        if (print) cout << "sigma: " << dual_cost.Sigma[j] <<endl;
    }

}
void ObjPricerDouble::updateDualCosts_customer(SCIP* scip, DualCostsCustomer & dual_cost, bool Farkas) {
    ///// RECUPERATION DES COUTS DUAUX

    int print = 0 ;
    int J = inst->getJ();
    int I = inst->getI();

    //cout << "solution duale :" << endl ;
    //couts duaux contrainte égalité en x
    for (int j = 0 ; j < J ; j++) {
        for (int i=1 ; i < I ; i++) {
            if (!Farkas) {
                dual_cost.Omega1[j*I + i] = SCIPgetDualsolLinear(scip, Master->eq_customer_facility_x.at(j*I + i) );
            }
            else{
                dual_cost.Omega1[j*I + i] = SCIPgetDualfarkasLinear(scip, Master->eq_customer_facility_x.at(j*I + i) );
            }
            if (print)
                cout << "omega1(" << j <<"," << i <<") = " << dual_cost.Omega1[j*I + i] <<endl;
        }
    }

    //cout << "solution duale :" << endl ;
    //couts duaux contrainte égalité en y
    for (int j = 0 ; j < J ; j++) {
        for (int i=1 ; i < I ; i++) {
            if (!Farkas) {
                dual_cost.Omega2[j*I + i] = SCIPgetDualsolLinear(scip, Master->eq_customer_facility_y.at(j*I + i) );
            }
            else{
                dual_cost.Omega2[j*I + i] = SCIPgetDualfarkasLinear(scip, Master->eq_customer_facility_y.at(j*I + i) );
            }

            if (print)
                cout << "omega2(" << j <<"," << i <<") = " << dual_cost.Omega2[j*I + i] <<endl;
        }
    }


    //couts duaux "customer convexity constraint"
    for (int i=1 ; i < I ; i++) {
        if (!Farkas) {
            dual_cost.Sigma.at(i) = SCIPgetDualsolLinear(scip, Master->conv_lambda_customer.at(i));
        }
        else{
            dual_cost.Sigma.at(i) = SCIPgetDualfarkasLinear(scip, Master->conv_lambda_customer.at(i));
        }
        if (print)
            cout << "sigma(" << i <<") = " << dual_cost.Sigma[i] <<endl;
    }

}

void ObjPricerDouble::pricingRCFLP( SCIP*              scip  , bool Farkas             /**< SCIP data structure */)
{
#ifdef OUTPUT_PRICER
    cout<<"**************PRICER************ "<< endl ;
    // SCIPprintBestSol(scip, NULL, FALSE);
#endif

    iteration++;

    // if (Master->nbIter == 1) {
    //     cout << "RMP value : " << endl;
    //     //SCIPprintSol(scip, NULL, NULL, FALSE);
    //     SCIPwriteLP(scip, "debug.lp");
    // }

    totalDualCost = 0;
    int J = inst->getJ();
    int I = inst->getI();

    int print = 1;

    facilityVarsToAdd.clear() ;
    customerVarsToAdd.clear() ;

    //int iteration_limit=5 ;
    //    /// PMR courant et sa solution
    // SCIPwriteTransProblem(scip, NULL, NULL, FALSE);

    //   // cout << "solution du PMR:" << endl ;
    //    SCIPprintSol(scip, NULL, NULL, FALSE);

    //    //cout << "solution réalisable:" << endl ;
    //    SCIPprintBestSol(scip, NULL, FALSE);

    ofstream convergence("convergence/" + std::to_string(Master->J) + "_" + std::to_string(Master->I) + ".csv", std::ofstream::out | std::ofstream::app);

    currentLowerBound = 0;

    double redcost = 0;
    double objvalue = 0;
    IloNumArray xPlan  ;
    IloNumArray yPlan  ;
    bool solutionFound ;


////////// MISE A JOUR DES COUTS DUAUX
    DualCostsCustomer dual_cost_customer = DualCostsCustomer(inst) ;
    updateDualCosts_customer(scip, dual_cost_customer, Farkas);
    DualCostsFacility dual_cost_facility = DualCostsFacility(inst) ;
    updateDualCosts_facility(scip, dual_cost_facility, Farkas);


    // Recherche par facility

    for (int j = 0 ; j < J ; j++) {

        if (print) cout << "facility "<< j << endl;

        if (!Param.DynProgFacility) {
            (AlgoCplex_facility[j])->updateObjCoefficients(inst, Param, dual_cost_facility, Farkas) ;
            xPlan = IloNumArray((AlgoCplex_facility[j])->env, I) ;
            if (Param.compactCapacityConstraints){
                yPlan = IloNumArray((AlgoCplex_facility[j])->env, 1) ;
            }
            solutionFound = (AlgoCplex_facility[j])->findImprovingSolution(inst, dual_cost_facility, objvalue) ;
        }

        else { // résolution par programmation dynamique
            (AlgoDynProg_facility[j])->updateObjCoefficients(inst, Param, dual_cost_facility, Farkas) ;
            xPlan = IloNumArray((AlgoDynProg_facility[j])->env, I) ;
            if (Param.compactCapacityConstraints){
                yPlan = IloNumArray((AlgoDynProg_facility[j])->env, 1) ;
            }
            solutionFound = (AlgoDynProg_facility[j])->findImprovingSolution(inst, dual_cost_facility, objvalue) ;
        }
        
        if (!solutionFound) {
            // Pricer detected an infeasibility : we should immediately stop pricing
            infeasibilityDetected = true ;
            break ;
        }

        // TODO : currentLowerBound += objvalue + dual_cost.Sigma[s] ;

        if (objvalue < - Param.Epsilon) {

            if (!Param.DynProgFacility) {
                (AlgoCplex_facility[j])->getSolution(inst, dual_cost_facility, xPlan, yPlan, Farkas);
            }
            else{
                (AlgoDynProg_facility[j])->getSolution(inst, dual_cost_facility, xPlan, yPlan, Farkas);
            }

            MasterFacility_Variable* lambda = new MasterFacility_Variable(j, xPlan, yPlan);
            if (print) cout << "Plan found for facility " << j << " with reduced cost = " << objvalue << " "  << endl ;

            totalDualCost += objvalue;
            facilityVarsToAdd.push_back(lambda) ;
        }
    }


    // Recherche par client

    for (int i=1 ; i < I ; i++) {

        if (print) cout << "client "<< i << endl;

        if (!Param.DynProgCustomer) {
            (AlgoCplex_customer[i])->updateObjCoefficients(inst, Param, dual_cost_customer, Farkas) ;
            xPlan = IloNumArray((AlgoCplex_customer[i])->env, J) ;
            yPlan = IloNumArray((AlgoCplex_customer[i])->env, J) ;
            solutionFound = (AlgoCplex_customer[i])->findImprovingSolution(inst, dual_cost_customer, objvalue) ;
        }

        else { // résolution par programmation dynamique
            (AlgoDynProg_customer[i])->updateObjCoefficients(inst, Param, dual_cost_customer, Farkas) ;
            xPlan = IloNumArray((AlgoDynProg_customer[i])->env, J) ;
            yPlan = IloNumArray((AlgoCplex_customer[i])->env, J) ;
            solutionFound = (AlgoDynProg_customer[i])->findImprovingSolution(inst, dual_cost_customer, objvalue) ;
            }

        if (!solutionFound) {
            // Pricer detected an infeasibility : we should immediately stop pricing
            infeasibilityDetected = true ;
            break ;
        }

        // TODO : currentLowerBound += objvalue + dual_cost_customer.Sigma[t] ;

        if (objvalue < - Param.Epsilon) {

            if (!Param.DynProgFacility) {
                (AlgoCplex_customer[i])->getSolution(inst, dual_cost_customer, xPlan, yPlan, Farkas);
            }
            else{
                (AlgoDynProg_customer[i])->getSolution(inst, dual_cost_customer, xPlan, yPlan, Farkas);
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

    if (!infeasibilityDetected){

        MasterFacility_Variable* lambdaFacility ;
        
        while (!facilityVarsToAdd.empty()){

            lambdaFacility = facilityVarsToAdd.front() ;

            //// CREATION D'UNE NOUVELLE VARIABLE DANS LE MASTER
            Master->initMasterFacilityVariable(scip, lambdaFacility) ;

            /* add new variable to the list of variables to price into LP (score: leave 1 here) */
            SCIP_RETCODE ajout = SCIPaddPricedVar(scip, lambdaFacility->ptr, 1.0);
            cout << "ajout var par facility: " << ajout << endl;

            ///// ADD COEFFICIENTS TO CONVEXITY and customer/SITE EQUALITY CONSTRAINTS
            Master->addCoefsToConstraints_facilityVar(scip, lambdaFacility) ;

            facilityColumns++;

            facilityVarsToAdd.pop_front() ;
        }

        MasterCustomer_Variable* lambdaCustomer ;

        while (!customerVarsToAdd.empty()){

            lambdaCustomer = customerVarsToAdd.front() ;

            //// CREATION D'UNE NOUVELLE VARIABLE
            Master->initMasterCustomerVariable(scip, lambdaCustomer) ;

            /* add new variable to the list of variables to price into LP (score: leave 1 here) */
            SCIP_RETCODE ajout = SCIPaddPricedVar(scip, lambdaCustomer->ptr, 1.0);
            cout << "ajout var par temps: " << ajout << endl;

            ///// ADD COEFFICIENTS TO DEMAND, POWER LIMITS and CONVEXITY CONSTRAINTS
            Master->addCoefsToConstraints_customerVar(scip, lambdaCustomer) ;

            customerColumns++;

            customerVarsToAdd.pop_front() ;
        }
    }

    infeasibilityDetected = false ;

    cout << "total: " << totalDualCost << endl;

    if (Param.nodeLimit == 1){
        Master->totalDualCostList.push_back(totalDualCost);
        convergence << iteration << "," << fmax(currentLowerBound,0) << "," << SCIPgetSolOrigObj(scip,NULL) << endl;
    }

#ifdef OUTPUT_PRICER
    SCIPwriteTransProblem(scip, "RCFLP.lp", "lp", FALSE);
    cout<<"************END PRICER******************"<<endl;
#endif

}

void ObjPricerDouble::addVarBound(SCIP_ConsData* consdata) {

        //A implémenter

}

void ObjPricerDouble::removeVarBound(SCIP_ConsData* consdata) {

        //A implémenter

}