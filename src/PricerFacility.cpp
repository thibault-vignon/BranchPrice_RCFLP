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
    Master=M ;
    AlgoCplex = vector<CplexPricingAlgoFacility*>(J, NULL) ;
    AlgoDynProg = vector<DynProgPricingAlgoFacility*>(J, NULL) ;

    J = inst->getJ();

    if (!Param.DynProg) {
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



void ObjPricerFacility::updateDualCosts(SCIP* scip, DualCosts & dual_cost, bool Farkas) {
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














    for (int s = 0 ; s < S ; s++) {

       //cout << "site "<< s << endl;

        ///// MISE A JOUR DES OBJECTIFS DES SOUS PROBLEMES
       // cout << "mise à jour des couts, farkas=" << Farkas << endl;
        if (!Param.DynProg) {
            (AlgoCplex[s])->updateObjCoefficients(inst, Param, dual_cost, Farkas) ;
        }
        else {
            //nothing to do
           // (AlgoDynProg[s])->updateObjCoefficients(inst, Param, dual_cost, Farkas) ;

        }

        //// CALCUL D'UN PLAN DE COUT REDUIT MINIMUM
        double objvalue = 0 ;
        IloNumArray upDownPlan  ;
        IloNumArray powerPlan  ;
        bool solutionFound ;


        if (!Param.DynProg) {
            upDownPlan = IloNumArray((AlgoCplex[s])->env, Param.nbUnits(s)*T) ;
            solutionFound = (AlgoCplex[s])->findUpDownPlan(inst, dual_cost, upDownPlan, objvalue) ;
            for (int index=0 ; index <Param.nbUnits(s)*T ; index++ ) {
                if (upDownPlan[index]>1-epsilon) {
                    upDownPlan[index]=1 ;
                }
                if (upDownPlan[index]< epsilon) {
                    upDownPlan[index]=0 ;
                }
            }
        }

        else { // résolution par programmation dynamique

            upDownPlan = IloNumArray((AlgoDynProg[s])->env, Param.nbUnits(s)*T) ;
            solutionFound = (AlgoDynProg.at(s))->findImprovingSolution(inst, dual_cost, objvalue);
            (AlgoDynProg.at(s))->getUpDownPlan(inst, upDownPlan) ;

            cout << "DP resolution done" << endl ;

        }

        // cout << "solution found: " << solutionFound << endl;
        if (!solutionFound) {
            // Pricer detected an infeasibility : we should immediately stop pricing
            infeasibilityDetected = true ;
            break ;
        }

        if (print) cout << "Minimum reduced cost plan: "<< objvalue << endl ;

        if (print) {
            for (int t=0 ; t < T ; t++)  {
                for (int i=0 ; i < Param.nbUnits(s) ; i++) {
                    cout << fabs(upDownPlan[i*T+t]) << " " ;
                }
                //cout << endl ;
            }
            cout << endl ;
        }

        cout << endl ;

        //if (SCIPisNegative(scip, objvalue)) {

        if (objvalue < -epsilon ) {

            Master_Variable* lambda = new Master_Variable(s, upDownPlan);
            cout << "Plan found for site " << s << " with reduced cost = " << objvalue << " "  << endl ;

            if (Param.powerPlanGivenByLambda && !Param.DynProg) {
                powerPlan = IloNumArray((AlgoCplex[s])->env, Param.nbUnits(s)*T) ;
                (AlgoCplex[s]->cplex).getValues(AlgoCplex[s]->p, powerPlan) ;
                cout << "power plan: " << powerPlan << endl;
                lambda->addPowerPlan(powerPlan);
            }

            siteVarsToAdd.push_back(lambda) ;
        }
    }

    // If we detected an infeasibility in one of the pricers, we add no variables to the Master
    // This will force termination of column generation, the master problem will remain infeasible
    // And this will cause the node to be pruned

    cout << "infeasibility detected: " << infeasibilityDetected << endl;
    if (!infeasibilityDetected){

        Master_Variable* lambdaSite ;
        
        while (!siteVarsToAdd.empty()){

            lambdaSite = siteVarsToAdd.front() ;

            //// CREATION D'UNE NOUVELLE VARIABLE DANS LE MASTER
            Master->initMasterVariable(scip, inst, lambdaSite) ;

            /* add new variable to the list of variables to price into LP (score: leave 1 here) */
            SCIP_RETCODE ajout = SCIPaddPricedVar(scip, lambdaSite->ptr, 1.0);
            cout << "ajout var par unité: " << ajout << endl;

            ///// ADD COEFFICIENTS TO CONVEXITY and TIME/SITE EQUALITY CONSTRAINTS
            Master->addCoefsToConstraints(scip, lambdaSite, inst) ;

            unitColumns++;

            siteVarsToAdd.pop_front() ;
        }
    }

    infeasibilityDetected = false ;

#ifdef OUTPUT_PRICER
    SCIPwriteTransProblem(scip, "RCFLP.lp", "lp", FALSE);
    cout<<"************END PRICER******************"<<endl;
#endif

}

void ObjPricerFacility::addVarBound(SCIP_ConsData* consdata) {

    cout << "Enter addVarBound:" << endl;

    if (!Param.DynProg) {

    }
    else {

        cout << "for unit " << consdata->site << ", at time " << consdata->time <<", bound set to " << consdata->bound << endl ;
        cout << "End addVarBound" << endl ;
    }
}

void ObjPricerFacility::removeVarBound(SCIP_ConsData* consdata) {

    if (!Param.DynProg) {

    }
    else {

    }
}


