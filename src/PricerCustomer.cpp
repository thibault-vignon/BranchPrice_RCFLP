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
ObjPricerCustomerRCFLP::ObjPricerCustomerRCFLP(
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

    if (!Param.DynProgTime) {
        for (int i=0 ; i < I ; i++) {
            AlgoCplex[i] = new CplexPricingAlgoCustomer(inst, param, i) ;
        }
    }
    else {
        for (int i=0 ; i < I ; i++) {
            AlgoDynProg.at(t) = new DynProgPricingAlgoCustomer(inst, param, t) ;
        }
    }

    cout << "ici fin constructeur" << endl ;
}


/** Destructs the pricer object. */
ObjPricerCustomerRCFLP::~ObjPricerCustomerRCFLP()
{
    cout<<"Destructeur du pricer"<<endl;
}

/** initialization method of variable pricer (called after problem was transformed)
 *
 *  Because SCIP transformes the original problem in preprocessing, we need to get the references to
 *  the variables and constraints in the transformed problem from the references in the original
 *  problem.
 */
SCIP_DECL_PRICERINIT(ObjPricerCustomerRCFLP::scip_init)
{

    cout<<"**************PRICER INIT************ "<<endl;
    int I = inst->getI() ;

    //convexity constraint
        for (int i=0 ; i < I ; i++) {
        SCIPgetTransformedCons( scip, Master->convexity_cstr.at(i), &(Master->convexity_cstr.at(i)) );
    }

    // TODO : autres contraintes

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

SCIP_RETCODE ObjPricerCustomerRCFLP::scip_redcost(SCIP* scip, SCIP_PRICER* pricer, SCIP_Real* lowerbound, SCIP_Bool* stopearly, SCIP_RESULT* result)
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

SCIP_RETCODE ObjPricerCustomerRCFLP::scip_farkas( SCIP* scip, SCIP_PRICER* pricer, SCIP_RESULT* result ){

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



void ObjPricerCustomerRCFLP::updateDualCosts(SCIP* scip, DualCostsTime & dual_cost, bool Farkas) {
    ///// RECUPERATION DES COUTS DUAUX

    int print = 0 ;
    int I = inst->getI() ;

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
            cout << "sigma(" << t <<") = " << dual_cost.Sigma[t] <<endl;
    }

    // TODO : autres contraintes
}

void ObjPricerCustomerRCFLP::pricingRCFLP( SCIP*              scip  , bool Farkas             /**< SCIP data structure */)
{
#ifdef OUTPUT_PRICER
   //cout<<"**************PRICER************ " << endl ;
    // SCIPprintBestSol(scip, NULL, FALSE);
#endif

    int print = 0 ;

    cout<<"**************PRICER************ " << endl ;

    
    iteration++;

  // SCIPwriteTransProblem(scip, NULL, NULL, FALSE);
    //SCIPprintSol(scip, NULL, NULL, FALSE);

//        /// PMR courant et sa solution
//        SCIPwriteTransProblem(scip, NULL, NULL, FALSE);

//        // cout << "solution du PMR:" << endl ;
//        SCIPprintSol(scip, NULL, NULL, FALSE);

//        //cout << "solution réalisable:" << endl ;
//        SCIPprintBestSol(scip, NULL, FALSE);

    int I = inst->getI() ;
    int J = inst->getJ() ;

    customerVarsToAdd.clear() ;

    // Cout duaux
    DualCostsCustomer dual_cost = DualCostsCustomer(inst) ;
    updateDualCosts(scip, dual_cost, Farkas);








    while (!oneImprovingSolution && cas < nb_cas) { // Si on n'a pas trouvé de colonne améliorante dans le cas 1, on passe au cas 2: on cherche une colonne pour tous les pas de temps

        cas ++ ;

        if (Param.OneTimeStepPerIter) {
            lastTimeStep = (lastTimeStep+1)%T;
            min=lastTimeStep ;
            max=min+1 ;
        }

        //cout << "cas: " << cas << endl ;

        for (int t=min ; t < max ; t++) {

            if (print) cout << "time "<< t << endl;

            if (TimeSolNotFound.at(t) < 2 || cas==2 ) { // si un plan pour t a été généré au cours des 10 dernières itérations

                ///// MISE A JOUR DES OBJECTIFS DES SOUS PROBLEMES
                // cout << "mise à jour des couts, farkas=" << Farkas << endl;
                if (!Param.DynProgTime) {
                    (AlgoCplex.at(t))->updateObjCoefficients(inst, Param, dual_cost, Farkas) ;
                }
                else {
                    (AlgoDynProg.at(t))->updateObjCoefficients(inst, Param, dual_cost, Farkas) ;
                }

                //// CALCUL D'UN PLAN DE COUT REDUIT MINIMUM
                double objvalue = 0 ;
                double temps ;
                bool solutionFound;

                if (!Param.DynProgTime) {
                    solutionFound = (AlgoCplex.at(t))->findImprovingSolution(inst, dual_cost, objvalue, temps, 1);
                }
                else {
                    solutionFound = (AlgoDynProg.at(t))->findImprovingSolution(inst, dual_cost, objvalue, temps, Param.heurPricingTime);
                }
                nbCallsToCplex++;
                Master->cumul_resolution_pricing += temps ;


                if (!solutionFound) {
                    // Pricer detected an infeasibility : we should immediately stop pricing
                    infeasibilityDetected = true ;
                    break ;
                }

                if (objvalue < -epsilon) {
                    oneImprovingSolution = true ;

                    TimeSolNotFound.at(t) = 0 ;

                    double realCost=0 ;
                    double totalProd=0 ;

                    IloNumArray upDownPlan ;
                    IloNumArray powerPlan ;
                    if (!Param.DynProgTime) {
                        upDownPlan = IloNumArray((AlgoCplex.at(t))->env, n) ;
                        powerPlan = IloNumArray((AlgoCplex.at(t))->env, n) ;
                        (AlgoCplex.at(t))->getUpDownPlan(inst, dual_cost, upDownPlan, powerPlan,realCost, totalProd, Farkas) ;
                    }
                    else {
                        upDownPlan = IloNumArray((AlgoDynProg.at(t))->env, n) ;
                        powerPlan = IloNumArray((AlgoDynProg.at(t))->env, n) ;
                        (AlgoDynProg.at(t))->getUpDownPlan(inst, dual_cost, upDownPlan, powerPlan, realCost, totalProd, Farkas) ;
                    }

                     //

                    //cout << "total prod: " << totalProd << endl ;

                    // if (print) {
                         cout << "Minimum reduced cost plan: "<< objvalue << "for time " << t << endl ;
                         for (int i=0 ; i < n ; i++) {
                             cout << fabs(upDownPlan[i]) << " " ;
                         }
                         cout << endl ;
                     //}




                    int kmax = t ;
                    int kmin=t ;

                    if (Param.AddColumnToOtherCustomerSteps) {
                        if (t < T-1) {
                            if (inst->getD(t) < inst->getD(t+1)) {
                                kmax = fmin(T-1, t+T/12) ;
                            }
                        }
                        if (t>0) {
                            if (inst->getD(t) < inst->getD(t-1)) {
                                kmin = fmax(0, t-T/12) ;
                            }
                        }
                    }

                    for (int k=kmin ; k <= kmax; k++) {

                        if (inst->getD(k) <= totalProd) {

                            //  cout << "ajout maitre" << endl;

                            /// AJOUT VARIABLE DANS LE MAITRE ////

                            MasterCustomer_Variable* lambda = new MasterCustomer_Variable(k, upDownPlan);
                            lambda->addPowerPlan(powerPlan);

                            timeVarsToAdd.push_back(lambda) ;
                        }

                    }

                }

                else {
                    if (Param.DontPriceAllTimeSteps) {
                        TimeSolNotFound.at(t) ++;
                    }
                }
            }
        }
    }

    // If we detected an infeasibility in one of the pricers, we add no variables to the Master
    // This will force termination of column generation, the master problem will remain infeasible
    // And this will cause the node to be pruned

    cout << "infeasibility detected: " << infeasibilityDetected << endl;
    if (!infeasibilityDetected){

        MasterCustomer_Variable* lambdaTime ;

        while (!timeVarsToAdd.empty()){

            lambdaTime = timeVarsToAdd.front() ;

            //// CREATION D'UNE NOUVELLE VARIABLE
            Master->initMasterCustomerVariable(scip, lambdaTime) ;

            /* add new variable to the list of variables to price into LP (score: leave 1 here) */
            SCIP_RETCODE ajout = SCIPaddPricedVar(scip, lambdaTime->ptr, 1.0);
            cout << "ajout var par temps: " << ajout << endl;

            ///// ADD COEFFICIENTS TO DEMAND, POWER LIMITS and CONVEXITY CONSTRAINTS
            Master->addCoefsToConstraints(scip, lambdaTime) ;

            timeColumns++;

            timeVarsToAdd.pop_front() ;
        }
    }

    infeasibilityDetected = false ;

  //  cout<<"************END PRICER******************"<<endl;
#ifdef OUTPUT_PRICER
    SCIPwriteTransProblem(scip, "RCFLP.lp", "lp", FALSE);
  //  cout<<"************END PRICER******************"<<endl;
#endif

}

void ObjPricerCustomerRCFLP::addVarBound(SCIP_ConsData* consdata) {

    int t = consdata->time ;
    int i = consdata->unit ;

    if (AlgoDynProg.at(t) != NULL) {
        if (consdata->bound == 0) {
            if (!Param.powerPlanGivenByMu) {
                (AlgoDynProg.at(t))->W -= inst->getPmax(i) ;
            }
            (AlgoDynProg.at(t))->init.at(i) = 0 ;
        }
        else {
            (AlgoDynProg.at(t))->init.at(i) = 1 ;
        }
    }
    else {
        //A implémenter
    }
}

void ObjPricerCustomerRCFLP::removeVarBound(SCIP_ConsData* consdata) {

    int t = consdata->time ;
    int i = consdata->unit ;

    if (AlgoDynProg.at(t) != NULL) {
        if (consdata->bound == 0) {
            (AlgoDynProg.at(t))->W += inst->getPmax(i) ;
        }

        (AlgoDynProg.at(t))->init.at(i) = -1 ;
    }
    else {
        //A implémenter
    }
}


