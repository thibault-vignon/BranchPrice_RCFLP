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
    ObjPricerRCFLP(scip, pp_name, instance, param), ParamMaster(param)
{
    Master=M ;

    J = inst->getJ();
    I = inst->getI();

    AlgoCplex_facility = vector<CplexPricingAlgo*>(J, NULL) ;
    AlgoDynProg_facility = vector<DynProgPricingAlgo*>(J, NULL) ;

    if (!Param.DynProg) {
        for (int j = 0 ; j < J ; j++) {
            AlgoCplex_facility.at(s) = new CplexPricingAlgoFacility(inst, param, s) ;
        }
    }
    else {
        for (int j = 0 ; j < J ; j++) {
            AlgoDynProg_facility.at(s) = new DynProgPricingAlgoFacility(inst, param, s) ;
        }
    }


    AlgoCplex_customer = vector<CplexPricingAlgoCustomer*>(I, NULL) ;
    AlgoDynProg_customer = vector<DynProgPricingAlgoCustomer*>(I, NULL) ;

    if (!Param.DynProgTime) {
        for (int i=0 ; i < I ; i++) {
            AlgoCplex_customer[i] = new CplexPricingAlgoCustomer(inst, param, i) ;
        }
    }
    else {
        for (int i=0 ; i < I ; i++) {
            AlgoDynProg_customer.at(t) = new DynProgPricingAlgoCustomer(inst, param, t) ;
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

    J = inst->getJ();
    I = inst->getI();

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



void ObjPricerDouble::updateDualCosts_site(SCIP* scip, DualCosts & dual_cost, bool Farkas) {
    ///// RECUPERATION DES COUTS DUAUX


    int print = 0 ;
    J = inst->getJ();
    I = inst->getI();

    //cout << "solution duale :" << endl ;
    //couts duaux contrainte égalité
    for (int j = 0 ; j < J ; j++) {
        for (int i=1 ; i < I ; i++) {
            if (!Farkas) {
                dual_cost.Omega1[j*I + i] = SCIPgetDualsolLinear(scip, Master->eq_customer_facility_x.at(j*I + i) );
            }
            else{
                dual_cost.Omega1[j*I + i] = SCIPgetDualfarkasLinear(scip, Master->eq_customer_facility_x.at(j*I + i) );
            }

            if (print)
                cout << "omega(" << i <<"," << t <<") = " << dual_cost.Omega[i*T+t] <<endl;
        }
    }


    //couts duaux contrainte convexité site
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
void ObjPricerDouble::updateDualCosts_time(SCIP* scip, DualCostsTime & dual_cost, bool Farkas) {
    ///// RECUPERATION DES COUTS DUAUX

    int print = 0 ;
    int n = inst->getn() ;
    int T = inst->getT() ;
    int S = Param.nbDecGpes ;

    //cout << "solution duale :" << endl ;
    //couts duaux contrainte égalité
    for (int i = 0; i < n; i++) {
        for (int t = 0 ; t < T ; t++) {
            if (!Farkas) {
                dual_cost.Omega[i*T+t] = SCIPgetDualsolLinear(scip, Master->eq_time_site.at(i*T+t) );
            }
            else{
                dual_cost.Omega[i*T+t] = SCIPgetDualfarkasLinear(scip, Master->eq_time_site.at(i*T+t) );
            }
            if (print)
                cout << "omega(" << i <<"," << t <<") = " << dual_cost.Omega[i*T+t] <<endl;
        }
    }


    //couts duaux "time convexity constraint"
    for (int t = 0 ; t < T ; t++) {
        if (!Farkas) {
            dual_cost.Sigma.at(t) = SCIPgetDualsolLinear(scip, Master->conv_lambda_time.at(t));
        }
        else{
            dual_cost.Sigma.at(t) = SCIPgetDualfarkasLinear(scip, Master->conv_lambda_time.at(t));
        }
        if (print)
            cout << "sigma(" << t <<") = " << dual_cost.Sigma[t] <<endl;
    }

    if (Param.minUpDownDouble) {
        for (int i = 0; i < n; i++) {
            for (int t = 1 ; t < T ; t++) {
                if (!Farkas) {
                    dual_cost.Mu.at(i*T+t) = SCIPgetDualsolLinear(scip, Master->logical.at(i*T+t));
                }
                else{
                    dual_cost.Mu.at(i*T+t) = SCIPgetDualfarkasLinear(scip, Master->logical.at(i*T+t));
                }
                if (print)
                    cout << "mu(" << i <<"," << t <<") = " << dual_cost.Mu[i*T+t] <<endl;
            }
        }

        //couts duaux "min-up constraint"
        for (int i = 0; i < n; i++) {
            int L = inst->getL(i) ;
            for (int t = L ; t < T ; t++) {
                if (!Farkas) {
                    dual_cost.Nu.at(i*T+t) = SCIPgetDualsolLinear(scip, Master->min_up.at(i*T+t));
                }
                else{
                    dual_cost.Nu.at(i*T+t) = SCIPgetDualfarkasLinear(scip, Master->min_up.at(i*T+t));
                }
                if (print)
                    cout << "nu(" << i <<"," << t <<") = " << dual_cost.Nu.at(i*T+t) <<endl;
            }
        }

        //couts duaux "min-down constraint"
        for (int i = 0; i < n; i++) {
            int l = inst->getl(i) ;
            for (int t = l ; t < T ; t++) {
                if (!Farkas) {
                    dual_cost.Xi.at(i*T+t) = SCIPgetDualsolLinear(scip, Master->min_down.at(i*T+t));
                }
                else{
                    dual_cost.Xi.at(i*T+t) = SCIPgetDualfarkasLinear(scip, Master->min_down.at(i*T+t));
                }
                if (print)
                    cout << "xi(" << i <<"," << t <<") = " << dual_cost.Xi.at(i*T+t) <<endl;
            }
        }
    }


}

void ObjPricerDouble::pricingRCFLP( SCIP*              scip  , bool Farkas             /**< SCIP data structure */)
{
#ifdef OUTPUT_PRICER
    cout<<"**************PRICER************ "<< endl ;
    // SCIPprintBestSol(scip, NULL, FALSE);
#endif

    Master->nbIter++;

    // if (Master->nbIter == 1) {
    //     cout << "RMP value : " << endl;
    //     //SCIPprintSol(scip, NULL, NULL, FALSE);
    //     SCIPwriteLP(scip, "debug.lp");
    // }

    totalDualCost = 0;
    int T = inst->getT() ;
    int n = inst->getn() ;

    int print = 1;
    iteration++;

    siteVarsToAdd.clear() ;
    timeVarsToAdd.clear() ;

    //int iteration_limit=5 ;
    //    /// PMR courant et sa solution
    // SCIPwriteTransProblem(scip, NULL, NULL, FALSE);

    //   // cout << "solution du PMR:" << endl ;
    //    SCIPprintSol(scip, NULL, NULL, FALSE);

    //    //cout << "solution réalisable:" << endl ;
    //    SCIPprintBestSol(scip, NULL, FALSE);

    ofstream convergence("convergence/" + std::to_string(Master->n) + "_" + std::to_string(Master->T) + + "_" + std::to_string((Master->inst)->id) + ".csv", std::ofstream::out | std::ofstream::app);

    currentLowerBound = 0;

    //// Cout duaux
    int S = Param.nbDecGpes ;

////////// MISE A JOUR DES COUTS DUAUX

// Cout duaux de time à mettre à jour avant ceux de Site, car la méthode computeObjCoef prend en arg dual_cost_time 
    DualCostsTime dual_cost_time = DualCostsTime(inst) ;
    updateDualCosts_time(scip, dual_cost_time, Farkas);
    DualCosts dual_cost = DualCosts(inst,Param) ;
    updateDualCosts_site(scip, dual_cost, Farkas);

    double epsilon= 0.0000001 ;
    double redcost = 0;

    if (!guidageTermine && (!Param.unitGEQTime || !Param.balanceCosts)) {
        guidageTermine = true ;
    }

    if (!guidageTermine){

        double facteur = 1 ;

        switch(Param.guidageRepartition){
            case 1:
                if (iteration > sqrt(n*T) ){
                    facteur = fmax( (3*sqrt(n*T) - iteration) / (2*sqrt(n*T)), 0) ;
                }
                break ;

            case 2:
                if (iteration > sqrt(n*T) ){
                    facteur = fmax( (2*sqrt(n*T) - iteration) / (2*sqrt(n*T)), 0) ;
                }
                break ;

            case 3:
                facteur = pow(0.5, ( (iteration - 1) / (2*sqrt(n*T)) ) ) ;
                break ;
        }

        if (facteur > 0.01){
            for (int i=0 ; i < n ; i++) {
                Param.costBalancingPricer.at(i) = Param.costBalancingPricer.at(i) * facteur + Param.costBalancingMaster.at(i) * (1 - facteur) ;
            }
        }
        else{
            guidageTermine = true ;
        }

    }

    if (!guidageTermine){

        dual_cost.computeObjCoef(inst,Param,Farkas, dual_cost_time);

        /////Recherche variable améliorante de type sites
        for (int s = 0 ; s < S ; s++) {

            //cout << "site "<< s << endl;

            ///// MISE A JOUR DES OBJECTIFS DES SOUS PROBLEMES
            // cout << "mise à jour des couts, farkas=" << Farkas << endl;
            if (!Param.DynProg) {
                (AlgoCplex_site[s])->updateObjCoefficients(inst, Param, dual_cost, Farkas) ;
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
                upDownPlan = IloNumArray((AlgoCplex_site[s])->env, Param.nbUnits(s)*T) ;
                solutionFound = (AlgoCplex_site[s])->findUpDownPlan(inst, dual_cost, upDownPlan, objvalue) ;
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
                cout << "début dynprogsite" << endl;
                cout << "s = " << s << endl;
                cout << AlgoDynProg_site[s]->Site << endl;
                upDownPlan = IloNumArray((AlgoDynProg_site[s])->env, Param.nbUnits(s)*T) ;
                cout << "updownplan initialisé" << endl;
                if (Param.DynProgSUSD) {
                    (AlgoDynProg_site.at(s))->findImprovingSolutionSUSD(inst, dual_cost, objvalue);
                    (AlgoDynProg_site.at(s))->getUpDownPlanSUSD(inst, upDownPlan) ;
                }
                else {
                    (AlgoDynProg_site.at(s))->findImprovingSolution(inst, dual_cost, objvalue);
                    (AlgoDynProg_site.at(s))->getUpDownPlan(inst, upDownPlan) ;
                }
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
            cout << endl ;
            }


            //if (SCIPisNegative(scip, objvalue)) {

            if (objvalue < -epsilon ) {

                Master_Variable* lambda = new Master_Variable(s, upDownPlan);
                if (print) cout << "Plan found for site " << s << " with reduced cost = " << objvalue << " "  << endl ;

                if (Param.powerPlanGivenByLambda && !Param.DynProg) {
                    powerPlan = IloNumArray((AlgoCplex_site[s])->env, Param.nbUnits(s)*T) ;
                    (AlgoCplex_site[s]->cplex).getValues(AlgoCplex_site[s]->p, powerPlan) ;
                    if (print) cout << "power plan: " << powerPlan << endl;
                    lambda->addPowerPlan(powerPlan);
                }

                //SCIPwriteOrigProblem(scip, "debug.lp", "lp", FALSE);

                dual_cost.computeObjCoef(inst,ParamMaster,Farkas, dual_cost_time);
                dual_cost.computeRedcost(inst, ParamMaster, lambda, redcost);
                cout << "cout reduit calculé: " << redcost << endl;

                if (redcost < -epsilon) {
                    totalDualCost += redcost;
                    siteVarsToAdd.push_back(lambda) ;
                }

                dual_cost.computeObjCoef(inst,Param,Farkas, dual_cost_time);

                if (Param.stopFirstSite){
                    break;
                }
            }
        }

        if (!Param.oneRoundTime || (totalDualCost > -Param.Epsilon) || (Master->nbIter == 2) ){

            if (print) cout << "RECHERCHE SUR LES PAS DE TEMPS" << endl ;


            /////Recherche variable améliorante de type time

            //// Cout duaux


            //cout << "cas: " << cas << endl ;
            int min=0 ;
            int max=T ;

            if (Param.OneTimeStepPerIter) {
                lastTimeStep = (lastTimeStep+1)%T;
                min=lastTimeStep ;
                max=min+1 ;
            }

            for (int t=min ; t < max ; t++) {

                if (print) cout << "time "<< t << endl;


                ///// MISE A JOUR DES OBJECTIFS DES SOUS PROBLEMES
                // cout << "mise à jour des couts, farkas=" << Farkas << endl;
                if (!Param.DynProgTime) {
                    cout << "début updateobjcoeff" << endl;
                    (AlgoCplex_time.at(t))->updateObjCoefficients(inst, Param, dual_cost_time, Farkas) ;
                    if(t == 0 && Master->nbIter == 2){
                        (AlgoCplex_time.at(t))->cplex.exportModel( ( std::to_string(Param.PminDifferentPmax) + "_bug.lp" ).c_str() );
                        for (int i=0 ; i<n ; i++) {
                            cout << "i:" << i << endl;
                            cout << dual_cost_time.Omega.at(i*T+t) << endl;
                        }
                    }
                }
                else {
                    cout << "début updateobjcoeff" << endl;
                    (AlgoDynProg_time.at(t))->updateObjCoefficients(inst, Param, dual_cost_time, Farkas) ;
                }

                //// CALCUL D'UN PLAN DE COUT REDUIT MINIMUM
                double objvalue = 0 ;
                double temps ;
                bool solutionFound;

                if (!Param.DynProgTime) {
                    cout << "début findimprovingsol" << endl;
                    solutionFound = (AlgoCplex_time.at(t))->findImprovingSolution(inst, dual_cost_time, objvalue, temps, 1);
                }
                else {
                    cout << "début findimprovingsol" << endl;
                    solutionFound = (AlgoDynProg_time.at(t))->findImprovingSolution(inst, dual_cost_time, objvalue, temps, Param.heurPricingTime);
                }
                Master->cumul_resolution_pricing += temps ;

                if (!solutionFound) {
                    // Pricer detected an infeasibility : we should immediately stop pricing
                    infeasibilityDetected = true ;
                    break ;
                }

                if (objvalue < -epsilon) {

                    timeStepColumns.at(t) += 1;

                    double realCost=0 ;
                    double totalProd=0 ;

                    IloNumArray upDownPlan  ;
                    IloNumArray powerPlan  ;

                    if (Param.DynProgTime){
                        upDownPlan = IloNumArray((AlgoDynProg_time.at(t))->env, n) ;
                        powerPlan = IloNumArray((AlgoDynProg_time.at(t))->env, n) ;

                        (AlgoDynProg_time.at(t))->getUpDownPlan(inst, dual_cost_time, upDownPlan, powerPlan, realCost, totalProd, Farkas) ;
                    }

                    else{
                        upDownPlan = IloNumArray((AlgoCplex_time.at(t))->env, n) ;
                        powerPlan = IloNumArray((AlgoCplex_time.at(t))->env, n) ;

                        (AlgoCplex_time.at(t))->getUpDownPlan(inst, dual_cost_time, upDownPlan, powerPlan, realCost, totalProd, Farkas) ;
                    }


                    if (print) {
                        cout << "Minimum reduced cost plan: "<< objvalue << "for time " << t << endl ;
                        for (int i=0 ; i < n ; i++) {
                            cout << fabs(upDownPlan[i]) << " " ;
                        }
                        cout << endl ;
                    }

                    //cout << "total prod: " << totalProd << endl ;

                    /// AJOUT VARIABLE DANS LE MAITRE ////

                    MasterTime_Variable* lambda = new MasterTime_Variable(t, upDownPlan);
                    // cout << "Plan found for time " << t << " with reduced cost = " << objvalue << " ";

                    if (Param.powerPlanGivenByMu) {
                        lambda->addPowerPlan(powerPlan);
                    }

                    timeVarsToAdd.push_back(lambda) ;
                    totalDualCost += objvalue ;

                    if (Param.stopFirstTime){
                        break;
                    }
                }
            }
        }

        if (totalDualCost > -epsilon){
            guidageTermine = true;
        }

    }

    if (guidageTermine) {
        cout << "début pricing test" << endl;

        dual_cost.computeObjCoef(inst,ParamMaster,Farkas, dual_cost_time);

        for (int s = 0 ; s < S ; s++) {

            double objvalue = 0;
            IloNumArray upDownPlan  ;
            IloNumArray powerPlan  ;
            bool solutionFound ;

            if (!Param.DynProg) {
                (AlgoCplex_site[s])->updateObjCoefficients(inst, ParamMaster, dual_cost, Farkas) ;
            }

            if (!Param.DynProg) {
                upDownPlan = IloNumArray((AlgoCplex_site[s])->env, Param.nbUnits(s)*T) ;
                solutionFound= (AlgoCplex_site[s])->findUpDownPlan(inst, dual_cost, upDownPlan, objvalue) ;
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
                cout << AlgoDynProg_site[s]->Site << endl;
                upDownPlan = IloNumArray((AlgoDynProg_site[s])->env, Param.nbUnits(s)*T) ;
                if (Param.DynProgSUSD) {
                    (AlgoDynProg_site.at(s))->findImprovingSolutionSUSD(inst, dual_cost, objvalue);
                    (AlgoDynProg_site.at(s))->getUpDownPlanSUSD(inst, upDownPlan) ;
                }
                else {
                    solutionFound = (AlgoDynProg_site.at(s))->findImprovingSolution(inst, dual_cost, objvalue);
                    (AlgoDynProg_site.at(s))->getUpDownPlan(inst, upDownPlan) ;
                }
            }

            if (!solutionFound) {
                // Pricer detected an infeasibility : we should immediately stop pricing
                infeasibilityDetected = true ;
                break ;
            }

            currentLowerBound += objvalue + dual_cost.Sigma[s] ;

            if (objvalue < -epsilon) {

                Master_Variable* lambda = new Master_Variable(s, upDownPlan);
                if (print) cout << "Plan found for site " << s << " with reduced cost = " << objvalue << " "  << endl ;

                if (Param.powerPlanGivenByLambda && !Param.DynProg) {
                    powerPlan = IloNumArray((AlgoCplex_site[s])->env, Param.nbUnits(s)*T) ;
                    (AlgoCplex_site[s]->cplex).getValues(AlgoCplex_site[s]->p, powerPlan) ;
                    if (print) cout << "power plan: " << powerPlan << endl;
                    lambda->addPowerPlan(powerPlan);
                }

                totalDualCost += objvalue;
                siteVarsToAdd.push_back(lambda) ;
            }
        }


        // RECHERCHE SACS A DOS PAS DE TEMPS

        int min=0 ;
        int max=T ;

        if (Param.OneTimeStepPerIter) {
            lastTimeStep = (lastTimeStep+1)%T;
            min=lastTimeStep ;
            max=min+1 ;
        }

        for (int t=min ; t < max ; t++) {

            if (print) cout << "time "<< t << endl;

            if (!Param.DynProgTime) {
                (AlgoCplex_time.at(t))->updateObjCoefficients(inst, ParamMaster, dual_cost_time, Farkas) ;
            }
            else {
                (AlgoDynProg_time.at(t))->updateObjCoefficients(inst, ParamMaster, dual_cost_time, Farkas) ;
            }

            double objvalue = 0 ;
            double temps ;
            bool solutionFound;

            if (!Param.DynProgTime) {
                solutionFound = (AlgoCplex_time.at(t))->findImprovingSolution(inst, dual_cost_time, objvalue, temps, 1);
            }
            else {
                solutionFound = (AlgoDynProg_time.at(t))->findImprovingSolution(inst, dual_cost_time, objvalue, temps, Param.heurPricingTime);
            }
            Master->cumul_resolution_pricing += temps ;

            if (!solutionFound) {
                // Pricer detected an infeasibility : we should immediately stop pricing
                infeasibilityDetected = true ;
                break ;
            }

            currentLowerBound += objvalue + dual_cost_time.Sigma[t] ;

            if (objvalue < -epsilon) {

                timeStepColumns.at(t) += 1;

                double realCost=0 ;
                double totalProd=0 ;

                IloNumArray upDownPlan  ;
                IloNumArray powerPlan  ;

                if (Param.DynProgTime){
                    upDownPlan = IloNumArray((AlgoDynProg_time.at(t))->env, n) ;
                    powerPlan = IloNumArray((AlgoDynProg_time.at(t))->env, n) ;

                    (AlgoDynProg_time.at(t))->getUpDownPlan(inst, dual_cost_time, upDownPlan, powerPlan, realCost, totalProd, Farkas) ;
                }

                else{
                    upDownPlan = IloNumArray((AlgoCplex_time.at(t))->env, n) ;
                    powerPlan = IloNumArray((AlgoCplex_time.at(t))->env, n) ;

                    (AlgoCplex_time.at(t))->getUpDownPlan(inst, dual_cost_time, upDownPlan, powerPlan, realCost, totalProd, Farkas) ;
                }

                MasterTime_Variable* lambda = new MasterTime_Variable(t, upDownPlan);
                cout << "Plan found for time " << t << " with reduced cost = " << objvalue << " ";

                if (Param.powerPlanGivenByMu) {
                    lambda->addPowerPlan(powerPlan);
                }

                timeVarsToAdd.push_back(lambda) ;
                totalDualCost += objvalue;
            }
        }
    }

    // Ajouts terme de droite lorsque les contraintes de demande sont dans le maître
    if (guidageTermine && Param.PminDifferentPmax && !Param.powerPlanGivenByMu){
        for (int t = 0 ; t < T ; t++) {
            currentLowerBound += +dual_cost.Mu[t] * inst->getD(t);
        }
    }

    // If we detected an infeasibility in one of the pricers, we add no variables to the Master
    // This will force termination of column generation, the master problem will remain infeasible
    // And this will cause the node to be pruned

    if (!infeasibilityDetected){

        Master_Variable* lambdaSite ;
        
        while (!siteVarsToAdd.empty()){

            lambdaSite = siteVarsToAdd.front() ;

            //// CREATION D'UNE NOUVELLE VARIABLE DANS LE MASTER
            Master->initMasterSiteVariable(scip, inst, lambdaSite) ;

            /* add new variable to the list of variables to price into LP (score: leave 1 here) */
            SCIP_RETCODE ajout = SCIPaddPricedVar(scip, lambdaSite->ptr, 1.0);
            cout << "ajout var par unité: " << ajout << endl;

            ///// ADD COEFFICIENTS TO CONVEXITY and TIME/SITE EQUALITY CONSTRAINTS
            Master->addCoefsToConstraints_siteVar(scip, lambdaSite, inst) ;

            unitColumns++;

            siteVarsToAdd.pop_front() ;
        }

        MasterTime_Variable* lambdaTime ;

        while (!timeVarsToAdd.empty()){

            lambdaTime = timeVarsToAdd.front() ;

            //// CREATION D'UNE NOUVELLE VARIABLE
            Master->initMasterTimeVariable(scip, lambdaTime) ;

            /* add new variable to the list of variables to price into LP (score: leave 1 here) */
            SCIP_RETCODE ajout = SCIPaddPricedVar(scip, lambdaTime->ptr, 1.0);
            cout << "ajout var par temps: " << ajout << endl;

            ///// ADD COEFFICIENTS TO DEMAND, POWER LIMITS and CONVEXITY CONSTRAINTS
            Master->addCoefsToConstraints_timeVar(scip, lambdaTime) ;

            timeColumns++;

            timeVarsToAdd.pop_front() ;
        }
    }

    infeasibilityDetected = false ;

    cout << "total: " << totalDualCost << endl;

    if (Param.nodeLimit == 1){
        Master->totalDualCostList.push_back(totalDualCost);
        convergence << Master->nbIter << "," << fmax(currentLowerBound,0) << "," << SCIPgetSolOrigObj(scip,NULL) << endl;
    }

#ifdef OUTPUT_PRICER
    SCIPwriteTransProblem(scip, "RCFLP.lp", "lp", FALSE);
    cout<<"************END PRICER******************"<<endl;
#endif

}

void ObjPricerDouble::addVarBound(SCIP_ConsData* consdata) {

    int t = consdata->time ;
    int i = consdata->unit ;

    cout << "Enter addVarBound:" << endl;

    if (!Param.DynProg) {
        cout << "branchement non supporté pour sous-problèmes par unité résolus par Cplex en double décomposition" << endl;
    }
    else {
        (AlgoDynProg_site[consdata->site])->branchingDecisions.at(t) = consdata->bound ;
        cout << "for unit " << consdata->site << ", at time " << t <<", bound set to " << consdata->bound << endl ;
        cout << "End addVarBound" << endl ;
    }

    if (AlgoDynProg_time.at(t) != NULL) {
        if (consdata->bound == 0) {
            if (!Param.powerPlanGivenByMu) {
                (AlgoDynProg_time.at(t))->W -= inst->getPmax(i) ;
            }
            (AlgoDynProg_time.at(t))->init.at(i) = 0 ; 
        }
        else {
            (AlgoDynProg_time.at(t))->init.at(i) = 1 ;
        }
    }
    else {
        //A implémenter
    }
    cout << "placed var bound on time pricer" << endl;

    guidageTermine = false ;
}

void ObjPricerDouble::removeVarBound(SCIP_ConsData* consdata) {

    int t = consdata->time ;
    int i = consdata->unit ;

    if (!Param.DynProg) {
        cout << "branchement non supporté pour sous-problèmes par unité résolus par Cplex en double décomposition" << endl;
    }
    else {
        (AlgoDynProg_site[consdata->site])->branchingDecisions.at(t) = 8 ;
    }

    if (AlgoDynProg_time.at(t) != NULL) {
        if (consdata->bound == 0 && !Param.powerPlanGivenByMu) {
            (AlgoDynProg_time.at(t))->W += inst->getPmax(i) ;
        }

        (AlgoDynProg_time.at(t))->init.at(i) = -1 ;
    }
    else {
        //A implémenter
    }
}


