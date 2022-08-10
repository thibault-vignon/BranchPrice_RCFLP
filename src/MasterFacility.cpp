#include "Master.h"

/* namespace usage */
using namespace std;
using namespace scip;

MasterFacility_Variable::MasterFacility_Variable(int f, IloNumArray x, IloNumArray y) {
    ptr = NULL ;
    facility = f ;
    cost = 0 ;
    x_plan = x ;
    y_plan = y ;
}

void MasterFacility_Variable::computeCost(InstanceRCFLP* inst, const Parameters & Param) {
    //compute cost of up/down plan lambda: fixed cost (including minimum power output cost) and start up cost
    cost=0 ;

    if (Param.doubleDecompo){
        for (int i=0 ; i < inst->getI() ; i++) {
            if (x_plan[i] > 1 - Param.Epsilon){
                cost += 0.5 * (2 - Param.balanceCostsX) * inst->geta(facility,i);
            }
            if (Param.compactCapacityConstraints && y_plan[0] > 1 - Param.Epsilon){
                cost += 0.5 * Param.balanceCostsY * inst->getb(facility);
            }
        }
    }

    else{
        for (int i=0 ; i < inst->getI() ; i++) {
            if (x_plan[i] > 1 - Param.Epsilon){
                cost += inst->geta(facility,i);
            }
        }
        if (y_plan[0] > 1 - Param.Epsilon){
            cost += inst->getb(facility);
        }
    }
}


void MasterFacility_Model::addCoefsToConstraints(SCIP* scip, MasterFacility_Variable* lambda) {

    int f = lambda->facility ;

    /* add coefficient to the convexity constraint for site s */
    SCIPaddCoefLinear(scip, convexity_cstr[f], lambda->ptr, 1.0) ;

    /* add coefficient to the assignment constraints*/
    for (int i = 0 ; i<I ; i++){
        if (lambda->x_plan[i] > 1 - Param.Epsilon){
            SCIPaddCoefLinear(scip, assignment_cstr[i], lambda->ptr, 1.0) ;
        }
    }

    /* add coefficients to the reliability constraints*/
    for (int i = 0 ; i<I ; i++){

        if (lambda->x_plan[i] > 1 - Param.Epsilon){
            SCIPaddCoefLinear(scip, reliability_cstr[f*I + i], lambda->ptr, - inst->getd(i) * inst->getK()) ;
        }

        for (int j = 0 ; j<J ; j++){
            for (int indice=0; indice<inst->getv(); indice++){
                if ( (inst->getV(j)[indice] == f) && (lambda->y_plan[0] > 1 - Param.Epsilon) ){
                    SCIPaddCoefLinear(scip, reliability_cstr[j*I + i], lambda->ptr, inst->getc(f)) ;
                }
            }
        }
    }
}


void MasterFacility_Model::initMasterFacilityVariable(SCIP* scip, MasterFacility_Variable* var) {
    char var_name[255];
    SCIPsnprintf(var_name, 255, "V_%d",L_var.size());
    SCIPdebugMsg(scip, "new variable <%s>\n", var_name);

    /* create the new variable: Use upper bound of infinity such that we do not have to care about
     * the reduced costs of the variable in the pricing. The upper bound of 1 is implicitly satisfied
     * due to the set partitioning constraints.
     */

    var->computeCost(inst, Param);
    double cost= var->cost;

    cout << "cost: " << cost << endl ;

    SCIP_Vartype type = SCIP_VARTYPE_INTEGER;
 
    SCIPcreateVar(scip, &(var->ptr), var_name,
                  0.0,                     // lower bound
                  SCIPinfinity(scip),      // upper bound
                  cost,                     // objective
                  type,    // variable type
                  false, false, NULL, NULL, NULL, NULL, NULL);

    //// Add new variable to the list
    L_var.push_back(var);

    cout << "Variable " << var_name << " added, with x plan: "  ;

    for (int i=0 ; i < inst->getI() ; i++) {
        cout << var->x_plan[i] << " "  ;
    }
    cout << endl ;

    cout << "and y plan: " << var->y_plan[0] ;

    cout << endl ;
}

MasterFacility_Model::MasterFacility_Model(InstanceRCFLP* inst, const Parameters & Parametres) : Master_Model(Parametres, inst) {
    convexity_cstr.resize(J, (SCIP_CONS*) NULL) ;
    reliability_cstr.resize(I*J, (SCIP_CONS*) NULL) ;
    assignment_cstr.resize(I, (SCIP_CONS*) NULL) ;
}

void  MasterFacility_Model::initScipMasterFacilityModel(SCIP* scip) {


    ////////////////////////////////////////////////////////////////
    /////////////   MASTER CONSTRAINT INITIALIZATION   /////////////
    ////////////////////////////////////////////////////////////////
    // Constraints form: lhs <= ax <= rhs

    cout << "convex cons" << endl ;

    ///// Convexity constraint ////
    char con_name_convex[255];
    for (int j = 0 ; j<J ; j++)
    {
        SCIP_CONS* con = NULL;
        (void) SCIPsnprintf(con_name_convex, 255, "Convexity(%d)", j); // nom de la contrainte
        SCIPcreateConsLinear( scip, &con, con_name_convex, 0, NULL, NULL,
                              1.0,   // lhs
                              1.0,   // rhs  SCIPinfinity(scip) if >=1
                              true,  /* initial */
                              false, /* separate */
                              true,  /* enforce */
                              true,  /* check */
                              true,  /* propagate */
                              false, /* local */
                              true,  /* modifiable */
                              false, /* dynamic */
                              false, /* removable */
                              false  /* stickingatnode */ );
        SCIPaddCons(scip, con);
        convexity_cstr[j] = con;
    }

    cout << "assignment cons" << endl ;

    ///// Assignment constraint ////
    char con_name_assignment[255];
    for (int i = 0 ; i<I ; i++)
    {
        SCIP_CONS* con = NULL;
        (void) SCIPsnprintf(con_name_assignment, 255, "Assignment(%d)", i); // nom de la contrainte
        SCIPcreateConsLinear( scip, &con, con_name_assignment, 0, NULL, NULL,
                              1.0,   // lhs
                              1.0,   // rhs  SCIPinfinity(scip) if >=1
                              true,  /* initial */
                              false, /* separate */
                              true,  /* enforce */
                              true,  /* check */
                              true,  /* propagate */
                              false, /* local */
                              true,  /* modifiable */
                              false, /* dynamic */
                              false, /* removable */
                              false  /* stickingatnode */ );
        SCIPaddCons(scip, con);
        assignment_cstr[i] = con;
    }

    ///// Reliability constraint ////
    char con_name_reliability[255];
    for (int j = 0 ; j<J ; j++){
        for (int i = 0 ; i<I ; i++)
        {
            SCIP_CONS* con = NULL;
            (void) SCIPsnprintf(con_name_reliability, 255, "Reliability(%d,%d)", j, i); // nom de la contrainte
            SCIPcreateConsLinear( scip, &con, con_name_reliability, 0, NULL, NULL,
                                0.0,   // lhs
                                SCIPinfinity(scip),   // rhs
                                true,  /* initial */
                                false, /* separate */
                                true,  /* enforce */
                                true,  /* check */
                                true,  /* propagate */
                                false, /* local */
                                true,  /* modifiable */
                                false, /* dynamic */
                                false, /* removable */
                                false  /* stickingatnode */ );
            SCIPaddCons(scip, con);
            reliability_cstr[j*I + i] = con;
        }
    }

    ///////////////////////////////////////////////////////////////
    //////////   MASTER LAMBDA VARIABLES INITIALIZATION   /////////
    ///////////////////////////////////////////////////////////////

    L_var.clear();

    if (!Param.Farkas){
        for (int j=0 ; j<J; j++)
        {
            IloNumArray x_plan = IloNumArray(env, I) ;
            IloNumArray y_plan = IloNumArray(env, 1) ;

            // TODO : Initialize plans corresponding to feasible solution

            // MasterFacility_Variable* lambda = new MasterFacility_Variable(j, x_plan, y_plan);

            // initMasterFacilityVariable(scip, lambda);

            // SCIPaddVar(scip, lambda->ptr);

            // addCoefsToConstraints(scip, lambda) ;
        }
    }
}


void MasterFacility_Model::computeFracSol(SCIP* scip) {
    list<MasterFacility_Variable*>::const_iterator itv;
    SCIP_Real frac_value;
    for (int ind=0 ; ind < I*J ; ind++) {
        x_frac[ind]=0;
    }
    for (int ind=0 ; ind < J ; ind++) {
        y_frac[ind]=0;
    }

    for (itv = L_var.begin(); itv!=L_var.end(); itv++) {

        frac_value = fabs(SCIPgetVarSol(scip,(*itv)->ptr));

        int f = (*itv)->facility ;
        
        if ((*itv)->y_plan[0] > 1 - Param.Epsilon) {
            y_frac[f] += frac_value ;
        }

        for (int i=0 ; i < I ; i++) {
            if ((*itv)->x_plan[i] > 1 - Param.Epsilon) {
                x_frac[f*I + i] += frac_value ;
            }
        }
    }
}


void MasterFacility_Model::discardVar(SCIP* scip, SCIP_ConsData* consdata) {

    /////On met à 0 les lambda incompatibles avec la contrainte

    consdata->L_var_bound.clear() ; // L_var_bound stocke les variables scip dont la borne a été effectivement changée (ie elle n'était pas déjà à 0)

    list<MasterFacility_Variable*>::const_iterator itv;

    /* TODO : modifier pour mettre en place le branchement

    for (itv = L_var.begin(); itv!=L_var.end(); itv++) {
        if () {
            if ( != consdata->bound ) {

                SCIP_Real old_bound =  SCIPgetVarUbAtIndex(scip, (*itv)->ptr, NULL, 0) ;

                ///  L_var_bound est mis à jour
                if (!SCIPisZero(scip,old_bound)) {
                    SCIPchgVarUbNode(scip, NULL, (*itv)->ptr, 0) ;
                    consdata->L_var_bound.push_back((*itv)->ptr) ;
                }
            }
        }
    }
    */
}

void MasterFacility_Model::restoreVar(SCIP* scip, SCIP_ConsData* consdata) {

    ////On remet à +inf les lambda qui étaient incompatibles avec la contrainte de branchement

    list<SCIP_VAR*>::const_iterator itv;

    for (itv = consdata->L_var_bound.begin(); itv!=consdata->L_var_bound.end(); itv++) {
        SCIPchgVarUbNode(scip, NULL, (*itv), SCIPinfinity(scip)) ;
    }
    consdata->L_var_bound.clear() ;

}


//////// Créé des variables lambda à partir d'une solution (x,y) ///////////
void MasterFacility_Model::createColumns(SCIP* scip, IloNumArray x, IloNumArray y) {

    for (int j = 0 ; j < J ; j++) {
        IloNumArray x_plan = IloNumArray(env, I) ;
        IloNumArray y_plan = IloNumArray(env, 1) ;

        y_plan[0] = y[j] ;
        for (int i=0 ; i < I ; i++) {
            x_plan[i] = x[j*I + i] ;
        }

        MasterFacility_Variable* lambda = new MasterFacility_Variable(j, x_plan, y_plan);
        initMasterFacilityVariable(scip, lambda);
        SCIPaddVar(scip, lambda->ptr);
        addCoefsToConstraints(scip, lambda) ;
    }
}
