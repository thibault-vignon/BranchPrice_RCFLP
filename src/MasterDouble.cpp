#include "Master.h"

/* namespace usage */
using namespace std;
using namespace scip;


void MasterDouble_Model::addCoefsToConstraints_facilityVar(SCIP* scip, MasterFacility_Variable* lambda) {

    int f = lambda->facility ;
    
    /* add coefficient to the convexity constraint for facility f */
    SCIPaddCoefLinear(scip, conv_lambda_facility.at(f), lambda->ptr, 1.0) ;

    // add coef to the equality x(facility) = x(customer)
    for (int i = 0 ; i < I ; i++) {
        if (lambda->x_plan[i] > 1 - Param.Epsilon){
            SCIPaddCoefLinear(scip, eq_customer_facility_x.at(f*I + i), lambda->ptr, 1.0) ;
        }
    }

    if (Param.compactCapacityConstraints && lambda->y_plan[0] > 1 - Param.Epsilon){
        // add coef to the equality y(facility) = y(customers)
        for (int i = 0 ; i < I ; i++) {
            SCIPaddCoefLinear(scip, eq_customer_facility_y.at(f*I + i), lambda->ptr, 1.0) ;
        }
    }
}

void MasterDouble_Model::addCoefsToConstraints_customerVar(SCIP* scip, MasterCustomer_Variable* lambda) {

    int c = lambda->customer ;

    /* add coefficient to the convexity constraint for site s */
    SCIPaddCoefLinear(scip, conv_lambda_customer.at(c), lambda->ptr, 1.0) ;

    // add coef to the equality x(facility) = x(customer)
    for (int j = 0 ; j < J ; j++) {
        if (lambda->x_plan[j] > 1 - Param.Epsilon){
            SCIPaddCoefLinear(scip, eq_customer_facility_x.at(j*I + c), lambda->ptr, -1.0) ;
        }
    }

    // add coef to the equality y(facility) = y(customers)
    if (Param.compactCapacityConstraints){
        for (int j = 0 ; j < J ; j++) {
            if (lambda->y_plan[j] > 1 - Param.Epsilon){
                SCIPaddCoefLinear(scip, eq_customer_facility_y.at(j*I + c), lambda->ptr, -1.0) ;
            }
        }
    }
    else{
        if (c > 0){
            for (int j = 0 ; j < J ; j++) {
                if (lambda->y_plan[j] > 1 - Param.Epsilon){
                    SCIPaddCoefLinear(scip, eq_customer_facility_y.at(j*I + c), lambda->ptr, -1.0) ;
                }
            }
        }
        else{
            // On utilise la copie du client 0 comme rhs
            for (int i = 1 ; i < I ; i++) {
                for (int j = 0 ; j < J ; j++) {
                    if (lambda->y_plan[j] > 1 - Param.Epsilon){
                        SCIPaddCoefLinear(scip, eq_customer_facility_y.at(j*I + i), lambda->ptr, 1.0) ;
                    }
                }
            }
        }
    }
}

void MasterDouble_Model::initMasterFacilityVariable(SCIP* scip, MasterFacility_Variable* var) {
    char var_name[255];
    SCIPsnprintf(var_name, 255, "V_%d",L_var_facility.size());
    SCIPdebugMsg(scip, "new variable <%s>\n", var_name);

    /* create the new variable: Use upper bound of infinity such that we do not have to care about
     * the reduced costs of the variable in the pricing. The upper bound of 1 is implicitly satisfied
     * due to the set partitioning constraints.
     */

    var->computeCost(inst, Param);
    double cost= var->cost;

    //cout << "cost: " << cost << endl ;

    SCIP_Vartype type = SCIP_VARTYPE_INTEGER;

    SCIPcreateVar(scip, &(var->ptr), var_name,
                  0.0,                     // lower bound
                  SCIPinfinity(scip),      // upper bound
                  cost,                     // objective
                  type,    // variable type
                  false, false, NULL, NULL, NULL, NULL, NULL);
    //// Add new variable to the list
    L_var_facility.push_back(var);
}

//////// Initialisation d'une variable lambda(time) /////////////
void MasterDouble_Model::initMasterCustomerVariable(SCIP* scip, MasterCustomer_Variable* var) {

    char var_name[255];
    SCIPsnprintf(var_name, 255, "W_%d",L_var_customer.size());
    SCIPdebugMsg(scip, "new variable <%s>\n", var_name);

    /* create the new variable: Use upper bound of infinity such that we do not have to care about
     * the reduced costs of the variable in the pricing. The upper bound of 1 is implicitly satisfied
     * due to the set partitioning constraints.
     */

    var->computeCost(inst, Param);
    double cost= var->cost;
    //cout << var_name << ", cost: " << cost << endl ;

    SCIP_Vartype type  = SCIP_VARTYPE_INTEGER;

    SCIPcreateVar(scip, &(var->ptr), var_name,
                  0.0,                     // lower bound
                  SCIPinfinity(scip),      // upper bound
                  cost,                     // objective
                  type,    // variable type
                  false, false, NULL, NULL, NULL, NULL, NULL);

    //// Add new variable to the list
    L_var_customer.push_back(var);
}

MasterDouble_Model::MasterDouble_Model(InstanceRCFLP* inst, const Parameters & Parametres) : Master_Model(Parametres, inst) {

    conv_lambda_facility.resize(J, (SCIP_CONS*) NULL) ;
    conv_lambda_customer.resize(I, (SCIP_CONS*) NULL) ;
    eq_customer_facility_x.resize(I*J, (SCIP_CONS*) NULL) ;
    eq_customer_facility_y.resize(I*J, (SCIP_CONS*) NULL) ;
}

void  MasterDouble_Model::initScipMasterDoubleModel(SCIP* scip) {

    ////////////////////////////////////////////////////////////////
    /////////////   MASTER CONSTRAINT INITIALIZATION   /////////////
    ////////////////////////////////////////////////////////////////
    // Constraints form: lhs <= ax <= rhs

    ///// Equality facility/client on x variables constraint /////
    char con_name_eq_customer_facility_x[255];
    SCIP_Real upper_bound =  0 ;

    for (int i = 0 ; i < I ; i++)
    {
        for (int j = 0; j < J; j++)
        {
            SCIP_CONS* con = NULL;
            (void) SCIPsnprintf(con_name_eq_customer_facility_x, 255, "EqCustomerFacilityX(%d,%d)", j, i); // nom de la contrainte
            SCIPcreateConsLinear( scip, &con, con_name_eq_customer_facility_x, 0, NULL, NULL,
                                  0.0,   // lhs
                                  upper_bound,   // rhs  SCIPinfinity(scip) if >=1
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
            eq_customer_facility_x.at(j*I + i) = con;
        }
    }

    ///// Equality facility/client on y variables constraint /////
    char con_name_eq_customer_facility_y[255];

    for (int i = 0 ; i < I ; i++)
    {
        for (int j = 0; j < J; j++)
        {
            SCIP_CONS* con = NULL;
            (void) SCIPsnprintf(con_name_eq_customer_facility_y, 255, "EqCustomerFacilityY(%d,%d)", j, i); // nom de la contrainte
            SCIPcreateConsLinear( scip, &con, con_name_eq_customer_facility_y, 0, NULL, NULL,
                                  0.0,   // lhs
                                  upper_bound,   // rhs  SCIPinfinity(scip) if >=1
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
            eq_customer_facility_y.at(j*I + i) = con;
        }
    }


    ///// Convexity constraint : site////
    char con_name_convex_facility[255];
    for (int j = 0 ; j < J ; j++)
    {
        SCIP_CONS* con = NULL;
        (void) SCIPsnprintf(con_name_convex_facility, 255, "ConvFacility(%d)", j); // nom de la contrainte
        SCIPcreateConsLinear( scip, &con, con_name_convex_facility, 0, NULL, NULL,
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
        conv_lambda_facility.at(j) = con;
    }



    ///// Convexity constraint : time ////
    char con_name_convex_customer[255];
    for (int i = 0 ; i < I ; i++)
    {
        SCIP_CONS* con = NULL;
        (void) SCIPsnprintf(con_name_convex_customer, 255, "ConvCustomer(%d)", i); // nom de la contrainte
        SCIPcreateConsLinear( scip, &con, con_name_convex_customer, 0, NULL, NULL,
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
        conv_lambda_customer.at(i) = con;
    }


    ///////////////////////////////////////////////////////////////
    //////////   SITE LAMBDA VARIABLES INITIALIZATION   /////////
    ///////////////////////////////////////////////////////////////

    L_var_facility.clear();

    if (!Param.Farkas){

    }

    ///////////////////////////////////////////////////////////////
    //////////   TIME LAMBDA VARIABLES INITIALIZATION   /////////
    ///////////////////////////////////////////////////////////////

    L_var_customer.clear();

    if (!Param.Farkas){
  
    }

}

void MasterDouble_Model::computeFracSol(SCIP* scip) {
    list<MasterFacility_Variable*>::const_iterator itv;
    SCIP_Real frac_value;
    for (int ind=0 ; ind < I*J ; ind++) {
        x_frac[ind]=0;
    }
    for (int ind=0 ; ind < J ; ind++) {
        y_frac[ind]=0;
    }

    for (itv = L_var_facility.begin(); itv!=L_var_facility.end(); itv++) {

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


void MasterDouble_Model::discardVar(SCIP* scip, SCIP_ConsData* consdata) {

//    /////On met à 0 les lambda incompatibles avec la contrainte

   consdata->L_var_bound.clear() ; // L_var_bound stocke les variables scip dont la borne a été effectivement changée (ie elle n'était pas déjà à 0)

   list<MasterFacility_Variable*>::const_iterator itv_facility;

   /* TODO : modifier pour mettre en place le branchement

   for (itv_facility = L_var_facility.begin(); itv_facility!=L_var_facility.end(); itv_facility++) {
       if ((*itv_facility)->Site == consdata->unit) {
           if ((*itv_facility)->UpDown_plan[consdata->time] != consdata->bound ) {

               SCIP_Real old_bound =  SCIPgetVarUbAtIndex(scip, (*itv_facility)->ptr, NULL, 0) ;

               ///  L_var_bound est mis à jour
               if (!SCIPisZero(scip,old_bound)) {
                   SCIPchgVarUbNode(scip, NULL, (*itv_facility)->ptr, 0) ;
                   consdata->L_var_bound.push_back((*itv_facility)->ptr) ;
               }
           }
       }
   }

   list<MasterCustomer_Variable*>::const_iterator itv_customer;

     for (itv_customer = L_var_customer.begin(); itv_customer!=L_var_customer.end(); itv_customer++) {
         if ((*itv_customer)->time == consdata->time) {
             if ((*itv_customer)->UpDown_plan[consdata->unit] != consdata->bound ) {

                 SCIP_Real old_bound =  SCIPgetVarUbAtIndex(scip, (*itv_customer)->ptr, NULL, 0) ;

                 ///  L_var_bound est mis à jour
                 if (!SCIPisZero(scip,old_bound)) {
                     SCIPchgVarUbNode(scip, NULL, (*itv_customer)->ptr, 0) ;
                     consdata->L_var_bound.push_back((*itv_customer)->ptr) ;
                 }
             }
         }
     }

    */
}

void MasterDouble_Model::restoreVar(SCIP* scip, SCIP_ConsData* consdata) {

   ////On remet à +inf les lambda qui étaient incompatibles avec la contrainte de branchement

   list<SCIP_VAR*>::const_iterator itv;

   for (itv = consdata->L_var_bound.begin(); itv!=consdata->L_var_bound.end(); itv++) {
       SCIPchgVarUbNode(scip, NULL, (*itv), SCIPinfinity(scip)) ;
   }
    consdata->L_var_bound.clear() ;

}


////////// Créé des variables lambda à partir d'une solution (x) ///////////
void MasterDouble_Model::createColumns(SCIP* scip, IloNumArray x, IloNumArray y) {

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
        addCoefsToConstraints_facilityVar(scip, lambda) ;
    }


    for (int i = 0 ; i < I ; i++) {
        IloNumArray x_plan = IloNumArray(env, J) ;
        IloNumArray y_plan = IloNumArray(env, J) ;

        for (int j = 0 ; j < J ; j++) {
            x_plan[j] = x[j*I + i] ;
            y_plan[j] = y[j] ;
        }

        MasterCustomer_Variable* lambda = new MasterCustomer_Variable(i, x_plan, y_plan);
        initMasterCustomerVariable(scip, lambda);
        SCIPaddVar(scip, lambda->ptr);
        addCoefsToConstraints_customerVar(scip, lambda) ;
    }
}
