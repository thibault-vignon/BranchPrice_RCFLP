#include "Master.h"

/* namespace usage */
using namespace std;
using namespace scip;


///////////////////////////////
////////// VARIABLES //////////
///////////////////////////////


////// Constructeur //////
MasterCustomer_Variable::MasterCustomer_Variable(int c, IloNumArray x, IloNumArray y) {
    ptr = NULL ;
    customer = c ;
    cost = 0 ;
    x_plan = x ;
    y_plan = y ;
    cout << y_plan.getSize() << endl;
}

void MasterCustomer_Variable::computeCost(InstanceRCFLP* inst, const Parameters & Param) {
    cost = 0;

    if (Param.doubleDecompo){
        for (int j = 0 ; j < inst->getJ() ; j++) {
            if (y_plan[j] > 1 - Param.Epsilon){
                cost += 0.5 * (2 - Param.balanceCostsY) * inst->getb(j) / inst->getI();
            }
            if (x_plan[j] > 1 - Param.Epsilon){
                cost += 0.5 * Param.balanceCostsX * inst->geta(j,customer);
            }
        }
    }
    else{
        for (int j = 0 ; j < inst->getJ() ; j++) {
            if (y_plan[j] > 1 - Param.Epsilon){
                cost += inst->getb(j) / inst->getI();
            }
            if (x_plan[j] > 1 - Param.Epsilon){
                cost += inst->geta(j,customer);
            }
        }
    }
}

/////// ajout des coefficients dans chaque contrainte pour variable lambda //////////
void MasterCustomer_Model::addCoefsToConstraints(SCIP* scip, MasterCustomer_Variable* lambda) {

    int c = lambda->customer ;

    /* add coefficient to the convexity constraint */
    SCIPaddCoefLinear(scip, convexity_cstr.at(c), lambda->ptr, 1.0) ;

    /* add coefficients to the capacity constraints */
    for (int j = 0 ; j<J ; j++){
        if (lambda->x_plan[j] > 1 - Param.Epsilon){
            SCIPaddCoefLinear(scip, capacity_cstr.at(j), lambda->ptr, - inst->getd(c)) ;
        }

        // il y a autant de copies de chaque y que de clients, donc il ne faut pas simplement additionner les capacités
        if (lambda->y_plan[j] > 1 - Param.Epsilon){
            SCIPaddCoefLinear(scip, capacity_cstr.at(j), lambda->ptr, inst->getc(j) / I) ;
        }
    }

    // add coefficient to the equality constraints between y variables
    if (c > 0){
        for (int j = 0 ; j < J ; j++) {
            if (lambda->y_plan[j] > 1 - Param.Epsilon){
                SCIPaddCoefLinear(scip, equality_cstr.at(j*(I-1) + c-1), lambda->ptr, -1.0) ;
            }
        }
    }
    else{
        // On utilise la copie du client 0 comme rhs
        for (int i = 0 ; i < I-1 ; i++) {
            for (int j = 0 ; j < J ; j++) {
                if (lambda->y_plan[j] > 1 - Param.Epsilon){
                    SCIPaddCoefLinear(scip, equality_cstr.at(j*(I-1) + i), lambda->ptr, 1.0) ;
                }
            }
        }
    }
}


//////// Initialisation d'une variable lambda /////////////
void MasterCustomer_Model::initMasterCustomerVariable(SCIP* scip, MasterCustomer_Variable* var) {

    char var_name[255];
    SCIPsnprintf(var_name, 255, "V_%d",L_var.size());
    SCIPdebugMsg(scip, "new variable <%s>\n", var_name);

    /* create the new variable: Use upper bound of infinity such that we do not have to care about
     * the reduced costs of the variable in the pricing. The upper bound of 1 is implicitly satisfied
     * due to the set partitioning constraints.
     */

    var->computeCost(inst, Param);
    double cost= var->cost;
    //cout << var_name << ", cost: " << cost << endl ;

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

    for (int j=0 ; j < inst->getJ() ; j++) {
        cout << var->x_plan[j] << " "  ;
    }
    cout << endl ;

    cout << "and y plan: " ;

    for (int j=0 ; j < inst->getJ() ; j++) {
        cout << var->y_plan[j] << " "  ;
    }
    cout << endl ;
}

//////// Créé des variables lambda à partir d'une solution (x,p) ///////////
void MasterCustomer_Model::createColumns(SCIP* scip, IloNumArray x, IloNumArray y) {

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
        addCoefsToConstraints(scip, lambda) ;
    }
}

///////////////////////////////////////////
////////// INITIALISATION MASTER //////////
///////////////////////////////////////////

MasterCustomer_Model::MasterCustomer_Model(InstanceRCFLP* instance, const Parameters & Parametres) : Master_Model(Parametres, instance) {
    convexity_cstr.resize(I, (SCIP_CONS*) NULL) ;
    capacity_cstr.resize(J, (SCIP_CONS*) NULL) ;
    equality_cstr.resize((I-1)*J, (SCIP_CONS*) NULL) ;
}

void  MasterCustomer_Model::initScipMasterCustomerModel(SCIP* scip) {

    ////////////////////////////////////////////////////////////////
    /////////////   MASTER CONSTRAINT INITIALIZATION   /////////////
    ////////////////////////////////////////////////////////////////
    // Constraints form: lhs <= ax <= rhs

    cout << "convex cons..." << endl;

    ///// Convexity constraint ////
    char con_name_convex[255];
    for (int i = 0 ; i<I ; i++)
    {
        SCIP_CONS* con = NULL;
        (void) SCIPsnprintf(con_name_convex, 255, "Convexity(%d)", i); // nom de la contrainte
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
        convexity_cstr.at(i) = con;
    }

    cout << "capacity cons..." << endl;

    ///// Capacity constraint ////
    char con_name_capacity[255];
    for (int j = 0 ; j<J ; j++)
    {
        SCIP_CONS* con = NULL;
        (void) SCIPsnprintf(con_name_capacity, 255, "Capacity(%d)", j); // nom de la contrainte
        SCIPcreateConsLinear( scip, &con, con_name_capacity, 0, NULL, NULL,
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
        capacity_cstr.at(j) = con;
    }

    cout << "equality cons..." << endl;

    //Equality between clients on y variables
    char con_name_equality[255];
    for (int i = 0 ; i < I-1 ; i++)
    {
        for (int j = 0 ; j<J ; j++){
            SCIP_CONS* con = NULL;
            (void) SCIPsnprintf(con_name_equality, 255, "Equality(%d,%d)", j, i); // nom de la contrainte
            SCIPcreateConsLinear( scip, &con, con_name_equality, 0, NULL, NULL,
                                0.0,   // lhs
                                0.0,   // rhs 
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
            equality_cstr.at(j*(I-1) + i) = con;
        }
    }

    cout << "constraints initialized" << endl;

    ///////////////////////////////////////////////////////////////
    //////////   MASTER LAMBDA VARIABLES INITIALIZATION   /////////
    ///////////////////////////////////////////////////////////////

    L_var.clear();

    if (!Param.Farkas){
        for (int i=0 ; i<I; i++)
        {
            IloNumArray x_plan = IloNumArray(env, J) ;
            IloNumArray y_plan = IloNumArray(env, J) ;

            // TODO : Initialize plans corresponding to feasible solution

            // MasterCustomer_Variable* lambda = new MasterCustomer_Variable(i, x_plan, y_plan);

            // initMasterCustomerVariable(scip, lambda);

            // SCIPaddVar(scip, lambda->ptr);

            // addCoefsToConstraints(scip, lambda) ;
        }
    }
}



//////////////////////////////////////
////////// METHODES VIRTUELLES ///////
//////////////////////////////////////

void MasterCustomer_Model::computeFracSol(SCIP* scip) {
    list<MasterCustomer_Variable*>::const_iterator itv;
    SCIP_Real frac_value;
    for (int ind=0 ; ind < I*J ; ind++) {
        x_frac[ind]=0;
    }
    for (int ind=0 ; ind < J ; ind++) {
        y_frac[ind]=0;
    }

    for (itv = L_var.begin(); itv!=L_var.end(); itv++) {

        frac_value = fabs(SCIPgetVarSol(scip,(*itv)->ptr));

        int c = (*itv)->customer ;

        for (int j=0 ; j < J ; j++) {
            if ((*itv)->x_plan[j] > 1 - Param.Epsilon) {
                x_frac[j*I + c] += frac_value ;
            }
        }

        // Variables y : splitting donc on divise par le nombre de copies

        for (int j=0 ; j < J ; j++) {
            if ((*itv)->y_plan[j] > 1 - Param.Epsilon) {
                y_frac[j] += frac_value / I;
            }
        }
    }
}


void MasterCustomer_Model::discardVar(SCIP* scip, SCIP_ConsData* consdata) {

    /////On met à 0 les lambda incompatibles avec la contrainte

     consdata->L_var_bound.clear() ; // L_var_bound stocke les variables scip dont la borne a été effectivement changée (ie elle n'était pas déjà à 0)

     list<MasterCustomer_Variable*>::const_iterator itv;

     /* TODO : modifier pour mettre en place le branchement

     for (itv = L_var.begin(); itv!=L_var.end(); itv++) {
         if ((*itv)->time == consdata->time) {
             if ((*itv)->UpDown_plan[consdata->unit] != consdata->bound ) {

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

void MasterCustomer_Model::restoreVar(SCIP* scip, SCIP_ConsData* consdata) {

    ////On remet à +inf les lambda qui étaient incompatibles avec la contrainte de branchement

    list<SCIP_VAR*>::const_iterator itv;

    for (itv = consdata->L_var_bound.begin(); itv!=consdata->L_var_bound.end(); itv++) {
        SCIPchgVarUbNode(scip, NULL, (*itv), SCIPinfinity(scip)) ;
    }
    consdata->L_var_bound.clear() ;
}
