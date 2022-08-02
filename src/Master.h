#ifndef MASTER
#define MASTER

#include <ilcplex/ilocplex.h>

#include <vector>
#include <list>

/* scip includes */
#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"

/* user defined includes */
#include "InstanceRCFLP.h"
#include "Process.h"


/* namespace usage */
using namespace std;
using namespace scip;



/////////////////////////////////////
////////// BRANCHING CONS DATA //////
/////////////////////////////////////


// Data associated to a constraint (each artificial constraint represents one branching constraint)
struct SCIP_ConsData {
    int VarX ; // =1 si la variable branchée est un x. Si branchement sur y: VarX=0
    int bound ; // variable fixée à 1 ou à 0
    int facility; 
    int customer ;
    list<SCIP_VAR*> L_var_bound;
};


////////////////////////////////////////////
////////// MASTER MODEL (virtual) //////////
////////////////////////////////////////////

class Master_Model {
public:
    int I ;
    int J ;

    const Parameters Param ;
    InstanceRCFLP* inst ;

    //solution fractionnaire et valeurs duales associes
    list<double> totalDualCostList;
    vector<double> x_frac ;
    vector<double> y_frac ;

    //compteur du nb d'appels au pricer
    double cumul_resolution_pricing ;

    int nbIter=0;

    Master_Model(const Parameters & Par, InstanceRCFLP* i) : Param(Par), inst(i) {
        I = inst->getI();
        J = inst->getJ() ;
        x_frac.resize(I*J,0);
        y_frac.resize(J,0);
    }
    virtual void computeFracSol(SCIP* scip) = 0;  // = 0 signifie "virtuelle pure"

    virtual void discardVar(SCIP* scip, SCIP_ConsData* consdata) = 0;
    virtual void restoreVar(SCIP* scip, SCIP_ConsData* consdata) = 0;

    virtual ~Master_Model() {}
};

///////////////////////////////
////////// VARIABLES //////////
///////////////////////////////


class MasterFacility_Variable{
public:

    /// Keep a pointer on every variable of the Master program
    SCIP_VAR* ptr;

    /// facility corresponding to variable ptr
    int facility ; 

    //// customer assignment plan corresponding to ptr
    IloNumArray x_plan ;

    //// facility opening plan corresponding to ptr
    IloNumArray y_plan ;

    double cost ;

    MasterFacility_Variable(int facility, IloNumArray x_plan, IloNumArray y_plan) ;
    void computeCost(InstanceRCFLP* inst, const Parameters & Param) ;

};



class MasterCustomer_Variable{
public:

    /// Keep a pointer on every variable of the Master program
    SCIP_VAR* ptr;

    /// customer corresponding to variable ptr
    int customer ;

    //// customer assignment plan corresponding to ptr
    IloNumArray x_plan ;

    //// facility opening plan corresponding to ptr
    IloNumArray y_plan ;

    double cost ;

    MasterCustomer_Variable(int customer, IloNumArray x_plan, IloNumArray y_plan) ;
    void computeCost(InstanceRCFLP* inst, const Parameters & Param) ;

};


////////////////////////////////////////////////////
////////// MASTER -- DECOMPOSITION PAR FACILITY ///////
////////////////////////////////////////////////////


class MasterFacility_Model : public Master_Model {
public:

    IloEnv env;

    // Keep a pointer on every constraint of the Master program
    vector<SCIP_CONS*> convexity_cstr;
    vector<SCIP_CONS*> assignment_cstr;
    vector<SCIP_CONS*> reliability_cstr;

    // Keep info on every variables of the Master program
    //NB: le fait d'utiliser une liste ne permet pas de supprimer des variables
    list<MasterFacility_Variable*> L_var;

    MasterFacility_Model(InstanceRCFLP* inst, const Parameters & Param) ;

    void addCoefsToConstraints(SCIP* scip, MasterFacility_Variable* lambda) ;
    void initScipMasterFacilityModel(SCIP* scip);
    void initMasterFacilityVariable(SCIP* scip, MasterFacility_Variable* lambda) ;
    void createColumns(SCIP* scip, IloNumArray x, IloNumArray y) ;

    void computeFracSol(SCIP* scip) ;

    void discardVar(SCIP* scip, SCIP_ConsData* consdata) ;
    void restoreVar(SCIP* scip, SCIP_ConsData* consdata) ;
};


////////////////////////////////////////////////////
////////// DECOMPOSITION PAR CLIENT //////////
////////////////////////////////////////////////////


class MasterCustomer_Model : public Master_Model {
public:

    IloEnv env;

    // Keep a pointer on every constraint of the MasterCustomer program 
    vector<SCIP_CONS*> convexity_cstr;
    vector<SCIP_CONS*> capacity_cstr;

    // Keep informations on every variables of the Master program
    //NB: le fait d'utiliser une liste ne permet pas de supprimer des variables
    list<MasterCustomer_Variable*> L_var;

    MasterCustomer_Model(InstanceRCFLP* inst, const Parameters & Param) ;

    void addCoefsToConstraints(SCIP* scip, MasterCustomer_Variable* lambda) ;

    void initScipMasterCustomerModel(SCIP* scip);
    void initMasterCustomerVariable(SCIP* scip, MasterCustomer_Variable* lambda) ;

    void createColumns(SCIP* scip, IloNumArray x, IloNumArray p) ;

    void computeFracSol(SCIP* scip) ;

    void discardVar(SCIP* scip, SCIP_ConsData* consdata) ;
    void restoreVar(SCIP* scip, SCIP_ConsData* consdata) ;

};


////////////////////////////////////////////////////
////////// MASTER -- DOUBLE DECOMPOSITION ///////
////////////////////////////////////////////////////


class MasterDouble_Model : public Master_Model {
public:

    IloEnv env;

    // Keep a pointer on every constraint of the Master program
    // convexity constraints
    vector<SCIP_CONS*> conv_lambda_facility;
    vector<SCIP_CONS*> conv_lambda_customer;
    // (in)equality between Customer and Facility solutions
    vector<SCIP_CONS*> eq_customer_facility_x;
    vector<SCIP_CONS*> eq_customer_facility_y;

    // Keep info on every variables of the Master program
    //NB: using a list does not allow to delete variables
    list<MasterFacility_Variable*> L_var_facility;
    list<MasterCustomer_Variable*> L_var_customer;

    MasterDouble_Model(InstanceRCFLP* inst, const Parameters & Param) ;

    void addCoefsToConstraints_facilityVar(SCIP* scip, MasterFacility_Variable* lambda) ;

    void addCoefsToConstraints_customerVar(SCIP* scip, MasterCustomer_Variable* lambda) ;

    void initScipMasterDoubleModel(SCIP* scip);
    void initMasterFacilityVariable(SCIP* scip, MasterFacility_Variable* lambda) ;
    void initMasterCustomerVariable(SCIP* scip, MasterCustomer_Variable* lambda) ;
    void createColumns(SCIP* scip, IloNumArray x, IloNumArray p) ;

    void computeFracSol(SCIP* scip) ;

    void discardVar(SCIP* scip, SCIP_ConsData* consdata) ;
    void restoreVar(SCIP* scip, SCIP_ConsData* consdata) ;
};
#endif /* MASTER INCLUDED */
