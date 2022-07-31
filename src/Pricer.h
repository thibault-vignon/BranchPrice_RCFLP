#ifndef PRICER
#define PRICER

#include "objscip/objscip.h"
#include "scip/pub_var.h"

#include <vector>
#include <list>

#include "InstanceRCFLP.h"
#include "Master.h"
#include "PricingAlgo.h"

using namespace std;
using namespace scip;




/** pricer class */
// hérité d'une classe de SCIP
class ObjPricerRCFLP : public ObjPricer {
public:
    InstanceRCFLP* inst ;
    Parameters Param ;

    int iteration;

    int facilityColumns;
    int customerColumns;

    double totalDualCost;

    ObjPricerRCFLP(
            SCIP*                               scip,
            const char*                         pp_name,
            InstanceRCFLP*                        i,
            const Parameters &                  p
            ) :
        ObjPricer(scip, pp_name, "Find plans with negative reduced costs for each subproblem", 0, TRUE),
        inst(i),
        Param(p)
    {
        iteration=0;
        facilityColumns=0;
        customerColumns=0;
        infeasibilityDetected = false ;
    }

    virtual void addVarBound(SCIP_ConsData* consdata) = 0;
    virtual void removeVarBound(SCIP_ConsData* consdata) = 0;

//    /** initialization method of variable pricer (called after problem was transformed) */
//    virtual SCIP_DECL_PRICERINIT(scip_init) {};

//    /** reduced cost pricing method of variable pricer for feasible LPs */
//   // virtual SCIP_DECL_PRICERREDCOST(scip_redcost);
//    virtual  SCIP_RETCODE scip_redcost(SCIP* scip, SCIP_PRICER* pricer, SCIP_Real* lowerbound, SCIP_Bool* stopearly, SCIP_RESULT* result) override {}


//    // recherche d'une variable pour la faisabilité du PMR et insertion si trouvée
//    virtual  SCIP_RETCODE scip_farkas(SCIP* scip, SCIP_PRICER* pricer, SCIP_RESULT* result) override {}
    virtual ~ObjPricerRCFLP() {}

    double currentLowerBound ;

    list<MasterFacility_Variable*> facilityVarsToAdd ;
    list<MasterCustomer_Variable*> customerVarsToAdd ;

    bool infeasibilityDetected ;
};


class ObjPricerFacility : public ObjPricerRCFLP {
public:


    MasterFacility_Model* Master ;
    vector<CplexPricingAlgoFacility*> AlgoCplex;
    vector<DynProgPricingAlgoFacility*> AlgoDynProg;

   /** Constructs the pricer object with the data needed */
   ObjPricerFacility(
      SCIP*                               scip,        /**< SCIP pointer */
      const char*                         p_name,       /**< name of pricer */
      MasterFacility_Model*                       M,
      InstanceRCFLP*                        inst,
      const Parameters &                  param
      );

   /** Destructs the pricer object. */
   virtual ~ObjPricerFacility();

   /** initialization method of variable pricer (called after problem was transformed) */
   virtual SCIP_DECL_PRICERINIT(scip_init);

   /** reduced cost pricing method of variable pricer for feasible LPs */
  // virtual SCIP_DECL_PRICERREDCOST(scip_redcost);
   virtual  SCIP_RETCODE scip_redcost(SCIP* scip, SCIP_PRICER* pricer, SCIP_Real* lowerbound, SCIP_Bool* stopearly, SCIP_RESULT* result) override;


   // recherche d'une variable pour la faisabilité du PMR et insertion si trouvée
   virtual  SCIP_RETCODE scip_farkas(SCIP* scip, SCIP_PRICER* pricer, SCIP_RESULT* result) override;

   /*put dual costs or farkas cost in vector dual_cost*/
   void updateDualCosts(SCIP* scip, DualCostsFacility & dual_cost, bool Farkas) ;

   void addVarBound(SCIP_ConsData* consdata) ;
   void removeVarBound(SCIP_ConsData* consdata) ;

   /** performs pricing */
   void pricingRCFLP(
      SCIP*              scip,               /**< SCIP data structure */
      bool               Farkas
      );

};


////////////////////////////////////////////////////
////////// DECOMPOSITION PAR PAS DE TEMPS //////////
////////////////////////////////////////////////////

class ObjPricerCustomerRCFLP : public ObjPricerRCFLP {
public:

    MasterCustomer_Model* Master ;
    vector<CplexPricingAlgoCustomer*> AlgoCplex;
    vector<DynProgPricingAlgoCustomer*> AlgoDynProg;

   /** Constructs the pricer object with the data needed */
   ObjPricerCustomerRCFLP(
      SCIP*                               scip,        /**< SCIP pointer */
      const char*                         p_name,       /**< name of pricer */
      MasterCustomer_Model*                       M,
      InstanceRCFLP*                        inst,
      const Parameters &                  param
      );

   /** Destructs the pricer object. */
   virtual ~ObjPricerCustomerRCFLP();

   /** initialization method of variable pricer (called after problem was transformed) */
   virtual SCIP_DECL_PRICERINIT(scip_init);

   /** reduced cost pricing method of variable pricer for feasible LPs */
  // virtual SCIP_DECL_PRICERREDCOST(scip_redcost);
   virtual  SCIP_RETCODE scip_redcost(SCIP* scip, SCIP_PRICER* pricer, SCIP_Real* lowerbound, SCIP_Bool* stopearly, SCIP_RESULT* result) override;


   // recherche d'une variable pour la faisabilité du PMR et insertion si trouvée
   virtual  SCIP_RETCODE scip_farkas(SCIP* scip, SCIP_PRICER* pricer, SCIP_RESULT* result) override;

   /*put dual costs or farkas cost in vector dual_cost*/
   void updateDualCosts(SCIP* scip, DualCostsCustomer & dual_cost, bool Farkas) ;

   void addVarBound(SCIP_ConsData* consdata) ;
   void removeVarBound(SCIP_ConsData* consdata) ;

   /** performs pricing */
   void pricingRCFLP(
      SCIP*              scip,               /**< SCIP data structure */
      bool               Farkas
      );

};



//////////////////////////////////////////
////////// DOUBLE DECOMPOSITION //////////
//////////////////////////////////////////

class ObjPricerDouble : public ObjPricerRCFLP {
public:


    MasterDouble_Model* Master ;

    vector<CplexPricingAlgoFacility*> AlgoCplex_facility;
    vector<DynProgPricingAlgoFacility*> AlgoDynProg_facility;

    vector<CplexPricingAlgoCustomer*> AlgoCplex_customer;
    vector<DynProgPricingAlgoCustomer*> AlgoDynProg_customer;

   /** Constructs the pricer object with the data needed */
   ObjPricerDouble(
      SCIP*                               scip,        /**< SCIP pointer */
      const char*                         p_name,       /**< name of pricer */
      MasterDouble_Model*                       M,
      InstanceRCFLP*                        inst,
      const Parameters &                  param
      );

   /** Destructs the pricer object. */
   virtual ~ObjPricerDouble();

   /** initialization method of variable pricer (called after problem was transformed) */
   virtual SCIP_DECL_PRICERINIT(scip_init);

   /** reduced cost pricing method of variable pricer for feasible LPs */
  // virtual SCIP_DECL_PRICERREDCOST(scip_redcost);
   virtual  SCIP_RETCODE scip_redcost(SCIP* scip, SCIP_PRICER* pricer, SCIP_Real* lowerbound, SCIP_Bool* stopearly, SCIP_RESULT* result) override;

   // recherche d'une variable pour la faisabilité du PMR et insertion si trouvée
   virtual  SCIP_RETCODE scip_farkas(SCIP* scip, SCIP_PRICER* pricer, SCIP_RESULT* result) override;

   /*put dual costs or farkas cost in vector dual_cost*/
   void updateDualCosts_facility(SCIP* scip, DualCostsFacility & dual_cost, bool Farkas) ;
   void updateDualCosts_customer(SCIP* scip, DualCostsCustomer & dual_cost, bool Farkas) ;

   void addVarBound(SCIP_ConsData* consdata) ;
   void removeVarBound(SCIP_ConsData* consdata) ;

   /** performs pricing */
   void pricingRCFLP(
      SCIP*              scip,               /**< SCIP data structure */
      bool               Farkas
      );

};


#endif
