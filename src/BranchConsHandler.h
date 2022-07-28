#ifndef BranchConsHandlerH
#define BranchConsHandlerH

#include "scip/scip.h"
#include "scip/cons_linear.h"

#include <ilcplex/ilocplex.h>

#include <scip/scipdefplugins.h>
#include "objscip/objscip.h"
#include "scip/cons_linear.h"

#include "Master.h"
#include "Pricer.h"


using namespace std;


#define SCIP_DEBUG




/*****
 * Create an "artificial" constraint that will be included locally in the node for containg the node data
 * in particular, creates cplex constraint "BranchConstraint" for ConsData structure
 *****/
void createBranchCstr(SCIP* scip_, int VarX, int bound, int unit, int time, int site, ObjPricerRCFLP* pricer, SCIP_CONS** cons);


///////////////////////////:::

class BranchConsHandler : public scip::ObjConshdlr {

public :

    ObjPricerRCFLP *Pricer;
    Master_Model* Master ;
    InstanceRCFLP* inst ;

    BranchConsHandler(SCIP* scip, Master_Model* m,  ObjPricerRCFLP* p) :
        scip::ObjConshdlr(
            scip,
            "BranchConsHandler",                    // const char *  	name,
            "Handler For Branching Constraints",   // const char *  	desc,
            2000000, -2000000, -2000000,           // int sepapriority, int enfopriority, int checkpriority,
            1, -1, 1, 0,                           // int sepafreq, int propfreq, int eagerfreq, int maxprerounds,
            FALSE, FALSE, FALSE,                   // delaysepa, delayprop, needscons,
            SCIP_PROPTIMING_BEFORELP,              // SCIP_PROPTIMING  	proptiming,
            SCIP_PRESOLTIMING_FAST                 // SCIP_PRESOLTIMING  	presoltiming
            )

    {
        Pricer = p ;
        Master = m ;
        inst = p->inst ;
    }


    /**
     *  Activation d'un noeud (appelée à chaque fois qu'on rentre dans un noeud)
     **/
    virtual SCIP_RETCODE scip_active(
            SCIP * 	scip,
            SCIP_CONSHDLR * conshdlr,
            SCIP_CONS * cons
            );

    /**
     * Désactivation d'un noeud (appelée à chaque fois qu'on quitte un noeud
     * sauf si on le quitte pour aller dans un noeud fils)
     **/
    virtual SCIP_RETCODE scip_deactive(
            SCIP * 	scip,
            SCIP_CONSHDLR * 	conshdlr,
            SCIP_CONS * 	cons
            );



    //////////////////////////:
    //////////////////////////

    /** transforms constraint data into data belonging to the transformed problem */
    virtual SCIP_RETCODE scip_trans(
            SCIP*              scip,               //**< SCIP data structure *
            SCIP_CONSHDLR*     conshdlr,           //**< the constraint handler itself *
            SCIP_CONS*         sourcecons,         //**< source constraint to transform *
            SCIP_CONS**        targetcons          //**< pointer to store created target constraint *
            );


    virtual SCIP_RETCODE scip_check(
            SCIP*              scip,               /**< SCIP data structure */
            SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
            SCIP_CONS**        conss,              /**< array of constraints to process */
            int                nconss,             /**< number of constraints to process */
            SCIP_SOL*          sol,                /**< the solution to check feasibility for */
            SCIP_Bool          checkintegrality,   /**< has integrality to be checked? */
            SCIP_Bool          checklprows,        /**< have current LP rows to be checked? */
            SCIP_Bool          printreason,        /**< should the reason for the violation be printed? */
            SCIP_Bool          completely,         /**< should all violations be checked? */
            SCIP_RESULT*       result              /**< pointer to store the result of the feasibility checking call */
            );

    virtual SCIP_RETCODE scip_enfolp(
            SCIP*              scip,               /**< SCIP data structure */
            SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
            SCIP_CONS**        conss,              /**< array of constraints to process */
            int                nconss,             /**< number of constraints to process */
            int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
            SCIP_Bool          solinfeasible,      /**< was the solution already declared infeasible by a constraint handler? */
            SCIP_RESULT*       result              /**< pointer to store the result of the enforcing call */
            );

    virtual SCIP_RETCODE scip_enfops(
            SCIP*              scip,               /**< SCIP data structure */
            SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
            SCIP_CONS**        conss,              /**< array of constraints to process */
            int                nconss,             /**< number of constraints to process */
            int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
            SCIP_Bool          solinfeasible,      /**< was the solution already declared infeasible by a constraint handler? */
            SCIP_Bool          objinfeasible,      /**< is the solution infeasible anyway due to violating lower objective bound? */
            SCIP_RESULT*       result              /**< pointer to store the result of the enforcing call */
            );

    virtual SCIP_RETCODE scip_lock(
            SCIP*              scip,               /**< SCIP data structure */
            SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
            SCIP_CONS*         cons,               /**< the constraint that should lock rounding of its variables, or NULL if the
                                *   constraint handler does not need constraints */
            SCIP_LOCKTYPE      locktype,
            int                nlockspos,          /**< no. of times, the roundings should be locked for the constraint */
            int                nlocksneg           /**< no. of times, the roundings should be locked for the constraint's negation */
            );

    virtual SCIP_RETCODE scip_sepalp(
            SCIP*              scip,               /**< SCIP data structure */
            SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
            SCIP_CONS**        conss,              /**< array of constraints to process */
            int                nconss,             /**< number of constraints to process */
            int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
            SCIP_RESULT*       result              /**< pointer to store the result of the separation call */
            );

    virtual SCIP_RETCODE scip_sepasol(
            SCIP*              scip,               /**< SCIP data structure */
            SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
            SCIP_CONS**        conss,              /**< array of constraints to process */
            int                nconss,             /**< number of constraints to process */
            int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
            SCIP_SOL*          sol,                /**< primal solution that should be separated */
            SCIP_RESULT*       result              /**< pointer to store the result of the separation call */
            );
};



#endif
