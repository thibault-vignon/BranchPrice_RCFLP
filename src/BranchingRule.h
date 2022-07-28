#ifndef BranchingRuleH
#define BranchingRuleH

#include <scip/scip.h>
#include <objscip/objscip.h>

#include "Master.h"
#include "InstanceRCFLP.h"
#include "Pricer.h"

using namespace std;


/**
 * Manage the branching decisions
 **/
class BranchingRule : public scip::ObjBranchrule {
public:

    InstanceRCFLP* inst ;
    Master_Model* master ;
    ObjPricerRCFLP* pricer ; // NULL dans le cas de la décomposition par pas de temps avec pricing résolu par prog dyn

    const Parameters Param ;


    BranchingRule(SCIP* scip, InstanceRCFLP* i, Master_Model* m, ObjPricerRCFLP* p, const Parameters & Par) :
        scip::ObjBranchrule(scip, "BranchingOn(x,u)", "", 2000000, -1, 1.0), // properties of the branching rule (see doc)
        Param(Par)
    {
        master = m;
        inst = i;
        pricer = p ;
    }

    virtual ~BranchingRule(){}

    /*
     * Exec the branching rule
     *  Possible return values for *result (if more than one applies, the first in the list should be used):
     *  - SCIP_CUTOFF     : the current node was detected to be infeasible
     *  - SCIP_CONSADDED  : an additional constraint (e.g. a conflict clause) was generated; this result code must not be
     *                      returned, if allowaddcons is FALSE
     *  - SCIP_REDUCEDDOM : a domain was reduced that rendered the current LP solution infeasible
     *  - SCIP_SEPARATED  : a cutting plane was generated
     *  - SCIP_BRANCHED   : branching was applied
     *  - SCIP_DIDNOTRUN  : the branching rule was skipped
     */
    SCIP_RETCODE scip_execlp(
            SCIP*              scip,               /** SCIP data structure */
            SCIP_BRANCHRULE*   branchrule,         /** the branching rule itself */
            SCIP_Bool          allowaddcons,       /** should adding constraints be allowed to avoid a branching? */
            SCIP_RESULT*       result              /** pointer to store the result of the branching call */
            );

};


#endif
