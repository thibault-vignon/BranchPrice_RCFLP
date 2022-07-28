#include "BranchingRule.h"
#include "BranchConsHandler.h"

#define OUTPUT_BRANCHRULE


using namespace std;

SCIP_RETCODE BranchingRule::scip_execlp(SCIP* scip, SCIP_BRANCHRULE* branchrule, SCIP_Bool allowaddcons, SCIP_RESULT* result) {

#ifdef OUTPUT_BRANCHRULE
    cout << " --------------------- Branching Rule EXECLP ---------------  \n";
    cout << "Nombre de noeuds actuel : " << SCIPgetNNodes(scip) << std::endl;
#endif



//    SCIP_NODE* node = SCIPgetCurrentNode(scip);
//    SCIP_ConsData *consdata;
//    if (node->conssetchg!=NULL) {
//        consdata=SCIPconsGetData(node->conssetchg->addedconss[0]);
//        //SCIP_CONSSETCHG* SCIP_Node::conssetchg --> constraint set changes at this node or NULL
//        //typedef struct SCIP_ConsSetChg SCIP_CONSSETCHG --> tracks additions and removals of the set of active constraints
//#ifdef OUTPUT_BRANCHRULE
//        cout<<"Consdata non null"<<endl;
//#endif
//    }
//    else {
//        consdata=NULL;
//#ifdef OUTPUT_BRANCHRULE
//        cout<<"Consdata null"<<endl;
//#endif
//    }

    int T = master->T ;

    double eps = master->Param.Epsilon;

    // Search for the "most fractional" unit
    master->computeFracSol(scip);

   cout << "solution x frac: " << endl;

//    for (int t=0 ; t < T ; t++) {
//        for (int i=0 ; i <master->n ; i++) {
//            cout << master->x_frac[i*T+t] << " " ;
//        }
//        cout << endl ;
//    }

    SCIP_Real bestfrac = 1;
    SCIP_Real tmp;
    int unit = -1 ;
    int time = -1 ;

    for (int i=0 ; i < master->n ; i++) {
        for (int t=0 ; t < T ; t++) {
            tmp = master->x_frac[i*T+t] ;
            if ( (tmp > eps ) && (tmp < 1-eps) && (fabs(tmp - 0.5) < fabs(bestfrac - 0.5) ) ) {
                bestfrac = tmp;
                unit = i ;
                time = t ;
            }
        }
    }



    if (bestfrac < 1 - eps) {

       // int unit = floor(bestfrac) + inst->getFirstG(group) ;

        int VarX=1 ;
        int Site = Param.getSiteOf(unit);


#ifdef OUTPUT_BRANCHRULE
        cout<<"Branch on var x(" << unit <<", " << time << ") ";
        cout<<" of value : "<< master->x_frac.at(unit*T+time) <<endl;
#endif

        SCIP_NODE *newnode;
        SCIP_CONS *newcons;


        // first node
        SCIPcreateChild(scip, &newnode, 1000.0, SCIPgetLocalTransEstimate(scip));
        createBranchCstr(scip, VarX, 0, unit, time, Site, pricer, &newcons );
        SCIPaddConsNode(scip, newnode, newcons, NULL);
        SCIPreleaseCons(scip, &newcons);

        // second node
        SCIPcreateChild(scip, &newnode, 1000.0, SCIPgetLocalTransEstimate(scip));
        createBranchCstr(scip, VarX, 1, unit, time, Site, pricer, &newcons );
        SCIPaddConsNode(scip, newnode, newcons, NULL);
        SCIPreleaseCons(scip, &newcons);

        *result = SCIP_BRANCHED;
    }
    else{
        cout<<"Every variable is integer!!!!!!"<<endl;
        *result = SCIP_CUTOFF;
    }

#ifdef DEBUG
    cout << "\n*****END OF Branching Rule EXECLP ****\n";
    cout << "****************************************\n";
#endif
    
    return SCIP_OKAY;

}
