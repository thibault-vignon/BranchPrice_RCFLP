#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "scip/dialog_default.h"

#include <ilcplex/ilocplex.h>


#include <iomanip>
/* scip includes */
#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"

/* user defined includes */
#include "InstanceRCFLP.h"
#include "Process.h"
#include "Master.h"
#include "Pricer.h"
#include "BranchConsHandler.h"
#include "BranchingRule.h"
#include "Checker.h"

/* namespace usage */
using namespace std;
using namespace scip;

//#define SCIP_OUTPUT

#define SCIP_DEBUG

int main(int argc, char** argv)
{
    ofstream fichier("result.txt", std::ofstream::out | std::ofstream::app);

    //////////////////////////////
    //////  INSTANCE DATA    /////
    //////////////////////////////


    int v ;
    string instanceSet;
    string id;
    int met = 0 ;

    v = atoi(argv[1]);
    instanceSet = argv[2];
    id = argv[3];
    met = atoi(argv[4]);

    InstanceProcessed Instance = InstanceProcessed(v, instanceSet, id) ;

    string nom = Instance.fileName() ;
    cout << nom << endl;
    const char* file = strdup(nom.c_str()) ;

    IloEnv env;
    InstanceRCFLP* inst = new InstanceRCFLP(env, file) ;

    ///////////////////////////
    //////  PARAMETERS    /////
    ///////////////////////////

    Parameters param = init_parameters(inst, met);
    cout << "met = " << met << endl;

    ////////////////////////////////////
    //////  SCIP INITIALIZATION    /////
    ////////////////////////////////////

    clock_t start;
    start = clock();

    // problem initialization
    SCIP *scip=NULL;
    SCIPcreate(&scip);

    SCIPprintVersion(scip, NULL);
    SCIPinfoMessage(scip, NULL, "\n");

    /* include default plugins */
   // SCIPincludeDefaultPlugins(scip);


    // include various SCIP features
    SCIPincludeConshdlrLinear(scip);
    SCIPincludeNodeselBfs(scip);
    SCIPincludeConshdlrIntegral(scip);
    SCIPincludeDispDefault(scip);
    //SCIPincludeDialogDefault(scip);
    SCIPincludeHeurActconsdiving(scip);
    SCIPincludeHeurClique(scip);
    SCIPincludeHeurCoefdiving(scip);
    SCIPincludeHeurCrossover(scip);
    SCIPincludeHeurDins(scip);
    SCIPincludeHeurFeaspump(scip);
    SCIPincludeHeurFixandinfer(scip);
    SCIPincludeHeurFracdiving(scip);
    SCIPincludeHeurGuideddiving(scip);
    SCIPincludeHeurIntdiving(scip);
    SCIPincludeHeurIntshifting(scip);
    SCIPincludeHeurLinesearchdiving(scip);
    SCIPincludeHeurLocalbranching(scip);
    SCIPincludeHeurMutation(scip);
    SCIPincludeHeurObjpscostdiving(scip);
    SCIPincludeHeurOctane(scip);
    SCIPincludeHeurOneopt(scip);
    SCIPincludeHeurPscostdiving(scip);
    SCIPincludeHeurRens(scip);
    SCIPincludeHeurRins(scip);
    SCIPincludeHeurShiftandpropagate(scip);
    SCIPincludeHeurShifting(scip);
    SCIPincludeHeurSimplerounding(scip);
    SCIPincludeHeurSubNlp(scip);
    SCIPincludeHeurTrivial(scip);
    SCIPincludeHeurTrySol(scip);
    SCIPincludeHeurTwoopt(scip);
    SCIPincludeHeurUndercover(scip);
    SCIPincludeHeurVbounds(scip);
    SCIPincludeHeurVeclendiving(scip);
    SCIPincludeHeurZirounding(scip);
    SCIPincludeHeurRootsoldiving(scip);
    SCIPincludeHeurRounding(scip);


    // resolution parameters (time and node limits)
    SCIPsetLongintParam(scip, "limits/nodes", param.nodeLimit);
    SCIPsetRealParam(scip, "limits/time", 3600);

    // Changements pour rendre le code compatible avec scip-8.0.0
    SCIPincludeDialogDefaultBasic(scip) ;
    SCIPincludeDialogDefaultSet(scip) ;
    SCIPincludeDialogDefaultFix(scip) ;

    /* set verbosity parameter */
    SCIPsetIntParam(scip, "display/verblevel", 5);
    //SCIPsetBoolParam(scip, "display/lpinfo", TRUE);

    /* create empty problem */
    SCIPcreateProb(scip, "UCP", 0, 0, 0, 0, 0, 0, 0);

    ////////////////////////////////////////////////////////
    //////  MASTER & PRICER PROBLEMS INITIALIZATION    /////
    ////////////////////////////////////////////////////////






    // frontal resolution with Cplex to verify results
    CplexChecker checker = CplexChecker(inst, param) ;
    
    cout << checker.getIntegerObjValue() << endl;

}

