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

    Parameters Param = init_parameters(inst, met);
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
    SCIPsetLongintParam(scip, "limits/nodes", Param.nodeLimit);
    SCIPsetRealParam(scip, "limits/time", 3600);

    // Changements pour rendre le code compatible avec scip-8.0.0
    SCIPincludeDialogDefaultBasic(scip) ;
    SCIPincludeDialogDefaultSet(scip) ;
    SCIPincludeDialogDefaultFix(scip) ;

    /* set verbosity parameter */
    SCIPsetIntParam(scip, "display/verblevel", 5);
    //SCIPsetBoolParam(scip, "display/lpinfo", TRUE);

    /* create empty problem */
    SCIPcreateProb(scip, "RCFLP", 0, 0, 0, 0, 0, 0, 0);

    ////////////////////////////////////////////////////////
    //////  MASTER & PRICER PROBLEMS INITIALIZATION    /////
    ////////////////////////////////////////////////////////


    // frontal resolution with Cplex to verify results
    CplexChecker checker = CplexChecker(inst, Param) ;

    // create master and pricing classes
    Master_Model* Master_ptr;
    ObjPricerRCFLP* Pricer = NULL ;

    static const char* PRICER_NAME = "Pricer_RCFLP";


    if (Param.doubleDecompo){

        Master_ptr = new MasterDouble_Model(inst, Param) ;
        MasterDouble_Model* MD ;
        MD = dynamic_cast<MasterDouble_Model*> (Master_ptr) ;

        if (MD != NULL) {

            ///Initialisation du master
            MD->initScipMasterDoubleModel(scip);

            /// Initialisation du pricer
            Pricer = new ObjPricerDouble(scip, PRICER_NAME, MD, inst, Param);
            SCIPincludeObjPricer(scip, Pricer, true);
            SCIPactivatePricer(scip, SCIPfindPricer(scip, PRICER_NAME));

        }

    }

    else if (Param.FacilityDecompo){

        Master_ptr = new MasterFacility_Model(inst, Param) ;
        MasterFacility_Model* MF ;
        MF = dynamic_cast<MasterFacility_Model*> (Master_ptr) ;

        if (MF != NULL) {

            ///Initialisation du master
            MF->initScipMasterFacilityModel(scip);

            /// Initialisation du pricer
            Pricer = new ObjPricerFacility(scip, PRICER_NAME, MF, inst, Param);
            SCIPincludeObjPricer(scip, Pricer, true);
            SCIPactivatePricer(scip, SCIPfindPricer(scip, PRICER_NAME));

        }
    }

    else {

        Master_ptr = new MasterCustomer_Model(inst, Param) ;
        MasterCustomer_Model* MC ;
        MC = dynamic_cast<MasterCustomer_Model*> (Master_ptr) ;

        if (MC != NULL) {

            ///Initialisation du master
            MC->initScipMasterCustomerModel(scip);

            /// Initialisation du pricer
            Pricer = new ObjPricerCustomer(scip, PRICER_NAME, MC, inst, Param);
            SCIPincludeObjPricer(scip, Pricer, true);
            SCIPactivatePricer(scip, SCIPfindPricer(scip, PRICER_NAME));

        }

    }



    if (Param.ColumnGeneration) {

        cout << "resolution..." << endl ;
        SCIPsolve(scip);

        fichier << met << " & " << inst->getJ() << " & " << inst->getI() << " & " << inst->getv()  << " & " << inst->getK() << " & " << " - ";

        fichier << " & " << Pricer->iteration ;
        fichier << " & " << SCIPgetNPricevarsFound(scip) ;
        fichier << " & " << Pricer->customerColumns ; 
        fichier << " & " << Pricer->facilityColumns ;

        double timeScip =  SCIPgetSolvingTime(scip) ;
        fichier << " &  " << timeScip;
        SCIP_PRICER ** scippricer = SCIPgetPricers(scip);
        fichier << " & " <<  timeScip - SCIPpricerGetTime(scippricer[0]);

        fichier << " &  " << SCIPgetGap(scip);
        fichier << " &  " << SCIPgetDualbound(scip) ;
        fichier << " &  " << SCIPgetPrimalbound(scip) ;
        fichier << " & " << checker.getLRValue() ; 
        fichier << " & " << checker.getLRCplex() ;
        fichier << " & " << checker.getIntegerObjValue() ;
        fichier << " & " << checker.cpuTime ;
    }


    else {

    }

    fichier <<" \\\\ " << endl ;

    Master_ptr->computeFracSol(scip);
    for (int j=0 ; j < inst->getJ() ; j++) {
        cout << "j=" << j << " : " ;
        for (int i=0 ; i < inst->getI() ; i++) {
            cout << Master_ptr->x_frac[j*inst->getI() + i] ;
            cout << "  " ;
        }
        cout << endl;
        cout << "y : " << Master_ptr->y_frac[j] << endl;
    }
    checker.checkSolution(Master_ptr->y_frac, Master_ptr->x_frac);

}

