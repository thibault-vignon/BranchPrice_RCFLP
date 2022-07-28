#ifndef CPLEXCHECK
#define CPLEXCHECK

#include <ilcplex/ilocplex.h>

#include <vector>
#include <list>

#include "Process.h"
#include "InstanceRCFLP.h"

using namespace std;


class CplexChecker {
 public:

  InstanceRCFLP* inst ;
  const Parameters Param ;
  IloEnv   env;
  IloModel model;
  IloCplex cplex;

  IloExpr cost;

  IloBoolVarArray y;
  IloBoolVarArray x;

  double PrimalBound;
  double DualBound ;  
  double PrimalBoundLowBound;
  double DualBoundLowBound ;
  double nbNodes;
  double nbNodesLowBound;
  double cpuTime ;
  double cpuTimeLowBound ;
  double gap ;


  double LRValue ;
  double LRCplexVal ;
  double valHeuristicCplex ;

  CplexChecker(InstanceRCFLP* inst, const Parameters & param) ;
  double getIntegerObjValue() ;
  double useLowBound(double lowbound) ;
  double getLRValue() ;
  double getLRCplex() ;
  void CplexPrimalHeuristic(IloNumArray solution_y, IloNumArray solution_x) ;

  double printSolution();
  void checkSolution(const vector<double> & y_frac, const vector<double> & x_frac);

};



#endif
