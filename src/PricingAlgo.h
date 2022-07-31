#ifndef CPLEXPRICER
#define CPLEXPRICER

#include <ilcplex/ilocplex.h>

#include <vector>
#include <list>

#include "InstanceRCFLP.h"
#include "Master.h"
#include "Process.h"

using namespace std;
using namespace scip;

class DualCostsCustomer {
public:

  vector<double> Sigma ; // contraintes convexité

  vector<double> Omega1 ; // contraintes d'égalité customer/facility en x
  vector<double> Omega2 ; // contraintes d'égalité customer/facility en y

  DualCostsCustomer(InstanceRCFLP* inst) ;
};


class DualCostsFacility {
public:

  vector<double> Sigma ; // contraintes convexité
    
  vector<double> Omega1 ; // contraintes d'égalités customer/facility en x
  vector<double> Omega2 ; // contraintes d'égalités customer/facility en y

  DualCostsFacility(InstanceRCFLP* inst) ;
};


class CplexPricingAlgoFacility {
 public:

  const Parameters Param ;
  int facility;
  IloEnv env;
  IloModel model;
  IloObjective obj;
  IloCplex cplex;

  IloBoolVarArray x;
  IloBoolVarArray y;

  vector<double> BaseObjCoefX ;
  vector<double> BaseObjCoefY ;

  double cpuTime ;

  CplexPricingAlgoFacility(InstanceRCFLP* inst, const Parameters & Param, int f);

  void updateObjCoefficients(InstanceRCFLP* inst, const Parameters & Param, const DualCostsFacility & Dual, bool Farkas);

  // Launch Cplex solver and get back an optimal solution
  bool findImprovingSolution(InstanceRCFLP* inst, const DualCostsFacility & Dual, double& objvalue) ;
  // returns true if an improving solution has been found. objvalue is updated in this case

  void getSolution(InstanceRCFLP* inst, const DualCostsFacility & Dual, IloNumArray xPlan, IloNumArray yPlan, bool Farkas) ; //updates plans and realCost

};


class CplexPricingAlgoCustomer {
 public:

  const Parameters Param ;
  int customer ;
  IloEnv   env;
  IloModel model;
  IloObjective obj;
  IloCplex cplex;

  IloBoolVarArray x;
  IloBoolVarArray y ;

  double cpuTime ;

  vector<double> BaseObjCoefX ;
  vector<double> BaseObjCoefY ;

  CplexPricingAlgoCustomer(InstanceRCFLP* inst, const Parameters & par, int c);

  void updateObjCoefficients(InstanceRCFLP* inst, const Parameters & Param, const DualCostsCustomer & Dual, bool Farkas);

  // Launch Cplex solver and get back an optimal up/down plan
  bool findImprovingSolution(InstanceRCFLP* inst, const DualCostsCustomer & Dual, double& objvalue) ;
  // returns true if an improving solution has been found. objvalue is updated in this case

  void getSolution(InstanceRCFLP* inst, const DualCostsCustomer & Dual, IloNumArray xPlan, IloNumArray yPlan, bool Farkas) ; //updates plans and realCost
};


class DynProgPricingAlgoFacility {
 public:

  Parameters Param ;
  int facility ;
  IloEnv env;

  vector<int> init_x ; // init[i]==0 si i fixé à 0, init[i]=1 si i fixé à 1, si ça reste à déterminer: -1
  vector<int> init_y ;

  vector<double> Table ;

  vector<double> BaseObjCoefX ;
  vector<double> ObjCoefX ;

  DynProgPricingAlgoFacility(InstanceRCFLP* inst, const Parameters & par, int f); // initialise les vecteurs

  void updateObjCoefficients(InstanceRCFLP* inst, const Parameters & Param, const DualCostsFacility & Dual, bool Farkas);

  bool findImprovingSolution(InstanceRCFLP* inst, const DualCostsFacility & Dual, double& objvalue) ;
  // returns true if an improving solution has been found. objvalue is updated in this case

  void getSolution(InstanceRCFLP* inst, const DualCostsFacility & Dual, IloNumArray xPlan, IloNumArray yPlan, bool Farkas) ;//updates plans and realCost
};


class DynProgPricingAlgoCustomer { // codé dans le cas Pmin=Pmax pour voir si c'est intéressant. Cas où D et les puissances sont entiers
 public:

  Parameters Param ;
  int customer ;
  IloEnv env;

  vector<int> init_x ; // init[i]==0 si i fixé à 0, init[i]=1 si i fixé à 1, si ça reste à déterminer: -1
  vector<int> init_y ;

  vector<double> Table ;

  vector<double> BaseObjCoefX ;
  vector<double> ObjCoefX ;

  DynProgPricingAlgoCustomer(InstanceRCFLP* inst, const Parameters & par, int c);

  void updateObjCoefficients(InstanceRCFLP* inst, const Parameters & Param, const DualCostsCustomer & Dual, bool Farkas);

  bool findImprovingSolution(InstanceRCFLP* inst, const DualCostsCustomer & Dual, double& objvalue) ;
  // computes Bellman table (vector Table)
  // returns true if an improving solution has been found. objvalue is updated in this case

  void getSolution(InstanceRCFLP* inst, const DualCostsCustomer & Dual, IloNumArray xPlan, IloNumArray yPlan, bool Farkas) ; //updates plans and realCost
};


#endif
