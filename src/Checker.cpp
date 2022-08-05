#include "Checker.h"
#include <iostream>

using namespace std;

CplexChecker::CplexChecker(InstanceRCFLP* instance, const Parameters & param) : Param(param) {

    inst = instance;

    int I = inst->getI();
    int J = inst->getJ() ;

    model = IloModel(env) ;

    valHeuristicCplex = -1 ;

    y = IloBoolVarArray(env, J) ;
    x = IloBoolVarArray(env, J*I) ;

    cost = IloExpr(env) ;

    // Objective Function: Minimize Cost

    for (int j=0 ; j < J ; j++) {
        cost += inst->getb(j) * y[j];
        for (int i=0; i<I; i++) {
            cost += inst->geta(j,i) * x[j*I + i];
        }
    }

    model.add(IloMinimize(env, cost));


    // Constraints

    // Assignment of each customer to a facility
    for (int i=0; i<I; i++) {
        IloExpr sum(env) ;
        for (int j=0 ; j < J ; j++) {
            sum += x[j*I + i];
        }
        model.add(sum == 1) ;
        sum.end() ;
    }

    // Capacity constraints
    for (int j=0 ; j < J ; j++) {
        IloExpr sum(env) ;
        for (int i=0; i<I; i++) {
            sum += x[j*I + i] * inst->getd(i);
        }
        model.add(sum <= y[j] * inst->getc(j)) ;
        sum.end() ;
    }

    // Reliability constraints
    for (int i=0; i<I; i++) {
        for (int j=0 ; j < J ; j++) {
            IloExpr sum(env) ;
            for (int indice=0; indice<inst->getv(); indice++){
                sum += y[inst->getV(j)[indice]] * inst->getc(inst->getV(j)[indice]);
            }
            model.add(sum >= x[j*I + i] * inst->getd(i) * inst->getK()) ;
        }
    }
}

double CplexChecker::getIntegerObjValue() {

    int print=1 ;


    clock_t start;
    start = clock();

    IloModel IntegerModel(env) ;
    IntegerModel.add(model) ;

    IloCplex IntegerObjCplex = IloCplex(IntegerModel) ; // ou juste valeur opt entière
    IntegerObjCplex.setParam(IloCplex::EpGap, Param.Epsilon) ;
    IntegerObjCplex.setParam(IloCplex::Param::ClockType, 1); //1 : CPU TIME
    //IntegerObjCplex.setParam(IloCplex::Param::TimeLimit, 30) ;


    IntegerObjCplex.solve() ;

    DualBound = IntegerObjCplex.getBestObjValue() ;
    PrimalBound = IntegerObjCplex.getObjValue() ;
    nbNodes = IntegerObjCplex.getNnodes() ;
    cpuTime = IntegerObjCplex.getCplexTime();
    gap = IntegerObjCplex.getMIPRelativeGap() ;

    cpuTime =  ( clock() - start ) / (double) CLOCKS_PER_SEC;


    int I = inst->getI();
    int J = inst->getJ() ;

    if (print) {

        IloNumArray solution_x = IloNumArray(env, J*I) ;
        IntegerObjCplex.getValues(solution_x, x) ;

        IloNumArray solution_y = IloNumArray(env, J) ;
        IntegerObjCplex.getValues(solution_y, y) ;

        cout.precision(6);
        cout << "X: " << endl ;
        for (int j=0 ; j < J ; j++) {
            cout << "j = " << j << " : " ;
            for (int i=0 ; i < I ; i++) {
                cout << fabs(solution_x[j*I + i]) << " " ;
            }
            cout << endl ;
        }
        cout << endl ;

        cout << "Y: " << endl ;
        for (int j=0 ; j < J ; j++) {
            cout << "j = " << j << " : " ;
            cout << fabs(solution_y[j]) << " " ;
            cout << endl ;
        }
        cout << endl ;
    }

    return PrimalBound;
}

double CplexChecker::useLowBound(double lowbound) {

    int print=0 ;


    clock_t start;
    start = clock();

    IloModel IntegerModel(env) ;
    IntegerModel.add(model) ;
    IntegerModel.add(cost >= lowbound);

    IloCplex useLowBoundCplex = IloCplex(IntegerModel) ; // ou juste valeur opt entière
    useLowBoundCplex.setParam(IloCplex::EpGap, Param.Epsilon) ;
    useLowBoundCplex.setParam(IloCplex::Param::ClockType, 1); //1 : CPU TIME
    //useLowBoundCplex.setParam(IloCplex::Param::TimeLimit, 30) ;



    useLowBoundCplex.solve() ;

    nbNodesLowBound = useLowBoundCplex.getNnodes() ;
    cpuTimeLowBound = useLowBoundCplex.getCplexTime();
    DualBoundLowBound = useLowBoundCplex.getBestObjValue() ;
    PrimalBoundLowBound = useLowBoundCplex.getObjValue() ;

    cpuTimeLowBound =  ( clock() - start ) / (double) CLOCKS_PER_SEC;

    return cpuTimeLowBound;
}

double CplexChecker::getLRValue() {

    //Modèle
    IloModel LRModel(env) ;
    LRModel.add(model) ;
    LRModel.add(IloConversion(env, x, IloNumVar::Float) ) ;
    LRModel.add(IloConversion(env, y, IloNumVar::Float) ) ;

    //Résolution
    IloCplex LRVal = IloCplex(LRModel) ;
    LRVal.setParam(IloCplex::EpGap, 0) ;
    LRVal.solve() ;

    LRValue = LRVal.getObjValue();
    return LRValue;
}

double CplexChecker::getLRCplex() {

    //Modèle
    IloModel LRCplexModel(env) ;
    LRCplexModel.add(model) ;

    //Résolution
    IloCplex LRCplex = IloCplex(LRCplexModel) ;
    LRCplex.setParam(IloCplex::EpGap, 0) ;
    LRCplex.setParam(IloCplex::Param::MIP::Limits::Nodes, 1) ;
    LRCplex.solve() ;

    LRCplexVal = LRCplex.getBestObjValue() ;
    return LRCplexVal ;
}

void CplexChecker::CplexPrimalHeuristic(IloNumArray solution_y, IloNumArray solution_x) {

    //Modèle
    IloModel LRCplexModel(env) ;
    LRCplexModel.add(model) ;

    //Résolution
    IloCplex LRCplex = IloCplex(LRCplexModel) ;
    LRCplex.setParam(IloCplex::EpGap, 0) ;
    LRCplex.setParam(IloCplex::Param::MIP::Limits::Nodes, 1) ;
    LRCplex.solve() ;

    valHeuristicCplex  = LRCplex.getObjValue() ;
    LRCplex.getValues(solution_x, x) ;
    LRCplex.getValues(solution_y, y) ;

}

//double CplexChecker::printSolution() {
//    int T = inst->getT() ;
//    int n = inst->getn();
//  /*  cplex.solve() ;
//    double objvalue = cplex.getObjValue() ;*/
//    IloNumArray solution = IloNumArray(env, n*T) ;
//    IloNumArray solution_u = IloNumArray(env, n*T) ;
//    IloNumArray solution_p = IloNumArray(env, n*T) ;
//    cplex.getValues(solution, x) ;
//    cplex.getValues(solution_u, u) ;
//    cplex.getValues(solution_p, pp) ;

//    cout.precision(6);
//    cout << "X: " << endl ;
//    for (int t=0 ; t < T ; t++) {
//        for (int i=0 ; i < n ; i++) {
//            cout << fabs(solution[i*T+t]) << " " ;
//        }
//        cout << endl ;
//    }
//    cout << endl ;

//    cout << "U: " << endl ;
//    for (int t=0 ; t < T ; t++) {
//        for (int i=0 ; i < n ; i++) {
//            cout << fabs(solution_u[i*T+t]) << " " ;
//        }
//        cout << endl ;
//    }
//    cout << endl ;



//    return 0 ;
//}



void CplexChecker::checkSolution(const vector<double> & y_frac, const vector<double> & x_frac) {

    double eps = 0.000001 ;

    cout << "start check..." << endl;

    IloModel CheckModel(env) ;

    int I = inst->getI();
    int J = inst->getJ() ;

    // Objective Function: Minimize Cost
    for (int j=0 ; j < J ; j++) {
        cost += inst->getb(j) * y_frac[j];
        for (int i=0; i<I; i++) {
            cost += inst->geta(j,i) * x_frac[j*I + i];
        }
    }
    CheckModel.add(IloMinimize(env, cost));


    // Assignment of each customer to a facility
    for (int i=0; i<I; i++) {
        double sum = 0 ;
        for (int j=0 ; j < J ; j++) {
            sum += x_frac[j*I + i];
        }
        if ((sum < 1 - eps) || (sum > 1 + eps)){
            cout << "assignment " << i << " non satisfaite" << endl ;
            cout << "somme des x : " << sum << endl;
        }
    }

    // Capacity constraints
    for (int j=0 ; j < J ; j++) {
        double sum = 0 ;
        for (int i=0; i<I; i++) {
            sum += x_frac[j*I + i] * inst->getd(i);
        }
        if(sum > y_frac[j] * inst->getc(j) + eps) {
            cout << "capacité " << j << " non satisfaite" << endl ;
            cout << "somme des demandes : " << sum << endl;
            cout << "capa : " << y_frac[j] * inst->getc(j) << endl;
        }
    }

    // Reliability constraints
    for (int i=0; i<I; i++) {
        for (int j=0 ; j < J ; j++) {
            double sum = 0 ;
            for (int indice=0; indice<inst->getv(); indice++){
                sum += y_frac[inst->getV(j)[indice]] * inst->getc(inst->getV(j)[indice]);
            }
            if(sum < x_frac[j*I + i] * inst->getd(i) * inst->getK() - eps) {
                cout << "fiabilité " << i << "," << j << " non satisfaite" << endl ;
                cout << "somme des capacités ouvertes : " << sum << endl;
                cout << "demande ouverte : " << x_frac[j*I + i] * inst->getd(i) * inst->getK() << endl;
            }
        }
    }



    IloCplex CheckCplex = IloCplex(CheckModel) ;

    cout << "Solve..." << endl ;
    CheckCplex.solve() ;
    cout << "end solve" << endl ;

    bool fea = CheckCplex.isPrimalFeasible();
    cout << "feasible: " << fea << endl ;
    double value = CheckCplex.getObjValue() ;
    cout << "value: " <<  value << endl ;
    
}

