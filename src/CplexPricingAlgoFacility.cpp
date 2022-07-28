#include "PricingAlgo.h"
#include <iostream>

using namespace std;

DualCostsFacility::DualCostsFacility(InstanceRCFLP* inst) {
    int I = inst->getI() ;
    int J = inst->getJ() ;

    Sigma.resize(J, 0) ;

    Omega1.resize(I*J, 0) ;
    Omega2.resize(I*J, 0) ;
}


void CplexPricingAlgo::AddSSBI(IloEnv env, IloModel model, IloBoolVarArray x, IloBoolVarArray u, int site, InstanceRCFLP* inst) {

    int allCoupleInequalities=0;

    int T = inst->getT() ;
    int firstOfSite = Param.firstUnit(site) ;
    int ns = Param.nbUnits(site) ;
    int k ;

    //inégalités symétries
    for (int i=0; i<ns; i++) {
        if (inst->getFirst(i+firstOfSite)) { // i est la première unité d'un groupe de symétrie
            int g = inst->getGroup(i+firstOfSite) ;

            int first = inst->getFirstG(g) - firstOfSite ;
            int last = inst->getLastG(g) - firstOfSite ;

            for (int i=first ; i <= last ; i++) {

                int l = inst->getl(i+firstOfSite) ;
                int L = inst->getL(i+firstOfSite) ;

                if (i < last) {
                    for (int t = l ; t < T ; t++) {

                        IloExpr rhs(env) ;
                        rhs += x[i*T+t] + x[i*T+t - l];

                        for (k=t -l + 1 ; k < t ; k++) {
                            rhs+= u[i*T+k] ;
                        }

                        int ub_j = i+1 ;
                        if (allCoupleInequalities) {
                            ub_j=last ;
                        }

                        for (int j=i+1 ; j <= ub_j ; j++) {
                            model.add(u[j*T + t] <= rhs) ;
                        }
                        rhs.end() ;
                    }
                }

                if (i>first) {

                    for (int t = L ; t < T ; t++) {

                        IloExpr rhs_w(env) ;
                        rhs_w += 2 - x[i*T+t] - x[i*T+t - L];

                        for (k=t -L + 1 ; k < t ; k++) {
                            rhs_w += x[i*T + k-1] - x[i*T+k] + u[i*T+k] ;
                        }


                        int lb_j = i-1 ;
                        if (allCoupleInequalities) {
                            lb_j=first ;
                        }
                        for (int j=i-1 ; j >= lb_j ; j--) {
                            model.add(x[j*T + t-1] - x[j*T+t] + u[j*T + t] <= rhs_w) ;
                        }

                        rhs_w.end() ;
                    }
                }
            }
        }
    }
}

CplexPricingAlgo::CplexPricingAlgo(InstanceRCFLP* inst, const Parameters & par, int site) : Param(par) {

    //env = IloEnv() ;
    Site=site ;

    int ns = Param.nbUnits(Site) ;
    int first = Param.firstUnit(Site) ;
    int last = first+ns-1;

    int T = inst->getT() ;
    model = IloModel(env) ;


    x = IloBoolVarArray(env, ns*T) ;
    u = IloBoolVarArray(env, ns*T) ;
    p = IloNumVarArray(env, ns*T, 0, 10000) ;

   cpuTime=0 ;


   if (Param.useSSBIinSubPb) {
    AddSSBI(env, model, x,u, site, inst) ;
   }

    obj = IloAdd(model, IloMinimize(env, 0.0));

    //    // Conditions initiales
    //    for (int i=0; i<ns; i++) {
    //        model.add(u[i*T] >= x[i*T] - 1 ) ;
    //    }

    //    for (int i=0; i<ns; i++) {
    //        IloExpr sum(env) ;
    //        for (int k= 0; k < inst->getl(first+i) ; k++) {
    //            sum += u[i*T + k] ;
    //        }
    //        model.add(sum <= 0 ) ;
    //        sum.end() ;
    //    }

    // Min up constraints
    for (int i=0; i<ns; i++) {
        for (int t=inst->getL(first+i) ; t < T ; t++) {
            IloExpr sum(env) ;
            for (int k= t - inst->getL(first+i) + 1; k <= t ; k++) {
                sum += u[i*T + k] ;
            }
            model.add(sum <= x[i*T + t]) ;
            sum.end() ;
        }
    }

    // Min down constraints
    for (int i=0; i<ns; i++) {
        for (int t=inst->getl(first+i) ; t < T ; t++) {
            IloExpr sum(env) ;
            for (int k= t - inst->getl(first+i) + 1; k <= t ; k++) {
                sum += u[i*T + k] ;
            }
            if (!Param.StartUpDecompo) {
                model.add(sum <= 1 - x[i*T + t - inst->getl(first+i)]) ;
            }
            else {
                model.add(sum == 1 - x[i*T + t - inst->getl(first+i)]) ;
            }
            sum.end() ;
        }
    }

    //Relation entre u et x
    for (int i=0; i<ns; i++) {
        for (int t=1 ; t < T ; t++) {
            model.add(x[i*T + t] - x[i*T + t-1] <= u[i*T + t]);
        }
    }

    //Contraintes intra-site
    if (Param.IntraSite && !Param.UnitDecompo) {
        for (int t=1 ; t < T ; t++) {
            IloExpr sum(env) ;
            for (int i=0 ; i <ns ; i++) {
                sum+=u[i*T+t] ;
            }
            model.add(sum <= 1);
            sum.end() ;
        }
    }

    if (Param.powerPlanGivenByLambda) {
        cout << "contrainte UB prise en compte" << endl ;
        //power limits
        for (int i=0; i<ns; i++) {
            for (int t=0 ; t < T ; t++) {
                model.add( (inst->getPmax(first+i) - inst->getPmin(first+i) )*x[i*T + t] >= p[i*T + t]);
            }
        }
    }

    if (Param.Ramp && Param.rampInSubPb) {
        for (int i = 0 ; i <ns ; i++) {

            double RU = (inst->getPmax(first+i) - inst->getPmin(first+i))/3 ;
            double RD = (inst->getPmax(first+i) - inst->getPmin(first+i))/2 ;

            for (int t = 1 ; t < T ; t++) {
                model.add(p[i*T + t] - p[i*T + t-1] <= RU*x[i*T + t-1] );
                model.add(p[i*T + t-1] - p[i*T + t] <= RD*x[i*T + t] );
            }
        }
    }

    if (Param.DemandeResiduelle) {

        double DR = 0 ;
        for (int i=0 ; i <inst->getn() ; i++) {
            if (i < first || i > last) {
                DR += inst->getPmax(i) ;
            }
        }

        //Demande
        for (int t=0; t < T ; t++) {
            IloExpr Prod(env) ;
            for (int i=0; i<ns; i++) {
                Prod += inst->getPmax(first+i)*x[i*T + t];
            }
            model.add(fmax(0,inst->getD(t) - DR) <= Prod);
            Prod.end() ;
        }
    }

    //Objectif pour u
    for (int i=0 ; i<ns ; i++) {
        for (int t=0 ; t < T ; t++) {
            obj.setLinearCoef(u[i*T +t],inst->getc0(first+i));
        }
    }

    cplex = IloCplex(model);
    cplex.setParam(IloCplex::EpGap, 0.00001) ;

}

void CplexPricingAlgo::updateObjCoefficients(InstanceRCFLP* inst, const Parameters & Param, const DualCosts & Dual, bool Farkas) {
    int ns = Param.nbUnits(Site) ;
    int T = inst->getT();
    int first = Param.firstUnit(Site) ;
    for (int i=0 ; i<ns ; i++) {
        for (int t=0 ; t < T ; t++) {

            obj.setLinearCoef(x[i*T +t], (Dual.ObjCoefX).at((first+i)*T+t)  );
            obj.setLinearCoef(u[i*T +t], (Dual.ObjCoefU).at((first+i)*T+t) ) ;
            if (Param.powerPlanGivenByLambda) {
                obj.setLinearCoef(p[i*T +t], (Dual.ObjCoefP).at((first+i)*T+t) ) ;
            }

        }
    }
}



bool CplexPricingAlgo::findUpDownPlan(InstanceRCFLP* inst, const DualCosts & Dual, IloNumArray UpDownPlan, double& objvalue) {
    //returns True if a solution has been found

    ofstream LogFile("LogFile.txt");
    cplex.setOut(LogFile);

    if ( !cplex.solve() ) {
        env.error() << "Failed to optimize Pricer with Cplex" << endl;
        exit(1);
    }

    if (cplex.getStatus()==CPX_STAT_INFEASIBLE){
        cout<<"NO SOLUTION TO PRICER"<<endl;
        cout<<endl<<" ************************* END PRICER with CPLEX"<<endl<<endl;

        return false;
    }
    else {

        cplex.getValues(x, UpDownPlan) ;

        /*cout << "for site " << Site << "; " << endl ;
       cout << "obj value without sigma: " << cplex.getObjValue() << endl;*/
        objvalue = cplex.getObjValue() - Dual.Sigma[Site] ;

        //cpuTime = cplex.getCpu
    }


    return true;
}

void DualCosts::computeRedcost(InstanceRCFLP* inst, const Parameters & Param, Master_Variable * lambda, double& redcost) {

    redcost = - Sigma[lambda->Site] ;
    int T= inst->getT() ;

    int first = Param.firstUnit(lambda->Site) ;
    int last = Param.firstUnit(lambda->Site) + Param.nbUnits(lambda->Site) - 1 ;

    for (int i = first ; i <= last ; i++) {
        int down_t_1 = 0 ;
        for (int t = 0 ; t < T ; t++) {
            if (down_t_1 && lambda->UpDown_plan[(i-first)*T+t] > 1 - Param.Epsilon) { // il y a un démarrage en t
                redcost+= ObjCoefU.at(i*T+t) ;
            }
            if (lambda->UpDown_plan[(i-first)*T+t] < Param.Epsilon ) {
                down_t_1=1;
            }
            else { // l'unité i est up en t
                redcost+= ObjCoefX.at(i*T+t) ;
                down_t_1=0;
            }
        }
    }

    if (Param.powerPlanGivenByLambda) {
        for (int i = first ; i <= last ; i++) {
            for (int t = 0 ; t < T ; t++) {
                redcost += lambda->Power_plan[(i-first)*T+t]*ObjCoefP.at(i*T+t) ;
            }
        }
    }

} 