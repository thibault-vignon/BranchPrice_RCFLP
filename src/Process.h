#ifndef PROCESSS
#define PROCESSS

#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <math.h>
#include <ctime>
#include <list>
#include <vector>

#include "InstanceRCFLP.h"
using namespace std ;

class Parameters
{
public:
    bool ColumnGeneration;
    int nodeLimit;
    bool PriceAndBranch ;

    double Epsilon ;

    // si on écrit les contraintes de capacité en une seule contrainte (false : en deux contraintes, la deuxieme etant x<=y)
    bool compactCapacityConstraints ; 

    bool heuristicInit ;
    bool Farkas ;

    bool DynProgFacility ;
    bool DynProgCustomer ;

    bool balanceCostsX;
    bool balanceCostsY;
    vector<double> costBalancingMaster;
    vector<double> costBalancingPricer;

    bool FacilityDecompo ;
    bool CustomerDecompo ;
    bool doubleDecompo ;

    //Options de pricing
    bool heurPricing;
    double heurPricingThreshold;
    bool useLowerBound;

    Parameters(
            bool ColumnGeneration, int nodeLimit, bool PriceAndBranch,
            double epsilon,
            bool compactCapacityConstraints,
            bool heuristicInit, bool Farkas,
            bool DynProgFacility, bool DynProgCustomer,
            bool balanceCostsX, bool balanceCostsY, 
            bool FacilityDecompo, bool CustomerDecompo, bool doubleDecompo, 
            bool heurPricing, double heurPricingThreshold,
            bool useLowerBound) ;

};

class InstanceProcessed {
public :

    //donnees
    int v ;
    double K;
    string instanceSet ;
    string id ;

    InstanceProcessed(int v_, double K_, string instanceSet_, string id_)
    {
        v = v_ ;
        K = K_;
        instanceSet = instanceSet_ ;
        id = id_ ;
    }
    
    string fileName() ;

    void test();

};



// initializes Parameters class using met indicator
Parameters init_parameters(InstanceRCFLP* inst, int met) ;

#endif 
