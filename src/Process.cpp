#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <math.h>
#include <ctime>
#include <list>

#include "InstanceRCFLP.h"
#include "Process.h"

using namespace std ;

string to_string2(int number){
    string number_string = "";
    char ones_char = '0';
    int ones = 0;
    while(true){
        ones = number % 10;
        switch(ones){
        case 0: ones_char = '0'; break;
        case 1: ones_char = '1'; break;
        case 2: ones_char = '2'; break;
        case 3: ones_char = '3'; break;
        case 4: ones_char = '4'; break;
        case 5: ones_char = '5'; break;
        case 6: ones_char = '6'; break;
        case 7: ones_char = '7'; break;
        case 8: ones_char = '8'; break;
        case 9: ones_char = '9'; break;
        default : ; //ErrorHandling("Trouble converting number to string.");
        }
        number -= ones;
        number_string = ones_char + number_string;
        if(number == 0){
            break;
        }
        number = number/10;
    }
    return number_string;
}

string InstanceProcessed::fileName() {
        return("Instances/Reliable/" + to_string2(v) + "_" + instanceSet + "_" + id + ".txt") ;
}

void InstanceProcessed::test(){
    cout << "ok" ;
}

Parameters::Parameters(bool ColumnGeneration, int nodeLimit, bool PriceAndBranch,
            double epsilon,
            bool compactCapacityConstraints,
            bool heuristicInit, bool Farkas,
            bool DynProgFacility, bool DynProgCustomer,
            bool balanceCosts, 
            bool FacilityDecompo, bool CustomerDecompo, bool doubleDecompo, 
            bool heurPricing, double heurPricingThreshold, 
            bool useLowerBound) :
    ColumnGeneration(ColumnGeneration),
    nodeLimit(nodeLimit),
    PriceAndBranch(PriceAndBranch),
    Epsilon(epsilon),
    compactCapacityConstraints(compactCapacityConstraints),
    heuristicInit(heuristicInit),
    Farkas(Farkas),
    DynProgFacility(DynProgFacility),
    DynProgCustomer(DynProgCustomer),
    balanceCosts(balanceCosts),
    FacilityDecompo(FacilityDecompo),
    CustomerDecompo(CustomerDecompo),
    doubleDecompo(doubleDecompo),
    heurPricing(heurPricing),
    heurPricingThreshold(heurPricingThreshold),
    useLowerBound(useLowerBound)
{



}



Parameters init_parameters(InstanceRCFLP* inst, int met) {

    bool ColumnGeneration = true ; // true: Branch & Price avec SCIP. false: résolution Cplex boîte noire

    if (met==-1) {
        ColumnGeneration= false ;
    }

    double eps = 0.0000001; // tolérance
    int node_limit = 1000000000;

    bool PriceAndBranch = false;

    bool compactCapacityConstraints = false;

    //// Paramètres: type de décomposition ////
    bool FacilityDecompo = false;
    bool CustomerDecompo = false;

    // Pricing
    bool DynProgFacility = false;
    bool DynProgCustomer = false;

    // double décompo
    bool doubleDecompo = false;
    bool balanceCosts = false;

    // Paramètres du pricing
    bool heurPricing = false;
    double heurPricingThreshold = 0;

    // Initialisation
    bool heuristicInit = false;
    bool Farkas = false ;

    // Interruption de la génération de colonnes par borne inf
    bool useLowerBound = false ;


    // Parse met value given as argument to infer parameters

    int arr[16];
    int indice = 0;
    int chiffre ;

    // On récupère tous les chiffres de met
    while(met != 0) {
        chiffre = met % 10 ;
        arr[indice] = chiffre ;
        indice++ ;
        met = met / 10 ;
    }

    // Puis on regarde les chiffres de gauche à droite pour fixer les paramètres

    // Type de décomposition
    switch (arr[indice - 1]) {
        case 1:
            // Par facility
            FacilityDecompo = true ;
            compactCapacityConstraints = true ;
            break ;

        case 2:
            // Par client
            CustomerDecompo = true ;
            compactCapacityConstraints = true ;
            break ;

        case 3:
            // Double contraintes de capacité compactes
            doubleDecompo = true ;
            compactCapacityConstraints = true ;
            break ;

        case 4:
            // Double contraintes de capacité alternatives
            doubleDecompo = true ;
            break ;
    }

    // Branchement ou non
    if (arr[indice - 2] == 0) {
        node_limit = 1 ;
    }

    // Paramètres généraux de génération de colonnes
    switch (arr[indice - 3]) {
        case 1:
            heuristicInit = true ;
            break ;
        case 2:
            useLowerBound = true ;
            break ;
    }

    // Réglages spécifiques double decomposition
    if (doubleDecompo){
        switch (arr[indice - 4]) {
            // Cas 0 : basique, sans repartition de couts

            case 1:
                balanceCosts = true ;
                break ;
        }
    }

    if (indice >= 5){
        DynProgFacility = false ;
    }

    if (indice >= 6){
        DynProgCustomer = false ;
    }

    Parameters param(
                ColumnGeneration, node_limit, PriceAndBranch, 
                eps, 
                compactCapacityConstraints,
                heuristicInit, Farkas,
                DynProgFacility, DynProgCustomer,
                balanceCosts, 
                FacilityDecompo, CustomerDecompo, doubleDecompo,
                heurPricing, heurPricingThreshold, 
                useLowerBound);


    return param;

}