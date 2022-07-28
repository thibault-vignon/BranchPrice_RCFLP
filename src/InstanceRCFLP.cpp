#include "InstanceRCFLP.h"

#include <ctime>
#include <stdlib.h> // use rand
#include <fstream> // lire fichier de données

using namespace std ;


InstanceRCFLP::InstanceRCFLP(IloEnv envir, const char* file) {
    env = envir ;
    Lecture(file) ;
    Initialise() ;
}




void InstanceRCFLP::Initialise() {

}

void InstanceRCFLP::Lecture(const char* file) {
    //Lecture de n et de T
    ifstream fichier(file, ios::in);

    string nom = "";
    fichier >> nom;

    while(nom!="="){
        fichier >> nom;
    }
    nom = "";
    fichier >> J;

    while(nom!="="){
        fichier >> nom;
    }
    nom = "";
    fichier >> I;

    while(nom!="="){
        fichier >> nom;
    }
    nom = "";
    fichier >> K;



    //Initialisation des vecteurs de taille I et J

    c = IloNumArray(env, J);
    b = IloNumArray(env, J);
    d = IloNumArray(env, I);
    a = IloNumArray2(env, I);
    V = IloIntArray2(env, J);

    //Lecture des données
    //Capacités
    while(nom!="="){
        fichier >> nom;
    }
    nom = "";
    fichier >> c;

    //Demandes
    while(nom!="="){
        fichier >> nom;
    }
    nom = "";
    fichier >> d;

    //Coûts fixes d'installation de site
    while(nom!="="){
        fichier >> nom;
    }
    nom = "";
    fichier >> b;

    //Coûts variables d'assignation
    while(nom!="="){
        fichier >> nom;
    }
    nom = "";
    fichier >> a;

    //Voisins
    while(nom!="="){
        fichier >> nom;
    }
    nom = "";
    fichier >> V;

    //fin lecture de fichier
}

// Accès aux données

IloInt InstanceRCFLP::getI() const {
    return I ;
}

IloInt InstanceRCFLP::getJ() const {
    return J ;
}

IloNum InstanceRCFLP::getK() const {
    return K ;
}

IloInt InstanceRCFLP::getv() const {
    return V[0].getSize() ;
}

IloNum InstanceRCFLP::getd(IloInt i) const {
    return d[i] ;
}

IloNum InstanceRCFLP::getc(IloInt j) const {
    return c[j] ;
}

IloNum InstanceRCFLP::getb(IloInt j) const {
    return b[j] ;
}

IloNum InstanceRCFLP::geta(IloInt j, IloInt i) const {
    return a[j][i] ;
}

IloIntArray InstanceRCFLP::getV(IloInt j) const {
    return V[j] ;
}