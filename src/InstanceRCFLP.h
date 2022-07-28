#ifndef INSTANCERCFLP
#define INSTANCERCFLP

#include <ilcplex/ilocplex.h>
#include <fstream>

class InstanceRCFLP {

private:
    IloEnv env ;

    ///// Recupéré du fichier instance /////

    IloInt I, J, v;
    IloNum K;
    IloNumArray c, d, b;
    IloNumArray2 a;
    IloIntArray2 V;

public:

    InstanceRCFLP(IloEnv envir, const char* file) ;
    ~InstanceRCFLP() {
        c.end();
        d.end();
        b.end();

        a.end();
        V.end();
    }

public:
    void Lecture(const char* file);
    void Initialise() ;

    IloEnv getenv() ;

    IloInt getI() const ;
    IloInt getJ() const ;
    IloInt getv() const ;
    IloNum getK() const ;

    IloNum getc(IloInt j) const;
    IloNum getd(IloInt i) const;
    IloNum getb(IloInt j) const;

    IloNum geta(IloInt j, IloInt i) const;
    IloIntArray getV(IloInt j) const;

};

#endif /* INSTANCERCFLP_INCLUDED */
