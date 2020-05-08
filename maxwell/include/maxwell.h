#ifndef MAXWELL_H
#define MAXWELL_H

#include <ode.h>

class Maxwell : public ODE {

    private:
        static const unsigned int U_EY = 0;
        static const unsigned int U_BZ = 1;

    protected:
//      virtual void applyBoundaries(bool intermediate);
        virtual void rhs(const Grid &grid, double **u, double **dtu);

    public:
        Maxwell(Domain &d, Solver &s);
        virtual ~Maxwell();

        virtual void initData();
};

#endif

