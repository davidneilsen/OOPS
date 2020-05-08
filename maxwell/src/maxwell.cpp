#include <maxwell.h>
#include <operators.h>
#include <maxwellparameters.h>
#include <iostream>
#include <cmath>

// Constructor
Maxwell::Maxwell(Domain &d, Solver &s) : ODE(2,0) {

    if (d.getGhostPoints() < 2) {
        std::cerr << "Warning: domain has fewer ghost points than expected." << std::endl;
    }
    domain = &d;
    solver= &s;

//    params = new MaxwellParameters();

    reallocateData();

}

// Destructor
Maxwell::~Maxwell() {
//    delete params;
}

// Initial data routine
void Maxwell::initData() {

    double x0 = 0.0;
    double sigma = 0.25;

    // Loop over grids
    for (auto it = data.begin(); it != data.end(); ++it) {
        const double *x = it->getGrid().getPoints();
        unsigned int nx = it->getGrid().getSize();
        double **u= it->getData();
        for (unsigned int i = 0; i < nx; i++) {
            double val = std::exp(-(x[i] - x0)*(x[i] - x0)/(sigma*sigma));
            u[U_EY][i] = val;
            u[U_BZ][i] = -val*std::sin(x[i]);
        }
    }

}

/*--------------------------------------------------------------------------
 *
 *    dt_Ey = - dx_Bz
 *    dt_Bz = - dx_Ey
 *
 *--------------------------------------------------------------------------*/
void Maxwell::rhs(const Grid &grid, double **u, double **dtu) {


    double dx = grid.getSpacing();
    int nx = grid.getSize();

    dtu[U_EY][0] = - (u[U_BZ][1] - u[U_BZ][0]) / (dx);
    dtu[U_BZ][0] = - (u[U_EY][1] - u[U_EY][0]) / (dx);

    for (unsigned int i = 1; i < nx-1 ; i++) {
        dtu[U_EY][i] = - (u[U_BZ][i+1] - u[U_BZ][i-1]) / (2.0*dx);
        dtu[U_BZ][i] = - (u[U_EY][i+1] - u[U_EY][i-1]) / (2.0*dx);
    }

    dtu[U_EY][nx-1] = - (u[U_BZ][nx-1] - u[U_BZ][nx-2]) / (dx);
    dtu[U_BZ][nx-1] = - (u[U_EY][nx-1] - u[U_EY][nx-2]) / (dx);

}

void Maxwell::applyBoundaries(bool intermediate)
{
    unsigned int nb = domain->getGhostPoints();

    auto left_it = data.begin();
    auto right_it = --data.end();

    double **left;
    double **right;

    if (!intermediate) {
        left = left_it->getData();
        right = right_it->getData();
    }
    else {
        left = left_it->getIntermediateData();
        right = right_it->getIntermediateData();
    }

    unsigned int nr = right_it->getGrid().getSize();
    double dx = right_it->getGrid().getSpacing();

    left[U_EY][nb] = (left[U_EY][nb+1] - left[U_EY][nb]) / dx;
    left[U_BZ][nb] = (left[U_BZ][nb+1] - left[U_BZ][nb]) / dx;

    right[U_EY][nr-1-nb] = - (right[U_EY][nr-1-nb] - right[U_EY][nr-1-nb-1])/dx;
    right[U_BZ][nr-1-nb] = - (right[U_BZ][nr-1-nb] - right[U_BZ][nr-1-nb-1])/dx;

}

