#include <domain.h>
#include <grid.h>
#include <rk4.h>
#include <cmath>
#include <cstdio>
#include <maxwell.h>
#include <maxwellparser.h>
#include <polynomialinterpolator.h>
#include <iostream>

int main(int argc, char *argv[]) 
{

    if (argc < 2) {
        std::cerr << "Usage:  " << argv[0] << " <parameter file>" << std::endl;
        exit(1);
    }
    MaxwellParameters pars;
    MaxwellParser parser;
    parser.updateParameters(argv[1], &pars); 

    char *fnames[2];
    fnames[0] = new char[16];
    fnames[1] = new char[16];
    sprintf(fnames[0], "Ey");
    sprintf(fnames[1], "Bz");


    Domain domain = Domain();
    double bounds[2];
    bounds[0] = pars.getxmin();
    bounds[1] = pars.getxmax();
    int N = pars.getnpoints();
    domain.setBounds(bounds);
    double cfl = pars.getcfl();
    double tmax = pars.gettmax();
    domain.setCFL(cfl);
    int output_frequency = pars.getoutput_frequency();

    std::cout << "Creating grid with " << N << " points and bounds [" << bounds[0] << ", " << bounds[1] << "]" << std::endl;

    domain.addGrid(bounds, N);

    RK4 rk4 = RK4();
    PolynomialInterpolator interp = PolynomialInterpolator(4);

    Maxwell ode = Maxwell(domain, rk4);
    ode.setInterpolator(&interp);

    ode.setParameters(&pars);
    ode.initData();

    double ti = 0.0;
    double dt = domain.getCFL()*(--domain.getGrids().end())->getSpacing();

    unsigned int M = (tmax - ti) / dt;
//  ode.output_frame(fnames[0], 0.0, 0);
//  ode.output_frame(fnames[1], 0.0, 1);
    int it = 0;
    ode.outputVTKScalar(fnames[0], 0.0, it, 0);
    ode.outputVTKScalar(fnames[1], 0.0, it, 1);

    std::cout << "Beginning evolution with " << M << " steps" << std::endl;

    double t = 0.0;
    for (unsigned int i = 1; i <= M; i++) {
        ode.evolveStep(dt);
        t += dt;
        if (i % output_frequency == 0) {
            std::cout << "Step " << i << " time " << t << std::endl;
//          ode.output_frame(fnames[0], t, 0);
//          ode.output_frame(fnames[1], t, 1);
            it++;
            ode.outputVTKScalar(fnames[0], t, it, 0);
            ode.outputVTKScalar(fnames[1], t, it, 1);
          
        }
    }

    delete [] fnames[0];
    delete [] fnames[1];

    return 0;
}
