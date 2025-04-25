#include <iostream>
#include <cmath>
#include <chrono>
#include <unistd.h>

#include "meshclass.hpp"
#include "quadrature.hpp"
#include "configclass.hpp"
#include "fem.hpp"


std::string ConfigClass::PATH = "../meshes/";

int main(int argc, char* argv[]){
    ConfigClass* conf;
    if (argc == 2){
        conf = new ConfigClass(argv[1]);
    }

    else {
        std::cout << "WARNING : No config file provided, resorting to default parameters (see doc)\n";
        conf = new ConfigClass;
    }

    Mesh mesh(conf->pathToMesh);

    mesh.buildConnectivity();
    mesh.domainSummary();

    FEMSolver solver(mesh, conf->functionF, conf->functionG, conf->integrationOrder);

        solver.matrixAssembly();
        solver.solveSystemGC();
        solver.exportSolution("test.vtk", conf->flabel);

    if (conf->solution != 0) {
        std::cout << "L2 err: " << solver.normL2(conf->solution, conf->integrationOrder) << std::endl;
        std::cout << "H1 err: " << solver.normH1(conf->solution, conf->dx_solution, conf->dy_solution, conf->integrationOrder) << std::endl;
    }

    return 0;
}