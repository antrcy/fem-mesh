#include "fem.hpp"
#include "bimap.hpp"

functionType operator*(functionType f, functionType g);

functionType P1TriangleBasis::functionPhi(int localDof) const {
    int a1 = (localDof + 1) % 3;
    int a2 = (localDof + 2) % 3;

    Node p1(mesh.getNodeFromElem(elementId, a1));
    Node p2(mesh.getNodeFromElem(elementId, a2));

    return functionType([&, p1, p2](double x, double y) {
        return (p1(0) * p2(1) - p2(0) * p1(1)
         + x * (p1(1) - p2(1))
         + y * (p2(0) - p1(0))) / 2.0 / surface;
    });
}

Eigen::Vector2d P1TriangleBasis::gradientPhi(int localDof) const {
    int a1 = (localDof + 1) % 3;
    int a2 = (localDof + 2) % 3;

    Node p1(mesh.getNodeFromElem(elementId, a1));
    Node p2(mesh.getNodeFromElem(elementId, a2));

    return Eigen::Vector2d((p1(1) - p2(1)) / 2.0 / surface,
                           (p2(0) - p1(0)) / 2.0 / surface);
}


double P1TriangleBasis::localL2dot(int order, std::array<double, 3> coeffF, std::array<double, 3> coeffG) const {
    MeshIntegration integrator(mesh);
    functionType phi0 = functionPhi(0);
    functionType phi1 = functionPhi(1);
    functionType phi2 = functionPhi(2);

    functionType funF([&, coeffF](double x, double y){
        return coeffF[0] * phi0(x, y) + coeffF[1] * phi1(x, y) + coeffF[2] * phi2(x, y);
    });

    functionType funG([&, coeffG](double x, double y){
        return coeffG[0] * phi0(x, y) + coeffG[1] * phi1(x, y) + coeffG[2] * phi2(x, y);
    });

    return integrator.integrateOverTriangle(funG * funF, order, elementId);
}

double P1TriangleBasis::localH1dot(int order, std::array<double, 3> coeffF, std::array<double, 3> coeffG) const {
    MeshIntegration integrator(mesh);

    double dot = localL2dot(order, coeffF, coeffG);

    Eigen::Vector2d gradF(coeffF[0] * gradientPhi(0) + coeffF[1] * gradientPhi(1) + coeffF[2] * gradientPhi(2));
    Eigen::Vector2d gradG(coeffG[0] * gradientPhi(0) + coeffG[1] * gradientPhi(1) + coeffG[2] * gradientPhi(2));

    return dot + integrator.integrateOverTriangle(gradF.dot(gradG), elementId);
}


functionType operator*(functionType f, functionType g) {
    return [&](double x, double y) {
        return f(x, y) * g(x, y);
    };
}

functionType operator*(const functionType& f, double scalar) {
    return [&, scalar](double x, double y) {
        return f(x, y) * scalar;
    };
}

functionType operator+(const functionType& f, const functionType& g) {
    return [&](double x, double y) {
        return f(x, y) + g(x, y);
    };
}

functionType operator-(const functionType& f, const functionType& g) {
    return [&](double x, double y) {
        return f(x, y) - g(x, y);
    };
}

void FEMSolver::submatrixAssembly(Eigen::Matrix3d& elemAk, Eigen::Vector3d& elemFk, int triangleId) const {
    P1TriangleBasis shape(domain, triangleId);
    MeshIntegration quadrature(domain);

    Eigen::Vector2d gradArray[] = {shape.gradientPhi(0),
                                   shape.gradientPhi(1),
                                   shape.gradientPhi(2)};
    
    for (int i = 0; i < 3; ++ i) {
        for (int j = i; j < 3; j ++) {
            double aij = gradArray[i].dot(gradArray[j]);
            elemAk(i, j) = elemAk(j, i) = quadrature.integrateOverTriangle(aij, triangleId);
        }
    }

    elemFk(0) = quadrature.integrateOverTriangle(functionF * shape.functionPhi(0), integrationOrder, triangleId);
    elemFk(1) = quadrature.integrateOverTriangle(functionF * shape.functionPhi(1), integrationOrder, triangleId);
    elemFk(2) = quadrature.integrateOverTriangle(functionF * shape.functionPhi(2), integrationOrder, triangleId);
}

void FEMSolver::matrixAssemblySym() {
    return;
}

void FEMSolver::matrixAssemblyAsym() {
    Eigen::Matrix3d Ak;
    Eigen::Vector3d Fk;

    std::vector<Eigen::Triplet<double>> vecInit;
    vecInit.reserve(dimension);

    for (int elemIndex = 0; elemIndex < domain.getNbElements(); elemIndex ++) {

        submatrixAssembly(Ak, Fk, elemIndex);
        std::array<int, 3> nodeIndex = domain.getNodeFromElem(elemIndex);

        for (int s = 0; s < 3; ++ s) {

            int i = nodeIndex[s];
            bool onBoundary = domain.isNodeOnBoundary(i);
        
            for (int t = 0; t < 3; ++ t) {
                int j = nodeIndex[t];

                if (onBoundary == false) {
                    vecInit.push_back(Eigen::Triplet<double>(i, j, Ak(s, t)));
                }
            }
            vectorF(i) = vectorF(i) + Fk(s);
        }
    }

    for (auto i : domain.getBoundary()) {
        vecInit.push_back(Eigen::Triplet<double>(i, i, 1.0));
        vectorF(i) = functionG(domain.getNode(i)(0),
                               domain.getNode(i)(1));
    }

    matrixA.setFromTriplets(vecInit.begin(), vecInit.end());

    for (int node = 0; node < domain.getNbNodes(); node ++) {

        if (!domain.isNodeOnBoundary(node)) {
            for (typename sparseType::InnerIterator it(matrixA, node); it; ++it) {
                double gi = vectorF(it.col());

                if (domain.isNodeOnBoundary(it.col())) {
                    vectorF(node) -= gi * it.valueRef();
                    it.valueRef() = 0.0;
                }
            }
        }
    }
}

void FEMSolver::solveSystemGC() {
    Eigen::ConjugateGradient<sparseType, Eigen::Lower | Eigen::Upper> gc;
    gc.compute(matrixA);

    if (gc.info() != Eigen::Success) {
        std::cout << "GC initialization failed" << std::endl;
    }

    solution.initialize(gc.solve(vectorF));

    if (gc.info() != Eigen::Success) {
        std::cout << "Failed to converge" << std::endl;
    }

    std::cout << "#iterations:     " << gc.iterations() << std::endl;
    std::cout << "estimated error: " << gc.error() << std::endl;
}

void FEMSolver::solveSystemLU() {
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    solver.analyzePattern(matrixA);
    solver.factorize(matrixA); 
    solution.initialize(solver.solve(vectorF));
}

double FEMSolver::normL2(functionType expr) const {
    FunctionSpace::FunctionElement realExpr = fspace.element();
    realExpr.evaluate(expr);

    double norm = 0.0;

    for (int elem = 0; elem < domain.getNbElements(); elem ++) {
        P1TriangleBasis shape(domain, elem);

        std::array<double, 3> coeffs = {realExpr.getValue(elem, 0) - solution.getValue(elem, 0),
                                        realExpr.getValue(elem, 1) - solution.getValue(elem, 1),
                                        realExpr.getValue(elem, 2) - solution.getValue(elem, 2)};

        norm += shape.localL2dot(integrationOrder, coeffs, coeffs);
    }

    return std::sqrt(norm);
}

double FEMSolver::normH1(functionType expr) const {
    FunctionSpace::FunctionElement realExpr = fspace.element();
    realExpr.evaluate(expr);

    double norm = 0.0;

    for (int elem = 0; elem < domain.getNbElements(); elem ++) {
        P1TriangleBasis shape(domain, elem);

        std::array<double, 3> coeffs = {realExpr.getValue(elem, 0) - solution.getValue(elem, 0),
                                        realExpr.getValue(elem, 1) - solution.getValue(elem, 1),
                                        realExpr.getValue(elem, 2) - solution.getValue(elem, 2)};

        norm += shape.localH1dot(integrationOrder, coeffs, coeffs);
    }

    return norm;
}



int FEMSolver::exportSolution(const std::string& path, const std::string& plotName) const {    
    std::ofstream ofile(path, std::ios::out);

    if (ofile) {
        int index = 0;
        ofile << "# vtk DataFile Version 2.0\n"
                << "plotToVTK in meshclass\n"
                << "ASCII\n"
                << "DATASET UNSTRUCTURED_GRID\n"
                << "POINTS " << domain.getNbNodes() << " float\n";


        for (int node = 0; node < domain.getNbNodes(); node ++) {
            const Node& p = domain.getNode(node);
            ofile << p(0) << ' ' << p(1) << ' ' << 0 << '\n';
        }

        const int nbTriangles = domain.getNbElements();

        ofile << "CELLS " << nbTriangles << ' ' << 4 * nbTriangles << '\n';
        for (int elem = 0; elem < domain.getNbElements(); elem ++) {
            std::array<int, 3> nodeIds = domain.getNodeFromElem(elem);
            ofile << 3 << ' ' << nodeIds[0]
                       << ' ' << nodeIds[1]
                       << ' ' << nodeIds[2] << '\n';
        }

        ofile << "CELL_TYPES " << nbTriangles;
        for (unsigned int i = 0; i < nbTriangles; i ++) {ofile << '\n' << 5;}

        ofile << "\nPOINT_DATA " << domain.getNbNodes() << '\n'
              << "SCALARS " << plotName << " float 1\n"
              << "LOOKUP_TABLE default" << '\n';
        
            
        for (int i = 0; i < domain.getNbNodes(); i ++) {
            ofile << solution.getValue(i) << '\n';
        }

        ofile.close();
        return 1;
    }
    
    ofile.close();
    return 0;
}