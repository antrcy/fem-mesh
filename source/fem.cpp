#include "fem.hpp"
#include "bimap.hpp"

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

    for (const auto& elemId : domain.idToIndexTriangles) {

        submatrixAssembly(Ak, Fk, elemId.first);
        std::array<int, 3> nodeIndex = domain.getElement(elemId.first).globalNodeIndex;
        std::array<int, 3> nodeId = domain.getNodeFromElem(elemId.first);

        for (int s = 0; s < 3; ++ s) {

            int i = nodeIndex[s];
            bool onBoundary = domain.isNodeOnBoundary(nodeId[s]);
        
            for (int t = 0; t < 3; ++ t) {
                int j = nodeIndex[t];

                if (onBoundary == false) {
                    vecInit.push_back(Eigen::Triplet<double>(i, j, Ak(s, t)));
                }
            }
            vectorF(i) = vectorF(i) + Fk(s);
        }
    }

    int k = 0;
    for (auto it : domain.idToIndexNodes) {
        if (domain.isNodeOnBoundary(it.first)) {
            k = it.second;
            vecInit.push_back(Eigen::Triplet<double>(k, k, 1.0));
            vectorF(k) = functionG(domain.getNode(it.first)(0),
                                   domain.getNode(it.first)(1));
        }
    }

    matrixA.setFromTriplets(vecInit.begin(), vecInit.end());

    std::cout << matrixA << vectorF << std::endl;

    for (auto it : domain.idToIndexNodes) {
        if (domain.isFacetOnBoundary(it.first)) {
            k = it.second;
            for (typename sparseType::InnerIterator it(matrixA, k); it; ++it) {
                double gi = vectorF(it.col());

                if (domain.isNodeOnBoundary(domain.getNodeId(it.col()))) {
                    vectorF(k) -= vectorF(it.col()) * it.valueRef();
                    it.valueRef() = 0.0;
                }
            }
        }
    }

    std::cout << matrixA << vectorF << std::endl;
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
    MeshIntegration quadrature(domain);
    realExpr.evaluate(expr);

    double norm = 0.0;

    for (auto elemId : domain.idToIndexTriangles) {
        P1TriangleBasis shape(domain, elemId.first);
        
        functionType squaredDiff([&](double x, double y) {
            return std::pow(shape.functionPhi(0)(x, y) * (realExpr.getValue(elemId.first, 0) - solution.getValue(elemId.first, 0))
                 + shape.functionPhi(1)(x, y) * (realExpr.getValue(elemId.first, 1) - solution.getValue(elemId.first, 1))
                 + shape.functionPhi(2)(x, y) * (realExpr.getValue(elemId.first, 2) - solution.getValue(elemId.first, 2)), 2);
        });
        norm += quadrature.integrateOverTriangle(squaredDiff, integrationOrder, elemId.first);
    }

    return std::sqrt(norm);
}



int FEMSolver::exportSolution(const std::string& path, const std::string& plotName) const {    
    std::ofstream ofile(path, std::ios::out);
    biMap<int, int> order;

    if (ofile) {
        int index = 0;
        ofile << "# vtk DataFile Version 2.0\n"
                << "plotToVTK in meshclass\n"
                << "ASCII\n"
                << "DATASET UNSTRUCTURED_GRID\n"
                << "POINTS " << domain.getNbNodes() << " float\n";


        for (auto id : domain.idToIndexNodes) {
            const Node& node = domain.getNode(id.first);
            ofile << node(0) << ' ' << node(1) << ' ' << 0 << '\n';

            order.insert({id.first, index});
            index ++;
        }

        const int nbTriangles = domain.getNbElements();

        ofile << "CELLS " << nbTriangles << ' ' << 4 * nbTriangles << '\n';
        for (auto elem : domain.idToIndexTriangles) {
            std::array<int, 3> nodeIds = domain.getNodeFromElem(elem.first);
            ofile << 3 << ' ' << order.at_first(nodeIds[0])
                       << ' ' << order.at_first(nodeIds[1])
                       << ' ' << order.at_first(nodeIds[2]) << '\n';
        }

        ofile << "CELL_TYPES " << nbTriangles;
        for (unsigned int i = 0; i < nbTriangles; i ++) {ofile << '\n' << 5;}

        ofile << "\nPOINT_DATA " << domain.getNbNodes() << '\n'
              << "SCALARS " << plotName << " float 1\n"
              << "LOOKUP_TABLE default" << '\n';
        
            
        for (int i = 0; i < domain.getNbNodes(); i ++) {
            ofile << solution.getValue(order.at_second(i)) << '\n';
        }

        ofile.close();
        return 1;
    }
    
    ofile.close();
    return 0;
}