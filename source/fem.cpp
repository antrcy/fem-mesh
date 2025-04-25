#include "fem.hpp"
#include "bimap.hpp"

using tranformation_t = Mesh::GeometricTransformation;

/**
 * P1 LAGRANGE BASIS METHODS
 */

Eigen::MatrixXd P1LagrangeBasis::evaluate(Eigen::Vector2d coord) const {
    Eigen::MatrixXd matrix(1, 3);
    matrix(0, 0) = 1 - coord(0) - coord(1);
    matrix(0, 1) = coord(0);
    matrix(0, 2) = coord(1);
    return matrix;
}

Eigen::MatrixXd P1LagrangeBasis::gradient(Eigen::Vector2d coord) const {
    return Eigen::MatrixXd({{-1, 1, 0},
                            {-1, 0, 1}});
}

/**
 * FEM SOLVER METHODS
 */

void FEMSolver::submatrixAssembly(Eigen::Matrix3d& elemAk, Eigen::Vector3d& elemFk, int triangleId) const {
    P1LagrangeBasis basis(domain);
    MeshIntegration quadrature(domain);

    tranformation_t T_K(domain, triangleId);
    double det = std::abs(T_K.jacobianDeterminant());

    Eigen::MatrixXd grad = T_K.jacobianInverse().transpose() * basis.gradient({0.0, 0.0});
    
    for (int i = 0; i < 3; ++ i) {
        for (int j = i; j < 3; j ++) {
            elemAk(i, j) = elemAk(j, i) = grad.col(i).dot(grad.col(j)) * det / 2.0; //gradArray[i].dot(gradArray[j]);
        }
    }

    functionType phi1 = [&](double x, double y) {
        Eigen::Vector2d node = T_K.map(Eigen::Vector2d({x, y}));
        return functionF(node(0), node(1)) * (1 - x - y) * det;
    }; elemFk(0) = quadrature.integrateOverRef(phi1, integrationOrder);

    functionType phi2 = [&](double x, double y) {
        Eigen::Vector2d node = T_K.map(Eigen::Vector2d({x, y}));
        return functionF(node(0), node(1)) * x * det;
    }; elemFk(1) = quadrature.integrateOverRef(phi2, integrationOrder);

    functionType phi3 = [&](double x, double y) {
        Eigen::Vector2d node = T_K.map(Eigen::Vector2d({x, y}));
        return functionF(node(0), node(1)) * y * det;
    }; elemFk(2) = quadrature.integrateOverRef(phi3, integrationOrder);
}

void FEMSolver::matrixAssembly() {
    if (assembled){return;}

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
    assembled = true;
}

void FEMSolver::solveSystemGC() {
    std::cout << std::endl;
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

double FEMSolver::normL2(functionType Solution, int order) const {
    double res = 0.0;

    MeshIntegration quadrature(domain);
    QuadratureRule quadRule = QuadratureRule::getQuadratureRule(order);

    for (int elemIndex = 0; elemIndex < domain.getNbElements(); elemIndex ++) {
        double u_0 = solution.getValue(elemIndex, 0);
        double u_1 = solution.getValue(elemIndex, 1);
        double u_2 = solution.getValue(elemIndex, 2);

        tranformation_t T_K(domain, elemIndex);
        double det = std::abs(T_K.jacobianDeterminant());

        functionType local_dot = [&](double x, double y) {
            Eigen::Vector2d node = T_K.map(Eigen::Vector2d({x, y}));

            double z = Solution(node(0), node(1)) - (1.0 - x - y) * u_0 - x * u_1 - y * u_2;
            
            return std::pow(z, 2) * det;
        };
        
        res += quadrature.integrateOverRef(local_dot, quadRule);
    }

    return std::sqrt(res);
}

double FEMSolver::normH1(functionType Solution, functionType grad_x, functionType grad_y, int order) const {
    double res = 0.0;

    P1LagrangeBasis basis(domain);
    MeshIntegration quadrature(domain);
    QuadratureRule quadRule = QuadratureRule::getQuadratureRule(order);

    for (int elemIndex = 0; elemIndex < domain.getNbElements(); elemIndex ++) {
        double u_0 = solution.getValue(elemIndex, 0);
        double u_1 = solution.getValue(elemIndex, 1);
        double u_2 = solution.getValue(elemIndex, 2);

        tranformation_t T_K(domain, elemIndex);
        double det = std::abs(T_K.jacobianDeterminant());

        Eigen::Vector2d grad = T_K.jacobianInverse().transpose() * basis.gradient({0.0, 0.0}) * Eigen::Vector3d({u_0, u_1, u_2});

        functionType local_dot = [&](double x, double y) {
            Eigen::Vector2d node = T_K.map(Eigen::Vector2d({x, y}));

            Eigen::Vector2d grad_diff = grad - Eigen::Vector2d({grad_x(node(0), node(1)),
                                                                grad_y(node(0), node(1))});
            
            return grad_diff.dot(grad_diff) * det;
        };

        res += quadrature.integrateOverRef(local_dot, quadRule);
    }

    double l2 = normL2(Solution, order);
    return std::sqrt(res + l2 * l2); 
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