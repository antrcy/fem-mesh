#ifndef FEM_H
#define FEM_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "meshclass.hpp"
#include "quadrature.hpp"

using functionType = std::function<double (double, double)>;
using funElem = FunctionSpace::FunctionElement;
using sparseType = Eigen::SparseMatrix<double, Eigen::RowMajor>;

class P1TriangleBasis {
private:
    const Mesh& mesh;
    int elementId;
    double surface;

    FunctionSpace fspace;

public:
    P1TriangleBasis(const Mesh& mesh, int elemId) : mesh(mesh), elementId(elemId), fspace(mesh) {
        surface = mesh.getTriangleAera(elementId);
    }

    functionType functionPhi(int localDof) const;

    Eigen::Vector2d gradientPhi(int localDof) const;
};

class FEMSolver {
private:
    bool assembled;
    const Mesh& domain;
    const functionType& functionF;
    const functionType& functionG;
    FunctionSpace fspace;

    int dimension;
    sparseType matrixA;
    Eigen::VectorXd vectorF;
    FunctionSpace::FunctionElement solution;

    int integrationOrder;

public:
    FEMSolver(Mesh& mesh, functionType argf, functionType argg, int order):
            assembled(false),
            domain(mesh),
            functionF(argf),
            functionG(argg),
            integrationOrder(order),
            fspace(mesh),
            solution(fspace.element())
    {
        dimension = mesh.getNbNodes();
        matrixA = sparseType(dimension, dimension);
        vectorF = Eigen::VectorXd(dimension);
        vectorF.setZero(); solution.setZero();
    }

    void submatrixAssembly(Eigen::Matrix3d& elemAk, Eigen::Vector3d& elemFk, int triangleId) const;
    void matrixAssemblySym();
    void matrixAssemblyAsym();

    void solveSystemGC();
    void solveSystemLU();

    double normL2(functionType expr) const;

    int exportSolution(const std::string& path, const std::string& plotName) const;

    void reset(){
        matrixA.resize(dimension, dimension);
        vectorF.setZero();
    }
};


#endif