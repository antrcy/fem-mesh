#ifndef FEM_H
#define FEM_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "meshclass.hpp"
#include "quadrature.hpp"

using functionType = std::function<double (double, double)>;
using funElem = FunctionSpace::FunctionElement;
using sparseType = Eigen::SparseMatrix<double, Eigen::RowMajor>;

/**
 * @brief Function basis for Lagrange elements in approximation space
 */
class BasisFunction {
protected:
    const Mesh& domain;

public:
    BasisFunction(const Mesh& mesh): domain(mesh) {}

    virtual Eigen::MatrixXd evaluate(Eigen::Vector2d coord) const = 0;
    virtual Eigen::MatrixXd gradient(Eigen::Vector2d coord) const = 0;
};

/**
 * @brief P1 Lagrange function basis implementation
 */
class P1LagrangeBasis: public BasisFunction {
public:
    P1LagrangeBasis(const Mesh& mesh): BasisFunction(mesh) {}

    Eigen::MatrixXd evaluate(Eigen::Vector2d coord) const;
    Eigen::MatrixXd gradient(Eigen::Vector2d coord) const;
};

/**
 * @brief Solver class for finite element system.
 */
class FEMSolver {
private:
    bool assembled;                // flag true if matrixAssembly step is done
    const Mesh& domain;            // discretized domain
    const functionType& functionF; // source
    const functionType& functionG; // boundary

    FunctionSpace fspace; // function space for the solution

    int dimension;           // equals to number of nodes
    sparseType matrixA;      // complete stiffness matrix
    Eigen::VectorXd vectorF; // complete stiffness vector
    FunctionSpace::FunctionElement solution; // function element of the solution

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

    /** @brief Assemble element matrix. */
    void submatrixAssembly(Eigen::Matrix3d& elemAk, Eigen::Vector3d& elemFk, int triangleId) const;
    /** @brief Assemble stiffness matrix. */
    void matrixAssembly();

    /** @brief Solve the system using conjugate gradient. */
    void solveSystemGC();
    /** @brief Solve the system with LU decomposition. */
    void solveSystemLU();

    /** @brief Computes the L2 distance between the solution and expr in P1 basis. */
    double normL2(functionType Solution, int order) const;
    /** @brief Computes the H1 distance between the solution and expr in P1 basis. */
    double normH1(functionType Solution, functionType grad_x, functionType grad_y, int order) const;

    /** @brief Builds a VTK file to visualize the solution with ParaView. */
    int exportSolution(const std::string& path, const std::string& plotName) const;

    /** @brief Reset matrix and vector content. */
    void reset(){
        matrixA.resize(dimension, dimension);
        vectorF.setZero();
    }
};


#endif