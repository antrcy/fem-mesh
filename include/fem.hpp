#ifndef FEM_H
#define FEM_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "meshclass.hpp"
#include "quadrature.hpp"

using functionType = std::function<double (double, double)>;

class P1TriangleBasis {
private:
    const Mesh& mesh;
    int elementId;
    double surface;

public:
    P1TriangleBasis(const Mesh& mesh, int elemId) : mesh(mesh), elementId(elemId) {
        surface = mesh.getTriangleAera(elementId);
    }

    functionType functionPhi(int localDof) const;

    Eigen::Vector2d gradientPhi(int localDof) const;
};


#endif