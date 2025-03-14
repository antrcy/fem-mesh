#include <cmath>
#include <Eigen/Dense>

#include "quadrature.hpp"

double MeshIntegration::integrateOverTriangle(const functionType& f, int order, int triangleId) const {
    QuadratureRule quadRule = QuadratureRule::getQuadratureRule(order);
    double surface = M_mesh.getTriangleAera(triangleId);
    double integral = 0.0;

    std::array<Node, 3> nodes = M_mesh.getNodeFromElem(triangleId);

    for (const auto& pair : quadRule.baryNodesAndWeights) {
        Node quadNode({pair.first[0] * nodes[0](0) +
                       pair.first[1] * nodes[1](0) +
                       pair.first[2] * nodes[2](0),

                       pair.first[0] * nodes[0](1) +
                       pair.first[1] * nodes[1](1) +
                       pair.first[2] * nodes[2](1)}, -1);
        integral += surface * pair.second * f(quadNode(0), quadNode(1));
    }

    return integral;
}

double MeshIntegration::integrateOverMesh(const functionType& f, int order) const {
        
    double integral = 0.0;

    for (const auto& triangle : M_mesh.idToIndexTriangles) {
        integral += integrateOverTriangle(f, order, triangle.first);
    }

    return integral;
}