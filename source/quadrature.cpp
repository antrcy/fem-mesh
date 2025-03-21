#include "quadrature.hpp"

double MeshIntegration::integrateOverTriangle(const functionType& f, const QuadratureRule& rule, int triangleId) const {
    double surface = domain.getTriangleAera(triangleId);
    double integral = 0.0;

    std::array<Node, 3> nodes = {domain.getNodeFromElem(triangleId, 0),
                                 domain.getNodeFromElem(triangleId, 1),
                                 domain.getNodeFromElem(triangleId, 2)};

    for (const auto& pair : rule.baryNodesAndWeights) {
        Node quadNode({pair.first[0] * nodes[0](0) +
                       pair.first[1] * nodes[1](0) +
                       pair.first[2] * nodes[2](0),

                       pair.first[0] * nodes[0](1) +
                       pair.first[1] * nodes[1](1) +
                       pair.first[2] * nodes[2](1)});

        integral += surface * pair.second * f(quadNode(0), quadNode(1));
    }

    return integral; 
}

double MeshIntegration::integrateOverTriangle(const functionType& f, int order, int triangleId) const {
    QuadratureRule quadRule = QuadratureRule::getQuadratureRule(order);
    return integrateOverTriangle(f, quadRule, triangleId);
}

double MeshIntegration::integrateOverTriangle(const double& value, int triangleId) const {
    return domain.getTriangleAera(triangleId) * value;
}


double MeshIntegration::integrateOverMesh(const functionType& f, int order) const {
    QuadratureRule quadRule = QuadratureRule::getQuadratureRule(order);
    double integral = 0.0;

    for (int triangleIndex = 0; triangleIndex < domain.getNbElements(); triangleIndex ++) {
        integral += integrateOverTriangle(f, quadRule, triangleIndex);
    }

    return integral;
}