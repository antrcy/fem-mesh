#ifndef QUADRATURE_H
#define QUADRATURE_H

#include <iostream>
#include <map>

#include "meshclass.hpp"

typedef std::pair<std::array<double, 3>, double> nodeWeightPair;

struct QuadratureRule {
    int integrationOrder;
    std::vector<nodeWeightPair> baryNodesAndWeights;

    QuadratureRule(int order): integrationOrder(order) {
        switch(integrationOrder) {
        case 2:
            baryNodesAndWeights.resize(3);
            baryNodesAndWeights.push_back({{2./3., 1./6., 1./6.}, 1./3.});
            baryNodesAndWeights.push_back({{1./6., 2./3., 1./6.}, 1./3.});
            baryNodesAndWeights.push_back({{1./6., 1./6., 2./3.}, 1./3.});
            break;
        
        case 3:
            baryNodesAndWeights.resize(4);
            baryNodesAndWeights.push_back({{1./3., 1./3., 1./3.}, -9./16.});
            baryNodesAndWeights.push_back({{3./5., 1./5., 1./5.}, 25./48.});
            baryNodesAndWeights.push_back({{1./5., 3./5., 1./5.}, 25./48.});
            baryNodesAndWeights.push_back({{1./5., 1./5., 3./5.}, 25./48.});
            break;
        
        default:
            baryNodesAndWeights.resize(1);
            baryNodesAndWeights.push_back({{1./3., 1./3., 1./3.}, 1.0});
            break;
        }
    }

    static QuadratureRule getQuadratureRule(int order){return QuadratureRule(order);}
};

class MeshIntegration {
public:
    MeshIntegration(const Mesh& mesh): M_mesh(mesh) {}

    double integrateOverTriangle(const functionType& f, int order, int triangleId) const;

    double integrateOverMesh(const functionType& f, int order) const;

private:
    const Mesh& M_mesh;
};

#endif