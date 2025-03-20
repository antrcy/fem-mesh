#ifndef QUADRATURE_H
#define QUADRATURE_H

#include <map>

#include "meshclass.hpp"

typedef std::pair<std::array<double, 3>, double> nodeWeightPair;

/**
 * @brief Container for integration points and wheights.
 */
struct QuadratureRule {
    // Attributes 
    int integrationOrder;                            // Order 1 to 3 come pre baked.
    std::vector<nodeWeightPair> baryNodesAndWeights; // Integration points in barycentric coordinates with their wiehgts

    // Constructor
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

/**
 * @brief Computes the integral of an expression over a triangle or the entire mesh.
 */
class MeshIntegration {
public:
    MeshIntegration(const Mesh& mesh): M_mesh(mesh) {}

    /** @brief Integrates an expression over a single triangle given a quadrature rule, ideal for repeated integration. */
    double integrateOverTriangle(const functionType& f, const QuadratureRule& rule, int triangleId) const;
    /** @brief Integrates an expression over a single triangle given a quadrature order. */
    double integrateOverTriangle(const functionType& f, int order, int triangleId) const;
    /** @brief Integrates a constant value over a single triangle. */
    double integrateOverTriangle(const double& value, int triangleId) const;
    /** @brief Integrates an expression over the entire mesh given a quadrature order. */
    double integrateOverMesh(const functionType& f, int order) const;

private:
    const Mesh& M_mesh;
};

#endif