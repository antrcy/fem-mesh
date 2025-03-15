#include "fem.hpp"

functionType P1TriangleBasis::functionPhi(int localDof) const {
    int a1 = (localDof + 1) % 3;
    int a2 = (localDof + 2) % 3;

    Node p1 = mesh.getNodeFromElem(elementId)[a1];
    Node p2 = mesh.getNodeFromElem(elementId)[a2];

    return functionType([=](double x, double y) {
        return (p1(0) * p2(1) - p2(0) * p1(1)
         + x * (p1(1) - p2(1))
         + y * (p2(0) - p1(0))) / 2.0 / surface;
    });
}

Eigen::Vector2d P1TriangleBasis::gradientPhi(int localDof) const {
    int a1 = (localDof + 1) % 3;
    int a2 = (localDof + 2) % 3;

    Node p1 = mesh.getNodeFromElem(elementId)[a1];
    Node p2 = mesh.getNodeFromElem(elementId)[a2];

    return Eigen::Vector2d((p1(1) - p2(1)) / 2.0 / surface,
                           (p2(0) - p1(0)) / 2.0 / surface);
}

functionType operator*(functionType f, functionType g) {
    return [=](double x, double y) {
        return f(x, y) * g(x, y);
    };
}