#include <iostream>
#include <cmath>
#include <chrono>
#include <unistd.h>

#include "meshclass.hpp"
#include "quadrature.hpp"

/*
static int allocations = 0;

void* operator new(std::size_t size) {
    allocations++;
    return malloc(size);
}
*/

inline double constant(double x, double y) {
    return 1.0;
}

inline double quadratic(double x, double y) {
    return x*x + y*y;
}


int main(int argc, char* argv[]){

    functionType fun(&quadratic);

    Mesh mesh("../meshes/square2d_M2.msh");

    mesh.domainSummary();
    mesh.buildConnectivity();

    /*
    mesh.printTriangles();
    mesh.printNodes();
    */

    std::cout << "Number of nodes: " << mesh.getNbNodes() << std::endl;
    std::cout << "Number of elements: " << mesh.getNbElements() << std::endl;
    std::cout << "Number of facets: " << mesh.getNbFacets() << std::endl;

    std::cout << "Number of segments: " << mesh.getNbSegments() << std::endl;

    mesh.getMarkedElements("\"Omega\"");
    mesh.exportToVTK("test.vtk", "name", fun);

    /*

    for (int elem : mesh.getElementsForNode(4)){
        std::cout << "Element " << elem << " is connected to node 4" << std::endl;
    }

    for (int facet : mesh.getFacetsForNode(0)){
        std::cout << "Facet " << facet << " is connected to node 0" << std::endl;
    }

    for (int elem : mesh.getElementsForFacets(5)){
        std::cout << "Element " << elem << " is connected to facet 5" << std::endl;
    }

    */

    //std::cout << "Allocations : " << allocations << std::endl;

    //allocations = 0;

    //std::cout << "Allocations : " << allocations << std::endl;
    std::cout << "Aera: " << mesh.meshAera() << std::endl;
    std::cout << "Perimeter: " << mesh.meshPerimeter() << std::endl;
    MeshIntegration quadMesh(mesh);
    std::cout << "Integral: " << quadMesh.integrateOverMesh(&quadratic, 2) << std::endl;

    return 0;
}