#include <iostream>
#include <cmath>
#include <chrono>
#include <unistd.h>

#include "meshclass.hpp"

/*
static int allocations = 0;

void* operator new(std::size_t size) {
    allocations++;
    return malloc(size);
}
*/

int main(int argc, char* argv[]){

    Mesh mesh("../meshes/square2d_perforated.msh");
    mesh.buildConnectivity();

    /*
    Mesh mesh("square2d_4elt.msh");

    std::cout << "Number of nodes: " << mesh.getNbNodes() << std::endl;
    std::cout << "Number of elements: " << mesh.getNbElements() << std::endl;
    std::cout << "Number of facets: " << mesh.getNbFacets() << std::endl;

    mesh.buildConnectivity();

    mesh.printFacets();
    mesh.printTriangles();
    mesh.printNodes();

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

    auto start = std::chrono::steady_clock::now();
    std::cout << "Perimeter: " << mesh.perimeter() << std::endl;
    auto end = std::chrono::steady_clock::now();
    std::cout << "Elapsed time (ms) : " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << std::endl;

    //std::cout << "Allocations : " << allocations << std::endl;

    testEigen();

    return 0;
}