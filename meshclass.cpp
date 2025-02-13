#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <set>
#include <cmath>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <chrono>
#include <unistd.h>

#include "meshclass.hpp"

static int allocations = 0;

void* operator new(std::size_t size) {
    allocations++;
    return malloc(size);
}

/**
 * NODE CLASS IMPLEMENTATION
*/

double Node::distance(const Node& other) const {
    return std::sqrt(std::pow(coefficients[0] - other.coefficients[0], 2) +
                     std::pow(coefficients[1] - other.coefficients[1], 2));
}

std::ostream& operator<<(std::ostream& os, const Node& node) {
    os << "(" << node.coefficients[0] << ", " << node.coefficients[1] << ")";
    return os;
}


/**
 * TRIANGLE ELEMENT CLASS IMPLEMENTATION
*/

std::ostream& operator<<(std::ostream& os, const TriangleElements& elem) {
    os << "A = " << (*elem.globalNodeTab)[elem.globalNodeId[0]]
       << ", B = " << (*elem.globalNodeTab)[elem.globalNodeId[1]]
       << ", C = " << (*elem.globalNodeTab)[elem.globalNodeId[2]];
    return os;
}


/**
 * FACET CLASS IMPLEMENTATION
*/

double Facet::length() const {
    return (*globalNodeTab)[globalNodeId[0]].distance((*globalNodeTab)[globalNodeId[1]]);
}

std::ostream& operator<<(std::ostream& os, const Facet& facet) {
    os << "Nodes on facet\n";
    os << "A = " << facet.globalNodeId[0]
       << ", B = " << facet.globalNodeId[1] << std::endl;

    os << "Elements on facet\n";
    os << "Elem1 - " << facet.globalElemId[0]
       << ", Elem2 - " << facet.globalElemId[1] << std::endl;

    os << "Is boundary facet ? " << facet.isBoundary << std::endl;
    return os;
}


/**
 * MESH CLASS IMPLEMENTATION
*/

Mesh::Mesh(std::string path) {

    std::ifstream meshFile(path, std::ios::in);

    nbFacets = 0;

    if (meshFile){
        std::string firstlines;
        firstlines.reserve(1048);
        
        // Extracting mesh nodes
        std::getline(meshFile, firstlines);
        while (firstlines != "$Nodes"){
            std::getline(meshFile, firstlines);
        }

        std::string lines;
        lines.reserve(1048);
        std::getline(meshFile, lines);
        std::istringstream stream(lines);

        stream >> nbNodes;
        tabNodes.resize(nbNodes);

        std::getline(meshFile, lines);

        while (lines != "$EndNodes"){
            stream.clear();
            stream.str(lines);

            unsigned int id; double x, y, z;            
            stream >> id >> x >> y >> z;
            id --;

            tabNodes[id] = Node(x, y);

            std::getline(meshFile, lines);
        }

        // Extracting elements (skipping segments)

        while (lines != "$Elements"){
            std::getline(meshFile, lines);
        }

        std::getline(meshFile, lines);
        stream.clear();
        stream.str(lines);

        int Ne;
        stream >> Ne;
        tabElements.reserve(Ne);
        nbElements = 0;

        std::getline(meshFile, lines);
        while (lines != "$EndElements"){
            stream.clear();
            stream.str(lines);
            unsigned int id; int marker; int tags;
            stream >> id >> marker >> tags >> tags >> tags;

            if (marker == 2){
                int p1, p2, p3;
                stream >> p1 >> p2 >> p3;
                const std::array <int, 3> &tab = {--p1, --p2, --p3};

                tabElements.push_back(TriangleElements(tab, tabNodes));
                nbElements ++;
            }

            std::getline(meshFile, lines);
        }

        meshFile.close();
    }

    else {
        std::cout << "Could not open " << path << std::endl;
    }
}

void Mesh::addNode(const Node& node) {
    tabNodes.push_back(node);
}

void Mesh::addElement(const TriangleElements& element) {
    tabElements.push_back(element);
}

int Mesh::getNbNodes() const {
    return nbNodes;
}

int Mesh::getNbElements() const {
    return nbElements;
}

int Mesh::getNbFacets() const {
    return nbFacets;
}

const Node& Mesh::getNode(int nodeId) const {
    return tabNodes[nodeId];
}

const TriangleElements& Mesh::getElement(int elementId) const {
    return tabElements[elementId];
}

const Facet& Mesh::getFacet(int facetId) const {
    return tabFacets[facetId];
}

void Mesh::buildConnectivity() {
    tabFacets.clear();
    nodeToFacets.clear();
    elementToFacets.clear();

    nodeToElements.reserve(nbNodes);
    nodeToFacets.reserve(nbNodes);
    elementToFacets.reserve(nbElements);

    nbFacets = 0;

    // Temporary container to handle duplicates
    std::map<std::array<int, 2>, Facet> facetMap;

    // Build facet map
    int facetId = 0;
    for (int element = 0; element < nbElements; ++element) {
        std::array<int, 3> globalFacet = {-1, -1, -1};

        for (int localFacet = 0; localFacet < 3; ++localFacet) {
            
            std::array<int, 2> nodeId(getFaceNodes(tabElements[element], localFacet));
            // Making sure the facet is oriented
            swapFaceNodes(nodeId);

            // Building nodeToElements connectivity
            nodeToElements[nodeId[0]].insert(element);
            nodeToElements[nodeId[1]].insert(element);

            if (facetMap.find(nodeId) == facetMap.end()) {
                facetMap[nodeId] = Facet(nodeId, {element, -1}, tabNodes, facetId);
                globalFacet[localFacet] = facetId;

                // Building nodeToFacets connectivity
                nodeToFacets[nodeId[0]].insert(facetId);
                nodeToFacets[nodeId[1]].insert(facetId);

                nbFacets ++; facetId ++;

            } else {
                facetMap[nodeId].setElemId(1, element);
                globalFacet[localFacet] = facetMap[nodeId].identifier;
            }
        }
        // Building elementToFacets connectivity
        elementToFacets[element] = globalFacet;
    }

    // Building facet array
    tabFacets.resize(facetMap.size());
    for (auto& facet : facetMap) {
        tabFacets[facet.second.identifier] = facet.second;
    }
}

std::vector<int> Mesh::getElementsForNode(int nodeId) const {
    return std::vector<int>(nodeToElements.at(nodeId).begin(), nodeToElements.at(nodeId).end());
}

std::vector<int> Mesh::getFacetsForNode(int nodeId) const {
    return std::vector<int>(nodeToFacets.at(nodeId).begin(), nodeToFacets.at(nodeId).end());
}

std::vector<int> Mesh::getElementsForFacets(int facetId) const {
    std::vector<int> elements;
    elements.reserve(2);

    if (tabFacets[facetId].globalElemId[0] != -1) {
        elements.push_back(tabFacets[facetId].globalElemId[0]);
    } if (tabFacets[facetId].globalElemId[1] != -1) {
        elements.push_back(tabFacets[facetId].globalElemId[1]);
    }

    return elements;
}

std::vector<int> Mesh::getNeighborTriangles(int triangleId) const {
    std::vector<int> neighbors;
    neighbors.reserve(3);

    for (int i = 0; i < 3; ++i) {
        int facetId = elementToFacets.at(triangleId)[i];

        // Copy to simplify
        std::array<int, 2> neighborElementsOfFacet = {tabFacets[facetId].globalElemId[0], tabFacets[facetId].globalElemId[1]};

        if (neighborElementsOfFacet[0] != -1 && neighborElementsOfFacet[0] != triangleId) {
            neighbors.push_back(neighborElementsOfFacet[0]);
        } 
        
        if (neighborElementsOfFacet[1] != -1 && neighborElementsOfFacet[1] != triangleId) {
            neighbors.push_back(neighborElementsOfFacet[1]);
        }
    }

    return neighbors;
}

bool Mesh::isFacetOnBoundary(int facetId) const {
    return tabFacets[facetId].isBoundary;
}

bool Mesh::isNodeOnBoundary(int nodeId) const {
    bool isOnBoundary = true;
    for (int facetId : nodeToFacets.at(nodeId)) {
        isOnBoundary = isOnBoundary && isFacetOnBoundary(facetId);
    }
    return !isOnBoundary;
}

bool Mesh::isTriangleOnBoundary(int triangleId) const {
    bool isOnBoundary = true;
    for (int facetId : elementToFacets.at(triangleId)) {
        isOnBoundary = isOnBoundary && isFacetOnBoundary(facetId);
    }
    return !isOnBoundary;
}

void Mesh::printNodes() const{
    for (auto& node : tabNodes){
        std::cout << node << std::endl;
    }
}

void Mesh::printTriangles() const{
    for (auto& triangle : tabElements){
        std::cout << triangle << std::endl;
    }   
}

void Mesh::printFacets() const{
    for (auto& facet : tabFacets){
        std::cout << facet << std::endl;
    }
}

double Mesh::perimeter() const{
    double perimeter = 0.0;

    for (auto facet : tabFacets){
        if (facet.isBoundary){
            perimeter += facet.length();
        }
    }

    return perimeter;
}

int main(int argc, char* argv[]){

    Mesh mesh("square2d_perforated.msh");
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

    std::cout << "Allocations : " << allocations << std::endl;

    allocations = 0;

    auto start = std::chrono::steady_clock::now();
    std::cout << "Perimeter: " << mesh.perimeter() << std::endl;
    auto end = std::chrono::steady_clock::now();
    std::cout << "Elapsed time (ms) : " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << std::endl;

    std::cout << "Allocations : " << allocations << std::endl;

    return 0;
}