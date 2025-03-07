#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <set>
#include <cmath>
#include <unordered_set>
#include <map>
#include <unordered_map>

#include <Eigen/Dense>
#include "meshclass.hpp"

/**
 * NODE CLASS IMPLEMENTATION
*/

double Node::distance(const Node& other) const {
    return std::sqrt(std::pow(coefficients[0] - other.coefficients[0], 2) +
                     std::pow(coefficients[1] - other.coefficients[1], 2));
}

double Node::operator()(int i) const{
    return coefficients[i];
}

std::ostream& operator<<(std::ostream& os, const Node& node) {
    os << "(" << node.coefficients[0] << ", " << node.coefficients[1] << ")";
    return os;
}

/**
 * TRIANGLE ELEMENT CLASS IMPLEMENTATION
*/

double TriangleElement::getAera() const {
    int a, b, c;
    Eigen::Matrix<double, 3, 3> M3;
    a = globalNodeIndex[0]; b = globalNodeIndex[1]; c = globalNodeIndex[2];
    
    M3 << 1.0, (*globalNodeTab)[a](0), (*globalNodeTab)[a](1),
          1.0, (*globalNodeTab)[b](0), (*globalNodeTab)[b](1),
          1.0, (*globalNodeTab)[c](0), (*globalNodeTab)[c](1);
    
    return M3.determinant() / 2.0;
}

double TriangleElement::getPerimeter() const{
    int a, b, c;
    a = globalNodeIndex[0]; b = globalNodeIndex[1]; c = globalNodeIndex[2];

    double perimeter = (*globalNodeTab)[a].distance((*globalNodeTab)[b])
                     + (*globalNodeTab)[b].distance((*globalNodeTab)[c])
                     + (*globalNodeTab)[c].distance((*globalNodeTab)[a]);

    return perimeter;
}

std::ostream& operator<<(std::ostream& os, const TriangleElement& elem) {
    os << "A = " << (*elem.globalNodeTab)[elem.globalNodeIndex[0]]
       << ", B = " << (*elem.globalNodeTab)[elem.globalNodeIndex[1]]
       << ", C = " << (*elem.globalNodeTab)[elem.globalNodeIndex[2]];
    return os;
}


/**
 * FACET CLASS IMPLEMENTATION
*/

double Facet::length() const {
    return (*globalNodeTab)[globalNodeIndex[0]].distance((*globalNodeTab)[globalNodeIndex[1]]);
}

std::ostream& operator<<(std::ostream& os, const Facet& facet) {
    os << "Nodes on facet\n";
    os << "A = " << facet.globalNodeIndex[0]
       << ", B = " << facet.globalNodeIndex[1] << std::endl
       << "Id = " << facet.identifier << std::endl;

    os << "Is boundary facet ? " << facet.isBoundary << std::endl;
    return os;
}


/**
 * MESH CLASS IMPLEMENTATION
*/

Mesh::Mesh(std::string path) {

    std::ifstream meshFile(path, std::ios::in);

    nbFacets = 0;
    nbTriangles = 0;

    if (meshFile){
        std::string firstlines;
        firstlines.reserve(1048);
        
        // Get the first lines out of the way
        std::getline(meshFile, firstlines);
        while (firstlines != "$PhysicalNames"){
            std::getline(meshFile, firstlines);
        }

        std::getline(meshFile, firstlines);

        std::string lines;
        lines.reserve(1048);
        std::getline(meshFile, lines);
        std::istringstream stream(lines);

        // Get the physical names
        int dimension; int physicalId; std::string physicalName;
        while (lines != "$EndPhysicalNames"){
            stream.clear();
            stream.str(lines);

            stream >> dimension >> physicalId >> physicalName;

            physicalMarkers.insert({physicalName, physicalId});
            markerDimension.insert({physicalId, dimension});
            std::getline(meshFile, lines);
        }

        while (lines != "$Nodes"){
            std::getline(meshFile, lines);
        }

        // Get the nodes
        std::getline(meshFile, lines);

        stream.clear();
        stream.str(lines);

        stream >> nbNodes;
        tabNodes.resize(nbNodes);
        idToIndexNodes.reserve(nbNodes);

        std::getline(meshFile, lines);

        int tabIndex = 0;
        while (lines != "$EndNodes"){
            stream.clear();
            stream.str(lines);

            int id; double x, y, z;            
            stream >> id >> x >> y >> z;
            
            tabNodes[tabIndex] = Node(x, y, id);
            idToIndexNodes[id] = tabIndex;

            std::getline(meshFile, lines);
            tabIndex ++;
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
        tabTriangle.reserve(Ne);
        tabFacets.reserve(Ne);
        nextFacetId = 0;

        std::getline(meshFile, lines);
        while (lines != "$EndElements"){
            stream.clear();
            stream.str(lines);

            unsigned int id; int marker; int nbTags; int physical;
            stream >> id >> marker >> nbTags >> physical; 

            for (int i = 1; i < nbTags; i ++){int trash; stream >> trash;}

            if (marker == 1){

                // The element is a segment
                int p1, p2;
                stream >> p1 >> p2;
                const std::array<int, 2> &tab = {idToIndexNodes[p1], 
                                                 idToIndexNodes[p2]};

                tabFacets.push_back(Facet(tab, tabNodes, id, true));

                markedFacets[id] = physical;
                idToIndexFacets[id] = nbFacets;

                nbFacets ++; nextFacetId = std::max<int>(id, nextFacetId);
            }

            if (marker == 2){

                // The element is a triangle
                int p1, p2, p3;
                stream >> p1 >> p2 >> p3;
                const std::array <int, 3> &tab = {idToIndexNodes[p1],
                                                  idToIndexNodes[p2],
                                                  idToIndexNodes[p3]};

                tabTriangle.push_back(TriangleElement(tab, tabNodes, id));

                // Build partial connectivity
                markedElements[id] = physical;
                idToIndexTriangles[id] = nbTriangles;

                nbTriangles ++;
            }

            std::getline(meshFile, lines);
        }

        tabTriangle.shrink_to_fit();
        tabFacets.shrink_to_fit();

        meshFile.close();
    }

    else {
        std::cout << "Could not open " << path << std::endl;
    }
}

void Mesh::addNode(const Node& node) {
    tabNodes.push_back(node);
}

void Mesh::addElement(const TriangleElement& element) {
    tabTriangle.push_back(element);
}

int Mesh::getNbNodes() const {
    return nbNodes;
}

int Mesh::getNbElements() const {
    return nbTriangles;
}

int Mesh::getNbFacets() const {
    return nbFacets;
}

int Mesh::getNbSegments() const {
    int res = 0;
    for (const auto& face : tabFacets) {
        res += face.isSegement;
    }
    return res;
}

const Node& Mesh::getNode(int nodeId) const {
    int index = idToIndexNodes.at(nodeId);
    return tabNodes[index];
}

const TriangleElement& Mesh::getElement(int elementId) const {
    int index = idToIndexTriangles.at(elementId);
    return tabTriangle[index];
}

const Facet& Mesh::getFacet(int facetId) const {
    int index = idToIndexFacets.at(facetId);
    return tabFacets[index];
}

void Mesh::buildConnectivity() {

    nodeToElements.clear();
    nodeToFacets.clear();
    facetToElements.clear();
    elementToFacets.clear();

    nodeToElements.reserve(nbNodes);
    nodeToFacets.reserve(nbNodes);
    elementToFacets.reserve(nbTriangles);

    // Temporary container to handle duplicates
    std::map<std::array<int, 2>, Facet> facetMap;
    std::map<int, std::array<int, 2>> facetToElementsTemp;

    for (int i = 0; i < nbFacets; i ++){
        std::array<int, 2> node = swapFaceNodes(tabFacets[i].globalNodeIndex);
        facetMap[node] = tabFacets[i];
        facetToElementsTemp[tabFacets[i].identifier] = {-1, -1};
    }

    // Build facet map
    int facetId = 0;
    nextFacetId ++;

    for (int t = 0; t < nbTriangles; ++ t) {
        int triangleId = tabTriangle[t].identifier;
        std::array<int, 3> globalFacet = {-1, -1, -1};

        for (int f = 0; f < 3; ++ f){
            std::array<int, 2> nodeId(getFaceNodes(tabTriangle[t], f));
            nodeId = swapFaceNodes(nodeId);

            if (facetMap.find(nodeId) == facetMap.end()){
                facetId = nextFacetId;
            } else {
                facetId = facetMap.at(nodeId).identifier;
            }

            globalFacet[f] = facetId;

            nodeToElements[nodeId[0]].insert(triangleId);
            nodeToElements[nodeId[1]].insert(triangleId);

            nodeToFacets[nodeId[0]].insert(facetId);
            nodeToFacets[nodeId[1]].insert(facetId);

            if (facetMap.find(nodeId) == facetMap.end()){
                facetMap[nodeId] = Facet(nodeId, tabNodes, facetId, false);
                facetToElementsTemp[facetId] = {triangleId, -1};

                nbFacets ++; nextFacetId ++;

            } else if (facetToElementsTemp[facetId][0] == -1){
                facetToElementsTemp[facetId][0] = triangleId;

            } else {
                facetToElementsTemp[facetId][1] = triangleId;
            }
        }

        elementToFacets[triangleId] = globalFacet;
    }

    // Setting boundary flag
    for (auto& facet : facetMap) {
        int currentId = facet.second.identifier;

        if (facetToElementsTemp.at(currentId)[0] == -1 ||
            facetToElementsTemp.at(currentId)[1] == -1) {
                facet.second.isBoundary = true;
        }
    }

    // Building facet array
    nbFacets = facetMap.size();
    tabFacets.resize(nbFacets);
    int index = 0;

    for (auto& facet : facetMap) {
        tabFacets[index] = facet.second;
        idToIndexFacets[facet.second.identifier] = index;
        index ++;
    }

    facetToElements.reserve(nbFacets);

    for (auto& facet : facetToElementsTemp) {
        facetToElements[facet.first] = facet.second;
    }
}

std::vector<int> Mesh::getElementsForNode(int nodeId) const {
    return std::vector<int>(nodeToElements.at(nodeId).begin(), nodeToElements.at(nodeId).end());
}

std::vector<int> Mesh::getFacetsForNode(int nodeId) const {
    return std::vector<int>(nodeToFacets.at(nodeId).begin(), nodeToFacets.at(nodeId).end());
}

std::array<int, 2> Mesh::getElementsForFacets(int facetId) const {
    if (facetToElements.find(facetId) != facetToElements.end()){
        return facetToElements.at(facetId);
    } else {
        return {-1, -1};
    }
}

std::array<int, 3> Mesh::getNeighborTriangles(int triangleId) const {
    std::array<int, 3> neighbors = {-1, -1, -1};

    for (int i = 0; i < 3; ++i) {
        int facetId = elementToFacets.at(triangleId)[i];

        std::array<int, 2> neighborElementsOfFacet = facetToElements.at(facetId);

        if (neighborElementsOfFacet[0] != -1 && neighborElementsOfFacet[0] != triangleId) {
            neighbors[i] = neighborElementsOfFacet[0];
        } 
        
        if (neighborElementsOfFacet[1] != -1 && neighborElementsOfFacet[1] != triangleId) {
            neighbors[i] = neighborElementsOfFacet[1];
        }
    }

    return neighbors;
}

bool Mesh::isFacetOnBoundary(int facetId) const {
    int index = idToIndexFacets.at(facetId);
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

void Mesh::domainSummary() const{
    std::cout << "$Physical Names" << std::endl;

    for (const auto& names : physicalMarkers.left) {
        std::cout << names.first << " - " << names.second;

        int counter = 0;
        for (const auto& elem : markedElements) {
            counter += (elem.second == names.second);
        } for (const auto& facet : markedFacets) {
            counter += (facet.second == names.second);
        }

        std::cout << " : " << counter << std::endl;
    }
}

std::vector<int> Mesh::getMarkedElements(const std::string& name) const {
    std::vector<int> elemID; elemID.reserve(markedElements.size());

    int physicalID = physicalMarkers.left.at(name);

    for (const auto& elem : markedElements){
        if (elem.second == physicalID) {
            elemID.push_back(elem.first);
        }
    }

    elemID.shrink_to_fit();
    return elemID;
}

std::vector<int> Mesh::getMarkedFacet(const std::string& name) const {
    std::vector<int> facID; facID.reserve(markedFacets.size());

    int physicalID = physicalMarkers.left.at(name);

    for (const auto& face : markedFacets){
        if (face.second == physicalID) {
            facID.push_back(face.first);
        }
    }

    facID.shrink_to_fit();
    return facID;
}

void Mesh::printNodes() const{
    for (auto& node : tabNodes){
        std::cout << node << std::endl;
    }
}

void Mesh::printTriangles() const{
    for (auto& triangle : tabTriangle){
        std::cout << triangle << std::endl;
    }   
}

void Mesh::printFacets() const{
    for (auto& facet : tabFacets){
        std::cout << facet << std::endl;
    }
}

double Mesh::meshPerimeter() const{
    double perimeter = 0.0;

    for (auto facet : tabFacets){
        if (facet.isBoundary){
            perimeter += facet.length();
        }
    }

    return perimeter;
}

double Mesh::meshAera() const{
    double aera = 0.0;

    for (auto elements : tabTriangle){
        aera += elements.getAera();
    }

    return aera;
}

int Mesh::exportToVTK(const std::string& path, const std::string& plotName, functionType function = 0){
    std::ofstream ofile(path, std::ios::out);

    if (ofile) {
        ofile << "# vtk DataFile Version 2.0\n"
                << "plotToVTK in meshclass\n"
                << "ASCII\n"
                << "DATASET UNSTRUCTURED_GRID\n"
                << "POINTS " << nbNodes << " float\n";
        
        for (unsigned int i = 0; i < nbNodes; i ++){
            std::array<double, 2> coeff = tabNodes[i].coefficients;
            ofile << coeff[0] << ' ' << coeff[1] << ' ' << 0 << '\n';
        }

        std::cout << "PING " << std::endl;

        ofile << "CELLS " << nbTriangles << ' ' << 4 * nbTriangles << '\n';
        for (unsigned int i = 0; i < nbTriangles; i ++){
            std::array<int, 3> nodeIndex = tabTriangle[i].globalNodeIndex;

            ofile << 3 << ' ' << nodeIndex[0]
                       << ' ' << nodeIndex[1]
                       << ' ' << nodeIndex[2] << '\n';
        }

        std::cout << "PING " << std::endl;
        ofile << "CELL_TYPES " << nbTriangles;
        for (unsigned int i = 0; i < nbTriangles; i ++){ofile << '\n' << 5;}

        if (function != 0){
            ofile << "\nPOINT_DATA " << nbNodes << '\n'
                  << "SCALARS " << plotName << " float 1\n"
                  << "LOOKUP_TABLE default" << '\n';

            for (int i = 0; i < nbNodes; i ++){
                ofile << function(tabNodes[i](0), tabNodes[i](1)) << '\n';
            }
        }

        return 1;
    } else {
        return 0;
    }

    ofile.close();
}
