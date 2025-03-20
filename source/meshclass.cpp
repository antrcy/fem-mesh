#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <set>
#include <cmath>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <stdexcept>

#include <Eigen/Dense>
#include "meshclass.hpp"
#include "bimap.hpp"

/**
 * NODE CLASS IMPLEMENTATION
*/

double Node::distance(const Node& other) const {
    return std::sqrt(std::pow(coefficients[0] - other.coefficients[0], 2) +
                     std::pow(coefficients[1] - other.coefficients[1], 2));
}

double Node::operator()(int localDof) const{
    return coefficients[localDof];
}

std::ostream& operator<<(std::ostream& os, const Node& node) {
    os << "(" << node.coefficients[0] << ", " << node.coefficients[1] << ")";
    return os;
}

bool Node::operator==(const Node& other) const{
    return (coefficients[0] - other.coefficients[0]) < 1e-6 && coefficients[1] == other.coefficients[1] < 1e-6;
}

/**
 * TRIANGLE ELEMENT METHODS
*/

double Mesh::getTriangleAera(int identifier) const {

    const TriangleElement& triangle = getElement(identifier);

    int a, b, c;
    Eigen::Matrix<double, 3, 3> M3;
    a = triangle[0]; b = triangle[1]; c = triangle[2];
    
    M3 << 1.0, (tabNodes)[a](0), (tabNodes)[a](1),
          1.0, (tabNodes)[b](0), (tabNodes)[b](1),
          1.0, (tabNodes)[c](0), (tabNodes)[c](1);
    
    return M3.determinant() / 2.0;
}

double Mesh::getTrianglePerimeter(int identifier) const {

    const TriangleElement& triangle = getElement(identifier);

    int a, b, c;
    a = triangle[0]; b = triangle[1]; c = triangle[2];

    double perimeter = (tabNodes)[a].distance((tabNodes)[b])
                     + (tabNodes)[b].distance((tabNodes)[c])
                     + (tabNodes)[c].distance((tabNodes)[a]);

    return perimeter;
}

/**
 * FACET METHODS
*/

double Mesh::getFacetLength(int identifier) const {

    const Facet& facet = getFacet(identifier);

    return tabNodes[facet[0]].distance(tabNodes[facet[1]]);
}

/**
 * MESH CLASS IMPLEMENTATION
*/

Mesh::Mesh(std::string path) {

    std::ifstream meshFile(path, std::ios::in);

    nbFacets = 0;
    nbTriangles = 0;
    nbNodes = 0;

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
        tabNodes.reserve(nbNodes);

        std::getline(meshFile, lines);

        while (lines != "$EndNodes"){
            stream.clear();
            stream.str(lines);

            int id; double x, y, z;            
            stream >> id >> x >> y >> z;
            
            tabNodes.push_back(Node(x, y));

            std::getline(meshFile, lines);
        }

        // Extracting elements
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

        int facetIndex = 0;
        int triangleIndex = 0;

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
                const std::array<int, 2> &tab = {--p1, --p2};

                tabFacets.push_back(Facet(tab, true));

                markedFacets[facetIndex] = physical;

                nbFacets ++; facetIndex ++;
            }

            if (marker == 2){

                // The element is a triangle
                int p1, p2, p3;
                stream >> p1 >> p2 >> p3;
                const std::array <int, 3> &tab = {--p1, --p2, --p3};

                tabTriangle.push_back(TriangleElement(tab));

                // Build partial connectivity
                markedElements[triangleIndex] = physical;

                nbTriangles ++; triangleIndex ++;
            }

            std::getline(meshFile, lines);
        }

        tabTriangle.shrink_to_fit();
        tabFacets.shrink_to_fit();

        meshFile.close();
    }

    else {
        throw std::runtime_error{"Could not open " + path};
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
        res += face.isSegment;
    }
    return res;
}

const Node& Mesh::getNode(int index) const {
    return tabNodes[index];
}

const TriangleElement& Mesh::getElement(int index) const {
    return tabTriangle[index];
}

const Facet& Mesh::getFacet(int index) const {
    return tabFacets[index];
}

std::array<int, 3> Mesh::getNodeFromElem(int index) const {
    return getElement(index).globalNodeIndex;
}

const Node& Mesh::getNodeFromElem(int identifier, int localDof) const {
    return tabNodes[getElement(identifier).globalNodeIndex[localDof]];
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
    std::map<std::array<int, 2>, std::pair<Facet, int>> facetMap;
    std::map<int, std::array<int, 2>> facetToElementsTemp;

    for (int i = 0; i < nbFacets; i ++){
        std::array<int, 2> nodeIndex = swapFaceNodes(tabFacets[i].globalNodeIndex);
        facetMap[nodeIndex] = {tabFacets[i], i};
        facetToElementsTemp[i] = {-1, -1};
    }

    // Build facet map
    int nextFacetIndex = nbFacets;

    for (int t = 0; t < nbTriangles; ++ t) {
        std::array<int, 3> globalFacet = {-1, -1, -1};

        for (int f = 0; f < 3; ++ f){
            int facetIndex;
            std::array<int, 2> nodeIndex(swapFaceNodes(getFaceNodes(tabTriangle[t], f)));

            if (facetMap.find(nodeIndex) == facetMap.end()){
                facetIndex = nextFacetIndex;
            } else {
                facetIndex = facetMap.at(nodeIndex).second;
            }

            globalFacet[f] = facetIndex;

            nodeToElements[nodeIndex[0]].insert(t);
            nodeToElements[nodeIndex[1]].insert(t);

            nodeToFacets[nodeIndex[0]].insert(facetIndex);
            nodeToFacets[nodeIndex[1]].insert(facetIndex);

            if (facetMap.find(nodeIndex) == facetMap.end()){
                facetMap[nodeIndex] = {Facet(nodeIndex, false), facetIndex};
                facetToElementsTemp[facetIndex] = {t, -1};

                nbFacets ++; nextFacetIndex ++;

            } else if (facetToElementsTemp.at(facetIndex)[0] == -1){
                facetToElementsTemp[facetIndex][0] = t;

            } else {
                facetToElementsTemp[facetIndex][1] = t;
            }
        }

        elementToFacets[t] = globalFacet;
    }


    // Setting boundary flag
    for (auto& facet : facetMap) {
        int currentIndex = facet.second.second;

        if (facetToElementsTemp.at(currentIndex)[0] == -1 ||
            facetToElementsTemp.at(currentIndex)[1] == -1) {
                facet.second.first.isBoundary = true;
        }
    }

    // Building facet array
    nbFacets = facetMap.size();
    tabFacets.resize(nbFacets);

    for (auto& facet : facetMap) {
        tabFacets[facet.second.second] = facet.second.first;
    }

    facetToElements.reserve(nbFacets);

    for (auto& facet : facetToElementsTemp) {
        facetToElements[facet.first] = facet.second;
    }

    for (int indexNodes = 0; indexNodes < nbNodes; indexNodes ++) {
        if (isNodeOnBoundary(indexNodes)) {
            boundaryNodes.insert(indexNodes);
        }
    }
}

const std::set<int>& Mesh::getBoundary() const {
    return boundaryNodes;
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
    return tabFacets[facetId].isBoundary;
}

bool Mesh::isNodeOnBoundary(int nodeId) const {
    int counter = 0;
    for (int facetId : nodeToFacets.at(nodeId)) {
        counter += isFacetOnBoundary(facetId);
    }
    return (counter >= 2);
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

    for (auto it = physicalMarkers.get_iterator(); !it.end(); ++it) {

        std::cout << it.get_left() << " - " << it.get_right();

        int counter = 0;
        for (const auto& elem : markedElements) {
            counter += (elem.second == it.get_right());
        } for (const auto& facet : markedFacets) {
            counter += (facet.second == it.get_right());
        }

        std::cout << " : " << counter << std::endl;
    }
}

std::vector<int> Mesh::getMarkedElements(const std::string& name) const {
    std::vector<int> elemID; elemID.reserve(markedElements.size());

    int physicalID = physicalMarkers.at_first(name);

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

    int physicalID = physicalMarkers.at_first(name);

    for (const auto& face : markedFacets){
        if (face.second == physicalID) {
            facID.push_back(face.first);
        }
    }

    facID.shrink_to_fit();
    return facID;
}

double Mesh::meshPerimeter() const{
    double perimeter = 0.0;

    for (int facetIndex = 0; facetIndex < nbFacets; facetIndex ++) {
        if (tabFacets[facetIndex].isBoundary){
            perimeter += getFacetLength(facetIndex);
        }
    }

    return perimeter;
}

double Mesh::meshAera() const{
    double aera = 0.0;

    for (int triIndex = 0; triIndex < nbTriangles; triIndex ++) {
        aera += getTriangleAera(triIndex);
    }

    return aera;
}

int Mesh::exportToVTK(const std::string& path, const std::string& plotName, const functionType& function) const {
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

        ofile << "CELLS " << nbTriangles << ' ' << 4 * nbTriangles << '\n';
        for (unsigned int i = 0; i < nbTriangles; i ++){
            std::array<int, 3> nodeIndex = tabTriangle[i].globalNodeIndex;

            ofile << 3 << ' ' << nodeIndex[0]
                       << ' ' << nodeIndex[1]
                       << ' ' << nodeIndex[2] << '\n';
        }

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

        ofile.close();
        return 1;
    }

    ofile.close();
    return 0;
}

/**
 * FUNCTIONSPACE CLASS IMPLEMENTATION
*/

unsigned int FunctionSpace::connectivity(int elementId, int localIndex) const {
    int globalIndex = (M_domain.getElement(elementId)).globalNodeIndex[localIndex];
    return globalIndex;
}

// FUNCTION ELEMENT

void FunctionSpace::FunctionElement::setZero() {
    int n = M_functionSpace.M_globalDof;
    for (int i = 0; i < n; i ++){
        M_data[i] = 0;
    }
}

void FunctionSpace::FunctionElement::initialize(const Eigen::VectorXd& dataVector) {
    if (dataVector.size() != M_functionSpace.M_globalDof) {
        throw std::invalid_argument{"Function dof does not match argument size.\n"};
    }

    for (int i = 0; i < dataVector.size(); i ++) {
        M_data[i] = dataVector(i);
    }
}

void FunctionSpace::FunctionElement::evaluate(const functionType& expression) {
    const Mesh& domain = M_functionSpace.M_domain;

    for (int nodeIndex = 0; nodeIndex < domain.getNbNodes(); nodeIndex ++){
        const Node& point = domain.getNode(nodeIndex);

        M_data[nodeIndex] = expression(point(0), point(1));
    }
}

void FunctionSpace::FunctionElement::setValue(int index, double value) {
    M_data[index] = value;
}

void FunctionSpace::FunctionElement::setValue(int elementId, int localIndex, double value) {
    M_data[M_functionSpace.connectivity(elementId, localIndex)] = value;
}

double FunctionSpace::FunctionElement::getValue(int id) const {
    return M_data[id];
}

double FunctionSpace::FunctionElement::getValue(int elemId, int localDof) const {
    return M_data[M_functionSpace.connectivity(elemId, localDof)];
}