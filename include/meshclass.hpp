#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <memory>
#include <functional>
#include <algorithm>
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <set>
#include <Eigen/Dense>

#include "bimap.hpp"

typedef std::function<double (double, double)> functionType;

/**
 * @brief Represents a node in the mesh.
 */
struct Node {

    // Attributes
    std::array<double, 2> coefficients; // x, y coordinates

    // Constructors
    Node(): coefficients({0.0, 0.0}) {}
    Node(double a, double b): coefficients({a, b}) {}
    Node(const std::array<double, 2>& c, int id): coefficients(c) {}

    // Methods
    /** @brief Returns the euclidean distance between two nodes.*/
    double distance(const Node& other) const; 

    // Operators
    double operator()(int) const;
    bool operator==(const Node&) const; // For debugging purposes only
};

// Output operator
std::ostream& operator<<(std::ostream& os, const Node& node);

/**
 * @brief Represents a triangle element in the mesh. Can only exist within a mesh
 */
struct TriangleElement {

    // Attributes
    std::array<int, 3> globalNodeIndex;  // global node ids

    // Constructor
    TriangleElement(const std::array<int, 3>& nodeId): globalNodeIndex(nodeId) {}
    // Operator
    int operator[](int localNodeId) const {
        return globalNodeIndex[localNodeId];
    }
};

/**
 * @brief Represents a facet in the mesh - elements edges
 */
struct Facet{

    // Attributes
    std::array<int, 2> globalNodeIndex; // global node index
    
    bool isBoundary;  // true if the facet is on the boundary of the mesh
    bool isSegment;   // true if the facet is an element segment in .msh

    // Constructors
    Facet(): globalNodeIndex({-1, -1}), isBoundary(false), isSegment(false) {}
    Facet(const std::array<int, 2>& node, bool flag) : 
                isBoundary(false),
                isSegment(flag),
                globalNodeIndex(node) {}

    int operator[](int localNodeId) const {
        return globalNodeIndex[localNodeId];
    }
};

/**
 * @brief Mesh class encapsulating finite element data structures
 */
class Mesh {
private:
    // Mesh properties
    unsigned int nbNodes;    // number of nodes
    unsigned int nbTriangles;// number of elements
    unsigned int nbFacets;   // number of facets

    // Containers (msh elements)
    std::vector<Node> tabNodes;                 // array of nodes
    std::vector<TriangleElement> tabTriangle;   // array of triangles
    std::vector<Facet> tabFacets;               // array of facets
    int nextFacetId;

    // Connectivity maps
    std::unordered_map<int, std::unordered_set<int>> nodeToElements; // nodeId -> set of elements
    std::unordered_map<int, std::unordered_set<int>> nodeToFacets;   // nodeId -> set of facets
    std::unordered_map<int, std::array<int, 2>> facetToElements;     // facetID -> set of elements
    std::unordered_map<int, std::array<int, 3>> elementToFacets;     // elementId -> set of facets

    // Region markers
    biMap<std::string, int> physicalMarkers;
    std::map<int, int> markerDimension;

    std::unordered_map<int, int> markedFacets;
    std::unordered_map<int, int> markedElements;

    std::set<int> boundaryNodes;

public:
    // Constructor
    Mesh(std::string);

    // Mapping utilities
    std::array<std::array<int, 2>, 3> localfacetToLocalNode = {{
        {1, 2}, 
        {2, 0},
        {0, 1}
    }};
    
    std::array<int, 2> getFaceNodes(const std::array<int, 3>& element, int localFacetId) {
        return {element[localfacetToLocalNode[localFacetId][0]], element[localfacetToLocalNode[localFacetId][1]]};
    }

    std::array<int, 2> getFaceNodes(const TriangleElement& element, int localFacetId) {
        return {element[localfacetToLocalNode[localFacetId][0]], element[localfacetToLocalNode[localFacetId][1]]};
    }

    std::array<int, 2> swapFaceNodes(const std::array<int, 2>& faceNodes) {
        std::array<int, 2> copy = faceNodes;
        if (faceNodes[0] > faceNodes[1]) {
            std::swap(copy[0], copy[1]);
        }
        return copy;
    }

        // ELEMENT RELATED METHODS

    /** @brief Returns the length of a facet given its identifier. */
    double getFacetLength(int identifier) const;
    /** @brief Returns the aera of a triangle given its identifier. */
    double getTriangleAera(int identifier) const;
    /** @brief Retruns the perimeter of a triangle given its identifier. */
    double getTrianglePerimeter(int identifier) const;

    const std::set<int>& getBoundary() const;

    /** @brief Adds a node to the mesh. */
    void addNode(const Node& node);
    /** @brief Adds an element to the mesh. */
    void addElement(const TriangleElement& element);
    /** @brief Returns the number of Nodes in the mesh. */
    int getNbNodes() const;
    /** @brief Returns the number of elements in the mesh. */
    int getNbElements() const;
    /** @brief Returns the number of facets in the mesh. */
    int getNbFacets() const;
    /** @brief Returns the number of segments marked facets in the mesh. */
    int getNbSegments() const;
    /** @brief Returns a reference to a node in the mesh given a vector index. */
    const Node& getNode(int identifier) const;
    /** @brief Returns a reference to an element in the mesh given a vector index. */
    const TriangleElement& getElement(int identifier) const;
    /** @brief Returns a reference to a facet in the mesh given a vector index. */
    const Facet& getFacet(int facetIndex) const;
    /** @brief Returns node Ids of a given element. */
    std::array<int, 3> getNodeFromElem(int identifier) const;
    /** @brief Returns node reference of a given element. */
    const Node& getNodeFromElem(int identifier, int localDof) const;

        // CONNECTIVITY RELATED

    int getNodeId(int nodeIndex) const;

    /** @brief Build connectivity maps for the mesh. */
    void buildConnectivity();
    /** @brief Returns a vector of element IDs that contain a given node. */
    std::vector<int> getElementsForNode(int nodeId) const;
    /** @brief Returns a vector of facet IDs that contain a given node. */
    std::vector<int> getFacetsForNode(int nodeId) const;
    /** @brief Returns a vector of element IDs that contain a given facet. */
    std::array<int, 2> getElementsForFacets(int facetId) const;
    /** @brief Returns a vector of element IDs that are neighbors of a given element. */
    std::array<int, 3> getNeighborTriangles(int triangleId) const;
    /** @brief Returns a boolean indicating whether a given facet is on the boundary. */
    bool isFacetOnBoundary(int facetId) const;
    /** @brief Returns a boolean indicating whether a given node is on the boundary. */
    bool isNodeOnBoundary(int nodeId) const;
    /** @brief Returns a boolean indicating whether a given element is on the boundary. */
    bool isTriangleOnBoundary(int triangleId) const;

        // FILTERS AND PHYSICAL NAMES

    /** @brief Prints the physical names and the number of elements associated. */
    void domainSummary() const;
    /** @brief Returns a vector of triangle elements marked with a physical name. */
    std::vector<int> getMarkedElements(const std::string&) const;
    /** @brief Returns a vector of facet marked with a physical name. */
    std::vector<int> getMarkedFacet(const std::string&) const;

    /** @brief Print nodes in (x, y) format. */
    //void printNodes() const;
    /** @brief Print triangle and their nodes. */
    //void printTriangles() const;
    /** @brief Print facets - triangle ids, node ids, and boundary or not. */
    //void printFacets() const;

    /** @brief Identifies boundary facets and sums their lengths. */
    double meshPerimeter() const;
    /** @brief Returns the mesh aera. */
    double meshAera() const;

    /** @brief Make a .vtk file of the mesh for plotting with optional functional. */
    int exportToVTK(const std::string& path, const std::string& plotName, const functionType& function = 0) const;
};

class FunctionSpace {
    public:
    
        FunctionSpace(const Mesh& mesh): M_domain(mesh) {
            M_globalDof = M_domain.getNbNodes();
        }
    
        unsigned int connectivity(int elementId, int localIndex) const;
        
        class FunctionElement {
        public:
            FunctionElement(const FunctionSpace& funcSpace) : M_functionSpace(funcSpace) {
                M_data.resize(funcSpace.M_globalDof);
            }
    
            /** @brief Sets all values to zero. */
            void setZero();
            /** @brief Initialize the elements data with an Eigen vector. */
            void initialize(const Eigen::VectorXd& dataVector);
            /** @brief Evaluate an expression over each dof. */
            void evaluate(const functionType& expression);
            /** @brief Sets a certain value given a node Id. */
            void setValue(int nodeId, double value);
            /** @brief Sets an element's node to a certain value. */
            void setValue(int elementId, int localIndex, double value);
            /** @brief Returns a value given a node Id. */
            double getValue(int nodeId) const;
            /** @brief Returns a value given an element Id and a local Dof. */
            double getValue(int elemId, int localDof) const;
    
        private:
            const FunctionSpace& M_functionSpace;
            std::vector<double> M_data;
        };
    
        FunctionElement element() const {return FunctionElement(*this);}
    
    private:
        const Mesh& M_domain;
        unsigned int M_globalDof;
    };

#endif