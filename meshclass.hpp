#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <memory>
#include <functional>
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>

/**
 * @brief Represents a node in the mesh.
 */
struct Node{

    // Attributes
    std::array<double, 2> coefficients; // x, y coordinates

    // Constrcutors
    Node(): coefficients({0.0, 0.0}) {}
    Node(double a, double b): coefficients({a, b}) {}
    Node(const std::array<double, 2>& c): coefficients(c) {}

    // Methods
    /** @brief Returns the euclidean distance between two nodes.*/
    double distance(const Node& other) const; 
};

// Output operator
std::ostream& operator<<(std::ostream& os, const Node& node);


/**
 * @brief Represents a triangle element in the mesh
 */
struct TriangleElements {

    // Attributes
    std::array<int, 3> globalNodeId;  // global node ids
    const std::vector<Node>* globalNodeTab; // pointer to the global node table

    // Constructor
    TriangleElements(const std::array<int, 3>& nodeId, const std::vector<Node>& nodeTab): globalNodeId(nodeId), globalNodeTab(&nodeTab) {}

    // Operator
    int operator[](int localNodeId) const {
        return globalNodeId[localNodeId];
    }

};

// Output operator
std::ostream& operator<<(std::ostream& os, const TriangleElements& elem);


/**
 * @brief Represents a facet in the mesh - elements edges
 */
struct Facet{

    // Attributes
    const std::vector<Node>* globalNodeTab; // pointer to mesh node tab
    std::array<int, 2> globalNodeId;        // global node ids
    std::array<int, 2> globalElemId;        // global element ids
    int identifier;  // unique identifier in the mesh
    bool isBoundary; // true if the facet is on the boundary of the mesh

    // Constructors
    Facet(): globalNodeId({-1, -1}), globalElemId({-1, -1}), identifier(-1), isBoundary(false), globalNodeTab(nullptr) {}
    Facet(const std::array<int, 2>& node, const std::array<int, 2>& elem, const std::vector<Node>& nodeTab, int id) : 
                identifier(id), 
                isBoundary(false),
                globalNodeTab(&nodeTab)

    {
        globalNodeId = node;
        globalElemId = elem;
        isBoundary = (globalElemId[1] == -1) || (globalElemId[0] == -1);
    }

    // Accessors and mutators
    void setElemId(int local, int global) {
        globalElemId[local] = global;
        isBoundary = (globalElemId[1] == -1) || (globalElemId[0] == -1);
    }

    // Methods
    double length() const;

};

// Output operator
std::ostream& operator<<(std::ostream& os, const Facet& facet);


/**
 * @brief Mesh class encapsulating finite element data structures
 */
class Mesh{
public:
    // Mesh properties
    unsigned int nbNodes;    // number of nodes
    unsigned int nbElements; // number of elements
    unsigned int nbFacets;   // number of facets

    // Containers
    std::vector<Node> tabNodes;                 // array of nodes
    std::vector<TriangleElements> tabElements;  // array of elements
    std::vector<Facet> tabFacets;               // array of facets

    // Connectivity maps
    std::unordered_map<int, std::unordered_set<int>> nodeToElements; // nodeId -> set of elements
    std::unordered_map<int, std::unordered_set<int>> nodeToFacets;   // nodeId -> set of facets
    std::unordered_map<int, std::array<int, 3>> elementToFacets;     // elementId -> set of facets

public:
    // Constructor
    Mesh(std::string);

    // Mapping utilities
    std::array<std::array<int, 2>, 3> facetToNode = {{
        {1, 2}, 
        {2, 0},
        {0, 1}
    }};
    
    std::array<int, 2> getFaceNodes(const std::array<int, 3>& element, int localFacetId) {
        return {element[facetToNode[localFacetId][0]], element[facetToNode[localFacetId][1]]};
    }

    std::array<int, 2> getFaceNodes(const TriangleElements& element, int localFacetId) {
        return {element[facetToNode[localFacetId][0]], element[facetToNode[localFacetId][1]]};
    }

    void swapFaceNodes(std::array<int, 2>& faceNodes) {
        if (faceNodes[0] > faceNodes[1]) {
            std::swap(faceNodes[0], faceNodes[1]);
        }
    }

    // Methods
    /** @brief Adds a node to the mesh. */
    void addNode(const Node& node);
    /** @brief Adds an element to the mesh. */
    void addElement(const TriangleElements& element);
    /** @brief Returns the number of Nodes in the mesh. */
    int getNbNodes() const;
    /** @brief Returns the number of elements in the mesh. */
    int getNbElements() const;
    /** @brief Returns the number of facets in the mesh. */
    int getNbFacets() const;
    /** @brief Returns a reference to a node in the mesh. */
    const Node& getNode(int nodeId) const;
    /** @brief Returns a reference to an element in the mesh. */
    const TriangleElements& getElement(int elementId) const;
    /** @brief Returns a reference to a facet in the mesh. */
    const Facet& getFacet(int facetId) const;
    /** @brief Build connectivity maps for the mesh. */
    void buildConnectivity();
    /** @brief Returns a vector of element IDs that contain a given node. */
    std::vector<int> getElementsForNode(int nodeId) const;
    /** @brief Returns a vector of facet IDs that contain a given node. */
    std::vector<int> getFacetsForNode(int nodeId) const;
    /** @brief Returns a vector of element IDs that contain a given facet. */
    std::vector<int> getElementsForFacets(int facetId) const;
    /** @brief Returns a vector of element IDs that are neighbors of a given element. */
    std::vector<int> getNeighborTriangles(int triangleId) const;
    /** @brief Returns a boolean indicating whether a given facet is on the boundary. */
    bool isFacetOnBoundary(int facetId) const;
    /** @brief Returns a boolean indicating whether a given node is on the boundary. */
    bool isNodeOnBoundary(int nodeId) const;
    /** @brief Returns a boolean indicating whether a given element is on the boundary. */
    bool isTriangleOnBoundary(int triangleId) const;
    /** @brief Print nodes in (x, y) format. */
    void printNodes() const;
    /** @brief Print triangle and their nodes. */
    void printTriangles() const;
    /** @brief Print facets - triangle ids, node ids, and boundary or not. */
    void printFacets() const;
    /** @brief Identifies boundary facets and sums their lengths. */
    double perimeter() const;
};