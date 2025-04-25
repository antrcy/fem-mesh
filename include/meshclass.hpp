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

    // Methods
    /** @brief Returns the euclidean distance between two nodes.*/
    double distance(const Node& other) const; 

    // Operators
    double operator()(int) const;
    bool operator==(const Node&) const; // For debugging purposes only
    friend std::ostream &operator<<(std::ostream&, const Node&); 
};

/**
 * @brief Represents a triangle element in the mesh. Can only exist within a mesh.
 */
struct TriangleElement {

    // Attributes
    std::array<int, 3> globalNodeIndex;  // global node ids

    // Constructor
    TriangleElement(const std::array<int, 3>& nodeIndex): globalNodeIndex(nodeIndex) {}
    // Operator
    int operator[](int localNodeIndex) const {
        return globalNodeIndex[localNodeIndex];
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

    int operator[](int localDof) const {
        return globalNodeIndex[localDof];
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

    // Connectivity maps
    std::unordered_map<int, std::unordered_set<int>> nodeToElements; // nodeIndex -> set of elements
    std::unordered_map<int, std::unordered_set<int>> nodeToFacets;   // nodeIndex -> set of facets
    std::unordered_map<int, std::array<int, 2>> facetToElements;     // facetIndex -> set of elements
    std::unordered_map<int, std::array<int, 3>> elementToFacets;     // elementIndex -> set of facets

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
    
    std::array<int, 2> getFaceNodes(const std::array<int, 3>& element, int localDof) {
        return {element[localfacetToLocalNode[localDof][0]], element[localfacetToLocalNode[localDof][1]]};
    }
    
    std::array<int, 2> getFaceNodes(const TriangleElement& element, int localFacetId) {
        return {element[localfacetToLocalNode[localFacetId][0]], element[localfacetToLocalNode[localFacetId][1]]};
    }
    /** @brief Swap two array coefficients so that they are ordered numerically. */
    std::array<int, 2> swapFaceNodes(const std::array<int, 2>& faceNodes) {
        std::array<int, 2> copy = faceNodes;
        if (faceNodes[0] > faceNodes[1]) {
            std::swap(copy[0], copy[1]);
        }
        return copy;
    }

        // GEOMETRIC TRANSFORMATION

    /**
     * @brief Implementation of the affine transformation from reference element to a real element.
     */
    class GeometricTransformation {
    public:
        GeometricTransformation(const Mesh& mesh): domain(mesh) {
            x0 = Eigen::Vector2d({0.0, 0.0});
            x1 = Eigen::Vector2d({1.0, 0.0});
            x2 = Eigen::Vector2d({0.0, 1.0});
        }

        GeometricTransformation(const Mesh& mesh, int elem): domain(mesh) {
            x0 = Eigen::Vector2d(domain.getNodeFromElem(elem, 0).coefficients.data());
            x1 = Eigen::Vector2d(domain.getNodeFromElem(elem, 1).coefficients.data());
            x2 = Eigen::Vector2d(domain.getNodeFromElem(elem, 2).coefficients.data());
        }

        /** @brief Switch which element the transformation points to. */
        void initialize(int elem) {
            x0 = Eigen::Vector2d(domain.getNodeFromElem(elem, 0).coefficients.data());
            x1 = Eigen::Vector2d(domain.getNodeFromElem(elem, 1).coefficients.data());
            x2 = Eigen::Vector2d(domain.getNodeFromElem(elem, 2).coefficients.data());
        }

        /** @brief Maps a given point in reference barycentring coordinates to its real element tranformation. */
        Eigen::Vector2d map(const Eigen::Vector2d& coord) const {
            return x0 * (1 - coord(0) - coord(1)) + x1 * (coord(0)) + x2 * (coord(1));
        }

        /** @brief Returns the jacobian of the affine transformation from reference to real. */
        Eigen::Matrix2d jacobian() const {
            return Eigen::Matrix2d({{-x0(0) + x1(0), -x0(0) + x2(0)},
                                    {-x0(1) + x1(1), -x0(1) + x2(1)}});
        }

        /** @brief Returns the inverse jacobian of the affine transformation from reference to real. */
        Eigen::Matrix2d jacobianInverse() const {
            return jacobian().inverse();
        }

        /** @brief Returns the determinant of the jacobian of the affine transformation from reference to real. */
        double jacobianDeterminant() const {
            return jacobian().determinant();
        }

    private:
        const Mesh& domain;
        Eigen::Vector2d x0, x1, x2;
    };

        // ELEMENT RELATED METHODS

    /** @brief Returns the length of a facet given its index. */
    double getFacetLength(int index) const;
    /** @brief Returns the aera of a triangle given its index. */
    double getTriangleAera(int index) const;
    /** @brief Retruns the perimeter of a triangle given its index. */
    double getTrianglePerimeter(int index) const;
    /** @brief Returns a reference to the set of boundary index. */
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
    const Node& getNode(int index) const;
    /** @brief Returns a reference to an element in the mesh given a vector index. */
    const TriangleElement& getElement(int index) const;
    /** @brief Returns a reference to a facet in the mesh given a vector index. */
    const Facet& getFacet(int facetIndex) const;
    /** @brief Returns node indexes of a given element. */
    std::array<int, 3> getNodeFromElem(int index) const;
    /** @brief Returns node reference of a given element. */
    const Node& getNodeFromElem(int index, int localDof) const;

        // CONNECTIVITY RELATED

    /** @brief Build connectivity maps for the mesh. */
    void buildConnectivity();
    /** @brief Returns a vector of element IDs that contain a given node. */
    std::vector<int> getElementsForNode(int nodeIndex) const;
    /** @brief Returns a vector of facet IDs that contain a given node. */
    std::vector<int> getFacetsForNode(int nodeIndex) const;
    /** @brief Returns a vector of element IDs that contain a given facet. */
    std::array<int, 2> getElementsForFacets(int facetIndex) const;
    /** @brief Returns a vector of element IDs that are neighbors of a given element. */
    std::array<int, 3> getNeighborTriangles(int triangleIndex) const;
    /** @brief Returns a boolean indicating whether a given facet is on the boundary. */
    bool isFacetOnBoundary(int facetIndex) const;
    /** @brief Returns a boolean indicating whether a given node is on the boundary. */
    bool isNodeOnBoundary(int nodeIndex) const;
    /** @brief Returns a boolean indicating whether a given element is on the boundary. */
    bool isTriangleOnBoundary(int triangleIndex) const;

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
    /** @brief Reference element testing */
    double meshAera_ref() const;

    /** @brief Make a .vtk file of the mesh for plotting with optional functional. */
    int exportToVTK(const std::string& path, const std::string& plotName, const functionType& function = 0) const;
};

/**
 * @brief FunctionSpace to handle data over a mesh.
 */
class FunctionSpace {
    public:
    
        FunctionSpace(const Mesh& mesh): domain(mesh) {
            globalDofs = domain.getNbNodes();
        }
        
        /** @brief Provides mapping between elements and global node indexes. */
        unsigned int connectivity(int elemIndex, int localDof) const;
        
        /**
         * @brief Function element within a function space
         */ 
        class FunctionElement {
        public:
            FunctionElement(const FunctionSpace& funcSpace) : fspace(funcSpace) {
                data.resize(funcSpace.globalDofs);
            }
    
            /** @brief Sets all values to zero. */
            void setZero();
            /** @brief Initialize the elements data with an Eigen vector. */
            void initialize(const Eigen::VectorXd& dataVector);
            /** @brief Evaluate an expression over each dof. */
            void evaluate(const functionType& expression);
            /** @brief Sets a certain value given a node Index. */
            void setValue(int nodeIndex, double value);
            /** @brief Sets an element's node to a certain value. */
            void setValue(int elemIndex, int localDof, double value);
            /** @brief Returns a value given a node Index. */
            double getValue(int nodeIndex) const;
            /** @brief Returns a value given an element Id and a local Dof. */
            double getValue(int elemIndex, int localDof) const;
    
        private:
            const FunctionSpace& fspace;
            std::vector<double> data;
        };
        
        /** Returns a function element within function space. */
        FunctionElement element() const {return FunctionElement(*this);}
    
    private:
        const Mesh& domain;
        unsigned int globalDofs;
    };

#endif