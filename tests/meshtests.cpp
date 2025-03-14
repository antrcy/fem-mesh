#include "meshclass.hpp"
#include <cmath>
#include <gtest/gtest.h>


class MeshTest : public ::testing::Test {
protected:
    Mesh* mesh_trivial;
    Mesh* mesh_complex;

    void SetUp() override {
        mesh_trivial = new Mesh("../meshes/square2d_4elt.msh");
        mesh_complex = new Mesh("../meshes/square2d_perforated.msh");
    }

    void TearDown() override {
        delete mesh_trivial;
        delete mesh_complex;
    }
};

TEST_F(MeshTest, LoadValidMshFile) {
    EXPECT_GT(mesh_trivial->getNbNodes(), 5);
    EXPECT_GT(mesh_trivial->getNbElements(), 4);
    EXPECT_GT(mesh_trivial->getNbFacets(), 4);
    EXPECT_GT(mesh_trivial->getNbSegments(), 4);

    EXPECT_GT(mesh_complex->getNbNodes(), 21900);
    EXPECT_GT(mesh_complex->getNbElements(), 42184);
    EXPECT_GT(mesh_complex->getNbFacets(), 1674);
    EXPECT_GT(mesh_complex->getNbSegments(), 1674);
}

TEST_F(MeshTest, CorrectMeshArea) {
    double expectedAreaTrivial = 1.0;
    EXPECT_NEAR(mesh_trivial->meshAera(), expectedAreaTrivial, 1e-14);

    double expectedAreaComplex = 0.837725;
    EXPECT_NEAR(mesh_complex->meshAera(), expectedAreaTrivial, 1e-6);
}

TEST_F(MeshTest, CorrectMeshPerimeter) {
    double expectedPerimeterTrivial = 4.0; // Replace with actual known value
    EXPECT_NEAR(mesh_trivial->meshPerimeter(), expectedPerimeterTrivial, 1e-14);

    double expectedPerimeterComplex = 11.5086;
    EXPECT_NEAR(mesh_complex->meshPerimeter(), expectedPerimeterComplex, 1e-6);
}

TEST_F(MeshTest, GetNode) {
    int testNodeId = 3;
    Node testNode = mesh_trivial->getNode(testNodeId);
    EXPECT_TRUE(Node({1, 1}, -1) == testNode);
}

// Test retrieving an element by ID
TEST_F(MeshTest, GetElement) {
    int testElementId = 0;
    EXPECT_NO_THROW(mesh->getElement(testElementId));
}

// Test retrieving a facet by ID
TEST_F(MeshTest, GetFacet) {
    int testFacetId = 0;
    EXPECT_NO_THROW(mesh->getFacet(testFacetId));
}

// Test individual triangle area calculation
TEST(TriangleElementTest, TriangleAreaCalculation) {
    std::vector<Node> nodes = {
        Node(0, 0, 0), Node(1, 0, 1), Node(0, 1, 2)
    };
    TriangleElement triangle({0, 1, 2}, nodes, 0);
    EXPECT_NEAR(triangle.getAera(), 0.5, 1e-6);
}

// Test individual triangle perimeter calculation
TEST(TriangleElementTest, TrianglePerimeterCalculation) {
    std::vector<Node> nodes = {
        Node(0, 0, 0), Node(1, 0, 1), Node(0, 1, 2)
    };
    TriangleElement triangle({0, 1, 2}, nodes, 0);
    EXPECT_NEAR(triangle.getPerimeter(), 2.0 + std::sqrt(2.0), 1e-6);
}

// Test connectivity map construction
TEST_F(MeshTest, BuildConnectivity) {
    EXPECT_NO_THROW(mesh_complex->buildConnectivity());
}

// Test exporting to VTK
TEST_F(MeshTest, ExportToVTK) {
    EXPECT_EQ(mesh_complex->exportToVTK("test_output.vtk", "testMesh"), 0);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
