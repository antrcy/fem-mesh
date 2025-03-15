#include "function.hpp"

// FUNCTION SPACE

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

void FunctionSpace::FunctionElement::evaluate(functionType expression) {
    const Mesh& domain = M_functionSpace.M_domain;

    for (const auto& index : domain.idToIndexNodes){
        const Node& point = domain.getNode(index.first);

        M_data[index.second] = expression(point(0), point(1));
    }
}

void FunctionSpace::FunctionElement::setValue(int id, double value) {
    M_data[M_functionSpace.M_domain.idToIndexNodes.at(id)] = value;
}

void FunctionSpace::FunctionElement::setValue(int elementId, int localIndex, double value) {
    M_data[M_functionSpace.connectivity(elementId, localIndex)] = value;
}

double FunctionSpace::FunctionElement::getValue(int id) const {
    return M_data[M_functionSpace.M_domain.idToIndexNodes.at(id)];
}