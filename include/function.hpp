#ifndef FUNC_H
#define FUNC_H

#include <iostream>
#include <vector>

#include "meshclass.hpp"

typedef std::function<double (double, double)> functionType;

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
        /** @brief Evaluate an expression over each dof. */
        void evaluate(functionType expression);
        /** @brief Sets a certain value given a node Id. */
        void setValue(int nodeId, double value);
        /** @brief Sets an element's node to a certain value. */
        void setValue(int elementId, int localIndex, double value);
        /** @brief Returns a value given a node Id. */
        double getValue(int nodeId) const;

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