#ifndef CONFIGCLASS_H
#define CONFIGCLASS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <unordered_map>
#include <functional>

# define M_PI 3.14159265358979323846

typedef std::function<double (double, double)> functionType;

// Implementations for g

inline double constantFun_g(double x, double y) {
    return 1.0;
}

inline double quadraticFun_g(double x, double y) {
    return x * x + y * y;
}

inline double linearFun_g(double x, double y) {
    return x + y;
}

inline double zeroFun(double x, double y) {
    return 0.0;
}

// Implementations for f

inline double quadraticFun_f(double x, double y) {
    return -4.0;
}

// Implementation for solution grad

inline double dx_quadratic(double x, double y) {
    return 2 * x;
}

inline double dy_quadratic(double x, double y) {
    return 2 * y;
}

inline double grad_linear(double x, double y) {
    return 1.0;
}

// Map of all avalible function
const std::unordered_map<std::string, functionType> FunctionMap = {
            {"constant", &constantFun_g},
            {"quadratic_g", &quadraticFun_g},
            {"linear_g", &linearFun_g},
            {"zero", &zeroFun},
            {"quadratic_f", &quadraticFun_f},
            {"quadratic_dx", &dx_quadratic},
            {"quadratic_dy", &dy_quadratic}
};

void lineFormat(std::unordered_map<std::string, std::string>&, std::string);

struct ConfigClass{

    static std::string PATH; // To be edited if the relative path from build repo to meshes is changed

    std::string pathToMesh;  // relative path to /meshes from /build
    std::string flabel;      // name to be displayed in ParaView
    
    int integrationOrder;    // order of quadrature
    functionType functionF;  // source function
    functionType functionG;  // boundary condition

    functionType solution;    // exact solution

    functionType dx_solution; // solution grad_x
    functionType dy_solution; // solution grad_y

    std::unordered_map<std::string, functionType> FunMap;   // FunMap : function map previously defined

    /* Main config init; default if 
        - f constant,
        - g zero,
        - order = 1,
        - minimalistic mesh 4elt (square comprised of 4 isosceles triangles) */

    ConfigClass(){
        // Default parameters
        std::string Path = "square2d_4elt.msh";

        flabel = "constant";
        std::string meshPath(PATH);
        meshPath += Path;
        pathToMesh = meshPath;
        FunMap = FunctionMap;
        integrationOrder = 1;

        functionF = &constantFun_g;
        functionG = &zeroFun;

        // No explicit solution to this problem.
        solution = 0; 
        dx_solution = 0;
        dy_solution = 0;
    }

    ConfigClass(std::string configFile): ConfigClass() {
        std::ifstream ifile(configFile);
        std::unordered_map<std::string, std::string> fileIn;

        if (ifile){
            std::string line;
            
            while(std::getline(ifile, line)){

                if (line.empty() || line[0] == '#'){
                    continue;
                }

                else {
                    lineFormat(fileIn, line);
                }
            }

            pathToMesh = PATH + fileIn["mesh_file"];
            flabel = fileIn["label"];
            integrationOrder = std::stoi(fileIn["integration_order"]);
            FunMap = FunctionMap;

            auto it = FunMap.find(fileIn["sol"]);
            if (it != FunMap.end()) {
                solution = it->second;
            }

            it = FunMap.find(fileIn["dx_sol"]);
            if (it != FunMap.end()) {
                dx_solution = it->second;
            }

            it = FunMap.find(fileIn["dy_sol"]);
            if (it != FunMap.end()) {
                dy_solution = it->second;
            }

            it = FunMap.find(fileIn["source"]);
            if (it != FunMap.end()) {
                functionF = it->second;
            }

            it = FunMap.find(fileIn["boundary"]);
            if (it != FunMap.end()) {
                functionG = it->second;
            }

            std::cout << "Loaded config:\n";
            
            for (auto it = fileIn.begin(); it != fileIn.end(); it ++){
                std::cout << " - " << it->first << ": " << it->second << std::endl;
            }
        }

        else {
            std::cout << "Loading config file failed, resorting to default parameters\n";
        }
    }
};

void lineFormat(std::unordered_map<std::string, std::string>& set, std::string line){
    std::stringstream stream(line);
    std::string key, val, dump;

    stream >> key;
    stream >> dump;
    stream >> val;

    set[key] = val;
}

#endif