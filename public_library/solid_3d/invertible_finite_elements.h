//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// ditto/public_library/solid_3d/invertible_finite_elements.h
// Copyright 2012, Chenfanfu Jiang
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

#ifndef DITTO_PUBLIC_LIBRARY_SOLID_3D_INVERTIBLE_FINITE_ELEMENTS_H
#define DITTO_PUBLIC_LIBRARY_SOLID_3D_INVERTIBLE_FINITE_ELEMENTS_H

#include <cmath>
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <string>

#include <ditto/public_library/algebra/linear_algebra.h>

namespace ditto { namespace solid_3d {

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Class: Invertible_Finite_Elements
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T, class MeshType>
class Invertible_Finite_Elements {
private:
    MeshType &mesh;
    
    //###################################################
    // elasticity_model:
    //                1:    neo-hookean
    //###################################################
    int elasticity_model;
    
    //###################################################
    // time_integration_scheme:
    //                1:    forward euler
    //###################################################    
    int time_integration_scheme;

    T dt;

public:
    Invertible_Finite_Elements(const std::string &input_elasticity_model, const std::string &input_time_integration_scheme, const MeshType& input_mesh, const T input_dt) {
        if (input_elasticity_model == "neo-hookean") 
            elasticity_model = 1;
        else { 
            std::cout << "ERROR: elasticity model not supported!\n";
            exit(0);
        }

        if (input_time_integration_scheme == "forward euler")
            time_integration_scheme = 1;
        else {
            std::cout << "ERROR: time integration scheme not supported!\n";
            exit(0);
        }

        mesh = input_mesh;
        dt = input_dt;
    }





};





} } // end namespaces
#endif
