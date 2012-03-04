//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// ditto/public_library/geometry/triangle_mesh_3d.h
// Copyright 2012, Chenfanfu Jiang
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

#ifndef DITTO_PUBLIC_LIBRARY_GEOMETRY_TRIANGLE_MESH_3D_H
#ifndef DITTO_PUBLIC_LIBRARY_GEOMETRY_TRIANGLE_MESH_3D_H

#include <cmath>
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <string>
#include <vector>

#include <ditto/public_library/algebra/linear_algebra.h>

namespace ditto { namespace geometry {

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Class: Triangle_Mesh_3d
// Description: Mainly used for 3d surfaces such as cloths.
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T>
class Triangle_Mesh_3d {
public:
    typedef ditto::algebra::VECTOR_3D<T> node_type;
    typedef typename std::vector<node_type> node_list_type;
    typedef ditto::algebra::VECTOR_3D<int> element_type;
    typedef typename std::vector<element_type> element_list_type;
    typedef std::vector<T> rho_list_type;
    typedef std::vector<T> area_list_type;

    int Nx;
    int Ny;
    node_list_type nodes;
    element_list_type elements;
    rho_list_type rho;
    area_list_type area;

    Triangle_Mesh_3d() {}
    void initialize_regular_mesh(const int m, const int n, const T input_rho = 50);
};

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Function: initialize_regular_mesh
// Description: Initialize a regular triangle mesh in x-y plane.
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T>
void Triangle_Mesh_3d<T>::initialize_regular_mesh(const int input_Nx, const int input_Ny, const T xmin, const T xmax, const T ymin, const T ymax, const T input_rho)
{
    Nx = input_Nx;
    Ny = input_Ny;

    for (int i=0; i<Nx; i++) {
        for (int j=0; j<Ny; j++) {

        }
    }

}




} } // end namespaces
