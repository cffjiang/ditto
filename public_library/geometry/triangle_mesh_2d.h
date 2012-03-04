//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// ditto/public_library/geometry/triangle_mesh_2d.h
// Copyright 2012, Chenfanfu Jiang
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

#ifndef DITTO_PUBLIC_LIBRARY_GEOMETRY_TRIANGLE_MESH_2D_H
#define DITTO_PUBLIC_LIBRARY_GEOMETRY_TRIANGLE_MESH_2D_H

#include <cmath>
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <string>
#include <vector>

#include <ditto/public_library/algebra/linear_algebra.h>
#include <ditto/public_library/geometry/simplex.h>

namespace ditto { namespace geometry {

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Class: Triangle_Mesh_2d
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T>
class Triangle_Mesh_2d {
public:
    typedef ditto::algebra::VECTOR_2D<T> node_type;
    typedef typename std::vector<node_type> node_list_type;
    typedef ditto::algebra::VECTOR_3D<int> element_type;
    typedef typename std::vector<element_type> element_list_type;
    typedef typename std::vector<T> rho_list_type;
    typedef typename std::vector<T> area_list_type;

    typedef ditto::geometry::Triangle_2d<T> simplex_type;

    int Nx;
    int Ny;
    T dx;
    T dy;
    T xmin;
    T xmax;
    T ymin;
    T ymax;

    node_list_type nodes;
    element_list_type elements;
    rho_list_type rho;
    area_list_type area;

    Triangle_Mesh_2d() {
        initialize_regular_mesh(3, 4, 0, 1, 0, 1);
    }
    void initialize_regular_mesh(const int m, const int n, const T input_xmin, const T input_xmax, const T input_ymin, const T input_ymax, const T input_rho = 50);
};

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Function: initialize_regular_mesh
// Description: Initialize a regular rectangle triangle mesh.
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T>
void Triangle_Mesh_2d<T>::initialize_regular_mesh(const int input_Nx, const int input_Ny, const T input_xmin, const T input_xmax, const T input_ymin, const T input_ymax, const T input_rho)
{
    std::cout << "Initializing regular rectangle triangle mesh...\n";

    Nx = input_Nx;
    Ny = input_Ny;
    xmin = input_xmin;
    xmax = input_xmax;
    ymin = input_ymin;
    ymax = input_ymax;
    dx = (xmax - xmin) / (Nx - 1);
    dy = (ymax - ymin) / (Ny - 1);

    // build nodes
    for (int j=0; j<Ny; j++) {
        for (int i=0; i<Nx; i++) {
            node_type current_node(xmin + dx*i, ymin + dy*j);
            nodes.push_back(current_node);
        }
    }

    // build elements
    for (int j=0; j<Ny-1; j++) {
        for (int i=0; i<Nx-1; i++) {
            element_type first_element(Nx*j+i, Nx*(j+1)+i+1, Nx*(j+1)+i);
            element_type second_element(Nx*j+i, Nx*j+i+1, Nx*(j+1)+i+1);
            elements.push_back(first_element);
            elements.push_back(second_element);
        }
    }

    // build rho and area
    for (unsigned int i=0; i<elements.size(); i++) {
        rho.push_back(input_rho);
        
        T positions[6] = { nodes[elements[i](0)](0), nodes[elements[i](0)](1), nodes[elements[i](1)](0), nodes[elements[i](1)](1), nodes[elements[i](2)](0), nodes[elements[i](2)](1) };
        simplex_type tri(positions);
        area.push_back(tri.get_area());
    }

    // debug
    for (unsigned int i=0; i<nodes.size(); i++) {
        std::cout << "\nnode: " << nodes[i](0) << " " << nodes[i](1) << std::endl;
    }

    for (unsigned int i=0; i<elements.size(); i++) {
        std::cout << "\nelement: " << elements[i](0) << " " << elements[i](1) << " " << elements[i](2) << std::endl;
        std::cout << "rho: " << rho[i] << std::endl;
        std::cout << "area: " << area[i] << std::endl;
    }
        

}


} } // end namespaces

#endif
