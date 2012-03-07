//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// ditto/public_library/geometry/triangle_mesh_3d.h
// Copyright 2012, Chenfanfu Jiang
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

#ifndef DITTO_PUBLIC_LIBRARY_GEOMETRY_TRIANGLE_MESH_3D_H
#define DITTO_PUBLIC_LIBRARY_GEOMETRY_TRIANGLE_MESH_3D_H

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
// Class: Triangle_Mesh_3d
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T>
class Triangle_Mesh_3d {
public:
    typedef ditto::algebra::VECTOR_3D<T> node_type;
    typedef typename std::vector<node_type> node_list_type;
    typedef ditto::algebra::VECTOR_3D<int> element_type;
    typedef typename std::vector<element_type> element_list_type;
    typedef typename std::vector<T> rho_list_type;
    typedef typename std::vector<T> area_list_type;
    typedef typename std::vector<T> mass_list_type;

    typedef ditto::geometry::Triangle_3d<T> simplex_type;

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
    mass_list_type mass;

    Triangle_Mesh_3d (const int m, const int n, const T input_xmin, const T input_xmax, const T input_ymin, const T input_ymax, const T input_rho = 1) {
        initialize_regular_mesh(m, n, input_xmin, input_xmax, input_ymin, input_ymax, input_rho); }

    void initialize_regular_mesh(const int m, const int n, const T input_xmin, const T input_xmax, const T input_ymin, const T input_ymax, const T input_rho = 1);
};

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Function: initialize_regular_mesh
// Description: Initialize a regular rectangle triangle mesh in x-y plane.
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T>
void Triangle_Mesh_3d<T>::initialize_regular_mesh(const int input_Nx, const int input_Ny, const T input_xmin, const T input_xmax, const T input_ymin, const T input_ymax, const T input_rho)
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
            node_type current_node(xmin + dx*i, ymin + dy*j, 0.0);
            nodes.push_back(current_node); } }

    // build elements
    for (int j=0; j<Ny-1; j++) {
        for (int i=0; i<Nx-1; i++) {
            element_type first_element(Nx*j+i, Nx*(j+1)+i+1, Nx*(j+1)+i);
            element_type second_element(Nx*j+i, Nx*j+i+1, Nx*(j+1)+i+1);
            elements.push_back(first_element);
            elements.push_back(second_element); }}

    // build rho and area
    for (unsigned int i=0; i<elements.size(); i++) {
        rho.push_back(input_rho);
        T positions[9] = { nodes[elements[i](0)](0), nodes[elements[i](0)](1), nodes[elements[i](0)](2),  
                           nodes[elements[i](1)](0), nodes[elements[i](1)](1), nodes[elements[i](1)](2), 
                           nodes[elements[i](2)](0), nodes[elements[i](2)](1), nodes[elements[i](2)](2) };
        simplex_type tri(positions);
        area.push_back(tri.get_area()); }

    // average rho to nodes -> get mass
    mass.resize(nodes.size(), 0.0);
    T total_mass = 0;
    for (unsigned int i=0; i<elements.size(); i++) {
        total_mass += rho[i]*area[i]; }
    for (unsigned int i=0; i<nodes.size(); i++) {
        mass[i] = total_mass / nodes.size(); }
    
    
    // debug code
    /*
    for (unsigned int i=0; i<nodes.size(); i++) {
    std::cout << "\nnode: " << nodes[i](0) << " " << nodes[i](1) << " " << node[i](2) << std::endl;
    }

    for (unsigned int i=0; i<elements.size(); i++) {
        std::cout << "\nelement: " << elements[i](0) << " " << elements[i](1) << " " << elements[i](2) << std::endl;
        std::cout << "rho: " << rho[i] << std::endl;
        std::cout << "area: " << area[i] << std::endl;
    }
    */    

}


} } // end namespaces

#endif
