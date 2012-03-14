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
    typedef ditto::algebra::VECTOR_2D<int> edge_type;
    typedef typename std::vector<edge_type> edge_list_type;
    
    typedef typename std::vector<T> rho_list_type;
    typedef typename std::vector<T> area_list_type;
    typedef typename std::vector<T> mass_list_type;

    typedef ditto::geometry::Triangle_3d<T> simplex_type;

    int Nx;
    int Ny;
    int Nz;
    T dx;
    T dy;
    T dz;
    T xmin;
    T xmax;
    T ymin;
    T ymax;
    T zmin;
    T zmax;

    node_list_type nodes;
    element_list_type elements;
    edge_list_type edges;
    rho_list_type rho;
    area_list_type area;
    mass_list_type mass;
    
    Triangle_Mesh_3d(){}

    void initialize_regular_mesh(const int input_Nx, const int input_Ny, const T input_xmin, const T input_xmax, const T input_ymin, const T input_ymax, const T input_rho = 1);

    void add_cloth_to_existing_mesh(const int input_Nx, const int input_Ny, const T input_xmin, const T input_xmax, const T input_ymin, const T input_ymax, const T input_rho = 1);

    void initialize_parellel_clothes(const int num_clothes, const int input_Nx, const int input_Ny, const T input_xmin, const T input_xmax, const T input_ymin, const T input_ymax, const T input_zmin, const T input_zmax, const T input_rho = 1);
};

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Function: initialize_regular_mesh
// Description: Initialize a regular rectangle triangle mesh in x-y plane.
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T>
void Triangle_Mesh_3d<T>::
initialize_regular_mesh(const int input_Nx, const int input_Ny, const T input_xmin, const T input_xmax, const T input_ymin, const T input_ymax, const T input_rho)
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

    // build edges
    for (int j=0; j<Ny; j++) {
        for (int i=0; i<Nx-1; i++) {
            edge_type current_edge(Nx*j+i, Nx*j+i+1);
            edges.push_back(current_edge); }}
    for (int i=0; i<Nx; i++) {
        for (int j=0; j<Ny-1; j++) {
            edge_type current_edge(Nx*j+i, Nx*(j+1)+i);
            edges.push_back(current_edge); }}
    for (int i=0; i<Nx-1; i++) {
        for (int j=0; j<Ny-1; j++) {
            edge_type current_edge(Nx*j+i, Nx*(j+1)+i+1);
            edges.push_back(current_edge); }}

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
    for (unsigned int i=0; i<elements.size(); i++) {
        int node1 = elements[i](0);
        int node2 = elements[i](1);
        int node3 = elements[i](2);
        mass[node1] += rho[i]*area[i] / 3.0;; 
        mass[node2] += rho[i]*area[i] / 3.0;; 
        mass[node3] += rho[i]*area[i] / 3.0;; 
    }

    // debug code
    /*
    for (unsigned int i=0; i<nodes.size(); i++) {
        std::cout << "\nnode: " << nodes[i](0) << " " << nodes[i](1) << " " << nodes[i](2) << std::endl;
    }
    
    for (unsigned int i=0; i<elements.size(); i++) {
        std::cout << "\nelement: " << elements[i](0) << " " << elements[i](1) << " " << elements[i](2) << std::endl;
    }

    for (unsigned int i=0; i<edges.size(); i++) {
        std::cout << "\nedge: " << edges[i](0) << " " << edges[i](1) << std::endl;
    }
    */    

}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Function: add_cloth_to_existing_mesh
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T>
void Triangle_Mesh_3d<T>::
add_cloth_to_existing_mesh(const int input_Nx, const int input_Ny, const T input_xmin, const T input_xmax, const T input_ymin, const T input_ymax, const T input_rho)
{

}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Function: initialize_parellel_clothes
// Description: Initialize multiple regular rectangle meshes, each one is on
//              x-y plane. 
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T>
void Triangle_Mesh_3d<T>::
initialize_parellel_clothes(const int num_clothes, const int input_Nx, const int input_Ny, const T input_xmin, const T input_xmax, const T input_ymin, const T input_ymax, const T input_zmin, const T input_zmax, const T input_rho)
{
    std::cout << "Initializing multiple clothes...\n";

    Nx = input_Nx;
    Ny = input_Ny;
    Nz = num_clothes;

    xmin = input_xmin;
    xmax = input_xmax;
    ymin = input_ymin;
    ymax = input_ymax;
    zmin = input_zmin;
    zmax = input_zmax;
    dx = (xmax - xmin) / (Nx - 1);
    dy = (ymax - ymin) / (Ny - 1);
    dz  =(zmax - zmin) / (Nz - 1);

    for (int cloth_index = 0; cloth_index < num_clothes; cloth_index++) {
        // compute z coordinate
        T z_position = zmin + dz*cloth_index;

        // build nodes
        for (int j=0; j<Ny; j++) {
            for (int i=0; i<Nx; i++) {
                node_type current_node(xmin + dx*i, ymin + dy*j, z_position);
                nodes.push_back(current_node); } }

        // build elements
        for (int j=0; j<Ny-1; j++) {
            for (int i=0; i<Nx-1; i++) {
                element_type first_element(Nx*Ny*cloth_index + Nx*j+i, Nx*Ny*cloth_index + Nx*(j+1)+i+1, Nx*Ny*cloth_index + Nx*(j+1)+i);
                element_type second_element(Nx*Ny*cloth_index + Nx*j+i, Nx*Ny*cloth_index + Nx*j+i+1, Nx*Ny*cloth_index + Nx*(j+1)+i+1);
                elements.push_back(first_element);
                elements.push_back(second_element); }}
     
        // build edges
        int shift = Nx*Ny*cloth_index;
        for (int j=0; j<Ny; j++) {
            for (int i=0; i<Nx-1; i++) {
                edge_type current_edge(shift+Nx*j+i, shift+Nx*j+i+1);
                edges.push_back(current_edge); }}
        for (int i=0; i<Nx; i++) {
            for (int j=0; j<Ny-1; j++) {
                edge_type current_edge(shift+Nx*j+i, shift+Nx*(j+1)+i);
                edges.push_back(current_edge); }}
        for (int i=0; i<Nx-1; i++) {
            for (int j=0; j<Ny-1; j++) {
                edge_type current_edge(shift+Nx*j+i, shift+Nx*(j+1)+i+1);
                edges.push_back(current_edge); }}
        
    }

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
    for (unsigned int i=0; i<elements.size(); i++) {
        int node1 = elements[i](0);
        int node2 = elements[i](1);
        int node3 = elements[i](2);
        mass[node1] += rho[i]*area[i] / 3.0;; 
        mass[node2] += rho[i]*area[i] / 3.0;; 
        mass[node3] += rho[i]*area[i] / 3.0;; 
    }
    
    // debug code
    /*
    for (unsigned int i=0; i<nodes.size(); i++) {
        std::cout << "\nnode: " << nodes[i](0) << " " << nodes[i](1) << " " << nodes[i](2) << std::endl;
    }
    
    for (unsigned int i=0; i<elements.size(); i++) {
        std::cout << "\nelement: " << elements[i](0) << " " << elements[i](1) << " " << elements[i](2) << std::endl;
    }

    for (unsigned int i=0; i<edges.size(); i++) {
        std::cout << "\nedge: " << edges[i](0) << " " << edges[i](1) << std::endl;
    }
    */

}

} } // end namespaces

#endif
