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
    typedef typename T rho_type;
    typedef std::vector<rho_type> rho_list_type; 

    
};






} } // end namespaces
