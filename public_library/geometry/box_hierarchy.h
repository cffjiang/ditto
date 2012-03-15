//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// ditto/public_library/geometry/box_hierarchy.h
// Copyright 2012, Chenfanfu Jiang
//
// Naive box hiearchy                        
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

#ifndef DITTO_PUBLIC_LIBRARY_GEOMETRY_BOX_HIERARCHY_H
#define DITTO_PUBLIC_LIBRARY_GEOMETRY_BOX_HIERARCHY_H

#include <vector>
#include <ditto/public_library/algebra/linear_algebra.h>
#include <ditto/public_library/geometry/box.h>

namespace ditto { namespace geometry {

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Class: Box_Hierarchy
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T>
class Box_Hierarchy {

public:
    typedef ditto::algebra::VECTOR_3D<T> PointType;
    typedef ditto::geometry::Box<T,PointType> BoxType; 
    typedef std::vector<BoxType> BoxListType;
    typedef int ChildType;
    typedef std::vector<ChildType> ChildrenType;
    typedef std::vector<ChildrenType> ChildrenListType;
    
    BoxListType boxes;
    ChildrenListType childrens;
    
    Box_Hierarchy() { }

    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // Function: build_tree
    // Buiding (and storing) the tree bottom up.
    // This function is only called once as initialization.
    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    template<typename TriangleListType, typename VertexPositionListType>
    void build_tree(TriangleListType &elements, VertexPositionListType &vertices)
    {
        int num_leafs = elements.size();
        int num_boxes = boxes.size();

        assert(num_boxes == 0);

        // compute how many boxes in total we have. something like 7+3+1 or 8+4+3+1
        for (int num_current_level = num_leafs; num_current_level >= 1; num_current_level/=2) {
            num_boxes += num_current_level; }
        
        boxes.clear();
        childrens.clear();
        boxes.resize(num_boxes);
        childrens.resize(num_boxes); 
        
        build_hierarchy_structure_based_on_mesh_connectivity(num_leafs, num_boxes);
        update_box_positions(elements, vertices);
    }

    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // Function: update_box_positions
    // This function is usually called in each time step
    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    template<typename TriangleListType, typename VertexPositionListType>
    void update_box_positions(TriangleListType &elements, VertexPositionListType &vertices)
    {
        //TODO
    }

    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // Function: build_hierarchy_structure_based_on_mesh_connectivity
    // This function assumes childrens has the correct size already.
    // This function is only called by build_tree function.
    // This function is tested!
    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    void build_hierarchy_structure_based_on_mesh_connectivity(int num_leafs, int num_boxes)
    {
        int start = num_leafs;
        int last_start = 0;
        for (int num_current_level = num_leafs/2; num_current_level >=1; num_current_level /= 2) {
            int next_start = start + num_current_level;

            for (int i=start; i<next_start; i++) {
                // push something to childrens[i]
                childrens[i].push_back(last_start++);
                childrens[i].push_back(last_start++); }
            if (last_start != start) {
                childrens[next_start-1].push_back(last_start++); }
            assert(last_start = start);
            start = next_start;
        }
    }
};




} } // end namespaces

#endif
