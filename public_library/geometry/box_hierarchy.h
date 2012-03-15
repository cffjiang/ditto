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
    typedef ditto::algebra::VECTOR_3D<PointType> TriangleType; 
    typedef ditto::geometry::Box<T,PointType> BoxType; 
    typedef std::vector<BoxType> BoxListType;
    typedef int ChildType;
    typedef std::vector<ChildType> ChildrenType;
    typedef std::vector<ChildrenType> ChildrenListType;
    
    BoxListType boxes;
    ChildrenListType childrens;
    int num_leafs;
    int num_boxes;
    T margin;
    
    Box_Hierarchy() { }

    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // Function: query_box
    // This function returns a list of indices.
    // They are indices of potential triangles that may intersect the input box
    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    template<typename InputBoxType, typename IndexListType>
    bool query_box(InputBoxType &B, IndexListType &intersection_list)
    {
        intersection_list.clear();
        int tree_root = boxes.size() - 1;
        query_box_recurser(tree_root, B, intersection_list);
        return (intersection_list.size() != 0);
    }

    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // Function: query_box_recurser
    // A recursive helper function for query_box. 
    // Recursively build intersection list.
    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    template<typename InputBoxType, typename IndexListType>
    void query_box_recurser(const int me, InputBoxType &B, IndexListType &intersection_list)
    {
        bool I_intersect_the_box = boxes[me].test_intersection_with_another_box(B);

        if (I_intersect_the_box) {
            if (me < num_leafs) { // I am a leaf
                intersection_list.push_back(me); }
            else {
                int how_many_children = childrens[me].size();
                for (int c=0; c<how_many_children; c++) {
                    query_box_recurser(childrens[me][c], B, intersection_list); } } }
    }

    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // Function: query_point
    // This function returns a list of indices.
    // They are indices of potential triangles that may intersect the input point
    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    template<typename InputPointType, typename IndexListType>
    bool query_point(const InputPointType &P, IndexListType &intersection_list)
    {
        intersection_list.clear();
        int tree_root = boxes.size() - 1;
        query_point_recurser(tree_root, P, intersection_list);
        return (intersection_list.size() != 0);
    }

    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // Function: query_point_recurser
    // A recursive helper function for query_point. 
    // Recursively build intersection list.
    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    template<typename InputPointType, typename IndexListType>
    void query_point_recurser(const int me, const InputPointType &P, IndexListType &intersection_list)
    {
        bool I_have_the_point = boxes[me].test_point_inside_box(P);
        
        if (I_have_the_point) {
            if (me < num_leafs) { // I am a leaf
                intersection_list.push_back(me); }
            else {
                int how_many_children = childrens[me].size();
                for (int c=0; c<how_many_children; c++) {
                    query_point_recurser(childrens[me][c], P, intersection_list); } } }
    }

    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // Function: build_tree
    // Buiding (and storing) the tree bottom up.
    // This function is only called once as initialization. Tested.
    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    template<typename TriangleListType, typename VertexPositionListType>
    void build_tree(TriangleListType &elements, VertexPositionListType &vertices, T input_margin)
    {
        margin = input_margin;

        num_leafs = elements.size();
        num_boxes = boxes.size();

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
    // This function is usually called in each time step. Tested.
    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    template<typename ElementListType, typename VertexPositionListType>
    void update_box_positions(ElementListType &elements, VertexPositionListType &vertices)
    {
        TriangleType tri;

        // build leaf boxes. can parallel.
        for (int i=0; i<num_leafs; i++) {
            tri(0) = vertices[elements[i](0)]; 
            tri(1) = vertices[elements[i](1)]; 
            tri(2) = vertices[elements[i](2)]; 
            boxes[i].build_box(i, tri, margin); }

        // build other levels by union boxes. can't parallel.
        for (int i=num_leafs; i<num_boxes; i++) {
            int children_size = childrens[i].size();
            if (children_size == 2) {
                int first_child = childrens[i][0];
                int second_child = childrens[i][1];
                boxes[i].build_union_box(i, boxes[first_child], boxes[second_child]); }
            else if (children_size == 3) {
                int first_child = childrens[i][0];
                int second_child = childrens[i][1];
                int third_child = childrens[i][2];
                boxes[i].build_union_box(i, boxes[first_child], boxes[second_child], boxes[third_child]); }
        }
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
