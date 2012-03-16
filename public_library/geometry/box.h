//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// ditto/public_library/geometry/box.h
// Copyright 2012, Chenfanfu Jiang
//
// Axis aligned box.
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

#ifndef DITTO_PUBLIC_LIBRARY_GEOMETRY_BOX_H
#define DITTO_PUBLIC_LIBRARY_GEOMETRY_BOX_H

#include <ditto/public_library/algebra/linear_algebra.h>
#include <algorithm>

namespace ditto { namespace geometry {

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Class: Box
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T, class PointType>
class Box {
public:
    int id;
    PointType Pmin;
    PointType Pmax;

    Box() {} // default constructor

    template<class TriangleType>
    void build_box(const int input_id, const TriangleType &tri, const T margin)
    {
        id = input_id;

        T x_min = 10000, x_max = -10000, y_min = 10000, y_max = -10000, z_min = 10000, z_max = -10000;
        for (int p=0; p<3; p++) {
            if (tri(p)(0) < x_min) x_min = tri(p)(0);
            if (tri(p)(0) > x_max) x_max = tri(p)(0);
            if (tri(p)(1) < y_min) y_min = tri(p)(1);
            if (tri(p)(1) > y_max) y_max = tri(p)(1);
            if (tri(p)(2) < z_min) z_min = tri(p)(2);
            if (tri(p)(2) > z_max) z_max = tri(p)(2); }

        Pmin(0) = x_min - margin;
        Pmin(1) = y_min - margin;
        Pmin(2) = z_min - margin;

        Pmax(0) = x_max + margin;
        Pmax(1) = y_max + margin;
        Pmax(2) = z_max + margin;
    }

    template<class TriangleType>
    void build_box(const int input_id, const TriangleType &triA, const TriangleType &triB, const T margin)
    {
        id = input_id;

        T x_min = 10000, x_max = -10000, y_min = 10000, y_max = -10000, z_min = 10000, z_max = -10000;
        for (int p=0; p<3; p++) {
            if (triA(p)(0) < x_min) x_min = triA(p)(0);
            if (triA(p)(0) > x_max) x_max = triA(p)(0);
            if (triA(p)(1) < y_min) y_min = triA(p)(1);
            if (triA(p)(1) > y_max) y_max = triA(p)(1);
            if (triA(p)(2) < z_min) z_min = triA(p)(2);
            if (triA(p)(2) > z_max) z_max = triA(p)(2); 
            
            if (triB(p)(0) < x_min) x_min = triB(p)(0);
            if (triB(p)(0) > x_max) x_max = triB(p)(0);
            if (triB(p)(1) < y_min) y_min = triB(p)(1);
            if (triB(p)(1) > y_max) y_max = triB(p)(1);
            if (triB(p)(2) < z_min) z_min = triB(p)(2);
            if (triB(p)(2) > z_max) z_max = triB(p)(2); 
        }

        Pmin(0) = x_min - margin;
        Pmin(1) = y_min - margin;
        Pmin(2) = z_min - margin;

        Pmax(0) = x_max + margin;
        Pmax(1) = y_max + margin;
        Pmax(2) = z_max + margin;        
    }

    template<class SegmentType>
    void build_box_from_segment(const int input_id, SegmentType &seg, const T margin)
    {
        id = input_id;

        T x_min = 10000, x_max = -10000, y_min = 10000, y_max = -10000, z_min = 10000, z_max = -10000;
        for (int p=0; p<2; p++) {
            if (seg(p)(0) < x_min) x_min = seg(p)(0);
            if (seg(p)(0) > x_max) x_max = seg(p)(0);
            if (seg(p)(1) < y_min) y_min = seg(p)(1);
            if (seg(p)(1) > y_max) y_max = seg(p)(1);
            if (seg(p)(2) < z_min) z_min = seg(p)(2);
            if (seg(p)(2) > z_max) z_max = seg(p)(2); }

        Pmin(0) = x_min - margin;
        Pmin(1) = y_min - margin;
        Pmin(2) = z_min - margin;

        Pmax(0) = x_max + margin;
        Pmax(1) = y_max + margin;
        Pmax(2) = z_max + margin;
    }

    void build_union_box(const int input_id, Box<T,PointType> &A, Box<T,PointType> &B)
    {
        id = input_id;

        for (int i=0; i<3; i++) {
            Pmin(i) = std::min(A.Pmin(i), B.Pmin(i));
            Pmax(i) = std::max(A.Pmax(i), B.Pmax(i)); }
    }

    void build_union_box(const int input_id, Box<T,PointType> &A, Box<T,PointType> &B, Box<T,PointType> &C)
    {
        id = input_id;

        for (int i=0; i<3; i++) {
            Pmin(i) = std::min ( std::min(A.Pmin(i), B.Pmin(i)) , C.Pmin(i) );
            Pmax(i) = std::max ( std::max(A.Pmax(i), B.Pmax(i)) , C.Pmax(i) ); 
        }
    }

    template<class TestPointType>
    bool test_point_inside_box(TestPointType &P)
    {
        if ( P(0) > Pmax(0) || P(0) < Pmin(0) || P(1) > Pmax(1) || P(1) < Pmin(1) || P(2) > Pmax(2) || P(2) < Pmin(2) ) {
            return false; }
        return true;
    }

    template<class TestBoxType>
    bool test_intersection_with_another_box(TestBoxType &B)
    {
        if (Pmin(0) > B.Pmax(0) || Pmax(0) < B.Pmin(0) ||
            Pmin(1) > B.Pmax(1) || Pmax(1) < B.Pmin(1) ||
            Pmin(2) > B.Pmax(2) || Pmax(2) < B.Pmin(2) ) {
            return false; }
        return true;
    }

};


}}

#endif
