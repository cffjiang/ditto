//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// ditto/public_library/geometry/simplex.h
// Copyright 2012, Chenfanfu Jiang
//
// Supporting shapes:
//     Triangle_2d
//     Triangle_3d
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

#ifndef DITTO_PUBLIC_LIBRARY_GEOMETRY_SIMPLEX_H
#define DITTO_PUBLIC_LIBRARY_GEOMETRY_SIMPLEX_H

#include <cmath>
#include <ditto/public_library/algebra/linear_algebra.h>

namespace ditto { namespace geometry {

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Class: Triangle_2d
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T>
class Triangle_2d {
public:
    T x[3][2];
    
    Triangle_2d(T input[]) {
        int k = 0;
        for (int i=0; i<3; i++) {
            for (int j=0; j<2; j++) {
                x[i][j] = input[k++];
            }
        }
    }

    T get_area() {
        return 0.5 * ( x[0][0]*x[1][1] + x[1][0]*x[2][1] + x[2][0]*x[0][1] - x[0][0]*x[2][1] - x[1][0]*x[0][1] - x[2][0]*x[1][1] );
    }
};

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Class: Triangle_3d
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T>
class Triangle_3d {
public:
    T x[3][3];
    
    Triangle_3d(T input[]) {
        int k = 0;
        for (int i=0; i<3; i++) {
            for (int j=0; j<3; j++) {
                x[i][j] = input[k++];
            }
        }
    }

    template<class Vec3>
    Triangle_3d(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2) {
        for (int i=0; i<3; i++) {
            x[0][i] = v0[i];
            x[1][i] = v1[i];
            x[2][i] = v2[i];
        }
    }

    T get_area() {
        T a = x[1][0] - x[0][0];
        T b = x[1][1] - x[0][1];
        T c = x[1][2] - x[0][2];
        T d = x[2][0] - x[0][0];
        T e = x[2][1] - x[0][1];
        T f = x[2][2] - x[0][2];
        T ci = b*f-c*e;
        T cj = c*d-a*f;
        T ck = a*e-b*d;
        return 0.5*std::sqrt(ci*ci + cj*cj + ck*ck); }

    template<class Vec3>
    void get_normal(Vec3 n) {
        ditto::algebra::VECTOR_3D<T> edge1(x[1][0]-x[0][0], x[1][1]-x[0][1], x[1][2]-x[0][2]);
        ditto::algebra::VECTOR_3D<T> edge2(x[2][0]-x[0][0], x[2][1]-x[0][1], x[2][2]-x[0][2]);
        ditto::algebra::VECTOR_3D<T> normal = (edge1.Cross_Product(edge2)).Normalize();
        for (int i=0; i<3; i++) n(i) = normal(i); }
};


} } // end namespaces

#endif
