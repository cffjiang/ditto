//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// ditto/public_library/geometry/simplex.h
// Copyright 2012, Chenfanfu Jiang
//
// Supporting shapes:
//     Segment_3d
//     Triangle_2d
//     Triangle_3d
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

#ifndef DITTO_PUBLIC_LIBRARY_GEOMETRY_SIMPLEX_H
#define DITTO_PUBLIC_LIBRARY_GEOMETRY_SIMPLEX_H

#include <cmath>
#include <ditto/public_library/algebra/linear_algebra.h>

namespace ditto { namespace geometry {

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Class: Segment_3d
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T>
class Segment_3d {
public:
    typedef ditto::algebra::VECTOR_3D<T> Point;
    typedef ditto::algebra::VECTOR_3D<T> Vec;
    Point A;
    Point B;

    Segment_3d(Point &iA, Point &iB) {
        for (int i=0; i<3; i++) {
            A(i) = iA(i);
            B(i) = iB(i); } }

    T get_length() {
        Vec AmB = A-B;
        return AmB.Magnitude(); }

    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // Function: find_closest_point
    // This function takes care of degeneracy case: Segment too short
    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    T find_closest_point(Point &P, Point &P_hat, T &lambda) { // P_hat = (1-lambda)A + lambda*B
        Vec PmA = P-A;
        Vec BmA = B-A;
        
        // degenerate case: A is too close to B
        T length_square = BmA.Dot(BmA);
        if (length_square < 1e-10) {
            lambda = 0;
            P_hat = A;
            return  (P-P_hat).Magnitude(); }

        lambda = PmA.Dot(BmA) / length_square;

        if (lambda < 0) lambda = 0;
        else if (lambda > 1) lambda = 1;

        P_hat = A*(1-lambda) + B*lambda;
        return (P-P_hat).Magnitude();
    }

};

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
    typedef ditto::algebra::VECTOR_3D<T> Point;
    typedef ditto::algebra::VECTOR_3D<T> Vec;

    T x[3][3];
    Point A;
    Point B;
    Point C;
    
    Triangle_3d(T input[]) {
        int k = 0;
        for (int i=0; i<3; i++) {
            for (int j=0; j<3; j++) {
                x[i][j] = input[k++];
            }
        }

        A.Set_Value(x[0][0], x[0][1], x[0][2]);
        B.Set_Value(x[1][0], x[1][1], x[1][2]);
        C.Set_Value(x[2][0], x[2][1], x[2][2]);
    }

    template<class Vec3>
    Triangle_3d(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2) {
        for (int i=0; i<3; i++) {
            x[0][i] = v0(i);
            x[1][i] = v1(i);
            x[2][i] = v2(i);
        }

        A.Set_Value(x[0][0], x[0][1], x[0][2]);
        B.Set_Value(x[1][0], x[1][1], x[1][2]);
        C.Set_Value(x[2][0], x[2][1], x[2][2]);
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

    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // Function: find_closest_point
    // Before calling this function, should make sure area > tol
    // in order to avoid divide-by-zero problem.
    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    T find_closest_point(Point &P, Point &P_hat, T &ksi1, T &ksi2) { // P_hat = (1-ksi1-ksi2)A + ksi1*B + ksi2*C
        Vec BmA = B - A;
        Vec CmA = C - A;
        Vec PmA = P - A;
        ditto::algebra::MATRIX_2X2<T> M(BmA.Dot(BmA), CmA.Dot(BmA), CmA.Dot(BmA), CmA.Dot(CmA));
        M.Invert();
        ditto::algebra::VECTOR_2D<T> rhs(PmA.Dot(BmA), PmA.Dot(CmA));
        ditto::algebra::VECTOR_2D<T> ksi = M*rhs;

        ksi1 = ksi(0);
        ksi2 = ksi(1);

        T d;

        Point P_on_plane = A*(1-ksi1-ksi2) + B*ksi1 + C*ksi2;

        if (ksi1>=0 && ksi1<=1 && ksi2>=0 && ksi2<=1 && (1-ksi1-ksi2)>=0 && (1-ksi1-ksi2)<=1) { 
            P_hat = P_on_plane;
            d = (P_hat-P).Magnitude(); }
        else {
            T dmin = 100000;
            T lambda1, lambda2, lambda3;
            int flag = -1;

            Segment_3d<T> AB(A, B);
            Point P_on_AB;
            T distance_to_AB = AB.find_closest_point(P_on_plane, P_on_AB, lambda1);
            if (distance_to_AB < dmin) {
                flag = 1; 
                dmin = distance_to_AB;}
            
            Segment_3d<T> BC(B, C);
            Point P_on_BC;
            T distance_to_BC = BC.find_closest_point(P_on_plane, P_on_BC, lambda2);
            if (distance_to_BC < dmin) {
                flag = 2; 
                dmin = distance_to_BC;}

            Segment_3d<T> CA(C, A);
            Point P_on_CA;
            T distance_to_CA = CA.find_closest_point(P_on_plane, P_on_CA, lambda3);
            if (distance_to_CA < dmin) {
                flag = 3; 
                dmin = distance_to_CA;}

            if (flag == 1) {
                P_hat = P_on_AB; 
                ksi2 = 0;
                ksi1 = lambda1;
            }
            else if (flag == 2) {
                P_hat = P_on_BC; 
                ksi1 = 1-lambda2;
                ksi2 = lambda2;
            }
            else if (flag == 3) {
                P_hat = P_on_CA; 
                ksi2 = 1-lambda3;
                ksi1 = 0;
            }
        }

        d = (P_hat - P).Magnitude();
        return d;

    }

};


} } // end namespaces

#endif
