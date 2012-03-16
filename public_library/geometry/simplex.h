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
#include <vector>
#include <cassert>
#include <algorithm>
#include <ditto/public_library/algebra/linear_algebra.h>
#include <ditto/public_library/algebra/cubic_equation_solve.h>
#include <ditto/public_library/algebra/Eigen3/Eigen/LU>
#include <ditto/public_library/algebra/Eigen3/Eigen/Dense>

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

    Point& operator()(const int i) {
        assert(i == 0 || i == 1);
        if (i == 0) return A;
        else if (i == 1) return B; }

    T get_length() {
        Vec AmB = A-B;
        return AmB.Magnitude(); }

    void get_box(T buffer, T &xmin, T &xmax, T &ymin, T &ymax, T &zmin, T &zmax) {
        std::vector<T> xs(2);
        xs[0] = A(0); xs[1] = B(0);
        xmin = *(std::min_element(xs.begin(), xs.end())) - buffer;
        xmax = *(std::max_element(xs.begin(), xs.end())) + buffer;
        
        std::vector<T> ys(2);
        ys[0] = A(1); ys[1] = B(1);
        ymin = *(std::min_element(ys.begin(), ys.end())) - buffer;
        ymax = *(std::max_element(ys.begin(), ys.end())) + buffer;
        
        std::vector<T> zs(2);
        zs[0] = A(2); zs[1] = B(2);
        zmin = *(std::min_element(zs.begin(), zs.end())) - buffer;
        zmax = *(std::max_element(zs.begin(), zs.end())) + buffer;
    }

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

    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // Function: find_closest_points_seg_seg
    // This function takes care of degeneracy cases
    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    T find_closest_points_seg_seg(Segment_3d &him, Point &P, Point &Q, T &s, T &t) { //P = (1-s)A + sB, Q = (1-t)C + tD, C = him.A, D = him.B
        Point C = him.A;
        Point D = him.B;
        Vec vBA = B - A;
        Vec vDC = D - C;
        Vec vCA = C - A;
        Vec cross = vBA.Cross_Product(vDC);
        T area_square = cross.Dot(cross);
        T tol = 1e-10;

        if (area_square > tol) { // two segments are not parallel
            ditto::algebra::MATRIX_2X2<T> M(vDC.Dot(vBA), vDC.Dot(vDC), -(vBA.Dot(vBA)), -(vBA.Dot(vDC)));
            M.Invert();
            ditto::algebra::VECTOR_2D<T> rhs(-(vCA.Dot(vBA)), -(vCA.Dot(vDC)));
            ditto::algebra::VECTOR_2D<T> ts = M*rhs;

            t = ts(0);
            s = ts(1);

            if (t < 0) t = 0;
            else if (t > 1) t = 1;

            if (s < 0) s = 0;
            else if (s > 1) s = 1;

            P = A*(1-s) + B*s;
            Q = C*(1-t) + D*t;
        } else { // AB || CD
            T dmin = 100000;
            int flag = -1;
            
            Point C_hat;
            T lambda_C_hat;
            T dC = find_closest_point(C, C_hat, lambda_C_hat); // C_hat = (1-lambda_C_hat) A + lambda_C_hat B;
            if (dC < dmin) {
                flag = 1;
                dmin = dC; }

            Point D_hat;
            T lambda_D_hat;
            T dD = find_closest_point(D, D_hat, lambda_D_hat); // D_hat = (1-lambda_D_hat) A + lambda_D_hat B;
            if (dD < dmin) {
                flag = 2;
                dmin = dD; }

            Point A_hat;
            T lambda_A_hat;
            T dA = him.find_closest_point(A, A_hat, lambda_A_hat); // A_hat = (1-lambda_A_hat) him.A + lambda_A_hat him.B;
            if (dA < dmin) {
                flag = 3;
                dmin = dA; }

            Point B_hat;
            T lambda_B_hat;
            T dB = him.find_closest_point(B, B_hat, lambda_B_hat); // B_hat = (1-lambda_B_hat) him.A + lambda_B_hat him.B;
            if (dB < dmin) {
                flag = 4;
                dmin = dB; }

            assert(flag != -1);

            if (flag == 1) {
                P = C_hat;
                Q = C;
                s = lambda_C_hat;
                t = 0; }
            else if (flag == 2) {
                P = D_hat;
                Q = D;
                s = lambda_D_hat;
                t = 1; }
            else if (flag == 3) {
                P = A;
                Q = A_hat;
                s = 0;
                t = lambda_A_hat; }
            else if (flag == 4) {
                P = B;
                Q = B_hat;
                s = 1;
                t = lambda_B_hat; }
        }
        
        return (P-Q).Magnitude();

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

    void get_box(T buffer, T &xmin, T &xmax, T &ymin, T &ymax, T &zmin, T &zmax) {
        std::vector<T> xs(3);
        xs[0] = A(0); xs[1] = B(0); xs[2] = C(0);
        xmin = *(std::min_element(xs.begin(), xs.end())) - buffer;
        xmax = *(std::max_element(xs.begin(), xs.end())) + buffer;

        std::vector<T> ys(3);
        ys[0] = A(1); ys[1] = B(1); ys[2] = C(1);
        ymin = *(std::min_element(ys.begin(), ys.end())) - buffer;
        ymax = *(std::max_element(ys.begin(), ys.end())) + buffer;

        std::vector<T> zs(3);
        zs[0] = A(2); zs[1] = B(2); zs[2] = C(2);
        zmin = *(std::min_element(zs.begin(), zs.end())) - buffer;
        zmax = *(std::max_element(zs.begin(), zs.end())) + buffer;
    }

    template<class Vec3>
    void get_normal(Vec3 n) {
        ditto::algebra::VECTOR_3D<T> edge1(x[1][0]-x[0][0], x[1][1]-x[0][1], x[1][2]-x[0][2]);
        ditto::algebra::VECTOR_3D<T> edge2(x[2][0]-x[0][0], x[2][1]-x[0][1], x[2][2]-x[0][2]);
        ditto::algebra::VECTOR_3D<T> normal = (edge1.Cross_Product(edge2)).Normalize();
        for (int i=0; i<3; i++) n(i) = normal(i); }

    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // Function: contain_point
    // This function only tests for point in the same plane
    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    template<class InputPointType>
    bool contain_point(InputPointType &P)
    {
        Vec BmA = B - A;
        Vec CmA = C - A;
        Vec PmA = P - A;
        ditto::algebra::MATRIX_2X2<T> M(BmA.Dot(BmA), CmA.Dot(BmA), CmA.Dot(BmA), CmA.Dot(CmA));
        M.Invert();
        ditto::algebra::VECTOR_2D<T> rhs(PmA.Dot(BmA), PmA.Dot(CmA));
        ditto::algebra::VECTOR_2D<T> ksi = M*rhs;

        T ksi1 = ksi(0);
        T ksi2 = ksi(1);
        
        if (ksi1 < 0 || ksi1 > 1 || ksi2 < 0 || ksi2 > 1 || 1-ksi1-ksi2 < 0 || 1-ksi1-ksi2 > 1) return false;
        else return true;
    }

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

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Class: Tetrahedron
// Orientation: ABC is counter clockwise at bottom, D is top
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T>
class Tetrahedron {
public:
    typedef ditto::algebra::VECTOR_3D<T> Point;
    typedef ditto::algebra::VECTOR_3D<T> Vec;

    Point A;
    Point B;
    Point C;
    Point D;

    template<class InputPointType>
    Tetrahedron(InputPointType &iA, InputPointType &iB, InputPointType &iC, InputPointType &iD)
    {
        A = iA; B = iB; C = iC; D = iD;
    }

    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // Function penetration_safety_test
    // This function is non-trivial.
    // It tests whether at some time between 0 and dt, D penetrates triangle ABC.
    // (Assuming every point moves along current velocity.)
    // True means safe, false means there will be penetration.
    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    template<class InputVelocityType>
    bool penetration_safety_test(InputVelocityType &vA, InputVelocityType &vB, InputVelocityType &vC, InputVelocityType &vD, T dt)
    {
        typedef ditto::algebra::VECTOR_2D<T> ComplexNumber;
        Eigen::Matrix4f M;

        M << 1, 1, 1, 1,
            A(0), B(0), C(0), D(0),
            A(1), B(1), C(1), D(1),
            A(2), B(2), C(2), D(2);
        T V0 = M.determinant()/6;

        M << 1, 1, 1, 1,
            A(0)+vA(0)*dt, B(0)+vB(0)*dt, C(0)+vC(0)*dt, D(0)+vD(0)*dt,
            A(1)+vA(1)*dt, B(1)+vB(1)*dt, C(1)+vC(1)*dt, D(1)+vD(1)*dt,
            A(2)+vA(2)*dt, B(2)+vB(2)*dt, C(2)+vC(2)*dt, D(2)+vD(2)*dt;
        T V1 = M.determinant()/6;

        M << 1, 1, 1, 1,
            A(0)-vA(0)*dt, B(0)-vB(0)*dt, C(0)-vC(0)*dt, D(0)-vD(0)*dt,
            A(1)-vA(1)*dt, B(1)-vB(1)*dt, C(1)-vC(1)*dt, D(1)-vD(1)*dt,
            A(2)-vA(2)*dt, B(2)-vB(2)*dt, C(2)-vC(2)*dt, D(2)-vD(2)*dt;
        T Vm1 = M.determinant()/6;

        M << 1, 1, 1, 1,
            A(0)+vA(0)*dt*2, B(0)+vB(0)*dt*2, C(0)+vC(0)*dt*2, D(0)+vD(0)*dt*2,
            A(1)+vA(1)*dt*2, B(1)+vB(1)*dt*2, C(1)+vC(1)*dt*2, D(1)+vD(1)*dt*2,
            A(2)+vA(2)*dt*2, B(2)+vB(2)*dt*2, C(2)+vC(2)*dt*2, D(2)+vD(2)*dt*2;
        T V2 = M.determinant()/6;

        M << 1, 1, 1, 1,
            A(0)-vA(0)*dt*2, B(0)-vB(0)*dt*2, C(0)-vC(0)*dt*2, D(0)-vD(0)*dt*2,
            A(1)-vA(1)*dt*2, B(1)-vB(1)*dt*2, C(1)-vC(1)*dt*2, D(1)-vD(1)*dt*2,
            A(2)-vA(2)*dt*2, B(2)-vB(2)*dt*2, C(2)-vC(2)*dt*2, D(2)-vD(2)*dt*2;
        T Vm2 = M.determinant()/6;
        
        T d = V0;
        T b = 0.5 * (V1 + Vm1 - 2*d);
        T a = (V2 - Vm2 - 2*V1 + 2*Vm1) / 12;
        T c = 0.5 * (V1 - Vm1 - 2*a);

        ComplexNumber root1, root2, root3;
        int root_nature = ditto::algebra::cubic_equation_solve(a, b, c, d, root1, root2, root3);

        T t,t1,t2,t3;
        if (root_nature == 1)
        {
                t = root1(0);
                if (t > 1 || t < 0) return true;
                ditto::geometry::Triangle_3d<T> key_tri(A+vA*t*dt, B+vB*t*dt, C+vC*t*dt);
                ditto::algebra::VECTOR_3D<T> key_point = D+vD*t+dt;
                if (key_tri.contain_point(key_point)) return false;
                else return true;
        }
        else if (root_nature == 2)
        {
                t1 = root1(0);
                t2 = root2(0);
                t3 = root3(0);
                if (t1 > 0 && t1 < 1) {
                    Triangle_3d<T> key_tri1(A+vA*t1*dt, B+vB*t1*dt, C+vC*t1*dt);
                    ditto::algebra::VECTOR_3D<T> key_point1 = D+vD*t1*dt;
                    if (key_tri1.contain_point(key_point1)) return false; }
                if (t2 > 0 && t2 < 1) {
                    Triangle_3d<T> key_tri2(A+vA*t2*dt, B+vB*t2*dt, C+vC*t2*dt);
                    ditto::algebra::VECTOR_3D<T> key_point2 = D+vD*t2*dt;
                    if (key_tri2.contain_point(key_point2)) return false; }
                if (t3 > 0 && t3 < 1) {
                    Triangle_3d<T> key_tri3(A+vA*t3*dt, B+vB*t3*dt, C+vC*t3*dt);
                    ditto::algebra::VECTOR_3D<T> key_point3 = D+vD*t3*dt;
                    if (key_tri3.contain_point(key_point3)) return false; }
                return true;
        }
        else if (root_nature == 3)
        {
                t = root1(0);
                if (t > 1 || t < 0) return true;
                Triangle_3d<T> key_tri0(A+vA*t*dt, B+vB*t*dt, C+vC*t*dt);
                ditto::algebra::VECTOR_3D<T> key_point0 = D+vD*t*dt;
                if (key_tri0.contain_point(key_point0)) return false;
                else return true;
        }

    }



};


} } // end namespaces

#endif
