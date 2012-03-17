//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// ditto/public_library/cloth_3d/neo_hookean_cloth_3d_fvm_explicit.h
// Copyright 2012, Chenfanfu Jiang
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

#ifndef DITTO_PUBLIC_LIBRARY_CLOTH_3D_NEO_HOOKEAN_CLOTH_3D_FVM_EXPLICIT_H
#define DITTO_PUBLIC_LIBRARY_CLOTH_3D_NEO_HOOKEAN_CLOTH_3D_FVM_EXPLICIT_H

#include <cmath>
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <ctime>

#include <ditto/public_library/geometry/simplex.h>
#include <ditto/public_library/geometry/box_hierarchy.h>
#include <ditto/public_library/geometry/box.h>
#include <ditto/public_library/algebra/linear_algebra.h>
#include <ditto/public_library/algebra/Eigen3/Eigen/Dense>
#include <ditto/public_library/algebra/Eigen3/Eigen/SVD>
#include <ditto/public_library/algebra/Eigen3/Eigen/Jacobi>
#include <ditto/public_library/visualization/vtk_writer.h>

namespace ditto { namespace cloth_3d {

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Class: Neo_Hookean_Cloth_3d_Fvm_Explicit
// MeshType: ditto::geometry::Triangle_Mesh_3d
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T, class MeshType>
class Neo_Hookean_Cloth_3d_Fvm_Explicit {
public:
    typedef ditto::algebra::VECTOR_3D<T> node_3d_type;
    typedef typename std::vector<node_3d_type> node_3d_list_type;
    typedef ditto::algebra::VECTOR_2D<T> node_2d_type;
    typedef typename std::vector<node_2d_type> node_2d_list_type;
    typedef ditto::algebra::VECTOR_3D<int> element_type;
    typedef typename std::vector<element_type> element_list_type;
    typedef typename std::vector<T> rho_list_type;
    typedef typename std::vector<T> area_list_type;
    typedef ditto::algebra::MATRIX_2X2<T> matrix_2x2_type;
    typedef ditto::algebra::MATRIX_3X3<T> matrix_3x3_type;
    typedef typename std::vector<matrix_2x2_type> matrix_2x2_list_type;
    typedef ditto::geometry::Box_Hierarchy<T> triangle_hierarchy_type;

    MeshType &mesh;
    T dt;

    node_3d_list_type X_3d;
    node_3d_list_type x;
    node_3d_list_type v;
    node_3d_list_type f;

    bool use_gravity;
    bool use_ball_collision;
    bool use_ground_collision;

    bool use_self_collision;
    int self_collision_repulsion_iters;
    int self_collision_collision_iters;
    T self_collision_distance_tolerance;

    triangle_hierarchy_type hierarchy;
    triangle_hierarchy_type hierarchy_predicting_future;

    node_3d_type ball_center;
    T ball_radius;

    matrix_2x2_list_type Dm_inverse_list;
    matrix_2x2_list_type Dm_list;
    
    std::vector<int> dirichlet_nodes;
    node_3d_list_type dirichlet_displacement;

    int num_of_clothes;

    T youngs_modulus;
    T poisson_ratio;
    T mu;
    T lambda;
    T gamma;
 
    Neo_Hookean_Cloth_3d_Fvm_Explicit(MeshType &input_mesh, const T input_dt, const T input_E, const T input_nu, const T input_gamma, const bool input_use_gravity)
        :mesh(input_mesh) {
        dt = input_dt;
        X_3d = mesh.nodes;
        x = mesh.nodes;
        v.resize(x.size(), T(0.0));
        youngs_modulus = input_E;
        poisson_ratio = input_nu;
        mu = youngs_modulus / (2*(1+poisson_ratio));
        lambda = youngs_modulus * poisson_ratio / ( (1+poisson_ratio)*(1-2*poisson_ratio) );
        gamma = input_gamma;
        use_gravity = input_use_gravity;
        pre_build_Dm_inverse();
        use_ball_collision = false;
        use_ground_collision = false;
        use_self_collision = false;
        num_of_clothes = 1;
    }

    void aware_of_num_of_clothes(int N);

    void pre_build_Dm_inverse();

    void compute_elasticity();

    void update_one_step();

    void add_gravity(T gx, T gy, T gz);

    void set_dirichlet_with_a_bounding_box(T x0, T xM, T y0, T yM, T z0, T zM, T xmove, T ymove, T zmove);

    void clear_previous_dirichlet_conditions();

    void add_ball_collision(T xb, T yb, T zb, T rb, T spring_constant);

    void add_ground_collision(T ground_level, T spring_constant, T friction_constant);
    
    void do_point_triangle_repulsion(int iterations, T d_tol);

    void do_segment_segment_repulsion(int iterations, T d_tol);

    void do_point_triangle_collision(int iterations);

    void initialize_hierarchy(T margin);

    void update_hierarchy();

    void switch_self_collision(bool s, int input_iter, T input_dtol);

    void write_output(int frame);

    void write_vtk(int frame);
};

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Function: aware_of_num_clothes
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T, class MeshType>
void Neo_Hookean_Cloth_3d_Fvm_Explicit<T, MeshType>::
aware_of_num_of_clothes(int N) {
    num_of_clothes = N;
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Function: pre_build_Dm_inverse
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T, class MeshType>
void Neo_Hookean_Cloth_3d_Fvm_Explicit<T, MeshType>::
pre_build_Dm_inverse() {
    for (int e = 0; e < mesh.elements.size(); e++) {
    
        // first map the 3d triangle to a 2d plane in x-y
        node_3d_type X1 = X_3d[mesh.elements[e](0)];
        node_3d_type X2 = X_3d[mesh.elements[e](1)];
        node_3d_type X3 = X_3d[mesh.elements[e](2)];
        node_3d_type X2mX1 = X2 - X1;
        node_3d_type X3mX1 = X3 - X1;
        node_3d_type edge_cross = X2mX1.Cross_Product(X3mX1);
        node_3d_type nvec = edge_cross*(1.0/edge_cross.Magnitude());
        node_3d_type nxX2mX1 = nvec.Cross_Product(X2mX1);
        node_3d_type wvec = X2mX1*(1.0/X2mX1.Magnitude());
        node_3d_type vvec = nxX2mX1*(1.0/nxX2mX1.Magnitude());
    
        node_2d_type X1_mapped(0.0, 0.0);
        node_2d_type X2_mapped(X2mX1.Dot(wvec), X2mX1.Dot(vvec));
        node_2d_type X3_mapped(X3mX1.Dot(wvec), X3mX1.Dot(vvec));

        // now build Dm_inverse
        matrix_2x2_type Dm(X2_mapped(0), X2_mapped(1), X3_mapped(0), X3_mapped(1));
        Dm_list.push_back(Dm);
        Dm.Invert();
        Dm_inverse_list.push_back(Dm);
    }
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Function: compute_elasticity
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T, class MeshType>
void Neo_Hookean_Cloth_3d_Fvm_Explicit<T, MeshType>::
compute_elasticity() {

    // build f
    f.clear();
    f.resize(mesh.nodes.size());

#pragma omp parallel for schedule(static)
    for (unsigned int i=0; i<f.size(); i++){
        f[i](0) = f[i](1) = f[i](2) = 0.0;}

#pragma omp parallel for schedule(static)
    for (int e = 0; e < mesh.elements.size(); e++) {
        // get x
        node_3d_type &x1 = x[mesh.elements[e](0)];
        node_3d_type &x2 = x[mesh.elements[e](1)];
        node_3d_type &x3 = x[mesh.elements[e](2)];

        // Ds = [x2-x1, x3-x1]
        ditto::algebra::MATRIX_MXN<T> Ds(3,2);
        Ds(0,0) = x2(0)-x1(0); Ds(0,1) = x3(0)-x1(0);
        Ds(1,0) = x2(1)-x1(1); Ds(1,1) = x3(1)-x1(1);
        Ds(2,0) = x2(2)-x1(2); Ds(2,1) = x3(2)-x1(2);

        // Get Dm_inverse
        ditto::algebra::MATRIX_MXN<T> Dm(2,2);
        for (int i=0; i<2; i++) {
            for (int j=0; j<2; j++) {
                Dm(i,j) = Dm_list[e](i,j); }}
        ditto::algebra::MATRIX_MXN<T> Dm_inverse(2,2);
        for (int i=0; i<2; i++) {
            for (int j=0; j<2; j++) {
                Dm_inverse(i,j) = Dm_inverse_list[e](i,j); }}

        // F = Ds*Dm_inverse
        ditto::algebra::MATRIX_MXN<T> F(3,2);
        ditto::algebra::Multiply(Ds, Dm_inverse, F);

        // SVD
        Eigen::MatrixXf A(3,2);
        for (int i=0; i<3; i++) {
            for (int j=0; j<2; j++) {
                A(i,j) = F(i,j); }}

        Eigen::JacobiSVD<Eigen::MatrixXf> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
        Eigen::MatrixXf U = svd.matrixU();
        Eigen::MatrixXf V = svd.matrixV();
        Eigen::MatrixXf D(2,2);
        D << svd.singularValues()(0), 0.0, 0.0, svd.singularValues()(1);
        T sigma1 = svd.singularValues()(0);
        T sigma2 = svd.singularValues()(1);

        // Compute P
        Eigen::MatrixXf dPhidSigma(2,2);
        dPhidSigma <<  2*mu*(sigma1-1) + lambda*(sigma1*sigma2-1)*sigma2, 0.0, 0.0,  2*mu*(sigma2-1) + lambda*(sigma1*sigma2-1)*sigma1;
        Eigen::MatrixXf P = U * dPhidSigma * V.transpose();

        // Compute Force matrix
        T my_area = 0.5 * std::abs(Dm(0,0)* Dm(1,1) -  Dm(0,1)* Dm(1,0));


        Eigen::MatrixXf Dminvt(2,2);
        for (int i=0; i<2; i++) {
            for (int j=0; j<2; j++) {
                Dminvt(i,j) = Dm_inverse_list[e](j,i); }}

        Eigen::MatrixXf Force(3,2);

        Force = -my_area * P * Dminvt;
        
        // Compute force contributions to nodes
        node_3d_type stress_node2(Force(0,0), Force(1,0), Force(2,0));
        node_3d_type stress_node3(Force(0,1), Force(1,1), Force(2,1));
        node_3d_type stress_node1 = (stress_node2 + stress_node3) * (-1.0);

        f[mesh.elements[e](0)] = f[mesh.elements[e](0)] + stress_node1;
        f[mesh.elements[e](1)] = f[mesh.elements[e](1)] + stress_node2;
        f[mesh.elements[e](2)] = f[mesh.elements[e](2)] + stress_node3;
    
        // Add damping forces
        node_3d_type v_bar = (v[mesh.elements[e](0)] + v[mesh.elements[e](1)] + v[mesh.elements[e](2)]) * (1.0/3.0);
        for (int i=0; i<3; i++) {
            f[mesh.elements[e](i)] =  f[mesh.elements[e](i)] + ( -gamma*youngs_modulus * (v[mesh.elements[e](i)] - v_bar) ); }
    }

}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Function: update_one_step
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T, class MeshType>
void Neo_Hookean_Cloth_3d_Fvm_Explicit<T, MeshType>::
update_one_step() {
    // update v
#pragma omp parallel for schedule(static)
    for (unsigned int i=0; i<mesh.nodes.size(); i++) {
        v[i] = v[i] + dt * f[i] * (1.0/mesh.mass[i]);}
    
    // modify v with self collision impulse
    if (use_self_collision) {
        update_hierarchy();

        do_point_triangle_repulsion(self_collision_repulsion_iters, self_collision_distance_tolerance);
        do_segment_segment_repulsion(self_collision_repulsion_iters, self_collision_distance_tolerance); 

        do_point_triangle_collision(self_collision_collision_iters);
    }

    // update x
#pragma omp parallel for schedule(static)
    for (unsigned int i=0; i<mesh.nodes.size(); i++) {
        x[i] = x[i] + dt * v[i]; }
    
    // fix x for dirichlet nodes
#pragma omp parallel for schedule(static)
    for (unsigned int i=0; i<dirichlet_nodes.size(); i++) {
        v[dirichlet_nodes[i]].Set_To_Zero();
        x[dirichlet_nodes[i]] = X_3d[dirichlet_nodes[i]] + dirichlet_displacement[i];
    }
}
    
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Function: add_ball_collision
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T, class MeshType>
void Neo_Hookean_Cloth_3d_Fvm_Explicit<T, MeshType>::
add_ball_collision(T xb, T yb, T zb, T rb, T spring_constant) {
    if (!use_ball_collision) {
        use_ball_collision = true; }
    ball_radius = rb;
    ball_center(0) = xb;
    ball_center(1) = yb;
    ball_center(2) = zb;

#pragma omp parallel for schedule(static)
    for (unsigned int i=0; i<mesh.nodes.size(); i++) {
        node_3d_type my_position(x[i](0), x[i](1), x[i](2));
        node_3d_type center2position = my_position - ball_center;
        T dist = center2position.Magnitude();
        if (dist < rb) {
            node_3d_type push_normal = center2position * (1.0/dist);
            f[i] = f[i] + push_normal*(rb - dist)*spring_constant; } }
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Function: add_ground_collision
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T, class MeshType>
void Neo_Hookean_Cloth_3d_Fvm_Explicit<T, MeshType>::
add_ground_collision(T ground_level, T spring_constant, T friction_constant)
{
    if (!use_ground_collision) {
        use_ground_collision = true; }

#pragma omp parallel for schedule(static)
    for (unsigned int i=0; i<mesh.nodes.size(); i++) {
        T my_y = x[i](1);
        if (my_y < ground_level) {
            T penetration_depth = ground_level - my_y;
            node_3d_type push_normal(0, 1, 0);
            f[i] = f[i] + push_normal*penetration_depth*spring_constant; 
            
            // friction
            node_3d_type vxz(v[i](0), 0, v[i](2));
            node_3d_type friction_normal = vxz*(-1.0);
            friction_normal.Normalize();
            node_3d_type friction = friction_normal * vxz.Magnitude() * friction_constant; 
            f[i] = f[i] + friction; } }
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Function: add_gravity
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T, class MeshType>
void Neo_Hookean_Cloth_3d_Fvm_Explicit<T, MeshType>::
add_gravity(T gx, T gy, T gz) {
    node_3d_type g(gx, gy, gz);

#pragma omp parallel for schedule(static)
    for (unsigned int i=0; i<mesh.nodes.size(); i++) {
        f[i] = f[i] + mesh.mass[i] * g; }
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Function: set_dirichlet_with_a_bounding_box
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T, class MeshType>
void Neo_Hookean_Cloth_3d_Fvm_Explicit<T, MeshType>::
set_dirichlet_with_a_bounding_box(T x0, T xM, T y0, T yM, T z0, T zM, T xmove, T ymove, T zmove) {
    for (unsigned int i=0; i<mesh.nodes.size(); i++) {
        if ( X_3d[i](0) >= x0 && X_3d[i](0) <= xM && X_3d[i](1) >= y0 && X_3d[i](1) <= yM && X_3d[i](2) >= z0 && X_3d[i](2) <= zM ) {
            dirichlet_nodes.push_back(i);
            node_3d_type fixed_movement(xmove, ymove, zmove);
            dirichlet_displacement.push_back(fixed_movement); }}
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Function: clear_previous_dirichlet_conditions
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T, class MeshType>
void Neo_Hookean_Cloth_3d_Fvm_Explicit<T, MeshType>::
clear_previous_dirichlet_conditions() {
    dirichlet_nodes.clear();
    dirichlet_displacement.clear();
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Function: do_point_triangle_repulsion
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T, class MeshType>
void Neo_Hookean_Cloth_3d_Fvm_Explicit<T, MeshType>::
do_point_triangle_repulsion(int iterations, T d_tol) {
    for (int iter_count = 0; iter_count < iterations; iter_count++) {

        for (unsigned int p = 0; p < mesh.nodes.size(); p++) {
            node_3d_type P = x[p];
            T mp = mesh.mass[p];

            // query hierarchy to get potential intersection pairs
            std::vector<int> intersection_list;
            hierarchy.query_point(P, intersection_list);
            
            for (unsigned int iliter = 0; iliter < intersection_list.size(); iliter++) {
                int t = intersection_list[iliter];

                int node0 = mesh.elements[t](0);
                int node1 = mesh.elements[t](1);
                int node2 = mesh.elements[t](2);
                if (node0 == p || node1 == p || node2 == p) { // don't need to check element that contains node p
                    continue; }
                node_3d_type A = x[node0];
                node_3d_type B = x[node1];
                node_3d_type C = x[node2];

                ditto::geometry::Triangle_3d<T> tri(A, B, C);
                if (tri.get_area() < 1e-10) { // if the triangle is tooooo small, ignore it
                    continue; }

                // a naive bounding box
                T box_minx, box_maxx, box_miny, box_maxy, box_minz, box_maxz;
                tri.get_box(d_tol, box_minx, box_maxx, box_miny, box_maxy, box_minz, box_maxz);
                if (!( P(0)>box_minx && P(0)<box_maxx && P(1)>box_miny && P(1)<box_maxy && P(2)>box_minz && P(2)<box_maxz)) {
                    continue;}
            
                node_3d_type P_hat;
                T ksi1;
                T ksi2;

                T d = tri.find_closest_point(P, P_hat, ksi1, ksi2);
                if (d < d_tol) {
                    node_3d_type n = (P - P_hat)*(1.0/d);
                    node_3d_type v_hat = v[node0]*(1-ksi1-ksi2) + v[node1]*ksi1 + v[node2]*ksi2;
                    node_3d_type v_rel = v[p] - v_hat;
                    T v_rel_dot_n = v_rel.Dot(n);
                    if (v_rel_dot_n < 0) {
                        T m0 = mesh.mass[node0];
                        T m1 = mesh.mass[node1];
                        T m2 = mesh.mass[node2];
                        T ksi0 = 1-ksi1-ksi2;
                        T impulse = -v_rel_dot_n / (  1.0/mp +  ksi0*ksi0/m0  +  ksi1*ksi1/m1  +  ksi2*ksi2/m2  );

                        v[p] = v[p] + n*(impulse/mp);
                        v[node0] = v[node0] - n*(ksi0*impulse/m0);
                        v[node1] = v[node1] - n*(ksi1*impulse/m1);
                        v[node2] = v[node2] - n*(ksi2*impulse/m2); }}}}}
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Function: do_segment_segment_repulsion
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T, class MeshType>
void Neo_Hookean_Cloth_3d_Fvm_Explicit<T, MeshType>::
do_segment_segment_repulsion(int iterations, T d_tol) {
    for (int iter_count = 0; iter_count < iterations; iter_count++) {

        for (unsigned int first_edge_index = 0; first_edge_index < mesh.edges.size(); first_edge_index++) {
            int node_A = mesh.edges[first_edge_index](0);
            int node_B = mesh.edges[first_edge_index](1);
            T mA = mesh.mass[node_A];
            T mB = mesh.mass[node_B];
            node_3d_type A = x[node_A];
            node_3d_type B = x[node_B];
            ditto::geometry::Segment_3d<T> edgeAB(A, B);
            
            ditto::geometry::Box<T,node_3d_type> first_edge_box;
            first_edge_box.build_box_from_segment(first_edge_index, edgeAB, hierarchy.margin);
            std::vector<int> intersection_list;
            hierarchy.query_box(first_edge_box, intersection_list);

            // for (unsigned int second_edge_index = 0; second_edge_index < mesh.edges.size(); second_edge_index++) {
            //     int node_C = mesh.edges[second_edge_index](0);
            //     int node_D = mesh.edges[second_edge_index](1);

            for (unsigned int iliter = 0; iliter < intersection_list.size(); iliter++) {
                int t = intersection_list[iliter];
                int node0 = mesh.elements[t](0);
                int node1 = mesh.elements[t](1);
                int node2 = mesh.elements[t](2);
                int edges_of_triangle[][2] = { {node0, node1}, {node1, node2}, {node2, node0} };
                for (unsigned int edgeIter = 0; edgeIter < 3; edgeIter++) {
                    int node_C = edges_of_triangle[edgeIter][0];
                    int node_D = edges_of_triangle[edgeIter][1];
                    
                    if (node_C == node_A || node_C == node_B || node_D == node_A || node_D == node_B) { // don't need to check edge pair that share node
                        continue; }
                    
                    node_3d_type C = x[node_C];
                    node_3d_type D = x[node_D];
                    ditto::geometry::Segment_3d<T> edgeCD(C, D);
                    
                    // a naive bounding box
                    // T ABminx, ABmaxx, ABminy, ABmaxy, ABminz, ABmaxz;
                    // T CDminx, CDmaxx, CDminy, CDmaxy, CDminz, CDmaxz;
                    // edgeAB.get_box(d_tol, ABminx, ABmaxx, ABminy, ABmaxy, ABminz, ABmaxz);
                    // edgeCD.get_box(d_tol, CDminx, CDmaxx, CDminy, CDmaxy, CDminz, CDmaxz);
                    // if (ABminx > CDmaxx || ABmaxx < CDminx || ABminy > CDmaxy || ABmaxy < CDminy || ABminz > CDmaxz || ABmaxz < CDminz) {
                    //     continue;}
                    
                    node_3d_type P;
                    node_3d_type Q;
                    T s,t;
                    T d = edgeAB.find_closest_points_seg_seg(edgeCD, P, Q, s, t);
                    if (d < d_tol) {
                        T ksiA = 1-s;
                        T ksiB = s;
                        T ksiC = 1-t;
                        T ksiD = t;
                        node_3d_type n = (P-Q)*(1.0/d);
                        node_3d_type vP = v[node_A]*ksiA + v[node_B]*ksiB;
                        node_3d_type vQ = v[node_C]*ksiC + v[node_D]*ksiD;
                        node_3d_type v_rel = vP - vQ;
                        T v_rel_dot_n = v_rel.Dot(n);
                        if (v_rel_dot_n < 0) {
                            T mC = mesh.mass[node_C];
                            T mD = mesh.mass[node_D];
                            T impulse = -v_rel_dot_n / ( ksiA*ksiA/mA  +  ksiB*ksiB/mB  +  ksiC*ksiC/mC  +  ksiD*ksiD/mD );
                            
                            v[node_A] = v[node_A] + n*(ksiA*impulse/mA);
                            v[node_B] = v[node_B] + n*(ksiB*impulse/mB);
                            v[node_C] = v[node_C] - n*(ksiC*impulse/mC);
                            v[node_D] = v[node_D] - n*(ksiD*impulse/mD); }}}}}}
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Function: do_point_triangle_collision
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T, class MeshType>
void Neo_Hookean_Cloth_3d_Fvm_Explicit<T, MeshType>::
do_point_triangle_collision(int iterations) {
    for (int iter_count = 0; iter_count < iterations; iter_count++) {

        for (unsigned int p = 0; p < mesh.nodes.size(); p++) {
            T mp = mesh.mass[p];
            node_3d_type P = x[p]; 
            node_3d_type P_future = x[p] + v[p]*dt;

            // query hierarchy predicting future to get potential future intersection pairs
            ditto::geometry::Segment_3d<T> PnP(P, P_future);
            ditto::geometry::Box<T,node_3d_type> PnP_box;
            PnP_box.build_box_from_segment(p, PnP, hierarchy_predicting_future.margin);
            std::vector<int> intersection_list;
            hierarchy_predicting_future.query_box(PnP_box, intersection_list);

            for (unsigned int iliter = 0; iliter < intersection_list.size(); iliter++) {
                int t = intersection_list[iliter];

                int node0 = mesh.elements[t](0);
                int node1 = mesh.elements[t](1);
                int node2 = mesh.elements[t](2);
                if (node0 == p || node1 == p || node2 == p) { // don't need to check element that contains node p
                    continue; }
                node_3d_type A = x[node0];
                node_3d_type B = x[node1];
                node_3d_type C = x[node2];

                ditto::geometry::Triangle_3d<T> tri(A, B, C);
                if (tri.get_area() < 1e-10) { // if the triangle is tooooo small, ignore it
                    continue; }

                ditto::geometry::Tetrahedron<T> tet(A, B, C, P);
                if (tet.penetration_safety_test(v[node0], v[node1], v[node2], v[p], dt) == false) { 
                    node_3d_type P_hat;
                    T ksi1;
                    T ksi2;
                    T d = tri.find_closest_point(P, P_hat, ksi1, ksi2);
                    node_3d_type n = (P - P_hat)*(1.0/d);
                    node_3d_type v_hat = v[node0]*(1-ksi1-ksi2) + v[node1]*ksi1 + v[node2]*ksi2;
                    node_3d_type v_rel = v[p] - v_hat;
                    T v_rel_dot_n = v_rel.Dot(n);

                    T m0 = mesh.mass[node0];
                    T m1 = mesh.mass[node1];
                    T m2 = mesh.mass[node2];
                    T ksi0 = 1-ksi1-ksi2;
                    T impulse = -v_rel_dot_n / (  1.0/mp +  ksi0*ksi0/m0  +  ksi1*ksi1/m1  +  ksi2*ksi2/m2  );

                    v[p] = v[p] + n*(impulse/mp);
                    v[node0] = v[node0] - n*(ksi0*impulse/m0);
                    v[node1] = v[node1] - n*(ksi1*impulse/m1);
                    v[node2] = v[node2] - n*(ksi2*impulse/m2); }}}}               
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Function: switch_self_collision
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T, class MeshType>
void Neo_Hookean_Cloth_3d_Fvm_Explicit<T, MeshType>::
switch_self_collision(bool s, int input_iter, T input_dtol) {
    use_self_collision = s;
    self_collision_repulsion_iters = input_iter;
    self_collision_collision_iters = input_iter;
    self_collision_distance_tolerance = input_dtol;

    if (use_self_collision == true) {
        initialize_hierarchy(input_dtol); }
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Function: initialize_hierarchy
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T, class MeshType>
void Neo_Hookean_Cloth_3d_Fvm_Explicit<T, MeshType>::
initialize_hierarchy(T margin) {
    hierarchy.build_tree(mesh.elements, x, margin);
    hierarchy_predicting_future.build_tree(mesh.elements, x, margin);
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Function: update_hierarchy
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T, class MeshType>
void Neo_Hookean_Cloth_3d_Fvm_Explicit<T, MeshType>::
update_hierarchy() {
    hierarchy.update_box_positions(mesh.elements, x);
    hierarchy_predicting_future.update_box_positions_predicting_future(mesh.elements, x, v, dt);
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Function: write_output
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T, class MeshType>
void Neo_Hookean_Cloth_3d_Fvm_Explicit<T, MeshType>::
write_output(int frame) {
       std::stringstream ss;
       ss << frame;
       std::string filename = "output/frame"+ss.str();
       FILE* file = fopen(filename.c_str(), "w");
       fprintf(file,"NODES");  
       fprintf(file," %d\n", mesh.nodes.size());
       for (size_t i=0; i < mesh.nodes.size(); i++) {
           fprintf(file,"%f %f %f\n", x[i](0), x[i](1), x[i](2) ); }
       fprintf(file,"TRIS");  
       fprintf(file," %d\n", mesh.elements.size());
       for (size_t i=0; i < mesh.elements.size(); i++) {
           fprintf(file,"%d %d %d\n", mesh.elements[i](0), mesh.elements[i](1), mesh.elements[i](2)); }

       if (use_ball_collision) {
            fprintf(file,"BALL");
            fprintf(file," %f %f %f %f\n", ball_center(0), ball_center(1), ball_center(2), ball_radius*0.98);
       }

       fprintf(file,"NUM_IDENTICAL_COMPONENTS");
       fprintf(file," %d\n", num_of_clothes);

       fclose(file);
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Function: write_vtk
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T, class MeshType>
void Neo_Hookean_Cloth_3d_Fvm_Explicit<T, MeshType>::
write_vtk(int frame) {
       std::stringstream ss;
       ss << frame;

       // write sphere
       std::string sphere_name = "output/sphere"+ss.str()+".vtk";
       if (use_ball_collision) {
           ditto::visualization::write_vtk_sphere(ball_center(0), ball_center(1), ball_center(2), ball_radius, sphere_name);
       }

       // write cloth(es)
       for (int cloth_count = 1; cloth_count <= num_of_clothes; cloth_count++) {
           std::stringstream ss2;
           ss2 << cloth_count;
           std::string cloth_name = "output/cloth"+ss2.str()+"_"+ss.str()+".vtk";
           std::fstream out(cloth_name.c_str(), std::ios_base::out);
           out << "# vtk DataFile Version 3.0\n"
               << "vtk output\n"
               << "ASCII\n"
               << "DATASET POLYDATA\n"
               << "POINTS " << mesh.nodes.size()/num_of_clothes << " float" << std::endl;

           int start_i =  (cloth_count-1)*mesh.nodes.size()/num_of_clothes;
           int next_start_i = cloth_count*mesh.nodes.size()/num_of_clothes;
           for(int i = start_i; i < next_start_i; i++) {
               out << x[i](0) << " " << x[i](1) << " " << x[i](2) << std::endl; }

           out << "POLYGONS " << mesh.elements.size()/num_of_clothes << " " << 4*mesh.elements.size()/num_of_clothes << std::endl;
           int start_e =  (cloth_count-1)*mesh.elements.size()/num_of_clothes; 
           int next_start_e = cloth_count*mesh.elements.size()/num_of_clothes;
           for(int e = start_e; e < next_start_e; e++) {
               out << "3 " << mesh.elements[e](0)-start_i << " "
                   << mesh.elements[e](1)-start_i << " "
                   << mesh.elements[e](2)-start_i << std::endl; }
       }
}


} } // end namespaces
#endif
