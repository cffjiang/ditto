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

    MeshType &mesh;
    T dt;

    node_3d_list_type X_3d;
    node_3d_list_type x;
    node_3d_list_type v;
    node_3d_list_type f;

    bool use_gravity;
    bool use_ball_collision;
    bool use_ground_collision;
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
        num_of_clothes = 1;
    }

    void aware_of_num_of_clothes(int N);

    void pre_build_Dm_inverse();

    void compute_elasticity();

    void update_one_step();

    void add_gravity();

    void set_dirichlet_with_a_bounding_box(T x0, T xM, T y0, T yM, T z0, T zM, T xmove, T ymove, T zmove);

    void add_ball_collision(T xb, T yb, T zb, T rb, T spring_constant);

    void add_ground_collision(T ground_level, T spring_constant);

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
    for (unsigned int i=0; i<f.size(); i++){
        f[i](0) = f[i](1) = f[i](2) = 0.0;}

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
    // forward euler
    for (unsigned int i=0; i<mesh.nodes.size(); i++) {
        v[i] = v[i] + dt * f[i] * (1.0/mesh.mass[i]);
        x[i] = x[i] + dt * v[i]; 
    }
    
    // fix positions for dirichlet nodes
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
add_ground_collision(T ground_level, T spring_constant)
{
    if (!use_ground_collision) {
        use_ground_collision = true; }
     for (unsigned int i=0; i<mesh.nodes.size(); i++) {
         T my_y = x[i](1);
         if (my_y < ground_level) {
             T penetration_depth = ground_level - my_y;
             node_3d_type push_normal(0, 1, 0);
             f[i] = f[i] + push_normal*penetration_depth*spring_constant; } }
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Function: add_gravity
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T, class MeshType>
void Neo_Hookean_Cloth_3d_Fvm_Explicit<T, MeshType>::
add_gravity() {
    node_3d_type g(0.0, -9.8, 0.0);
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
           ditto::visualization::write_vtk_sphere(ball_center(0), ball_center(1), ball_center(2), ball_radius*0.98, sphere_name);
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
