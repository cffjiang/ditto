//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// ditto/public_library/cloth_3d/neo_hookean_cloth_3d_fvm_explicit.h
// Copyright 2012, Chenfanfu Jiang
//
// Require NEWMAT library compiled
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
#include <ditto/public_library/algebra/newmat11/newmatap.h>
#include <ditto/public_library/algebra/newmat11/newmat.h>
#include <ditto/public_library/algebra/newmat11/newmatio.h>

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

    matrix_2x2_list_type Dm_inverse_list;
    
    std::vector<int> dirichlet_nodes;
    node_3d_list_type dirichlet_displacement;

    T youngs_modulus;
    T poisson_ratio;
    T mu;
    T lambda;
 
    Neo_Hookean_Cloth_3d_Fvm_Explicit(MeshType &input_mesh, const T input_dt, const T input_E, const T input_nu, const bool input_use_gravity)
        :mesh(input_mesh) {
        dt = input_dt;
        X_3d = mesh.nodes;
        x = mesh.nodes;
        v.resize(x.size(), T(0.0));
        youngs_modulus = input_E;
        poisson_ratio = input_nu;
        mu = youngs_modulus / (2*(1+poisson_ratio));
        lambda = youngs_modulus * poisson_ratio / ( (1+poisson_ratio)*(1-2*poisson_ratio) );
        use_gravity = input_use_gravity;
        
        pre_build_Dm_inverse();}

    void pre_build_Dm_inverse();

    void update();

    void add_gravity();

    void set_dirichlet_with_a_bounding_box(T x0, T xM, T y0, T yM, T z0, T zM, T xmove, T ymove, T zmove);

    void write_output(int frame);
};

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
        node_3d_type edge_cross = X2mX1.Cross_Product(X3-X1);
        node_3d_type nvec = edge_cross*(1.0/edge_cross.Magnitude());
        node_3d_type nxX2mX1 = nvec.Cross_Product(X2mX1);
        node_3d_type wvec = X2mX1*(1.0/X2mX1.Magnitude());
        node_3d_type vvec = nxX2mX1*(1.0/nxX2mX1.Magnitude());

        matrix_3x3_type R;
        R(0,0) = wvec(0); R(0,1) = wvec(1); R(0,2) = wvec(2); 
        R(1,0) = vvec(0); R(1,1) = vvec(1); R(1,2) = vvec(2); 
        R(2,0) = nvec(0); R(2,1) = nvec(1); R(2,2) = nvec(2); 

        node_3d_type X1_mapped_3d = R*(X1-X1);
        node_3d_type X2_mapped_3d = R*(X2-X1);
        node_3d_type X3_mapped_3d = R*(X3-X1);
        
        // debug code
        T eps = 1e-10;
        assert(std::abs(X1_mapped_3d(2)) < eps  && std::abs(X2_mapped_3d(2)) < eps  && std::abs(X3_mapped_3d(2)) < eps);
    
        node_2d_type X1_mapped(X1_mapped_3d(0), X1_mapped_3d(1));
        node_2d_type X2_mapped(X2_mapped_3d(0), X2_mapped_3d(1));
        node_2d_type X3_mapped(X3_mapped_3d(0), X3_mapped_3d(1));

        // now build Dm_inverse
        matrix_2x2_type Dm_inverse(X2_mapped-X1_mapped, X3_mapped-X1_mapped);
        Dm_inverse.Invert();

        Dm_inverse_list.push_back(Dm_inverse);
    }
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Function: update
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T, class MeshType>
void Neo_Hookean_Cloth_3d_Fvm_Explicit<T, MeshType>::
update() {

    static int nframe = 1;
    std::cout << "---- Frame " << nframe << " -----------------------------------\n";

    // build f
    f.clear();
    f.resize(mesh.nodes.size());
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
        ditto::algebra::MATRIX_MXN<T> Dm_inverse(2,2);
        for (int i=0; i<2; i++) {
            for (int j=0; j<2; j++) {
                Dm_inverse(i,j) = Dm_inverse_list[e](i,j); }}

        // F = Ds*Dm_inverse
        ditto::algebra::MATRIX_MXN<T> F(3,2);
        ditto::algebra::Multiply(Ds, Dm_inverse, F);

        // SVD
        NEWMAT::Matrix A(3,2); A = 0.0;
        NEWMAT::Matrix U(3,2); U = 0.0;
        NEWMAT::Matrix V(2,2); V = 0.0;
        NEWMAT::DiagonalMatrix D(2);
        for (int i=0; i<3; i++) {
            for (int j=0; j<2; j++) {
                A.element(i,j) = F(i,j); }}
        NEWMAT::SVD(A, D, U, V);

        // Compute P
        NEWMAT::DiagonalMatrix dPhidSigma(2);
        dPhidSigma.element(0) = 2*mu*(D.element(0)-1) + lambda*(D.element(0)*D.element(1)-1)*D.element(1);
        dPhidSigma.element(1) = 2*mu*(D.element(1)-1) + lambda*(D.element(0)*D.element(1)-1)*D.element(0);
        NEWMAT::Matrix P(3,2); P = 0.0;
        P = U * dPhidSigma * (V.t());

        // Compute Force matrix
        T my_area = mesh.area[e];
        //T my_area = 0.5 * (Dm_inverse(0,0)* Dm_inverse(1,1) -  Dm_inverse(0,1)* Dm_inverse(1,0));
        NEWMAT::Matrix Dminvt(2,2); Dminvt = 0.0;
        for (int i=0; i<2; i++) {
            for (int j=0; j<2; j++) {
                Dminvt.element(i,j) = Dm_inverse_list[e](j,i); }}
        NEWMAT::Matrix Force(3,2); Force = 0.0;
        Force = -my_area * P * Dminvt;
        
        // Add force contributions to nodes
        for (int i=0; i<3; i++) {
            f[mesh.elements[e](1)](i) += Force.element(i,0); 
            f[mesh.elements[e](2)](i) += Force.element(i,1); }
        f[mesh.elements[e](0)] = f[mesh.elements[e](0)] + (f[mesh.elements[e](1)] +  f[mesh.elements[e](2)]) * (-1.0);
    }

    // add gravity
    if (use_gravity) add_gravity();

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

    nframe++;
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Function: add_gravity
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T, class MeshType>
void Neo_Hookean_Cloth_3d_Fvm_Explicit<T, MeshType>::
add_gravity() {
    node_3d_type g(0.0, 0.0, -9.8);
    for (unsigned int i=0; i<mesh.nodes.size(); i++) {
        f[i] = f[i] + mesh.mass[i] * g; }
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Function: write_output
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
       fclose(file);
}


} } // end namespaces
#endif
