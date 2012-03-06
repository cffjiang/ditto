//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// ditto/public_library/solid_2d/neo_hookean_2d_fem_implicit.h
// Copyright 2012, Chenfanfu Jiang
//
// Cloth in 3d is a special case because its material space is 2d.
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

#ifndef DITTO_PUBLIC_LIBRARY_SOLID_2D_NEO_HOOKEAN_2D_FEM_IMPLICIT_H
#define DITTO_PUBLIC_LIBRARY_SOLID_2D_NEO_HOOKEAN_2D_FEM_IMPLICIT_H

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

namespace ditto { namespace solid_2d {

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Class: Neo_Hookean_2d_Fem_Implicit
// MeshType: ditto::geometry::Triangle_Mesh_2d
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T, class MeshType>
class Neo_Hookean_2d_Fem_Implicit {
public:
    MeshType &mesh;
    T dt;
    std::vector<T> u;
    std::vector<T> u0;
    std::vector<T> v;
    std::vector<T> v0;
    std::vector<int> dirichlet_nodes;
    std::vector<T> dirichlet_u;
    T youngs_modulus;
    T poisson_ratio;
    T mu;
    T lambda;
 
    Neo_Hookean_2d_Fem_Implicit(MeshType &input_mesh, const T input_dt, const T input_E, const T input_nu)
        :mesh(input_mesh) {
        dt = input_dt;
        u.resize(mesh.nodes.size() * 2, T(0.0));
        v.resize(mesh.nodes.size() * 2, T(0.0));
        u0 = u;
        v0 = v;
        youngs_modulus = input_E;
        poisson_ratio = input_nu;
        mu = youngs_modulus / (2*(1+poisson_ratio));
        lambda = youngs_modulus * poisson_ratio / ( (1+poisson_ratio)*(1-2*poisson_ratio) );}

    void update(T newton_tol, int newton_max_iter, T minres_tol, int minres_max_iter);

    template<class SparseMatrixType, class VectorType>
    void build_linear_system(SparseMatrixType &dq_du, VectorType &delta_u, VectorType &minus_q, T minres_tol, int minres_max_iter);
    
    void write_output(int frame);
};


//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Function: update
// MeshType: ditto::geometry::Triangle_Mesh_2d
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T, class MeshType>
void Neo_Hookean_2d_Fem_Implicit<T, MeshType>::
update(T newton_tol, int newton_max_iter, T minres_tol, int minres_max_iter) {
    int system_dimension = u.size();
    static int nframe = 1;
    std::cout << "---- Frame " << nframe << " -----------------------------------\n";
    //################## Newton Iterations ############################################
    int N_iters = 0;
    int max_nlinear = 0;
    for (unsigned int newton_iter = 1; newton_iter <= newton_max_iter; newton_iter++) {
        ditto::algebra::SPARSE_MATRIX<T> dq_du(system_dimension, system_dimension);
        ditto::algebra::VECTOR<T> minus_q(system_dimension);
        ditto::algebra::VECTOR<T> delta_u(system_dimension);
        build_linear_system(dq_du, delta_u, minus_q, minres_tol, minres_max_iter);
        ditto::algebra::MINRES<T> minres_solver(dq_du, delta_u, minus_q, minres_max_iter);
        minres_solver.Set_Tolerance(minres_tol);
        int nlinear = minres_solver.Solve();
        if (nlinear > max_nlinear) {
            max_nlinear = nlinear; }
        N_iters++;
        for (unsigned int i=0; i<u.size(); i++) {
            u[i] += delta_u(i); }
        if (std::abs(delta_u.Max()) < newton_tol && std::abs(delta_u.Min()) < newton_tol) {
            break; } }
    //########################################################### End Newton Iterations
    std::cout << " Newton iterations: " << N_iters << std::endl;
    std::cout << " Linear iterations: " << max_nlinear << std::endl;
    for (unsigned int i=0; i<u.size(); i++) {
        v[i] = (u[i] - u0[i])/dt; }
    u0 = u;
    v0 = v;
    nframe++;
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Function: build_linear_system
// MeshType: ditto::geometry::Triangle_Mesh_2d
// SparseMatrixType: ditto::algebra::SPARSE_MATRIX
// VectorType: ditto::algebra::VECTOR
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T, class MeshType>
template<class SparseMatrixType, class VectorType>
void Neo_Hookean_2d_Fem_Implicit<T, MeshType>::
build_linear_system(SparseMatrixType &dq_du, VectorType &delta_u, VectorType &minus_q, T minres_tol, int minres_max_iter) {
    for (unsigned int alpha = 0; alpha < mesh.elements.size(); alpha++) {
        ditto::algebra::MATRIX_MXN<T> dN_dX(3,2);
        ditto::algebra::MATRIX_3X3<T> dN_dX_first;
        for(int j=0; j<3; j++) {
            dN_dX_first(0,j) = mesh.nodes[mesh.elements[alpha](j)](0);
            dN_dX_first(1,j) = mesh.nodes[mesh.elements[alpha](j)](1);
            dN_dX_first(2,j) = 1.0;}
        dN_dX_first.Invert();
        dN_dX(0,0) = dN_dX_first(0,0);
        dN_dX(0,1) = dN_dX_first(0,1);
        dN_dX(1,0) = dN_dX_first(1,0);
        dN_dX(1,1) = dN_dX_first(1,1);
        dN_dX(2,0) = dN_dX_first(2,0);
        dN_dX(2,1) = dN_dX_first(2,1);
        ditto::algebra::MATRIX_MXN<T> dN_dX_transpose(2,3);
        dN_dX.Transpose(dN_dX_transpose);
        ditto::algebra::MATRIX_MXN<T> F(2,2);
        ditto::algebra::MATRIX_MXN<T> u_tri(3,2);
        ditto::algebra::MATRIX_MXN<T> u_tri_transpose(2,3);
        for(int j=0; j<2; j++) {
            u_tri(0,j) = u[2*mesh.elements[alpha](0)+j];
            u_tri(1,j) = u[2*mesh.elements[alpha](1)+j];
            u_tri(2,j) = u[2*mesh.elements[alpha](2)+j];}
        u_tri.Transpose(u_tri_transpose);
        ditto::algebra::Multiply(u_tri_transpose, dN_dX, F);
        for(int j=0; j<2; j++) {
            F(j,j) += 1.0;}
        ditto::algebra::MATRIX_MXN<T> P(2,2);
        ditto::algebra::MATRIX_MXN<T> P_transpose(2,2);
        T J;
        ditto::algebra::MATRIX_MXN<T> dJ_dF(2,2);
        T rJ;
        T drJ;
        T d2rJ;
        T c1;
        T c2;
        J = F(1,1)*F(0,0) - F(0,1)*F(1,0);
        dJ_dF(0,0) = F(1,1);
        dJ_dF(0,1) = -F(1,0);
        dJ_dF(1,0) = -F(0,1);
        dJ_dF(1,1) = F(0,0);
        rJ = (J-1.0) - (J-1.0)*(J-1.0)/2.0 + (J-1.0)*(J-1.0)*(J-1.0)/3.0 - (J-1.0)*(J-1.0)*(J-1.0)*(J-1.0)/4.0;
        drJ = 1.0- (J-1.0) + (J-1.0)*(J-1.0) - (J-1.0)*(J-1.0)*(J-1.0);
        d2rJ = -1.0 + 2.0*(J-1.0) - 3.0*(J-1.0)*(J-1.0);
        c1 = lambda*drJ*drJ + (lambda*rJ-mu)*d2rJ;
        c2 = (lambda*rJ-mu)*drJ;
        for(int i=0; i<2; i++) {
            for(int j=0; j<2; j++) {
                P(i,j) = mu*F(i,j) + (lambda*rJ-mu)*drJ*dJ_dF(i,j);}}
        P.Transpose(P_transpose);
        ditto::algebra::MATRIX_MXN<T> dN_dX_Ptrans(3,2);
        ditto::algebra::Multiply(dN_dX, P_transpose, dN_dX_Ptrans);

        // build minus_q
        for (int a=0; a<3; a++) {
            int node_a = mesh.elements[alpha](a);
            for (int b=0; b<3; b++) {
                int node_b = mesh.elements[alpha](b);
                T int_Na_Nb = (a == b)? (1.0/6.0 * mesh.area[alpha]) : (1.0/12.0 * mesh.area[alpha]);
                for (int i=0; i<2; i++) {
                    minus_q(node_a*2+i) += int_Na_Nb * mesh.rho[alpha] * (u0[node_b*2+i] + dt*v0[node_b*2+i]);
                    minus_q(node_a*2+i) -= int_Na_Nb * mesh.rho[alpha] * u[node_b*2+i];
                    minus_q(node_a*2+i) -= dt * dt * mesh.area[alpha] * dN_dX_Ptrans(a,i); }}}

        // build dq_du
        for (int i=0; i<2; i++) {
            for (int j=0; j<2; j++) {
                ditto::algebra::MATRIX_MXN<T> dPik_dFjm(2,2);
                dPik_dFjm(0,0) = c1*dJ_dF(i,0)*dJ_dF(0,j);
                dPik_dFjm(0,1) = c1*dJ_dF(i,0)*dJ_dF(1,j);
                dPik_dFjm(1,0) = c1*dJ_dF(i,1)*dJ_dF(0,j);
                dPik_dFjm(1,1) = c1*dJ_dF(i,1)*dJ_dF(1,j);
                if (i == j) dPik_dFjm(i,j) += mu;
                if (i != j) {
                    dPik_dFjm(i,j) += c2;
                    dPik_dFjm(j,i) -= c2; }
                ditto::algebra::MATRIX_MXN<T> local_dq_du(3,3);
                ditto::algebra::MATRIX_MXN<T> local_dq_du_temp(3,2);
                ditto::algebra::Multiply(dN_dX,dPik_dFjm,local_dq_du_temp);				
                ditto::algebra::Multiply(local_dq_du_temp,dN_dX_transpose,local_dq_du);
                local_dq_du.Scale_By(dt*dt);
                if(i==j) {
                    ditto::algebra::MATRIX_MXN<T> scale_matrix(3,3);
                    scale_matrix.Set_To_Value(1.0);
                    for(int iii=0 ; iii<3; iii++)
                        scale_matrix(iii,iii) = 2.0;										
                    scale_matrix.Scale_By(mesh.rho[alpha]/12.0);
                    for(int a=0; a<3; a++) {
                        for(int b=0; b<3; b++) {												
                            local_dq_du(a,b) += scale_matrix(a,b);}}}					
                local_dq_du.Scale_By(mesh.area[alpha]);
                for(int a=0; a<3; a++) {
                    int node_a = mesh.elements[alpha](a);
                    for(int b=0; b<3; b++){
                        int node_b = mesh.elements[alpha](b);
                        int r = 2*node_a + i;
                        int c = 2*node_b + j;
                        dq_du.Add_To_Or_Create_Entry(r, c, local_dq_du(a,b)); }}}}}
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Function: write_output
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T, class MeshType>
void Neo_Hookean_2d_Fem_Implicit<T, MeshType>::
write_output(int frame) {
    std::stringstream ss;
    ss << frame;
    std::string filename = "output/frame"+ss.str();
    FILE* file = fopen(filename.c_str(), "w");
     fprintf(file,"NODES");  
    fprintf(file," %d\n", mesh.nodes.size());
     for (size_t i=0; i < mesh.nodes.size(); i++) {
        fprintf(file,"%f %f %f\n", mesh.nodes[i](0)+u[2*i+0], mesh.nodes[i](1)+u[2*i+1], 0.0); }
     fprintf(file,"TRIS");  
    fprintf(file," %d\n", mesh.elements.size());
     for (size_t i=0; i < mesh.elements.size(); i++) {
        fprintf(file,"%d %d %d\n", mesh.elements[i](0), mesh.elements[i](1), mesh.elements[i](2)); }
     fclose(file);
}


} } // end namespaces
#endif
