/*
 *  DEFORMABLE_OBJECT_3D.h
 *
 *  Created by Joseph Teran on 12/27/10.
 *  Copyright 2010 UCLA. All rights reserved.
 *
 */
#ifndef _deformable_object_
#define _deformable_object_

#include "../source/GEOMETRY.h"

using namespace GEOMETRY;
using namespace ALGEBRA;


template <class T>
class HYPERELASTICITY_CONSTITUTIVE_MODEL_3D{
public:
    virtual MATRIX_3X3<T> P(const MATRIX_3X3<T>& F){return MATRIX_3X3<T>();}
    virtual MATRIX_3X2<T> P(const MATRIX_3X2<T>& F){return MATRIX_3X2<T>();}
    virtual MATRIX_3X3<T> dP(const MATRIX_3X3<T>& F,const MATRIX_3X3<T>& dF){return MATRIX_3X3<T>();}
    virtual MATRIX_3X2<T> dP(const MATRIX_3X2<T>& F,const MATRIX_3X2<T>& dF){return MATRIX_3X2<T>();}
    virtual void Element_Stifness_Matrix(const int element,const MATRIX_3X3<T>& Dm_Inverse,MATRIX_MXN<T>& element_stiffness){}
    virtual void Update_Position_Based_State(){}
    int Index(const int node, const int component){
        return 3*node+component;
    }
    VECTOR_3D<T> Natural_Interpolating_Function_Gradient(const int node_number){
        assert(node_number>=0);
        assert(node_number<4);
		
        if(node_number==0) return VECTOR_3D<T>(-1,-1,-1);
        else if(node_number==1) return VECTOR_3D<T>(1,0,0);
        else if(node_number==2) return VECTOR_3D<T>(0,1,0);
        else if(node_number==3) return VECTOR_3D<T>(0,0,1);
        else return VECTOR_3D<T>(0,0,0);		
    }
};

template <class T>
class LINEAR_ELASTICITY_3D:public HYPERELASTICITY_CONSTITUTIVE_MODEL_3D<T>{
    T lambda;
    T mu;
	
  public:
    LINEAR_ELASTICITY_3D(const T youngs_modulus,const T poisson_ratio){
        lambda=youngs_modulus*poisson_ratio/(((T)1+poisson_ratio)*((T)1-(T)2*poisson_ratio));
        mu=youngs_modulus/((T)2*((T)1+poisson_ratio));
    }
	
    MATRIX_3X3<T> P(const MATRIX_3X3<T>& F){
        MATRIX_3X3<T> epsilon=(T).5*(F+F.Transposed())-MATRIX_3X3<T>::Identity();
        return 2*mu*epsilon+lambda*epsilon.Trace()*MATRIX_3X3<T>::Identity();
    }
	
    MATRIX_3X3<T> dP(const MATRIX_3X3<T>& F, const MATRIX_3X3<T>& dF){
        MATRIX_3X3<T> d_epsilon=(T).5*(dF+dF.Transposed());
        return 2*mu*d_epsilon+lambda*d_epsilon.Trace()*MATRIX_3X3<T>::Identity();
    }
	
    void Element_Stifness_Matrix(const int element,const MATRIX_3X3<T>& Dm_Inverse,MATRIX_MXN<T>& element_stiffness){
        MATRIX_3X3<T> M=Dm_Inverse*Dm_Inverse.Transposed();
		
        T entry_ab;
        T jacobian_determinant=(T)1/Dm_Inverse.Determinant();
		
        for(int a=0;a<4;a++){
            VECTOR_3D<T> grad_Na=HYPERELASTICITY_CONSTITUTIVE_MODEL_3D<T>::Natural_Interpolating_Function_Gradient(a);
            for(int b=0;b<4;b++){
                VECTOR_3D<T> grad_Nb=HYPERELASTICITY_CONSTITUTIVE_MODEL_3D<T>::Natural_Interpolating_Function_Gradient(b);
                VECTOR_3D<T> image=M*grad_Nb;
                entry_ab=-(mu/(T)6)*jacobian_determinant*VECTOR_3D<T>::Dot_Product(grad_Na,image);
                MATRIX_3X3<T> diag=entry_ab*MATRIX_3X3<T>::Identity();
                for(int i=0;i<3;i++){
                    for(int j=0;j<3;j++){
                        element_stiffness(HYPERELASTICITY_CONSTITUTIVE_MODEL_3D<T>::Index(a,i),HYPERELASTICITY_CONSTITUTIVE_MODEL_3D<T>::Index(b,j))=diag(i,j);}}}}
	
        for(int a=0;a<4;a++){
            VECTOR_3D<T> grad_Na=Dm_Inverse.Transposed()*HYPERELASTICITY_CONSTITUTIVE_MODEL_3D<T>::Natural_Interpolating_Function_Gradient(a);
            for(int i=0;i<3;i++){
                for(int b=0;b<4;b++){
                    VECTOR_3D<T> grad_Nb=Dm_Inverse.Transposed()*HYPERELASTICITY_CONSTITUTIVE_MODEL_3D<T>::Natural_Interpolating_Function_Gradient(b);
                    for(int j=0;j<3;j++){
                        element_stiffness(HYPERELASTICITY_CONSTITUTIVE_MODEL_3D<T>::Index(a,i),HYPERELASTICITY_CONSTITUTIVE_MODEL_3D<T>::Index(b,j))-=((mu+lambda)*jacobian_determinant/(T)6)*grad_Na(i)*grad_Nb(j);}}}}
    }
};

template <class T>
class LAGRANGIAN_FORCES_3D{
    VECTOR<T> forces;
    SPARSE_MATRIX<T> stiffness_matrix;//this is -df/dx(x_k)
    bool stiffness_matrix_initialized;
    VECTOR<T> delta_x;//for linearizing the forces around the current configutaion: f_lin(x+dx)=df/dx(x)dx+f(x)
    int number_nodes;
public:
    LAGRANGIAN_FORCES_3D(const int number_nodes_input):forces(3*number_nodes_input),
                                                      stiffness_matrix(3*number_nodes_input,3*number_nodes_input),
                                                      stiffness_matrix_initialized(false),
                                                      delta_x(3*number_nodes_input),number_nodes(number_nodes_input){}
    virtual void Update_Position_Based_State(){}
    virtual void Compute_Position_Based_Forces(){}
	
    virtual bool Stiffness_Matrix_Initialized() const {return stiffness_matrix_initialized;}
    virtual void Initialize_Stiffness_Matrix(){stiffness_matrix_initialized=true;}
	
    VECTOR<T>& Forces(){return forces;}
    SPARSE_MATRIX<T>& Stiffness_Matrix(){return stiffness_matrix;}
    VECTOR<T>& Delta_X(){return delta_x;}
	
    virtual int Index(const int node, const int component){
        return 3*node+component;
    }
	
    void Zero_Out_Forces(){
        forces.Set_To_Zero();
    }
	
    void Add_Force(const int particle_index,const VECTOR_3D<T>& force){
        forces(3*particle_index)+=force.x();
        forces(3*particle_index+1)+=force.y();
        forces(3*particle_index+2)+=force.z();
    }
};

template <class T>
class DEFORMABLE_OBJECT_3D{
    //State variables/geom
    //Currently can only have either a tet mesh or a tri mesh.
    TETRAHEDRON_MESH tet_mesh;//mesh type = 0 (see constructor)
    TRIANGLE_MESH tri_mesh;//mesh type = 1 (see constructor)
    VECTOR<T> positions;
    VECTOR<T> velocities;
    int number_nodes;
	
public:
    DEFORMABLE_OBJECT_3D(const int num_elements_input,const int number_nodes_input,const int mesh_type_input):
	tet_mesh(mesh_type_input==0?num_elements_input:0,number_nodes_input),tri_mesh(mesh_type_input==1?num_elements_input:0,number_nodes_input),positions(3*number_nodes_input),velocities(3*number_nodes_input),number_nodes(number_nodes_input)
    {}
	
    DEFORMABLE_OBJECT_3D(GRID_3D<T>& grid):
	tet_mesh(grid.M(),grid.N(),grid.MN()),tri_mesh(0,grid.M()*grid.N()*grid.MN()),positions(3*grid.M()*grid.N()*grid.MN()),velocities(3*grid.M()*grid.N()*grid.MN()),number_nodes(grid.M()*grid.N()*grid.MN())
    {}
	
    VECTOR<T>& Positions(){return positions;}
    VECTOR<T>& Velocities(){return velocities;}
	
    TETRAHEDRON_MESH& Tetrahedron_Mesh(){return tet_mesh;}
    TRIANGLE_MESH& Triangle_Mesh(){return tri_mesh;}
	
    void Set_Position(const int particle_index,const VECTOR_3D<T>& x){
        positions(3*particle_index)=x.x();
        positions(3*particle_index+1)=x.y();
        positions(3*particle_index+2)=x.z();
    }
	
    void Set_Velocity(const int particle_index,const VECTOR_3D<T>& v){
        velocities(3*particle_index)=v.x();
        velocities(3*particle_index+1)=v.y();
        velocities(3*particle_index+2)=v.z();
    }
	
    VECTOR_3D<T> X(const int particle_index){return VECTOR_3D<T>(positions(3*particle_index),positions(3*particle_index+1),positions(3*particle_index+2));}
    VECTOR_3D<T> V(const int particle_index){return VECTOR_3D<T>(velocities(3*particle_index),velocities(3*particle_index+1),velocities(3*particle_index+2));}
	
};

template <class T>
class BACKWARD_EULER_TIME_STEPPING_3D{
    T dt,final_time,start_time,time;
    DEFORMABLE_OBJECT_3D<T>& deformable_object;
    LAGRANGIAN_FORCES_3D<T>* elastic_forces;
    CONJUGATE_GRADIENT<T>* cg;
    SPARSE_MATRIX<T> be_matrix;
    VECTOR<T> be_rhs;
    VECTOR<T> mass;
    T rho;
    bool matrix_intialized;
    bool use_gravity;
public:
    BACKWARD_EULER_TIME_STEPPING_3D(const T dt_input,const T final_time_input, const T start_time_input,DEFORMABLE_OBJECT_3D<T>& object_input):
	dt(dt_input),final_time(final_time_input),start_time(start_time_input),deformable_object(object_input),elastic_forces(0),cg(0),
       be_matrix(deformable_object.Positions().Size(),deformable_object.Positions().Size()),
       mass(deformable_object.Positions().Size()),be_rhs(deformable_object.Positions().Size()),
       rho((T)1000),matrix_intialized(false),use_gravity(true){
        time=start_time;
    }
	
    ~BACKWARD_EULER_TIME_STEPPING_3D(){if(cg) delete cg;}
	
    SPARSE_MATRIX<T>& BE_Matrix(){
        return be_matrix;
    }
	
    void Initialize_BE_Matrix(const VECTOR<T>& nodal_volume){
        assert(elastic_forces);
		
        if(!elastic_forces->Stiffness_Matrix_Initialized()) elastic_forces->Initialize_Stiffness_Matrix();
		
        SPARSE_MATRIX<T>& elasticity_stiffness_matrix=elastic_forces->Stiffness_Matrix();
		
        for(int i=0;i<elasticity_stiffness_matrix.M();i++){
            SPARSE_ROW<T>& k_row=elasticity_stiffness_matrix.Row(i);
            SPARSE_ROW<T>& be_row=be_matrix.Row(i);
            for(int j=0;j<k_row.Number_Nonzero();j++){
                T kij=k_row.Value_At_Sparse_Index(j);
                be_row.Add_Entry(k_row.Index(j),-dt*dt*kij);}}
		
        for(int i=0;i<nodal_volume.Size();i++) mass(i)=rho*nodal_volume(i);
		
        for(int i=0;i<be_matrix.M();i++){
            SPARSE_ROW<T>& row=be_matrix.Row(i);
            row(i)+=mass(i);}
		
        matrix_intialized=true;
    }
	
    void Initialize_CG(){
        cg=new CONJUGATE_GRADIENT<T>(be_matrix,elastic_forces->Delta_X(),be_rhs,200);
    }
	
    void Set_Elastic_Forces(LAGRANGIAN_FORCES_3D<T>& forces_input){
        elastic_forces=&forces_input;
    }
	
    void Update_BE_RHS_And_System_Matrix(){
        elastic_forces->Update_Position_Based_State();//this updates the elasticity stiffness matrix and computes elastic forces at the current configutation
        VECTOR<T>& f_xn=elastic_forces->Forces();
        VECTOR<T>& vn=deformable_object.Velocities();
		
        //first set the RHS
        for(int i=0;i<be_rhs.Size();i++) be_rhs(i)=-dt*dt*f_xn(i)+dt*mass(i)*vn(i);
        if(use_gravity) for(int i=0;i<be_rhs.Size()/3;i++) be_rhs(3*i+1)-=dt*dt*mass(3*i+1)*9.8;
		
        //now update the BE matrix
        assert(matrix_intialized);
        be_matrix.Zero_Out_Without_Changing_Sparsity();
		
        SPARSE_MATRIX<T>& elasticity_stiffness_matrix=elastic_forces->Stiffness_Matrix();
		
        for(int i=0;i<elasticity_stiffness_matrix.M();i++){
            SPARSE_ROW<T>& k_row=elasticity_stiffness_matrix.Row(i);
            SPARSE_ROW<T>& be_row=be_matrix.Row(i);
            for(int j=0;j<k_row.Number_Nonzero();j++){
                T kij=k_row.Value_At_Sparse_Index(j);
                be_row(k_row.Index(j))-=dt*dt*kij;}}
		
        for(int i=0;i<be_matrix.M();i++){
            SPARSE_ROW<T>& be_row=be_matrix.Row(i);
            be_row(i)+=mass(i);}
    }
	
    void Advance_One_Time_Step(){
        //assuming the user set the boundary conditions already
        assert(elastic_forces);

        Update_BE_RHS_And_System_Matrix();
        //the unknowns are the positions
        Solve_Linearized_System();//this just does one step of Newton
        for(int p=0;p<deformable_object.Positions().Size();p++){
            deformable_object.Positions()(p)+=elastic_forces->Delta_X()(p);
            deformable_object.Velocities()(p)=elastic_forces->Delta_X()(p)/dt;}
		
        time+=dt;
    }
	
    T Time(){return time;}
	
    void Solve_Linearized_System(){if(cg) cg->Solve();}
	
    void Set_Boundary_Conditions(VECTOR<int>& constrained_nodes,VECTOR<T>& constrained_node_locations){
        if(cg) cg->Set_Dirichlet_Dofs(constrained_nodes);
        for(int i=0;i<constrained_nodes.Size();i++){
            deformable_object.Positions()(constrained_nodes(i))=constrained_node_locations(i);
            elastic_forces->Delta_X()(constrained_nodes(i))=(T)0;}
    }
};

template <class T>
class QUASISTATIC_TIME_STEPPING_3D{
    T dt,final_time,start_time,time;
    DEFORMABLE_OBJECT_3D<T>& deformable_object;
    LAGRANGIAN_FORCES_3D<T>* elastic_forces;
    CONJUGATE_GRADIENT<T>* cg;
public:
    QUASISTATIC_TIME_STEPPING_3D(const T dt_input,const T final_time_input, const T start_time_input,DEFORMABLE_OBJECT_3D<T>& object_input):
	dt(dt_input),final_time(final_time_input),start_time(start_time_input),deformable_object(object_input),elastic_forces(0),cg(0){
        time=start_time;
    }
	
    ~QUASISTATIC_TIME_STEPPING_3D(){if(cg) delete cg;}
	
    void Initialize_CG(SPARSE_MATRIX<T>& system_matrix,VECTOR<T>&x,VECTOR<T>&b){
        cg=new CONJUGATE_GRADIENT<T>(system_matrix,x,b,10*b.Size());
    }
	
    void Set_Elastic_Forces(LAGRANGIAN_FORCES_3D<T>& forces_input){
        elastic_forces=&forces_input;
    }
	
    void Advance_One_Time_Step(){
        //assuming the user set the boundary conditions already
        assert(elastic_forces);
        elastic_forces->Update_Position_Based_State();//this updates the matrix and computes forces at the current configutation
        Solve_For_Equilibrium_Configuration();//this just does one step of Newton
        for(int p=0;p<deformable_object.Positions().Size();p++) deformable_object.Positions()(p)+=elastic_forces->Delta_X()(p);
        time+=dt;
    }
	
    T Time(){return time;}
	
    void Solve_For_Equilibrium_Configuration(){if(cg) cg->Solve();}
	
    void Set_Boundary_Conditions(VECTOR<int>& constrained_nodes,VECTOR<T>& constrained_node_locations){
        if(cg) cg->Set_Dirichlet_Dofs(constrained_nodes);
        for(int i=0;i<constrained_nodes.Size();i++){
            deformable_object.Positions()(constrained_nodes(i))=constrained_node_locations(i);
            elastic_forces->Delta_X()(constrained_nodes(i))=(T)0;}
    }
};

#endif

