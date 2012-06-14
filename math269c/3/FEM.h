#include "../poisson_1d_math_270a/algebra.h"
#include "GEOMETRY.h"

template<class T>
class FEM_Poisson{
	VECTOR<VECTOR_2D<T> >& nodes;
	VECTOR<VECTOR<int>*>& elements;
	TRIANGLE_MESH mesh;
	VECTOR<MATRIX_2X2<T> > element_inverse_Jacobians;
	VECTOR<T> u,F,f_cell_centered;
	int number_nodes,number_elements;
	SPARSE_MATRIX<T> A;
	VECTOR<int>* dirichlet_nodes;
	VECTOR<VECTOR_2D<T> > grad_Ni_hats;

public:
	FEM_Poisson(VECTOR<VECTOR_2D<T> >& nodes_input,VECTOR<VECTOR<int>*>& elements_input):nodes(nodes_input),elements(elements_input),mesh(elements_input,nodes_input.Size()),element_inverse_Jacobians(elements_input.Size()),
	f_cell_centered(elements_input.Size()),dirichlet_nodes(0),number_nodes(nodes_input.Size()),number_elements(elements_input.Size()),u(nodes_input.Size()),F(nodes_input.Size()),A(nodes_input.Size(),nodes_input.Size()),grad_Ni_hats(3){
		grad_Ni_hats(0)=VECTOR_2D<T>((double)-1,(double)-1);
		grad_Ni_hats(1)=VECTOR_2D<T>((double)1,(double)0);
		grad_Ni_hats(2)=VECTOR_2D<T>((double)0,(double)1);
		mesh.Intialize_Oriented_Boundary_Segments();
		Initialize_System_Matrix();
	}
	
	TRIANGLE_MESH& Mesh(){return mesh;}
	
	VECTOR<int>& Nodes_Of_Element(const int element){return (*(elements(element)));}
	
	int Number_Of_Nodes(){return number_nodes;}
	
	int Number_Of_Elements(){return number_elements;}
	
	VECTOR<T>& RHS(){return F;}
	
	VECTOR<T>& Solution(){return u;}
	
	SPARSE_MATRIX<T>& System_Matrix(){return A;}
	
	VECTOR_2D<T> X_At_Node(const int i){return nodes(i);}
	
	VECTOR_2D<T> X_At_Cell_Center(const int element){
		VECTOR<int>& nodes_in_element=(*(elements(element)));
		VECTOR_2D<T> v1=nodes(nodes_in_element(0));
		VECTOR_2D<T> v2=nodes(nodes_in_element(1));
		VECTOR_2D<T> v3=nodes(nodes_in_element(2));
		T one_third=(T)1/(T)3;
		return one_third*(v1+v2+v3);}
	
	void Initialize_System_Matrix(){
		for(int e=0;e<number_elements;e++){
			VECTOR<int>& nodes_in_element=(*(elements(e)));
			VECTOR_2D<T> v1=nodes(nodes_in_element(0));
			VECTOR_2D<T> v2=nodes(nodes_in_element(1));
			VECTOR_2D<T> v3=nodes(nodes_in_element(2));
			element_inverse_Jacobians(e)=MATRIX_2X2<T>(v2-v1,v3-v1);
			T J=element_inverse_Jacobians(e).Determinant();
			element_inverse_Jacobians(e).Invert();
			for(int i_e=0;i_e<3;i_e++){
				for(int j_e=0;j_e<3;j_e++){
					VECTOR_2D<T> grad_Ni_transpose=element_inverse_Jacobians(e).Transpose()*grad_Ni_hats(i_e);
					VECTOR_2D<T> grad_Nj_transpose=element_inverse_Jacobians(e).Transpose()*grad_Ni_hats(j_e);
					T A_e_i_e_j_e=(double).5*J*grad_Ni_transpose.Dot(grad_Nj_transpose);
					int i=nodes_in_element(i_e);
					int j=nodes_in_element(j_e);
					SPARSE_ROW<T>& Ai=A.Row(i);
					if(Ai.Value_Exists_At_Entry(j))
						Ai(j)+=A_e_i_e_j_e;
					else
						Ai.Add_Entry(j,A_e_i_e_j_e);}}}
	}
	
	void Print_Mesh(){
		for(int e=0;e<number_elements;e++){
			VECTOR<int>& nodes_in_element=(*(elements(e)));
			std::cout<<"Element " << e << " : "<< nodes_in_element(0) << " , " << nodes_in_element(1) << " , " << nodes_in_element(2) <<std::endl;}
	}
	
	void Print_System_Matrix(){
		std::cout << "FEM system matrix" << std::endl;
		A.Print();
	}
	
	void Write_Dat_Files(){
		FILE* fpointer;
		std::string mesh_file("mesh.dat");
		fpointer=fopen(mesh_file.c_str(),"w");
		for(int e=0;e<number_elements;e++){
			VECTOR<int>& nodes_in_element=Nodes_Of_Element(e);
			for(int i=0;i<3;i++){
				fprintf(fpointer,"%i ",nodes_in_element(i)+1);}
			fprintf(fpointer,"\n");}		
		
		fclose(fpointer);
		
		std::string verts_file("vertices.dat");
		fpointer=fopen(verts_file.c_str(),"w");
		for(int i=0;i<number_nodes;i++){
			VECTOR_2D<T>& x=nodes(i);
			fprintf(fpointer,"%g ",x.x_copy());
			fprintf(fpointer,"%g ",x.y_copy());
			fprintf(fpointer,"\n");}		
		
		fclose(fpointer);
						
		u.Write_DAT_File(std::string("u.dat"));	
		
		mesh.Write_Boundary_Mesh_To_Dat();
	}			

};