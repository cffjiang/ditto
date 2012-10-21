#ifndef _geometry_
#define _geometry_

#include "ALGEBRA.h"
#include "GRID.h"
#include <algorithm>
#include <map>

using namespace GRIDS;
using namespace ALGEBRA;

namespace GEOMETRY{
class TRIANGLE_MESH{
    VECTOR<VECTOR_3D<int> > mesh;//assuming these are right hand oriented
    VECTOR<VECTOR_2D<int> >* oriented_boundary_segments;
    int number_nodes;//assuming all nodes are in at least one triangle
    int number_tris;
    int number_boundary_segments;
public:
    TRIANGLE_MESH(const int num_tris_input,const int number_nodes_input):mesh(num_tris_input),number_tris(num_tris_input),number_nodes(number_nodes_input),
                                                                        number_boundary_segments(0),oriented_boundary_segments(0)
    {}
	
    ~TRIANGLE_MESH(){if(oriented_boundary_segments) delete oriented_boundary_segments;}
	
    int Number_Of_Triangles(){return number_tris;}
	
    int Number_Boundary_Segments(){return number_boundary_segments;}
	
    VECTOR_2D<int> Nodes_Of_Boundary_Segment(const int s){
        assert(oriented_boundary_segments);
        assert(number_boundary_segments==oriented_boundary_segments->Size());
        assert(s<number_boundary_segments && s>=0);
        return (*oriented_boundary_segments)(s);}
	
    void Intialize_Oriented_Boundary_Segments(){
        if(oriented_boundary_segments) delete oriented_boundary_segments;
		
        //assumes triangle are right hand oriented
        SPARSE_MATRIX<int> edges(number_nodes,number_nodes);//e_ij says edge ij is counter clockwise oriented in some triangle
        for(int t=0;t<number_tris;t++){
            VECTOR_3D<int>& nodes=Nodes_Of_Element(t);
            //edge 01
            SPARSE_ROW<int>& r0=edges.Row(nodes(0));
            r0.Add_Entry(nodes(1),1);
            //edge 12
            SPARSE_ROW<int>& r1=edges.Row(nodes(1));
            r1.Add_Entry(nodes(2),1);
            //edge 20
            SPARSE_ROW<int>& r2=edges.Row(nodes(2));
            r2.Add_Entry(nodes(0),1);}
		
        number_boundary_segments=0;
        for(int i=0;i<number_nodes;i++){
            SPARSE_ROW<int>& ri=edges.Row(i);
            for(int j_sp=0;j_sp<ri.Number_Nonzero();j_sp++){
                int j=ri.Index(j_sp);//egde ij is counterclockwise oriented in some triangle
                SPARSE_ROW<int>& rj=edges.Row(j);
                //check if edge ji is also in some triangle, if not then edge ij is a boundary edge
                if(!(rj.Value_Exists_At_Entry(i))) number_boundary_segments++;}}
		
        oriented_boundary_segments=new VECTOR<VECTOR_2D<int> >(number_boundary_segments);
        int count=0;
        for(int i=0;i<number_nodes;i++){
            SPARSE_ROW<int>& ri=edges.Row(i);
            for(int j_sp=0;j_sp<ri.Number_Nonzero();j_sp++){
                int j=ri.Index(j_sp);
                SPARSE_ROW<int>& rj=edges.Row(j);
                if(!(rj.Value_Exists_At_Entry(i))) (*oriented_boundary_segments)(count++)=VECTOR_2D<int>(i,j);}}
        assert(count==number_boundary_segments);
    }
	
    VECTOR_3D<int>& Nodes_Of_Element(const int element){return mesh(element);}
	
    void Write_Boundary_Mesh_To_Dat(){
        FILE* fpointer;
        std::string mesh_file("boundary_segment_mesh.dat");
        fpointer=fopen(mesh_file.c_str(),"w");
        for(int s=0;s<number_boundary_segments;s++){
            VECTOR_2D<int> nodes_in_segment=(*oriented_boundary_segments)(s);
            fprintf(fpointer,"%i ",nodes_in_segment.x_copy()+1);
            fprintf(fpointer,"%i ",nodes_in_segment.y_copy()+1);
            fprintf(fpointer,"\n");}		
			
        fclose(fpointer);
    }
};

class TETRAHEDRON_MESH{
    VECTOR<VECTOR_4D<int> > mesh;//assuming these are right hand oriented, i.e. face i,j,k has a normal that points at l
    int number_nodes;//assuming all nodes are in at least one tet
    int number_tets;
public:
    VECTOR<VECTOR_3D<int> >* boundary_triangle_mesh;
	
    TETRAHEDRON_MESH(const int nodes_in_x_dimension,const int nodes_in_y_dimension,const int nodes_in_z_dimension): 
	mesh(5*(nodes_in_x_dimension-1)*(nodes_in_y_dimension-1)*(nodes_in_z_dimension-1)),
       number_tets(5*(nodes_in_x_dimension-1)*(nodes_in_y_dimension-1)*(nodes_in_z_dimension-1)),
       number_nodes(nodes_in_x_dimension*nodes_in_y_dimension*nodes_in_z_dimension),
       boundary_triangle_mesh(0){
        //generate tets from cubes
        //*note* this functions assumes the particles have been created with the ordering as in the variable index below
		
        int count=0;
        for(int k=0;k<nodes_in_z_dimension-1;k++)
            for(int j=0;j<nodes_in_y_dimension-1;j++)
                for(int i=0;i<nodes_in_x_dimension-1;i++){
                    //*note* we assume the particles have been ordered as follows:
                    int index=i+nodes_in_x_dimension*j+nodes_in_x_dimension*nodes_in_y_dimension*k;
                    //m->nodes_in_x_dimension, n->nodes_in_y_dimension, l->nodes_in_z_dimension                    
                    //0->i,j,k ->index	                            //		  4----6              
                    //1->i+1,j,k ->index+1						    //		 /|   /|               
                    //2->i+1,j+1,k ->index+m					    //		5----7 | 
                    //3->i,j+1,k ->index+m+1						//		| 0--|-2
                    //4->i,j,k+1 ->index+m*n						//		|/   |/
                    //5->i+1,j,k+1 ->index+1+m*n					//		1----3
                    //6->i+1,j+1,k+1 ->index+m+m*n
                    //7->i,j+1,k+1 ->index+m+1+m*n
					
                    //Even tets:
                    //(0,1,2,4);
                    //(1,3,2,7);
                    //(1,4,5,7);
                    //(2,4,7,6);
                    //(1,2,4,7);
                    if((i+j+k)%2 == 0){
                        mesh(count++)=VECTOR_4D<int>(index, index+1,index+nodes_in_x_dimension,index+nodes_in_x_dimension*nodes_in_y_dimension);
                        mesh(count++)=VECTOR_4D<int> (index+1, index+nodes_in_x_dimension+1,index+nodes_in_x_dimension,index+1+nodes_in_x_dimension+nodes_in_x_dimension*nodes_in_y_dimension);
                        mesh(count++)=VECTOR_4D<int> (index+1, index+nodes_in_x_dimension*nodes_in_y_dimension,index+1+nodes_in_x_dimension*nodes_in_y_dimension,index+1+nodes_in_x_dimension+nodes_in_x_dimension*nodes_in_y_dimension);
                        mesh(count++)=VECTOR_4D<int> (index+nodes_in_x_dimension, index+nodes_in_x_dimension*nodes_in_y_dimension,index+1+nodes_in_x_dimension+nodes_in_x_dimension*nodes_in_y_dimension,index+nodes_in_x_dimension+nodes_in_x_dimension*nodes_in_y_dimension);
                        mesh(count++)=VECTOR_4D<int> (index+1, index+nodes_in_x_dimension,index+nodes_in_x_dimension*nodes_in_y_dimension,index+1+nodes_in_x_dimension+nodes_in_x_dimension*nodes_in_y_dimension);}
					
                    //Odd tets:
                    //(5,4,6,0);
                    //(5,6,7,3);
                    //(0,1,3,5);
                    //(0,3,2,6);
                    //(0,6,5,3);
                    else{
                        mesh(count++)=VECTOR_4D<int>(index+1+nodes_in_x_dimension*nodes_in_y_dimension,index+nodes_in_x_dimension*nodes_in_y_dimension,index+nodes_in_x_dimension+nodes_in_x_dimension*nodes_in_y_dimension,index);
                        mesh(count++)=VECTOR_4D<int>(index+1+nodes_in_x_dimension*nodes_in_y_dimension,index+nodes_in_x_dimension+nodes_in_x_dimension*nodes_in_y_dimension,index+1+nodes_in_x_dimension+nodes_in_x_dimension*nodes_in_y_dimension,index+1+nodes_in_x_dimension);
                        mesh(count++)=VECTOR_4D<int>(index,index+1,index+1+nodes_in_x_dimension,index+1+nodes_in_x_dimension*nodes_in_y_dimension);
                        mesh(count++)=VECTOR_4D<int>(index,index+nodes_in_x_dimension+1,index+nodes_in_x_dimension,index+nodes_in_x_dimension*(nodes_in_y_dimension+1));
                        mesh(count++)=VECTOR_4D<int>(index,index+nodes_in_x_dimension*(nodes_in_y_dimension+1),index+1+nodes_in_x_dimension*nodes_in_y_dimension,index+1+nodes_in_x_dimension);}}
    }
	
    TETRAHEDRON_MESH(const int num_tets_input,const int number_nodes_input):mesh(num_tets_input),number_tets(num_tets_input),number_nodes(number_nodes_input),boundary_triangle_mesh(0)
    {}
	
    void Initialize_Oriented_Boundary_Triangles(){
        if(boundary_triangle_mesh) delete boundary_triangle_mesh;
		
        VECTOR<VECTOR_4D<VECTOR_3D<int> > > triangles_per_tet(number_tets);
		
        for(int t=0;t<number_tets;t++){
            triangles_per_tet(t)(0)=VECTOR_3D<int>(mesh(t)(0),mesh(t)(2),mesh(t)(1));
            triangles_per_tet(t)(1)=VECTOR_3D<int>(mesh(t)(0),mesh(t)(3),mesh(t)(2));
            triangles_per_tet(t)(2)=VECTOR_3D<int>(mesh(t)(0),mesh(t)(1),mesh(t)(3));
            triangles_per_tet(t)(3)=VECTOR_3D<int>(mesh(t)(1),mesh(t)(2),mesh(t)(3));}
		
        std::map<std::string,int> map_of_triangles;
        std::map<std::string,VECTOR_3D<int> > unsorted_triangle_indices;
		
        for(int t=0;t<number_tets;t++){
            for(int k=0;k<4;k++){
                VECTOR_3D<int> u=triangles_per_tet(t)(k).Sorted();
                char ss[30];
                sprintf(ss,"%d,%d,%d",u(0),u(1),u(2));
                std::string s=ss;
                if (map_of_triangles.find(s)!=map_of_triangles.end())
                    map_of_triangles[s]=map_of_triangles[s]+1;
                else{
                    unsorted_triangle_indices[s]=triangles_per_tet(t)(k);
                    map_of_triangles[s]=1;}}}

        int count=0;
        for(std::map<std::string,int>::iterator it1=map_of_triangles.begin();it1!=map_of_triangles.end();it1++){
            if(it1->second==1) count++;}
		
        boundary_triangle_mesh=new VECTOR<VECTOR_3D<int> >(count);
		
        count=0;
        for(std::map<std::string,int>::iterator it2=map_of_triangles.begin();it2!=map_of_triangles.end();it2++){
            if(it2->second==1) (*boundary_triangle_mesh)(count++)=unsorted_triangle_indices[it2->first];}
    }
	
    int Number_Of_Tetrahedra(){return number_tets;}
	
    VECTOR_4D<int>& Nodes_Of_Element(const int element){return mesh(element);}
	
    void Write_DAT_File(std::string mesh_file){
        FILE* fpointer;
        fpointer=fopen(mesh_file.c_str(),"w");
        for(int t=0;t<number_tets;t++){
            VECTOR_4D<int> nodes=Nodes_Of_Element(t);
            fprintf(fpointer,"%i\n",nodes(0));
            fprintf(fpointer,"%i\n",nodes(1));
            fprintf(fpointer,"%i\n",nodes(2));
            fprintf(fpointer,"%i\n",nodes(3));}		
        fclose(fpointer);
    }
};
	
}
#endif
