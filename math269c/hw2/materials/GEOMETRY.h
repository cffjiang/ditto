
class TRIANGLE_MESH{
	VECTOR<VECTOR<int>*>& mesh;
	VECTOR<VECTOR_2D<int> >* oriented_boundary_segments;
	int number_nodes;
	int number_tris;
	int number_boundary_segments;
public:
	TRIANGLE_MESH(VECTOR<VECTOR<int>*>& mesh_input,const int number_nodes_input):mesh(mesh_input),number_tris(mesh_input.Size()),number_nodes(number_nodes_input),
	number_boundary_segments(0),oriented_boundary_segments(0)
	{}
	
	int Number_Boundary_Segments(){return number_boundary_segments;}
	
	VECTOR_2D<int> Nodes_Of_Boundary_Segment(const int s){
		assert(oriented_boundary_segments);
		assert(number_boundary_segments==oriented_boundary_segments->Size());
		assert(s<number_boundary_segments && s>=0);
		return (*oriented_boundary_segments)(s);}
	
	void Intialize_Oriented_Boundary_Segments(){
		if(oriented_boundary_segments) delete oriented_boundary_segments;
		
		//assumes triangle are right hand oriented
		SPARSE_MATRIX<int> edges(number_nodes,number_nodes);
		for(int t=0;t<number_tris;t++){
			VECTOR<int>& nodes=Nodes_Of_Element(t);
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
				int j=ri.Index(j_sp);
				SPARSE_ROW<int>& rj=edges.Row(j);
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
	
	VECTOR<int>& Nodes_Of_Element(const int element){return (*(mesh(element)));}
	
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