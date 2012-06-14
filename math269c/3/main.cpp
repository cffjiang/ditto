#include <iostream>
#include <math.h>
#include "FEM.h"
#include <fstream>
#include <string>
#include <vector>

#define PI 3.14159265

void Initialize_Square_Mesh(int number_nodes_per_dimension,VECTOR<VECTOR_2D<double> >& nodes,VECTOR<VECTOR<int>*>& elements){
	double dx=(double)1/(double)(number_nodes_per_dimension-1);
	for(int i=0;i<number_nodes_per_dimension;i++){
		for(int j=0;j<number_nodes_per_dimension;j++){
			nodes(i+j*number_nodes_per_dimension)=VECTOR_2D<double>(dx*(double)i,dx*(double)j);}}
	
	for(int e=0;e<elements.Size();e++) elements(e)=new VECTOR<int>(3);
	
	int count=0;
	
	for(int ei=0;ei<number_nodes_per_dimension-1;ei++){
		for(int ej=0;ej<number_nodes_per_dimension-1;ej++){
			int square=ei+ej*(number_nodes_per_dimension-1);
			count++;
			(*(elements(2*square)))(0)=ei+ej*number_nodes_per_dimension;
			(*(elements(2*square)))(1)=ei+1+(ej+1)*number_nodes_per_dimension;
			(*(elements(2*square)))(2)=ei+(ej+1)*number_nodes_per_dimension;
			count++;
			(*(elements(2*square+1)))(0)=ei+ej*number_nodes_per_dimension;
			(*(elements(2*square+1)))(1)=ei+1+ej*number_nodes_per_dimension;
			(*(elements(2*square+1)))(2)=ei+1+(ej+1)*number_nodes_per_dimension;}}
	assert(count==elements.Size());
}

double Source(const VECTOR_2D<double>& X){
	//return (double)8*PI*PI*sin((double)2*PI*X.x_copy())*cos((double)2*PI*X.y_copy());
	return -(double)4;
}

double Boundary_Flux(const VECTOR_2D<double>& X,const VECTOR_2D<double>& N){
	
	//VECTOR_2D<double> grad_u((double)2*PI*cos((double)2*PI*X.x_copy())*(double)2*PI*cos((double)2*PI*X.y_copy()),
	//						 -(double)2*PI*sin((double)2*PI*X.x_copy())*(double)2*PI*sin((double)2*PI*X.y_copy()));
	//return VECTOR_2D<double>::Dot_Product(grad_u,N);

	VECTOR_2D<double> grad_u((double)2*X.x_copy(),(double)2*X.y_copy());
	return VECTOR_2D<double>::Dot_Product(grad_u,N);
}

void Initialize_RHS(FEM_Poisson<double>& femp){
	VECTOR<double>& rhs=femp.RHS();
	rhs.Set_To_Zero();
	for(int e=0;e<femp.Number_Of_Elements();e++){
		VECTOR<int>& nodes_in_element=femp.Nodes_Of_Element(e);
		VECTOR_2D<double> cell_center=femp.X_At_Cell_Center(e);
		VECTOR_2D<double> x1=femp.X_At_Node(nodes_in_element(0));
		VECTOR_2D<double> x2=femp.X_At_Node(nodes_in_element(1));
		VECTOR_2D<double> x3=femp.X_At_Node(nodes_in_element(2));
		double area=VECTOR_2D<double>::Signed_Triangle_Area(x2-x1,x3-x1);
		assert(area>0);
		double one_third=(double)1/(double)3;
		for(int i=0;i<3;i++)
			rhs(nodes_in_element(i))+=one_third*area*Source(cell_center);}
	
	for(int s=0;s<femp.Mesh().Number_Boundary_Segments();s++){
		VECTOR_2D<int> nodes_in_segment=femp.Mesh().Nodes_Of_Boundary_Segment(s);
		VECTOR_2D<double> X_mid=(double).5*(femp.X_At_Node(nodes_in_segment.x_copy())+femp.X_At_Node(nodes_in_segment.y_copy()));
		VECTOR_2D<double> T=femp.X_At_Node(nodes_in_segment.x_copy())-femp.X_At_Node(nodes_in_segment.y_copy());
		VECTOR_2D<double> N=T.Right_Handed_Perp_Vector();
		rhs(nodes_in_segment.x_copy())+=(double).5*Boundary_Flux(X_mid,N);
		rhs(nodes_in_segment.y_copy())+=(double).5*Boundary_Flux(X_mid,N);}
	
}

bool Parse_CSV_Line(std::ifstream& in,std::vector<std::string>& entries){
	char line_char[256];
	bool fail=!in.getline(line_char,256);
	if(fail) return fail;
	
	std::string line(line_char);	
	std::string delimiter(",");	
	std::string::size_type lastPos=line.find_first_not_of(delimiter, 0);
	// Find first "non-delimiter".
	std::string::size_type pos=line.find_first_of(delimiter, lastPos);
	
	while (std::string::npos != pos || std::string::npos != lastPos)
	{
		// Found a token, add it to the vector.
		entries.push_back(line.substr(lastPos, pos - lastPos));
		// Skip delimiter.  Note the "not_of"
		lastPos = line.find_first_not_of(delimiter, pos);
		// Find next "non-delimiter"
		pos = line.find_first_of(delimiter, lastPos);
	}
	return true;
}

void Load_Mesh_From_File(VECTOR<VECTOR<int>*>*& elements,VECTOR<VECTOR_2D<double > >*& nodes){
	std::ifstream in("mesh_data.dat");
    if (in.fail())  { std::cout << "File not found" <<std::endl;}
	
	std::vector<std::string> entries;
	bool worked=Parse_CSV_Line(in,entries);
	assert(worked);
	int number_nodes=atoi(entries[0].c_str());
	nodes=new VECTOR<VECTOR_2D<double> >(number_nodes);
	int number_tris=atoi(entries[1].c_str());
	elements=new VECTOR<VECTOR<int>*>(number_tris);
	for(int t=0;t<number_tris;t++) (*elements)(t)=new VECTOR<int>(3);
	
	in.close();
	
	std::ifstream mesh_in("mesh_with_holes.dat");
    if (mesh_in.fail())  { std::cout << "File not found" <<std::endl;}
	
	bool no_problem=true;
	int count=0;
	while(mesh_in.good() )
    {
		std::vector<std::string> entries;
		no_problem=Parse_CSV_Line(mesh_in,entries);
		assert(no_problem);
		if(no_problem && entries.size()==3){
			int i=atoi(entries[0].c_str())-1;
			int j=atoi(entries[1].c_str())-1;
			int k=atoi(entries[2].c_str())-1;
			(*((*elements)(count)))(0)=i;
			(*((*elements)(count)))(1)=j;
			(*((*elements)(count)))(2)=k;
			assert(count<number_tris);
			count++;}
	}
		mesh_in.close();
	
	std::ifstream nodes_in("nodes.dat");
    if (nodes_in.fail())  { std::cout << "File not found" <<std::endl;}
	
	no_problem=true;
	count=0;
	while(nodes_in.good() )
    {
		std::vector<std::string> entries;
		no_problem=Parse_CSV_Line(nodes_in,entries);
		assert(no_problem);
		if(no_problem && entries.size()==2){
			double x=(double)atof(entries[0].c_str());
			double y=(double)atof(entries[1].c_str());
			(*nodes)(count).x()=x;
			(*nodes)(count).y()=y;
			count++;}
	}
	nodes_in.close();
}

int main (int argc, char * const argv[]) {
	
	//int number_nodes_per_dimension=25;
	//int number_nodes=number_nodes_per_dimension*number_nodes_per_dimension;
	//int number_elements=2*(number_nodes_per_dimension-1)*(number_nodes_per_dimension-1);
	//int max_it=number_nodes_per_dimension*number_nodes_per_dimension;
    
	//VECTOR<VECTOR_2D<double> > nodes(number_nodes);
	//VECTOR<VECTOR<int>*> elements(number_elements);

	//Initialize_Square_Mesh(number_nodes_per_dimension,nodes,elements);
	
	VECTOR<VECTOR_2D<double> >* nodes;
	VECTOR<VECTOR<int>*>* elements;
	
	Load_Mesh_From_File(elements,nodes);
	FEM_Poisson<double> femp(*nodes,*elements);

	int max_it=nodes->Size();
	
	Conjugate_Gradient<double> cg(femp.System_Matrix(),femp.Solution(),femp.RHS(),max_it);
	
	Initialize_RHS(femp);
	
	cg.Solve(true);
	femp.Write_Dat_Files();
	
	for(int e=0;e<elements->Size();e++) delete (*elements)(e);
	delete elements;	
	delete nodes;
	
    return 0;
}
