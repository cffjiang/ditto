#ifndef _grids_
#define _grids_

#include <cassert>
#include <fstream>
#include "ALGEBRA.h"

using namespace ALGEBRA;

//All grids are node based. That is, dofs are assumed to live at the nodes of the grid.

namespace GRIDS{
template<class T>
class GRID_2D{
    int number_cells_per_x_dimension,number_cells_per_y_dimension;
    int number_nodes_per_x_dimension,number_nodes_per_y_dimension;
    T dx;//will always have the same dx in each dimension
    T xmin;
    T xmax;
    T ymin;
    T ymax;
    T epsilon;
    bool periodic;
	
public:
    GRID_2D(const int m_input,const T dx_input,const T xmin_input,const T ymin_input):
	number_cells_per_x_dimension(m_input),
       number_cells_per_y_dimension(m_input),
       dx(dx_input),xmin(xmin_input),ymin(ymin_input){
        xmax=xmin+(T)number_cells_per_x_dimension*dx;
        ymax=ymin+(T)number_cells_per_y_dimension*dx;
        Set_Epsilon((T).0000001);
        Set_Periodic(true);//this will initialize the number nodes per dimension
    }
	
    GRID_2D(const int m_input,const int n_input,const T dx_input,const T xmin_input,const T ymin_input):
	number_cells_per_x_dimension(m_input),
       number_cells_per_y_dimension(n_input),
       dx(dx_input),xmin(xmin_input),ymin(ymin_input){
        xmax=xmin+(T)number_cells_per_x_dimension*dx;
        ymax=ymin+(T)number_cells_per_y_dimension*dx;
        Set_Epsilon((T).0000001);
        Set_Periodic(false);//this will initialize the number nodes per dimension
    }
	
    void Resize(const int m_input,const int n_input,const T dx_input,const T xmin_input,const T ymin_input){
        number_cells_per_x_dimension=m_input;
        number_cells_per_y_dimension=n_input;
        dx=dx_input;
        xmin=xmin_input;
        ymin=ymin_input;
        Set_Periodic(periodic);
    }
	
    T Dx(){return dx;}
	
    T X_Max(){return xmax;}
    T X_Min(){return xmin;}
    T Y_Max(){return ymax;}
    T Y_Min(){return ymin;}
	
    int M(){return number_cells_per_x_dimension;}
    int N(){return number_cells_per_y_dimension;}
	
    int Number_Nodes_Per_X_Dimension(){return number_nodes_per_x_dimension;}
    int Number_Nodes_Per_Y_Dimension(){return number_nodes_per_y_dimension;}
	
    void Set_Periodic(const bool input){
        periodic=input;
        if(periodic){
            assert(number_cells_per_x_dimension==number_cells_per_y_dimension);
            number_nodes_per_x_dimension=number_cells_per_x_dimension;
            number_nodes_per_y_dimension=number_cells_per_y_dimension;} 
        else{ 
            number_nodes_per_x_dimension=number_cells_per_x_dimension+1;
            number_nodes_per_y_dimension=number_cells_per_y_dimension+1;}
    }
		
    void Set_Epsilon(T epsilon_input){epsilon=epsilon_input;}
		
    T X(const INDEX_2D& index){if(periodic) return (T)index.i_Periodic()*dx+xmin;else return (T)index.I()*dx+xmin;}
    T Y(const INDEX_2D& index){if(periodic) return (T)index.j_Periodic()*dx+ymin;else return (T)index.J()*dx+ymin;}
	
    void Cell_Containing_Point(T x,T y,T& lambda_x,T& lambda_y,INDEX_2D& index){
        /*
          This returns the cell in a grid that covers all of the plane. I.e. the lower left cell in the grid is (0,0), but you could get
          (i,j) with i and j being any integers, negative, larger than number_cells_per_dimension etc. The user can decide how to handle cells not on the grid
        */
		
        int i_min,j_min;
        T normalized=(x-xmin)/dx;
        if(normalized<0){
            i_min=(int)normalized-1;}
        else{
            i_min=(int)normalized;}
        normalized=(y-ymin)/dx;
        if(normalized<0){
            j_min=(int)normalized-1;}
        else{
            j_min=(int)normalized;}
		
        index.I()=i_min;
        index.J()=j_min;
        index.M()=number_nodes_per_x_dimension;
        index.N()=number_nodes_per_y_dimension;
        lambda_x=(x-(dx*(T)i_min+xmin))/dx;
        lambda_y=(y-(dx*(T)j_min+ymin))/dx;
    }
	
    void Clamp_To_Grid(T& x,T& y){
        if(x<xmin+epsilon) 
            x=xmin+epsilon;
        else if(x>xmax-epsilon) 
            x=xmax-epsilon;
        if(y<ymin+epsilon) 
            y=ymin+epsilon;
        else if(y>ymax-epsilon) 
            y=ymax-epsilon;
    }
	
    T Interpolate(VECTOR<T>& field,T x, T y){
        if(periodic){
            T lambda_x,lambda_y;
            INDEX_2D index,ip1,jp1,ip1jp1;
            Cell_Containing_Point(x,y,lambda_x,lambda_y,index);
            index.M()=number_nodes_per_x_dimension;
            index.N()=number_nodes_per_y_dimension;
            ip1.I()=index.I()+1;ip1.J()=index.J();ip1.M()=index.M();ip1.N()=index.N();
            jp1.I()=index.I();jp1.J()=index.J()+1;jp1.M()=index.M();jp1.N()=index.N();
            ip1jp1.I()=index.I()+1;ip1jp1.J()=index.J()+1;ip1jp1.M()=index.M();ip1jp1.N()=index.N();
            T w1=((T)1-lambda_x)*((T)1-lambda_y);
            T w2=lambda_x*((T)1-lambda_y);
            T w3=((T)1-lambda_x)*lambda_y;
            T w4=lambda_x*lambda_y;
            return w1*field(index.Index_Periodic())+w2*field(ip1.Index_Periodic())+w3*field(jp1.Index_Periodic())+w4*field(ip1jp1.Index_Periodic());}
        else{
            T lambda_x,lambda_y;
            INDEX_2D index,ip1,jp1,ip1jp1;
			
            //clamp to grid
            if(x<xmin+epsilon) 
                x=xmin+epsilon;
            else if(x>xmax-epsilon) 
                x=xmax-epsilon;
            if(y<ymin+epsilon) 
                y=ymin+epsilon;
            else if(y>ymax-epsilon) 
                y=ymax-epsilon;
			
            Cell_Containing_Point(x,y,lambda_x,lambda_y,index);//associates the cell with its lower left node
            ip1.I()=index.I()+1;ip1.J()=index.J();ip1.M()=index.M();ip1.N()=index.N();
            jp1.I()=index.I();jp1.J()=index.J()+1;jp1.M()=index.M();jp1.N()=index.N();
            ip1jp1.I()=index.I()+1;ip1jp1.J()=index.J()+1;ip1jp1.M()=index.M();ip1jp1.N()=index.N();
            T w1=((T)1-lambda_x)*((T)1-lambda_y);
            T w2=lambda_x*((T)1-lambda_y);
            T w3=((T)1-lambda_x)*lambda_y;
            T w4=lambda_x*lambda_y;
            return w1*field(index.Index())+w2*field(ip1.Index())+w3*field(jp1.Index())+w4*field(ip1jp1.Index());}
    }
	
    void Write_Domain_DAT_File(const std::string& output_dir){
        FILE* fpointer;
        std::string filename0("/Grid_domain_info.dat");
        fpointer=fopen((output_dir+filename0).c_str(),"w");
        fprintf(fpointer,"%i\n",number_nodes_per_x_dimension);
        fprintf(fpointer,"%i\n",number_nodes_per_y_dimension);
        fprintf(fpointer,"%g\n",xmin);
        fprintf(fpointer,"%g\n",ymin);
        fprintf(fpointer,"%g\n",xmax);
        fprintf(fpointer,"%g\n",ymax);
        fclose(fpointer);
		
        std::string filename("/Grid_domain_node_locations.dat");
        fpointer=fopen((output_dir+filename).c_str(),"w");
        for(int i=0;i<number_nodes_per_x_dimension;i++){
            for(int j=0;j<number_nodes_per_y_dimension;j++){
                INDEX_2D index(i,j,number_nodes_per_x_dimension,number_nodes_per_y_dimension);
                fprintf(fpointer,"%g ",X(index));fprintf(fpointer,"%g ",Y(index));}
            fprintf(fpointer,"\n");}			
		
        fclose(fpointer);
    }
};
	
template<class T>
class GRID_3D{
    int number_cells_per_x_dimension,number_cells_per_y_dimension,number_cells_per_z_dimension;
    int number_nodes_per_x_dimension,number_nodes_per_y_dimension,number_nodes_per_z_dimension;
    T dx;//will always have the same dx in each dimension
    T xmin;
    T xmax;
    T ymin;
    T ymax;
    T zmin;
    T zmax;
    T epsilon;
public:
    typedef enum {X_MAC,Y_MAC,Z_MAC,NODE} GRID_TYPE;
private:
    GRID_TYPE type;
		
public:
    GRID_3D(const int m_input,const T dx_input,const T xmin_input,const T ymin_input,const T zmin_input):
	number_cells_per_x_dimension(m_input),
       number_cells_per_y_dimension(m_input),
       number_cells_per_z_dimension(m_input),
       number_nodes_per_x_dimension(m_input+1),
       number_nodes_per_y_dimension(m_input+1),
       number_nodes_per_z_dimension(m_input+1),
       dx(dx_input),xmin(xmin_input),ymin(ymin_input),zmin(zmin_input),type(NODE){
        xmax=xmin+(T)number_cells_per_x_dimension*dx;
        ymax=ymin+(T)number_cells_per_y_dimension*dx;
        zmax=zmin+(T)number_cells_per_z_dimension*dx;
        Set_Epsilon((T).0000001);
    }
	
    GRID_3D(const int m_input,const int n_input,const int mn_input,const T dx_input,const T xmin_input,const T ymin_input,const T zmin_input):
	number_cells_per_x_dimension(m_input),
       number_cells_per_y_dimension(n_input),
       number_cells_per_z_dimension(mn_input),
       number_nodes_per_x_dimension(m_input+1),
       number_nodes_per_y_dimension(n_input+1),
       number_nodes_per_z_dimension(mn_input+1),
       dx(dx_input),xmin(xmin_input),ymin(ymin_input),zmin(zmin_input),type(NODE){
        xmax=xmin+(T)number_cells_per_x_dimension*dx;
        ymax=ymin+(T)number_cells_per_y_dimension*dx;
        zmax=zmin+(T)number_cells_per_z_dimension*dx;
        Set_Epsilon((T).0000001);
    }
	
    GRID_TYPE Grid_Type(){return type;}
	
    void Set_Grid_Type(GRID_TYPE type_input){type=type_input;}
	
    void Get_Four_Surrounding_X_MAC_Indices(INDEX_3D& x_mac_index_0,INDEX_3D& x_mac_index_1,INDEX_3D& x_mac_index_2,INDEX_3D& x_mac_index_3,const INDEX_3D& mac_index){
        /*This is used to interpolate staggered x MAC velocities to a y or z face
         */
		
        assert(type==Y_MAC || type==Z_MAC);
        if(type==Y_MAC){
            assert(mac_index.J()>0 && mac_index.J()<mac_index.N()-1);
            x_mac_index_0=INDEX_3D(mac_index.I(),mac_index.J()-1,mac_index.K(),mac_index.M()+1,mac_index.N()-1,mac_index.MN());
            x_mac_index_1=INDEX_3D(mac_index.I()+1,mac_index.J()-1,mac_index.K(),mac_index.M()+1,mac_index.N()-1,mac_index.MN());
            x_mac_index_2=INDEX_3D(mac_index.I()+1,mac_index.J(),mac_index.K(),mac_index.M()+1,mac_index.N()-1,mac_index.MN());
            x_mac_index_3=INDEX_3D(mac_index.I(),mac_index.J(),mac_index.K(),mac_index.M()+1,mac_index.N()-1,mac_index.MN());}
        else if(type==Z_MAC){
            assert(mac_index.K()>0 && mac_index.K()<mac_index.MN()-1);
            x_mac_index_0=INDEX_3D(mac_index.I(),mac_index.J(),mac_index.K()-1,mac_index.M()+1,mac_index.N(),mac_index.MN()-1);
            x_mac_index_1=INDEX_3D(mac_index.I(),mac_index.J(),mac_index.K(),mac_index.M()+1,mac_index.N(),mac_index.MN()-1);
            x_mac_index_2=INDEX_3D(mac_index.I()+1,mac_index.J(),mac_index.K(),mac_index.M()+1,mac_index.N(),mac_index.MN()-1);
            x_mac_index_3=INDEX_3D(mac_index.I()+1,mac_index.J(),mac_index.K()-1,mac_index.M()+1,mac_index.N(),mac_index.MN()-1);}
    }
	
    void Get_Four_Surrounding_Y_MAC_Indices(INDEX_3D& y_mac_index_0,INDEX_3D& y_mac_index_1,INDEX_3D& y_mac_index_2,INDEX_3D& y_mac_index_3,const INDEX_3D& mac_index){
        /*This is used to interpolate staggered y MAC velocities to an x or z face
         */
		
        assert(type==X_MAC || type==Z_MAC);
        if(type==X_MAC){
            assert(mac_index.I()>0 && mac_index.I()<mac_index.M()-1);//this doesn't make sense if we are at the left or right boundary of the x mac grid
            y_mac_index_0=INDEX_3D(mac_index.I()-1,mac_index.J(),mac_index.K(),mac_index.M()-1,mac_index.N()+1,mac_index.MN());
            y_mac_index_1=INDEX_3D(mac_index.I(),mac_index.J(),mac_index.K(),mac_index.M()-1,mac_index.N()+1,mac_index.MN());
            y_mac_index_2=INDEX_3D(mac_index.I(),mac_index.J()+1,mac_index.K(),mac_index.M()-1,mac_index.N()+1,mac_index.MN());
            y_mac_index_3=INDEX_3D(mac_index.I()-1,mac_index.J()+1,mac_index.K(),mac_index.M()-1,mac_index.N()+1,mac_index.MN());}
        else if(type==Z_MAC){
            assert(mac_index.K()>0 && mac_index.K()<mac_index.MN()-1);//this doesn't make sense if we are at the top or bottom boundary of the z mac grid
            y_mac_index_0=INDEX_3D(mac_index.I(),mac_index.J(),mac_index.K()-1,mac_index.M(),mac_index.N()+1,mac_index.MN()-1);
            y_mac_index_1=INDEX_3D(mac_index.I(),mac_index.J()+1,mac_index.K()-1,mac_index.M(),mac_index.N()+1,mac_index.MN()-1);
            y_mac_index_2=INDEX_3D(mac_index.I(),mac_index.J()+1,mac_index.K(),mac_index.M(),mac_index.N()+1,mac_index.MN()-1);
            y_mac_index_3=INDEX_3D(mac_index.I(),mac_index.J(),mac_index.K(),mac_index.M(),mac_index.N()+1,mac_index.MN()-1);}
    }
	
    void Get_Four_Surrounding_Z_MAC_Indices(INDEX_3D& z_mac_index_0,INDEX_3D& z_mac_index_1,INDEX_3D& z_mac_index_2,INDEX_3D& z_mac_index_3,const INDEX_3D& mac_index){
        /*This is used to interpolate staggered z MAC velocities to an x or y face
         */
		
        assert(type==X_MAC || type==Y_MAC);
        if(type==X_MAC){
            assert(mac_index.I()>0 && mac_index.I()<mac_index.M()-1);//this doesn't make sense if we are at the left or right boundary of the x mac grid
            z_mac_index_0=INDEX_3D(mac_index.I()-1,mac_index.J(),mac_index.K(),mac_index.M()-1,mac_index.N(),mac_index.MN()+1);
            z_mac_index_1=INDEX_3D(mac_index.I()-1,mac_index.J(),mac_index.K()+1,mac_index.M()-1,mac_index.N(),mac_index.MN()+1);
            z_mac_index_2=INDEX_3D(mac_index.I(),mac_index.J(),mac_index.K()+1,mac_index.M()-1,mac_index.N(),mac_index.MN()+1);
            z_mac_index_3=INDEX_3D(mac_index.I(),mac_index.J(),mac_index.K(),mac_index.M()-1,mac_index.N(),mac_index.MN()+1);}
        else if(type==Y_MAC){
            assert(mac_index.J()>0 && mac_index.J()<mac_index.N()-1);
            z_mac_index_0=INDEX_3D(mac_index.I(),mac_index.J()-1,mac_index.K(),mac_index.M(),mac_index.N()-1,mac_index.MN()+1);
            z_mac_index_1=INDEX_3D(mac_index.I(),mac_index.J(),mac_index.K(),mac_index.M(),mac_index.N()-1,mac_index.MN()+1);
            z_mac_index_2=INDEX_3D(mac_index.I(),mac_index.J(),mac_index.K()+1,mac_index.M(),mac_index.N()-1,mac_index.MN()+1);
            z_mac_index_3=INDEX_3D(mac_index.I(),mac_index.J()-1,mac_index.K()+1,mac_index.M(),mac_index.N()-1,mac_index.MN()+1);}
    }
			
		
    T X_Max(){return xmax;}
    T X_Min(){return xmin;}
    T Y_Max(){return ymax;}
    T Y_Min(){return ymin;}
    T Z_Max(){return zmax;}
    T Z_Min(){return zmin;}
	
    int M(){return number_cells_per_x_dimension;}
    int N(){return number_cells_per_y_dimension;}
    int MN(){return number_cells_per_z_dimension;}
	
    int Number_Nodes_Per_X_Dimension(){return number_nodes_per_x_dimension;}
    int Number_Nodes_Per_Y_Dimension(){return number_nodes_per_y_dimension;}
    int Number_Nodes_Per_Z_Dimension(){return number_nodes_per_z_dimension;}
	
    void Set_Epsilon(T epsilon_input){epsilon=epsilon_input;}
	
    T X(const INDEX_3D& index){return (T)index.I()*dx+xmin;}
    T Y(const INDEX_3D& index){return (T)index.J()*dx+ymin;}
    T Z(const INDEX_3D& index){return (T)index.K()*dx+zmin;}
	
    void Cell_Containing_Point(T x,T y,T z,T& lambda_x,T& lambda_y,T& lambda_z,INDEX_3D& index){
        /*
          This returns the cell in a grid that covers all of space. I.e. the lower left cell in the grid is (0,0,0), but you could get
          (i,j,k) with i, j and k being any integers, negative, larger than number_cells_per_dimension etc. The user can decide how to handle cells not on the grid
        */
		
        int i_min,j_min,k_min;
        T normalized=(x-xmin)/dx;
        if(normalized<0){
            i_min=(int)normalized-1;}
        else{
            i_min=(int)normalized;}
        normalized=(y-ymin)/dx;
        if(normalized<0){
            j_min=(int)normalized-1;}
        else{
            j_min=(int)normalized;}
        normalized=(z-zmin)/dx;
        if(normalized<0){
            k_min=(int)normalized-1;}
        else{
            k_min=(int)normalized;}
		
        index.I()=i_min;
        index.J()=j_min;
        index.K()=k_min;
        index.M()=number_nodes_per_x_dimension;
        index.N()=number_nodes_per_y_dimension;
        index.MN()=number_nodes_per_z_dimension;
        lambda_x=(x-(dx*(T)i_min+xmin))/dx;
        lambda_y=(y-(dx*(T)j_min+ymin))/dx;
        lambda_z=(z-(dx*(T)k_min+zmin))/dx;
    }
	
    T Interpolate(VECTOR<T>& field,T x, T y, T z){
        T lambda_x,lambda_y,lambda_z;
        INDEX_3D index;
			
        //clamp to grid
        if(x<xmin+epsilon) 
            x=xmin+epsilon;
        else if(x>xmax-epsilon) 
            x=xmax-epsilon;
        if(y<ymin+epsilon) 
            y=ymin+epsilon;
        else if(y>ymax-epsilon) 
            y=ymax-epsilon;
        if(z<zmin+epsilon) 
            z=zmin+epsilon;
        else if(z>zmax-epsilon) 
            z=zmax-epsilon;
			
        Cell_Containing_Point(x,y,z,lambda_x,lambda_y,lambda_z,index);//associates the cell with its lower left node
        VECTOR_2D<T> x_edge_weights((T)1-lambda_x,lambda_x);
        VECTOR_2D<T> y_edge_weights((T)1-lambda_y,lambda_y);
        VECTOR_2D<T> z_edge_weights((T)1-lambda_z,lambda_z);
		
        T result=(T)0;
        for(int i=0;i<2;i++){
            for(int j=0;j<2;j++){
                for(int k=0;k<2;k++){
                    INDEX_3D node_index(index.I()+i,index.J()+j,index.K()+k,index.M(),index.N(),index.MN());
                    result+=x_edge_weights(i)*y_edge_weights(j)*z_edge_weights(k)*field(node_index.Index());}}}
		
        return result;
    }	
};

template<class T>
class GRID_2D_MAC_X{
    int number_cells_per_x_dimension,number_cells_per_y_dimension;//cells are between four nodes, dofs are at the nodes
    int number_nodes_per_x_dimension,number_nodes_per_y_dimension;
    T dx;
    T xmin;
    T xmax;
    T ymin;
    T ymax;
    T epsilon;
	
public:
	
    //Note that for both constructors, m and n are taken to be the number of cells in the MAC grid
	
    GRID_2D_MAC_X(const int m_input,const T dx_input,const T xmin_input,const T ymin_input):number_cells_per_x_dimension(m_input),number_cells_per_y_dimension(m_input-1),
                                                                                           dx(dx_input),xmin(xmin_input),ymin(ymin_input),number_nodes_per_x_dimension(m_input+1),number_nodes_per_y_dimension(m_input){
        xmax=xmin+(T)number_cells_per_x_dimension*dx;
        ymax=ymin+(T)number_cells_per_y_dimension*dx;
        Set_Epsilon((T).0000001);
    }
	
    GRID_2D_MAC_X(const int m_input,const int n_input,const T dx_input,const T xmin_input,const T ymin_input):
	number_cells_per_x_dimension(m_input),
       number_cells_per_y_dimension(n_input-1),
       dx(dx_input),xmin(xmin_input),ymin(ymin_input),
       number_nodes_per_x_dimension(m_input+1),
       number_nodes_per_y_dimension(n_input){
        xmax=xmin+(T)number_cells_per_x_dimension*dx;
        ymax=ymin+(T)number_cells_per_y_dimension*dx;
        Set_Epsilon((T).0000001);
    }
	
    void Get_Two_Adjacent_Pressure_Indices(INDEX_2D& left_pressure,INDEX_2D& right_pressure,const INDEX_2D& x_mac_index){
        assert(x_mac_index.M()==number_nodes_per_x_dimension);
        assert(x_mac_index.N()==number_nodes_per_y_dimension);
        left_pressure=INDEX_2D(x_mac_index.I()-1,x_mac_index.J(),x_mac_index.M()-1,x_mac_index.N());
        right_pressure=INDEX_2D(x_mac_index.I(),x_mac_index.J(),x_mac_index.M()-1,x_mac_index.N());
    }
	
    int Number_Nodes_Per_X_Dimension(){return number_nodes_per_x_dimension;}
	
    int Number_Nodes_Per_Y_Dimension(){return number_nodes_per_y_dimension;}
	
    void Set_Epsilon(T epsilon_input){epsilon=epsilon_input;}
	
    T X(const INDEX_2D& index){return (T)index.I()*dx+xmin;}
    T Y(const INDEX_2D& index){return (T)index.J()*dx+ymin;}
	
    T X_Max(){return xmax;}
    T X_Min(){return xmin;}
    T Y_Max(){return ymax;}
    T Y_Min(){return ymin;}
	
    void Get_Four_Surrounding_Y_MAC_Indices(INDEX_2D& index_ll,INDEX_2D& index_lr,INDEX_2D& index_ur,INDEX_2D& index_ul,INDEX_2D& x_mac_index){
        assert(x_mac_index.I()>0 && x_mac_index.I()<x_mac_index.M()-1);//this doesn't make sense if we are at the left or right boundary of the x mac grid
        index_ll=INDEX_2D(x_mac_index.I()-1,x_mac_index.J(),x_mac_index.M()-1,x_mac_index.N()+1);
        index_lr=INDEX_2D(x_mac_index.I(),x_mac_index.J(),x_mac_index.M()-1,x_mac_index.N()+1);
        index_ur=INDEX_2D(x_mac_index.I(),x_mac_index.J()+1,x_mac_index.M()-1,x_mac_index.N()+1);
        index_ul=INDEX_2D(x_mac_index.I()-1,x_mac_index.J()+1,x_mac_index.M()-1,x_mac_index.N()+1);
    }
	
    void Cell_Containing_Point(T x,T y,T& lambda_x,T& lambda_y,INDEX_2D& index){
        /*
          This returns the cell in a grid that covers all of the plane. I.e. the lower left cell in the grid is (0,0,number_cells_per_dimension), but you could get
          (i,j,number_cells_per_dimension) with i and j being any integers, negative, larger than number_cells_per_dimension etc. The user can decide how to handle cells not on the grid
        */
		
        int i_min,j_min;
        T normalized=(x-xmin)/dx;
        if(normalized<0){
            i_min=(int)normalized-1;}
        else{
            i_min=(int)normalized;}
        normalized=(y-ymin)/dx;
        if(normalized<0){
            j_min=(int)normalized-1;}
        else{
            j_min=(int)normalized;}
		
        index.I()=i_min;
        index.J()=j_min;
        index.M()=number_nodes_per_x_dimension;
        index.N()=number_nodes_per_y_dimension;
        lambda_x=(x-(dx*(T)i_min+xmin))/dx;
        lambda_y=(y-(dx*(T)j_min+ymin))/dx;
    }
	
    void Clamp_To_Grid(T& x,T& y){
        if(x<xmin+epsilon) 
            x=xmin+epsilon;
        else if(x>xmax-epsilon) 
            x=xmax-epsilon;
        if(y<ymin+epsilon) 
            y=ymin+epsilon;
        else if(y>ymax-epsilon) 
            y=ymax-epsilon;
    }
		
    T Interpolate(VECTOR<T>& field,T x, T y){
        /*
          Clamps to the boundaries of the grid, then uses the special indexing that would be seen taking (i,j,number_cells_per_dimension) on a U velocity MAC grid where (0,0,number_cells_per_dimension)
          is the index at the lower left U velocity, (i.e. at spatial location (0,dx/2) )
        */
		
        T lambda_x,lambda_y;
        INDEX_2D index,ip1,jp1,ip1jp1;
        //clamp to grid
        if(x<xmin+epsilon) 
            x=xmin+epsilon;
        else if(x>xmax-epsilon) 
            x=xmax-epsilon;
        if(y<ymin+epsilon) 
            y=ymin+epsilon;
        else if(y>ymax-epsilon) 
            y=ymax-epsilon;
        Cell_Containing_Point(x,y,lambda_x,lambda_y,index);//associates the cell with its lower left node
        assert(index.J()<number_cells_per_y_dimension && index.I()<number_cells_per_x_dimension);
        assert(index.J()>=0 && index.I()>=0);
        ip1.I()=index.I()+1;ip1.J()=index.J();ip1.M()=index.M();ip1.N()=index.N();
        jp1.I()=index.I();jp1.J()=index.J()+1;jp1.M()=index.M();jp1.N()=index.N();
        ip1jp1.I()=index.I()+1;ip1jp1.J()=index.J()+1;ip1jp1.M()=index.M();ip1jp1.N()=index.N();
        T w1=((T)1-lambda_x)*((T)1-lambda_y);
        T w2=lambda_x*((T)1-lambda_y);
        T w3=((T)1-lambda_x)*lambda_y;
        T w4=lambda_x*lambda_y;
        return w1*field(index.Index())+w2*field(ip1.Index())+w3*field(jp1.Index())+w4*field(ip1jp1.Index());
    }
	
    void Write_Domain_DAT_File(const std::string& output_dir){
        FILE* fpointer;
        std::string filename0("/MAC_X_domain_info.dat");
        fpointer=fopen((output_dir+filename0).c_str(),"w");
        fprintf(fpointer,"%i\n",number_nodes_per_x_dimension);
        fprintf(fpointer,"%i\n",number_nodes_per_y_dimension);
        fprintf(fpointer,"%g\n",xmin);
        fprintf(fpointer,"%g\n",ymin);
        fprintf(fpointer,"%g\n",xmax);
        fprintf(fpointer,"%g\n",ymax);
        fclose(fpointer);
		
        std::string filename("/MAC_X_domain_node_locations.dat");
        fpointer=fopen((output_dir+filename).c_str(),"w");
        for(int i=0;i<number_nodes_per_x_dimension;i++){
            for(int j=0;j<number_nodes_per_y_dimension;j++){
                INDEX_2D index(i,j,number_cells_per_x_dimension,number_cells_per_y_dimension);
                fprintf(fpointer,"%g ",X(index));fprintf(fpointer,"%g ",Y(index));}
            fprintf(fpointer,"\n");}			
		
        fclose(fpointer);
    }
	
};

template<class T>
class GRID_2D_MAC_Y{
    int number_cells_per_x_dimension,number_cells_per_y_dimension;//cells are between four nodes, dofs are at the nodes
    int number_nodes_per_x_dimension,number_nodes_per_y_dimension;
    T dx;
    T xmin;
    T xmax;
    T ymin;
    T ymax;
    T epsilon;
	
public:
	
    //Note that for both constructors, m and n are taken to be the number of cells in the MAC grid
	
    GRID_2D_MAC_Y(const int m_input,const T dx_input,const T xmin_input,const T ymin_input):
	number_cells_per_x_dimension(m_input-1),
       number_cells_per_y_dimension(m_input),
       dx(dx_input),xmin(xmin_input),ymin(ymin_input),
       number_nodes_per_x_dimension(m_input),
       number_nodes_per_y_dimension(m_input+1){
        xmax=xmin+(T)number_cells_per_x_dimension*dx;
        ymax=ymin+(T)number_cells_per_y_dimension*dx;
        Set_Epsilon((T).0000001);
    }
	
    GRID_2D_MAC_Y(const int m_input,const int n_input,const T dx_input,const T xmin_input,const T ymin_input):
	number_cells_per_x_dimension(m_input-1),
       number_cells_per_y_dimension(n_input),
       dx(dx_input),xmin(xmin_input),ymin(ymin_input),
       number_nodes_per_x_dimension(m_input),
       number_nodes_per_y_dimension(n_input+1){
        xmax=xmin+(T)number_cells_per_x_dimension*dx;
        ymax=ymin+(T)number_cells_per_y_dimension*dx;
        Set_Epsilon((T).0000001);
    }
	
    int Number_Nodes_Per_X_Dimension(){return number_nodes_per_x_dimension;}
	
    int Number_Nodes_Per_Y_Dimension(){return number_nodes_per_y_dimension;}

    void Set_Epsilon(T epsilon_input){epsilon=epsilon_input;}
	
    T X(const INDEX_2D& index){return (T)index.I()*dx+xmin;}
    T Y(const INDEX_2D& index){return (T)index.J()*dx+ymin;}
	
    T X_Max(){return xmax;}
    T X_Min(){return xmin;}
    T Y_Max(){return ymax;}
    T Y_Min(){return ymin;}
	
    void Get_Two_Adjacent_Pressure_Indices(INDEX_2D& below_pressure,INDEX_2D& above_pressure,const INDEX_2D& y_mac_index){
        assert(y_mac_index.M()==number_nodes_per_x_dimension);
        assert(y_mac_index.N()==number_nodes_per_y_dimension);
        below_pressure=INDEX_2D(y_mac_index.I(),y_mac_index.J()-1,y_mac_index.M(),y_mac_index.N()-1);
        above_pressure=INDEX_2D(y_mac_index.I(),y_mac_index.J(),y_mac_index.M(),y_mac_index.N()-1);
    }
	
    void Get_Four_Surrounding_X_MAC_Indices(INDEX_2D& index_ll,INDEX_2D& index_lr,INDEX_2D& index_ur,INDEX_2D& index_ul,INDEX_2D& y_mac_index){
        assert(y_mac_index.J()>0 && y_mac_index.J()<y_mac_index.N()-1);//this doesn't make sense if we are at the bottom or top boundary of the y mac grid
        index_ll=INDEX_2D(y_mac_index.I(),y_mac_index.J()-1,y_mac_index.M()+1,y_mac_index.N()-1);
        index_lr=INDEX_2D(y_mac_index.I()+1,y_mac_index.J()-1,y_mac_index.M()+1,y_mac_index.N()-1);
        index_ur=INDEX_2D(y_mac_index.I()+1,y_mac_index.J(),y_mac_index.M()+1,y_mac_index.N()-1);
        index_ul=INDEX_2D(y_mac_index.I(),y_mac_index.J(),y_mac_index.M()+1,y_mac_index.N()-1);
    }
	
    void Cell_Containing_Point(T x,T y,T& lambda_x,T& lambda_y,INDEX_2D& index){
        /*
          This returns the cell in a grid that covers all of the plane. I.e. the lower left cell in the grid is (0,0), but you could get
          (i,j) with i and j being any integers, negative, larger than number_cells_per_dimension etc. The user can decide how to handle cells not on the grid
        */
		
        int i_min,j_min;
        T normalized=(x-xmin)/dx;
        if(normalized<0){
            i_min=(int)normalized-1;}
        else{
            i_min=(int)normalized;}
        normalized=(y-ymin)/dx;
        if(normalized<0){
            j_min=(int)normalized-1;}
        else{
            j_min=(int)normalized;}
		
        index.I()=i_min;
        index.J()=j_min;
        index.M()=number_nodes_per_x_dimension;
        index.N()=number_nodes_per_y_dimension;
        lambda_x=(x-(dx*(T)i_min+xmin))/dx;
        lambda_y=(y-(dx*(T)j_min+ymin))/dx;
    }
	
    void Clamp_To_Grid(T& x,T& y){
        if(x<xmin+epsilon) 
            x=xmin+epsilon;
        else if(x>xmax-epsilon) 
            x=xmax-epsilon;
        if(y<ymin+epsilon) 
            y=ymin+epsilon;
        else if(y>ymax-epsilon) 
            y=ymax-epsilon;
    }
	
    T Interpolate(VECTOR<T>& field,T x, T y){
        /*
          Clamps to the boundaries of the grid, then uses the special indexing that would be seen taking (i,j) on a V velocity MAC grid where (0,0)
          is the index at the lower left V velocity, (i.e. at spatial location (dx/2,0) )
        */
		
        T lambda_x,lambda_y;
        INDEX_2D index,ip1,jp1,ip1jp1;
        //clamp to grid
        if(x<xmin+epsilon) 
            x=xmin+epsilon;
        else if(x>xmax-epsilon) 
            x=xmax-epsilon;
        if(y<ymin+epsilon) 
            y=ymin+epsilon;
        else if(y>ymax-epsilon) 
            y=ymax-epsilon;
		
        Cell_Containing_Point(x,y,lambda_x,lambda_y,index);//associates the cell with its lower left node
        assert(index.J()<number_cells_per_y_dimension && index.I()<number_cells_per_x_dimension);
        assert(index.J()>=0 && index.I()>=0);
        ip1.I()=index.I()+1;ip1.J()=index.J();ip1.M()=index.M();ip1.N()=index.N();
        jp1.I()=index.I();jp1.J()=index.J()+1;jp1.M()=index.M();jp1.N()=index.N();
        ip1jp1.I()=index.I()+1;ip1jp1.J()=index.J()+1;ip1jp1.M()=index.M();ip1jp1.N()=index.N();
        T w1=((T)1-lambda_x)*((T)1-lambda_y);
        T w2=lambda_x*((T)1-lambda_y);
        T w3=((T)1-lambda_x)*lambda_y;
        T w4=lambda_x*lambda_y;
        return w1*field(index.Index())+w2*field(ip1.Index())+w3*field(jp1.Index())+w4*field(ip1jp1.Index());
    }
	
    void Write_Domain_DAT_File(const std::string& output_dir){
        FILE* fpointer;
        std::string filename0("/MAC_Y_domain_info.dat");
        fpointer=fopen((output_dir+filename0).c_str(),"w");
        fprintf(fpointer,"%i\n",number_nodes_per_x_dimension);
        fprintf(fpointer,"%i\n",number_nodes_per_y_dimension);
        fprintf(fpointer,"%g\n",xmin);
        fprintf(fpointer,"%g\n",ymin);
        fprintf(fpointer,"%g\n",xmax);
        fprintf(fpointer,"%g\n",ymax);
        fclose(fpointer);
		
        std::string filename("/MAC_Y_domain_node_locations.dat");
        fpointer=fopen((output_dir+filename).c_str(),"w");
        for(int i=0;i<number_nodes_per_x_dimension;i++){
            for(int j=0;j<number_nodes_per_y_dimension;j++){
                INDEX_2D index(i,j,number_nodes_per_x_dimension,number_nodes_per_y_dimension);
                fprintf(fpointer,"%g ",X(index));fprintf(fpointer,"%g ",Y(index));}
            fprintf(fpointer,"\n");}			
		
        fclose(fpointer);
    }
	
};
}
#endif	
