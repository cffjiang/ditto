#include <cassert>

template<class T>
class GRID_2D{
	int m;
	T dx;
	T xmin;
	T xmax;
	T ymin;
	T ymax;
	
public:
	GRID_2D(const int m_input,const T dx_input,const T xmin_input,const T ymin_input):m(m_input),dx(dx_input),xmin(xmin_input),ymin(ymin_input){
		xmax=xmin+(T)m*dx;
		ymax=ymin+(T)m*dx;
	}
		
	T X(INDEX_2D index){return (T)index.i_Periodic()*dx+xmin;}
	T Y(INDEX_2D index){return (T)index.j_Periodic()*dx+ymin;}
	
	void Cell_Containing_Point(T x,T y,T& lambda_x,T& lambda_y,INDEX_2D& index){
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
		
		index.i=i_min;
		index.j=j_min;
		index.m=m;
		lambda_x=(x-(dx*(double)i_min+xmin))/dx;
		lambda_y=(y-(dx*(double)j_min+ymin))/dx;
	}
	
	T Interpolate(VECTOR<T>& field,T x, T y){
		T lambda_x,lambda_y;
		INDEX_2D index,ip1,jp1,ip1jp1;
		Cell_Containing_Point(x,y,lambda_x,lambda_y,index);
		ip1.i=index.i+1;ip1.j=index.j;ip1.m=index.m;
		jp1.i=index.i;jp1.j=index.j+1;jp1.m=index.m;
		ip1jp1.i=index.i+1;ip1jp1.j=index.j+1;ip1jp1.m=index.m;
		T w1=((T)1-lambda_x)*((T)1-lambda_y);
		T w2=lambda_x*((T)1-lambda_y);
		T w3=((T)1-lambda_x)*lambda_y;
		T w4=lambda_x*lambda_y;
		return w1*field(index)+w2*field(ip1)+w3*field(jp1)+w4*field(ip1jp1);
	}
};
	
