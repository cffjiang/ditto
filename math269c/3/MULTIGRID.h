#include "algebra.h"
#include <stdio.h>

using namespace std;

template<class T>
class MULTIGRID_POISSON{
	int levels,n_fine,dimension,num_v_cycles;
	T dx_fine;
	VECTOR<T>** v_hierarchy;
	VECTOR<T>** temp_hierarchy;
	VECTOR<T>** rhs_hierarchy;
	SPARSE_MATRIX<T>** A_hierarchy;
	SPARSE_MATRIX<T>** R_hierarchy;
	SPARSE_MATRIX<T>** P_hierarchy;
	
public:
	MULTIGRID_POISSON(const int levels_input,const int n_fine_input,const int dimension_input,const int num_v_cycles_input):levels(levels_input),n_fine(n_fine_input),
	dimension(dimension_input),num_v_cycles(num_v_cycles_input)
	{
		assert(dimension==1 || dimension==2);
		v_hierarchy=new VECTOR<T>*[levels];
		temp_hierarchy=new VECTOR<T>*[levels];
		rhs_hierarchy=new VECTOR<T>*[levels];
		A_hierarchy=new SPARSE_MATRIX<T>*[levels];
		R_hierarchy=new SPARSE_MATRIX<T>*[levels-1];
		P_hierarchy=new SPARSE_MATRIX<T>*[levels-1];

		dx_fine=(T)1/(T)n_fine;
		
		//build hierarchy
		int factor=1;
		for(int l=0;l<levels;l++){
			int n=n_fine/factor;
			std::cout << "n = " << n << " at level " << l << std::endl;
			T dx=1./(T)n;
			int number_dofs;
			if(dimension==1) number_dofs=n;
			else number_dofs=n*n;
			A_hierarchy[l]=new SPARSE_MATRIX<T>(number_dofs,number_dofs);
			v_hierarchy[l]=new VECTOR<T>(number_dofs);
			temp_hierarchy[l]=new VECTOR<T>(number_dofs);
			rhs_hierarchy[l]=new VECTOR<T>(number_dofs);
			for(int i=0;i<number_dofs;i++) (*(rhs_hierarchy[l]))(i)=(T)0;
			Build_Periodic_Poisson_Matrix(*(A_hierarchy[l]),n,dx);
			factor=2*factor;}
		
		//build transfer operators
		factor=1;
		for(int l=0;l<levels-1;l++){
			int n=n_fine/factor;
			int n_over_2=n/2;
			int number_dofs;int number_coarse_dofs;
			if(dimension==1){number_dofs=n;number_coarse_dofs=n_over_2;}
			else{number_dofs=n*n;number_coarse_dofs=n_over_2*n_over_2;}
			R_hierarchy[l]=new SPARSE_MATRIX<T>(number_coarse_dofs,number_dofs);
			Build_Restriction((*(R_hierarchy[l])),n);
			P_hierarchy[l]=new SPARSE_MATRIX<T>(number_dofs,number_coarse_dofs);
			R_hierarchy[l]->Transpose((*(P_hierarchy[l])));
			P_hierarchy[l]->Scale_Rows((T)2*(T)dimension);
			factor=2*factor;}
	}
	
	SPARSE_MATRIX<T>& Prolongation(const int level){assert(level<levels && level>=0);return (*(P_hierarchy[level]));}
	SPARSE_MATRIX<T>& Restriction(const int level){assert(level<levels && level>=0);return (*(R_hierarchy[level]));}
	SPARSE_MATRIX<T>& System_Matrix(const int level){assert(level<levels && level>=0);return (*(A_hierarchy[level]));}
	
	void Write_DAT_File(){
		VECTOR<T>& v=(*(v_hierarchy[0]));
	
		FILE* fpointer;
		fpointer=fopen("solution.dat","w");
		if(dimension==1){
			for(int i=0;i<n_fine;i++)
				fprintf(fpointer,"%g\n",v(i));}
		else if(dimension==2){
			for(int i=0;i<n_fine;i++){
				for(int j=0;j<n_fine;j++){
					INDEX_2D index={i,j,n_fine};
					fprintf(fpointer,"%g ",v(index));}
				fprintf(fpointer,"\n");}}		
			
		fclose(fpointer);}
	
	void Write_DAT_File(std::string filename){
		VECTOR<T>& v=(*(v_hierarchy[0]));
		
		FILE* fpointer;
		fpointer=fopen(filename.c_str(),"w");
		if(dimension==1){
			for(int i=0;i<n_fine;i++)
				fprintf(fpointer,"%g\n",v(i));}
		else if(dimension==2){
			for(int i=0;i<n_fine;i++){
				for(int j=0;j<n_fine;j++){
					INDEX_2D index={i,j,n_fine};
					fprintf(fpointer,"%g ",v(index));}
				fprintf(fpointer,"\n");}}		
		
		fclose(fpointer);}
	
	VECTOR<T>& Fine_Level_RHS(){return (*(rhs_hierarchy[0]));}
	SPARSE_MATRIX<T>& Fine_Level_A(){return (*(A_hierarchy[0]));}
	VECTOR<T>& Fine_Level_V(){return (*(v_hierarchy[0]));}
	
	void Galerkin_Check(){
		for(int l=1;l<levels;l++){
			SPARSE_MATRIX<T>&A_fine=System_Matrix(l-1);
			SPARSE_MATRIX<T>&A_coarse=System_Matrix(l);
			SPARSE_MATRIX<T>&R=Restriction(l-1);
			SPARSE_MATRIX<T>&P=Prolongation(l-1);
			SPARSE_MATRIX<T> RA(R.M(),A_fine.N());
			SPARSE_MATRIX<T> RAP(R.M(),P.N());
			R.Right_Multiply(A_fine,RA);
			RA.Right_Multiply(P,RAP);
			std::cout << "Galerkin product = " << std::endl;
			RAP.Print();
			std::cout << "A coarse = " << std::endl;
			A_coarse.Print();}}
	
	T Dx(){return dx_fine;}
	
	~MULTIGRID_POISSON(){
		//clean up
		for(int l=0;l<levels;l++){
			delete A_hierarchy[l];
			delete v_hierarchy[l];
			delete temp_hierarchy[l];
			delete rhs_hierarchy[l];}
		
		delete[] A_hierarchy;delete[] v_hierarchy;delete[] rhs_hierarchy;delete[] temp_hierarchy;
	}

	void Build_Periodic_Poisson_Matrix(SPARSE_MATRIX<T>& A,const int n,const T dx)
	{
		if(dimension==1){
			for(int i=0;i<A.M();i++){
				SPARSE_ROW<T>& Ai=A.Row(i);
				int ip1=i+1;if(ip1==n) ip1=0;
				int im1=i-1;if(im1==-1) im1=n-1;
				Ai.Add_Entry(im1,-1./(dx*dx));
				Ai.Add_Entry(ip1,-1./(dx*dx));
				Ai.Add_Entry(i,2./(dx*dx));}}
		else if(dimension==2){
			for(int i=0;i<n;i++){
				for(int j=0;j<n;j++){
					INDEX_2D index={i,j,n};
					SPARSE_ROW<T>& Ai=A.Row(index.Index());
					INDEX_2D index_ip1={i+1,j,n};
					INDEX_2D index_im1={i-1,j,n};
					INDEX_2D index_jp1={i,j+1,n};
					INDEX_2D index_jm1={i,j-1,n};
					Ai.Add_Entry(index_jm1.Index(),-1./(dx*dx));
					Ai.Add_Entry(index_jp1.Index(),-1./(dx*dx));
					Ai.Add_Entry(index_im1.Index(),-1./(dx*dx));
					Ai.Add_Entry(index_ip1.Index(),-1./(dx*dx));
					Ai.Add_Entry(index.Index(),4./(dx*dx));}}}
	}	

	void Build_Restriction(SPARSE_MATRIX<T>& R,const int n)
	{	
		if(dimension==1){
			for(int i_coarse=0;i_coarse<R.M();i_coarse++){
				int i=2*i_coarse;
				int ip1=i+1;if(ip1==n) ip1=0;
				int im1=i-1;if(im1==-1) im1=n-1;
				SPARSE_ROW<T>& Ri=R.Row(i_coarse);
				Ri.Add_Entry(i,(T).5);
				Ri.Add_Entry(im1,(T).25);
				Ri.Add_Entry(ip1,(T).25);}}
		else if(dimension==2){
			for(int i_coarse=0;i_coarse<n/2;i_coarse++){
				for(int j_coarse=0;j_coarse<n/2;j_coarse++){
					INDEX_2D coarse_index={i_coarse,j_coarse,n/2};
					SPARSE_ROW<T>& Ri=R.Row(coarse_index.Index());
					for(int i=-1;i<=1;i++){
						for(int j=-1;j<=1;j++){
							INDEX_2D coarse_index_ij={2*coarse_index.i+i,2*coarse_index.j+j,n};
							T wi=(T).5-(T)abs(i)*(T).25;
							T wj=(T).5-(T)abs(j)*(T).25;
							T weight=wi*wj;
							Ri.Add_Entry(coarse_index_ij.Index(),weight);}}}}}
	}
	
	void Gauss_Seidel_Smooth(SPARSE_MATRIX<T>& A,VECTOR<T>& rhs,VECTOR<T>& v,const int number_steps,const int level)
	{
		//VECTOR<T> residual(A.M());
		VECTOR<T>& residual=(*(temp_hierarchy[level]));		
		
		for(int iteration=0;iteration<number_steps;iteration++){
			for(int i=0;i<A.M();i++){
				A.Residual(rhs,v,residual);
				SPARSE_ROW<T>& Ai=A.Row(i);
				T Aii=Ai(i);
				T delta=residual(i)/Aii;
				v(i)=v(i)+delta;}}
	}
	
	void Jacobi_Smooth(SPARSE_MATRIX<T>& A,VECTOR<T>& rhs,VECTOR<T>& v,const int number_steps,const int level)
	{
		VECTOR<T>& residual=(*(temp_hierarchy[level]));		
		
		for(int iteration=0;iteration<number_steps;iteration++){
			A.Residual(rhs,v,residual);
			for(int i=0;i<A.M();i++){
				SPARSE_ROW<T>& Ai=A.Row(i);
				T Aii=Ai(i);
				T delta=residual(i)/Aii;
				v(i)=v(i)+((T)2/(T)3)*delta;}}
	}
	
	void V_cycle()
	{
		for(int l=0;l<levels-1;l++){
			VECTOR<T>& v=(*(v_hierarchy[l]));
			VECTOR<T>& rhs=(*(rhs_hierarchy[l]));
			SPARSE_MATRIX<T>& R=Restriction(l);
			SPARSE_MATRIX<T>& A=System_Matrix(l);
			Jacobi_Smooth(A,rhs,v,1,l);
			VECTOR<T>& r=(*(temp_hierarchy[l]));
			A.Residual(rhs,v,r);
			VECTOR<T>& rhs_2h=(*(rhs_hierarchy[l+1]));
			R.Multiply(r,rhs_2h);
			v_hierarchy[l+1]->Set_To_Zero();
			//int n=A.M();
			//for(int i=0;i<=n/2-1;i++) (*v_hierarchy[l+1])(i)=0;
		
		}
		
		VECTOR<T>& v=(*(v_hierarchy[levels-1]));
		VECTOR<T>& rhs=(*(rhs_hierarchy[levels-1]));
		SPARSE_MATRIX<T>& A=System_Matrix(levels-1);
		Jacobi_Smooth(A,rhs,v,1000,levels-1);
		
		for(int l=levels-2;l>=0;l--){
			VECTOR<T>& v=(*(v_hierarchy[l]));
			VECTOR<T>& v_2h=(*(v_hierarchy[l+1]));
			VECTOR<T>& rhs=(*(rhs_hierarchy[l]));
			SPARSE_MATRIX<T>& P=Prolongation(l);
			SPARSE_MATRIX<T>& A=System_Matrix(l);
			VECTOR<T>& c=(*(temp_hierarchy[l]));
			P.Multiply(v_2h,c);
			v+=c;
			Jacobi_Smooth(A,rhs,v,1,l);}
	}
	
	void Solve(){
		for(int it=0;it<num_v_cycles;it++){
			V_cycle();
			VECTOR<T>& v=(*(v_hierarchy[0]));
			VECTOR<T>& rhs=(*(rhs_hierarchy[0]));
			SPARSE_MATRIX<T>& A=System_Matrix(0);
			VECTOR<T>& r=(*(temp_hierarchy[0]));
			A.Residual(rhs,v,r);
			//std::cout<<"Residual at iteration " <<it<< " = " << r.L_inf() << std::endl;
		}
	}
	
	void Recursive_V_Cycle(SPARSE_MATRIX<T>& A,VECTOR<T>& rhs,VECTOR<T>& v,const int level)
	{
		if(level<levels-1){
			VECTOR<T> temp(A.M());//going to be used for the residual initially
			Gauss_Seidel_Smooth(A,rhs,v,1);
			A.Residual(rhs,v,temp);
			VECTOR<T>& rhs_2h=(*(rhs_hierarchy[level+1]));
			SPARSE_MATRIX<T>& R=(*(R_hierarchy[level]));
			R.Multiply(temp,rhs_2h);
			Solve((*(A_hierarchy[level+1])),rhs_2h,(*(v_hierarchy[level+1])),level+1);
			SPARSE_MATRIX<T>& P=(*(P_hierarchy[level]));
			P.Multiply((*(v_hierarchy[level+1])),temp);//now make temp the correction
			v+=temp;
			Gauss_Seidel_Smooth(A,rhs,v,1);
			//v.Enforce_Zero_Sum();
		}
		else{
			Gauss_Seidel_Smooth(A,rhs,v,1000);}
	}
};