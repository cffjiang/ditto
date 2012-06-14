#include <cassert>
#include <iostream>
#include <string>

struct INDEX_2D{
	int i,j,m;
	
	void Print(){
		std::cout<<"Index 2D = {" << i << " , " << j << " , " << m << "}" << std::endl;}
	
	int Index_X_MAC(){
		assert(i>=0 && i<m+1 && j>=0 && j<m);		
		return i+(m+1)*j;}
	
	int Index_Y_MAC(){
		assert(i>=0 && i<m && j>=0 && j<m+1);
		return i+m*j;}
	
	int Index_Non_Periodic(){
		assert(i>=0 && i <=m);
		assert(j>=0 && j <=m);
		return j*m+i;
	}
	
	int Index_Periodic(){
		int i_periodic=i_Periodic();
		int j_periodic=j_Periodic();
		
		assert(j_periodic*m+i_periodic>=0 && j_periodic*m+i_periodic<m*m);
		return j_periodic*m+i_periodic;}
	
	int Index(){
		assert(j*m+i>=0 && j*m+i<m*m);
		return j*m+i;}
	
	int i_Periodic(){
		int i_periodic;
		if(i<0)
			i_periodic=((((-i)/m)+1)*m+i)%m;
		else
			i_periodic=i%m;
		assert(i_periodic>=0 && i_periodic<m);
		return i_periodic;}
	
	int j_Periodic(){
		int j_periodic;
		if(j<0)
			j_periodic=((((-j)/m)+1)*m+j)%m;
		else
			j_periodic=j%m;
		assert(j_periodic>=0 && j_periodic<m);
		return j_periodic;}
};

template<class T>
class VECTOR_2D{
    T v1,v2;
public:
	VECTOR_2D(const T input):v1(input),v2(input){}
	VECTOR_2D():v1((T)0),v2((T)0){}
	VECTOR_2D(T v1_input,T v2_input):v1(v1_input),v2(v2_input){}
	VECTOR_2D(const VECTOR_2D<T>& v_input):v1(v_input.v1),v2(v_input.v2){}
	
	VECTOR_2D<T>& operator=(const VECTOR_2D<T>& input){v1=input.v1;v2=input.v2;return *this;}
	VECTOR_2D<T> operator-(const VECTOR_2D<T>& input){return VECTOR_2D<T>(v1-input.v1,v2-input.v2);}
	VECTOR_2D<T> operator+(const VECTOR_2D<T>& input){return VECTOR_2D<T>(v1+input.v1,v2+input.v2);}
	VECTOR_2D<T> operator*(const T scale){return VECTOR_2D<T>(scale*v1,scale*v2);}
	
	VECTOR_2D<T> Right_Handed_Perp_Vector(){return VECTOR_2D<T>(-v2,v1);}
	
	static T Dot_Product(const VECTOR_2D<T>& u,const VECTOR_2D<T>& v){return u.Dot(v);}
	static T Signed_Triangle_Area(const VECTOR_2D<T>& u,const VECTOR_2D<T>& v){return (T).5*(u.x_copy()*v.y_copy()-v.x_copy()*u.y_copy());}
	static VECTOR_2D<T> ei(const int i){
		assert(i>=0 && i<2);
		if(i==0)
			return VECTOR_2D<T>((T)1,0);
		else
			return VECTOR_2D<T>(0,(T)1);}
	
	T& x(){return v1;} 
	T& y(){return v2;}
	
	T x_copy()const{return v1;}
	T y_copy()const{return v2;}
	
	T Dot(const VECTOR_2D<T>& v)const {return v.v1*v1+v.v2*v2;}
};

static VECTOR_2D<double> operator*(const double scale,const VECTOR_2D<double>& input){return VECTOR_2D<double>(scale*input.x_copy(),scale*input.y_copy());}
static VECTOR_2D<float> operator*(const float scale,const VECTOR_2D<float>& input){return VECTOR_2D<float>(scale*input.x_copy(),scale*input.y_copy());}

template<class T>
class MATRIX_2X2{
    T a11,a21,a12,a22;//column major
public:
	MATRIX_2X2(const T input):a11(input),a21(input),a12(input),a22(input){}
	MATRIX_2X2(const MATRIX_2X2<T>& A_input):a11(A_input.a11),a21(A_input.a21),a12(A_input.a12),a22(A_input.a22){}
	MATRIX_2X2():a11((T)0),a21((T)0),a12((T)0),a22((T)0){}
	MATRIX_2X2(const T a11_input,const T a21_input,const T a12_input,const T a22_input):a11(a11_input),a21(a21_input),a12(a12_input),a22(a22_input){}
	MATRIX_2X2(const VECTOR_2D<T>& c1,const VECTOR_2D<T>& c2):a11(c1.x_copy()),a21(c1.y_copy()),a12(c2.x_copy()),a22(c2.y_copy()){}
	
	MATRIX_2X2<T>& operator=(const MATRIX_2X2<T>& A_input){a11=A_input.a11;a21=A_input.a21;a12=A_input.a12;a22=A_input.a22;return *this;}
	MATRIX_2X2<T> operator+(const MATRIX_2X2<T>& A_input){return MATRIX_2X2<T>(a11+A_input.a11,a21+A_input.a21,a12+A_input.a12,a22+A_input.a22);}
	MATRIX_2X2<T> operator*(const T scale){return MATRIX_2X2<T>(scale*a11,scale*a21,scale*a12,scale*a22);}
	
	T Determinant(){return a11*a22-a21*a12;}
	
	MATRIX_2X2<T> Transpose(){return MATRIX_2X2<T>(a11,a12,a21,a22);}
	
	T Trace(){return a11+a22;}
	
	static MATRIX_2X2<double> Outer_Product(const VECTOR_2D<T>& u,const VECTOR_2D<T>& v){return MATRIX_2X2(u.x_copy()*v.x_copy(),u.y_copy()*v.x_copy(),u.x_copy()*v.y_copy(),u.y_copy()*v.y_copy());}
	static T Contract(const MATRIX_2X2<T>& A,const MATRIX_2X2<T>& B){return A.a11*B.a11+A.a21*B.a21+A.a12*B.a12+A.a22*B.a22;}
	static MATRIX_2X2<double> Identity(){return MATRIX_2X2((T)1,(T)0,(T)0,(T)1);}
	
	void Invert(){
		T d=Determinant();
		T a11_old=a11;
		T a22_old=a22;
		a22=a11_old/d;
		a11=a22_old/d;
		a12=-a12/d;
		a21=-a21/d;
	}
	
	T operator()(const int i,const int j) const {
		if(i==0 && j==0) return a11;
		if(i==1 && j==0) return a21;
		if(i==0 && j==1) return a12;
		if(i==1 && j==1) return a22;
		assert(false);
		return (T)0;
	};
	
	VECTOR_2D<T> operator*(VECTOR_2D<T>& v_input){return VECTOR_2D<T>(a11*v_input.x()+a12*v_input.y(),a21*v_input.x()+a22*v_input.y());}
	
};

static MATRIX_2X2<double> operator*(const double scale,const MATRIX_2X2<double>& input){return MATRIX_2X2<double>(scale*input(0,0),scale*input(1,0),scale*input(0,1),scale*input(1,1));}
static MATRIX_2X2<float> operator*(const float scale,const MATRIX_2X2<float>& input){return MATRIX_2X2<float>(scale*input(0,0),scale*input(1,0),scale*input(0,1),scale*input(1,1));}

template<class T>
class VECTOR{
    const int n;
    T* values;
public:
    VECTOR(const int n_input):n(n_input) {
        values=new T[n];
		for(int i=0;i<=n-1;i++) values[i]=T(0);}

    ~VECTOR() {delete[] values;}

	T& operator()(INDEX_2D index){
		assert(index.Index()<n && index.Index()>=0);
		return values[index.Index()];}
	
	VECTOR<T>& operator=(const VECTOR<T>& input){
		assert(input.Size()==this->Size());
		for(int i=0;i<n;i++) values[i]=input.values[i];
		return *this;}
	
	T Dot(VECTOR<T>& x){
		T result=(T)0;
		for(int i=0;i<n;i++) result+=values[i]*x(i);
		return result;
	}
	
    T& operator()(const int i) {assert(0<=i && i<=n-1);return values[i];}
	
	void Set_To_Zero(){for(int i=0;i<n;i++) values[i]=0;}
	
	T L_inf(){
		T max_norm=(T)0;
		for(int i=0;i<n;i++) if(fabs(values[i])>max_norm) max_norm=fabs(values[i]);
		return max_norm;}
	
	T Sum(){T sum=(T)0;for(int i=0;i<n;i++) sum+=values[i];return sum;}
	
	int Size() const {return n;}
	
	void Enforce_Zero_Sum(){
		T sum=(T)0;
		for(int i=0;i<n;i++) sum+=values[i];
		for(int i=0;i<n;i++) values[i]-=sum/((T)n);}
	
	void operator+=(VECTOR<T>& x) {assert(n==x.Size());for(int i=0;i<n;i++) values[i]=values[i]+x(i);}
	
	void operator-=(VECTOR<T>& x) {assert(n==x.Size());for(int i=0;i<n;i++) values[i]=values[i]-x(i);}
	
	void Write_DAT_File(std::string file){
		FILE* fpointer;
		fpointer=fopen(file.c_str(),"w");
		for(int i=0;i<n;i++)
			fprintf(fpointer,"%g\n",values[i]);
		fclose(fpointer);}
	
	void Print(){
		std::cout << "Vector =";
		for(int i=0;i<n;i++) std::cout << " " << values[i] << " , ";
		std::cout << "." << std::endl;}
};

template<class T>
class SPARSE_ROW{
    const int n;
    int size;
    int* indices;
    T* values;
public:
    SPARSE_ROW(const int n_input):n(n_input),size(0),indices(0),values(0) {}

    ~SPARSE_ROW() {delete[] indices;delete[] values;}
	
	bool Is_Non_Zero(const int index){for(int i=0;i<size;i++) if(indices[i]==index) return true;return false;}

    T& operator()(const int i){
        assert(0<=i && i<=n-1);
        for(int j=0;j<=size-1;j++) if(indices[j]==i) return values[j];
		assert(false);
		return values[0];}
	
	void Normalize_Row_Sum(){
		T sum=(T)0;
		for(int i=0;i<=size-1;i++) sum+=values[i];
		assert(sum!=0);
		for(int i=0;i<=size-1;i++) values[i]=values[i]/sum;
	}
	
	void Scale(T scale){
		for(int i=0;i<=size-1;i++) values[i]*=scale;}
	
	void Fill_Vector(VECTOR<T>& v){
		assert(v.Size()==n);
		for(int i=0;i<n;i++) v(i)=(T)0;
		for(int i=0;i<size;i++) v(indices[i])=values[i];}
	
	void Print(){
		std::cout << "Sparse Row =";
		for(int i=0;i<n;i++){
			bool found=false;
			for(int j=0;j<=size-1;j++){
				if(indices[j]==i){
					std::cout << " " << values[j] << " , ";
					found=true;}}
			if(!found)
				std::cout << " " << 0 << " , ";}
		std::cout << "." << std::endl;}

	int Number_Nonzero(){return size;}
	
	bool Value_Exists_At_Entry(const int index){
		for(int i=0;i<size;i++) 
			if(indices[i]==index) return true;
		return false;
	}
	
	int Index(const int i_hat){assert(i_hat<size);return indices[i_hat];}
	
	T Value_At_Sparse_Index(const int i_hat){assert(i_hat<size);return values[i_hat];}
	
    T Dot_Product(VECTOR<T>& v){
		assert(v.Size()==n);
        T result=0;for(int i=0;i<=size-1;i++) result+=values[i]*v(indices[i]);
        return result;}

    void Add_Entry(const int index,const T value){
		bool found=false;int entry=0;
		for(int i=0;i<size;i++) if(indices[i]==index){found=true;entry=i;}
		if(found){
			int non_zero_index=indices[entry];
			values[non_zero_index]=value;
			return;}
        size++;int* new_indices=new int[size];T* new_values=new T[size];
        for(int i=0;i<=size-2;i++){
            new_indices[i]=indices[i];new_values[i]=values[i];}
        new_indices[size-1]=index;delete[] indices;indices=new_indices;
        new_values[size-1]=value;delete[] values;values=new_values;}
};

template<class T>
class SPARSE_MATRIX{
    const int m,n;
    SPARSE_ROW<T>** rows;
public:
    SPARSE_MATRIX(const int m_input,const int n_input):m(m_input),n(n_input){
        rows=new SPARSE_ROW<T>*[m];
        for(int i=0;i<=m-1;i++) rows[i]=new SPARSE_ROW<T>(n);}

    ~SPARSE_MATRIX() {for(int i=0;i<=m-1;i++) delete rows[i];delete[] rows;}

    T& operator()(const int i,const int j){
        assert(0<=i && i<=m-1);return (*rows[i])(j);}
	
	void Column(const int j,VECTOR<T>& c){
		assert(j>=0 && j<n && c.Size()==m);
		for(int i=0;i<m;i++){
			SPARSE_ROW<T>& Ai=Row(i);
			if(Ai.Is_Non_Zero(j)) c(i)=Ai(j);
			else c(i)=(T)0;}}
	
	void Normalize_Row_Sums(){
		for(int i=0;i<m;i++) Row(i).Normalize_Row_Sum();
	}
	
	void Right_Multiply(SPARSE_MATRIX<T>& B,SPARSE_MATRIX<T>&AB){
		assert(n==B.M());
		VECTOR<T> column(B.M());
		for(int i=0;i<AB.M();i++){
			for(int j=0;j<B.N();j++){
				B.Column(j,column);
				SPARSE_ROW<T>& Ai=Row(i);
				SPARSE_ROW<T>& ABi=AB.Row(i);
				ABi(j)=Ai.Dot_Product(column);}}} 
	
	void Scale_Rows(T scale){
		for(int i=0;i<=m-1;i++) rows[i]->Scale(scale);}
	
	void Print(){
		std::cout << "Sparse Matrix = " << std::endl;
		for(int i=0;i<m;i++)
			Row(i).Print();}

    SPARSE_ROW<T>& Row(const int i){
        assert(0<=i && i<=m-1);return *rows[i];}
	
	void Residual(VECTOR<T>& rhs,VECTOR<T>& x,VECTOR<T>& r)
	{for(int i=0;i<m;i++){
		SPARSE_ROW<T>& Ai=Row(i);
		r(i)=rhs(i)-Ai.Dot_Product(x);}}
	
	void Multiply(VECTOR<T>& x,VECTOR<T>& b){//as in b=Ax
		assert(b.Size()==m);
		for(int i=0;i<m;i++){
			SPARSE_ROW<T>& Ai=Row(i);
			b(i)=Ai.Dot_Product(x);}}
	
	T A_Norm_Squared(const VECTOR<T>& x){
		T result;
		assert(x.Size()==m);
		for(int i=0;i<m;i++){
			SPARSE_ROW<T>& Ai=Row(i);
			result+=Ai.Dot_Product(x)*x(i);}
		return result;
	}
	
	int M(){return m;}
	
	int N(){return n;}
	
	void Multiply(const VECTOR<T>& x,VECTOR<T>& b)
	{for(int i=0;i<m;i++){
		SPARSE_ROW<T>& ri=Row(i);
		b(i)=ri.Dot_Product(x);}}
	
	void Transpose(SPARSE_MATRIX<T>& transpose)
	{
		assert(transpose.m==n && transpose.n==m);
		for(int i=0;i<m;i++){
			SPARSE_ROW<T>& Ri=Row(i);
			for(int j_hat=0;j_hat<Ri.Number_Nonzero();j_hat++){
				int j=Ri.Index(j_hat);
				SPARSE_ROW<T>& Pj=transpose.Row(j);
				Pj.Add_Entry(i,Ri.Value_At_Sparse_Index(j_hat));}}
	}
};

template<class T>
class Conjugate_Gradient{
	SPARSE_MATRIX<T>& A;
	VECTOR<T>& x;
	VECTOR<T>& b;
	VECTOR<T> r,p,q;
	VECTOR<int>* dirichlet_dofs;
	int max_iterations;
	T tolerance;
	
public:
	Conjugate_Gradient(SPARSE_MATRIX<T>& A_input,VECTOR<T>& x_input,VECTOR<T>& b_input,const int max_it):A(A_input),x(x_input),b(b_input),max_iterations(max_it),r(x_input.Size()),
	p(x_input.Size()),q(x_input.Size()),dirichlet_dofs(0){
		Set_Tolerance((T)1e-6);
	}
	
	void Set_Tolerance(const T tol_input){tolerance=tol_input;}
	
	void Set_Dirichlet_Dofs(VECTOR<int>& dirichlet){
		dirichlet_dofs=&dirichlet;
	}
	
	void Zero_Dirichlet_Residual(){
		if(dirichlet_dofs){
			for(int i=0;i<dirichlet_dofs->Size();i++) r((*dirichlet_dofs)(i))=(T)0;}
	}
	
	int Solve(const bool verbose=false){
		A.Residual(b,x,r);
		Zero_Dirichlet_Residual();
		if(r.L_inf()<tolerance) return 0;
		p=r;
		A.Multiply(p,q);
		T r_dot_r=r.Dot(r);
		T alpha=r_dot_r/p.Dot(q);
		
		for(int it=0;it<max_iterations;it++){
			for(int i=0;i<r.Size();i++){
				x(i)+=alpha*p(i);
				r(i)-=alpha*q(i);}
			
			Zero_Dirichlet_Residual();
			
			if(verbose) std::cout << "Residual at iteration it+1 = " << r.L_inf() << std::endl;
			if(r.L_inf()<tolerance) return it+1;
			
			T r_dot_r_new=r.Dot(r);
			T beta=r_dot_r_new/r_dot_r;
			r_dot_r=r_dot_r_new;
			
			for(int i=0;i<p.Size();i++) p(i)=beta*p(i)+r(i);
			
			A.Multiply(p,q);
			alpha=r_dot_r/p.Dot(q);}
		
		return max_iterations;
	}
	
};
