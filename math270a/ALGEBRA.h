#ifndef _algebra_
#define _algebra_

#include <cassert>
#include <iostream>
#include <string>
#include <cmath>
#include "float.h"
#include <math.h>

namespace ALGEBRA{
class INDEX_2D{
    /* This assumes that we have data at 0 <= i < n, 0 <=j< m. This then maps dofs at these coords to the integer 0 <= i*n+j < m*n
     */
	
    int i,j,m,n;
	
public:
    INDEX_2D(){i=0;j=0;m=0;n=0;}
    INDEX_2D(const int i_input){i=i_input;j=0;m=0;n=0;}
    //INDEX_2D(const int i_input,const int j_input,const int m_input){i=i_input,j=j_input;m=m_input;n=m_input;}
    INDEX_2D(const int i_input,const int j_input,const int m_input,const int n_input){
        i=i_input;
        j=j_input;
        m=m_input;
        n=n_input;}
    INDEX_2D(INDEX_2D& input){
        i=input.i;
        j=input.j;
        m=input.m;
        n=input.n;
    }
    INDEX_2D(const INDEX_2D& input){
        i=input.i;
        j=input.j;
        m=input.m;
        n=input.n;
    }

    INDEX_2D& operator=(const INDEX_2D& index){
        i=index.i;
        j=index.j;
        m=index.m;
        n=index.n;
        return *this;
    }
	
    bool Is_Valid(){
        bool valid=true;
        if(i<0 || i>=m) valid=false;
        if(j<0 || j>=n) valid=false;
        return valid;
    }
	
    static int Index(const int i,const int j,const int m,const int n){
        //assert(j*m+i>=0 && j*m+i<m*n);
        //return j*m+i;
        assert(i*n+j>=0 && i*n+j<m*n);
        return i*n+j;
    }
	
    int& I(){return i;}
    int& J(){return j;}
    int& M(){return m;}
    int& N(){return n;}
	
    int I() const {return i;}
    int J() const {return j;}
    int M() const {return m;}
    int N() const {return n;}
	
    void Print(){
        std::cout<<"Index 2D = {" << i << " , " << j << " , " << m << " , " << n << "}" << std::endl;}
	
/*	int Index_X_MAC(){
        assert(i>=0 && i<m+1 && j>=0 && j<n);		
        return i+(m+1)*j;}
	
	int Index_Y_MAC(){
        assert(i>=0 && i<m && j>=0 && j<n+1);
        return i+m*j;}
	
	int Index_Non_Periodic(){
        assert(i>=0 && i <=m);
        assert(j>=0 && j <=n);
        return j*m+i;
	}
*/
	
    int Index_Periodic(){
        int i_periodic=i_Periodic();
        int j_periodic=j_Periodic();
		
        assert(m==n);
		
        assert(j_periodic*m+i_periodic>=0 && j_periodic*m+i_periodic<m*m);
        return j_periodic*m+i_periodic;}
	
    int Index(){
        //assert(j*m+i>=0 && j*m+i<m*n);
        //return j*m+i;
        assert(i*n+j>=0 && i*n+j<m*n);
        return i*n+j;
    }
	
    int i_Periodic() const {
        assert(m==n);
        int i_periodic;
        if(i<0)
            i_periodic=((((-i)/m)+1)*m+i)%m;
        else
            i_periodic=i%m;
        assert(i_periodic>=0 && i_periodic<m);
        return i_periodic;}
	
    int j_Periodic() const {
        int j_periodic;
        if(j<0)
            j_periodic=((((-j)/m)+1)*m+j)%m;
        else
            j_periodic=j%m;
        assert(j_periodic>=0 && j_periodic<m);
        return j_periodic;}
};

class INDEX_3D{
    /* This assumes that we have data at 0 <= i < n, 0 <= j < m, 0 <= k < mn. This then maps dofs at these coords to the integer 0 <= j*m + i < m*n
     */
    int i,j,k,m,n,mn;
			
public:
    INDEX_3D(){i=0;j=0;k=0;m=0;n=0;mn=0;}
    INDEX_3D(const int i_input){i=i_input;j=0;k=0;m=0;n=0;mn=0;}
    INDEX_3D(const int i_input,const int j_input,const int k_input,const int m_input,const int n_input,const int mn_input){
        i=i_input;
        j=j_input;
        k=k_input;
        m=m_input;
        n=n_input;
        mn=mn_input;
    }
    INDEX_3D(INDEX_3D& input){
        i=input.i;
        j=input.j;
        k=input.k;
        m=input.m;
        n=input.n;
        mn=input.mn;
    }
    INDEX_3D(const INDEX_3D& input){
        i=input.i;
        j=input.j;
        k=input.k;
        m=input.m;
        n=input.n;
        mn=input.mn;
    }
			
    INDEX_3D& operator=(const INDEX_3D& index){
        i=index.i;
        j=index.j;
        k=index.k;
        m=index.m;
        n=index.n;
        mn=index.mn;
        return *this;
    }
			
    static int Index(const int i,const int j,const int k,const int m,const int n,const int mn){
        assert(i*mn*n+j*mn+k>=0 && i*mn*n+j*mn+k<m*n*mn);
        return i*mn*n+j*mn+k;
    }
			
    int& I(){return i;}
    int& J(){return j;}
    int& K(){return k;}
    int& M(){return m;}
    int& N(){return n;}
    int& MN(){return mn;}
			
    int I() const {return i;}
    int J() const {return j;}
    int K() const {return k;}
    int M() const {return m;}
    int N() const {return n;}
    int MN() const {return mn;}
			
    void Print(){
        std::cout<<"Index 3D = {" << i << " , " << j << " , " << k << " , " << m << " , " << n << " , " << mn << "}" << std::endl;}
			
    int Index(){
        assert(i*mn*n+j*mn+k>=0 && i*mn*n+j*mn+k<m*n*mn);
        return i*mn*n+j*mn+k;
    }
};	

template<class T>
class VECTOR_2D{
    T v1,v2;
public:
    VECTOR_2D(const T input):v1(input),v2(input){}
    VECTOR_2D():v1((T)0),v2((T)0){}
    VECTOR_2D(T v1_input,T v2_input):v1(v1_input),v2(v2_input){}
    VECTOR_2D(const VECTOR_2D<T>& v_input):v1(v_input.v1),v2(v_input.v2){}
    VECTOR_2D(const INDEX_2D& indices):v1(indices.I()),v2(indices.J()){}
	
    VECTOR_2D<T>& operator=(const VECTOR_2D<T>& input){v1=input.v1;v2=input.v2;return *this;}
    VECTOR_2D<T> operator-(const VECTOR_2D<T>& input)const{return VECTOR_2D<T>(v1-input.v1,v2-input.v2);}
    VECTOR_2D<T> operator+(const VECTOR_2D<T>& input){return VECTOR_2D<T>(v1+input.v1,v2+input.v2);}
    VECTOR_2D<T> operator*(const T scale){return VECTOR_2D<T>(scale*v1,scale*v2);}
	
    VECTOR_2D<T> operator-(){return VECTOR_2D<T>(-v1,-v2);}
	
    T operator()(const int component)const {
        assert(component==0 || component==1);
        return component?v2:v1;
    }
	
    T& operator()(const int component) {
        assert(component==0 || component==1);
        return component?v2:v1;
    }
	
    void Sort()//from smallest to largest
    {
        if(v1>v2){
            T temp=v1;v1=v2;v2=temp;}
    }
	
    VECTOR_2D<T> Right_Handed_Perp_Vector(){return VECTOR_2D<T>(-v2,v1);}
	
    static T Dot_Product(const VECTOR_2D<T>& u,const VECTOR_2D<T>& v){return u.Dot(v);}
    static T Signed_Triangle_Area(const VECTOR_2D<T>& u,const VECTOR_2D<T>& v){return (T).5*(u.x_copy()*v.y_copy()-v.x_copy()*u.y_copy());}
    static VECTOR_2D<T> ei(const int i){
        assert(i>=0 && i<2);
        if(i==0)
            return VECTOR_2D<T>((T)1,0);
        else
            return VECTOR_2D<T>(0,(T)1);}
	
    T& x() {return v1;} 
    T& y() {return v2;}
    const T& x() const {return v1;} 
    const T& y() const {return v2;}
	
    T x_copy()const{return v1;}
    T y_copy()const{return v2;}
	
    T Dot(const VECTOR_2D<T>& v)const {return v.v1*v1+v.v2*v2;}
	
    T Magnitude() const{return sqrt(v1*v1+v2*v2);}
	
    void Normalize(){
        T n=sqrt(v1*v1+v2*v2);
        if(n!=(T)0){
            v1=v1/n;
            v2=v2/n;}}
	
    VECTOR_2D<T> Normalized() const{
        VECTOR_2D<T> copy(*this);
        T n=sqrt(copy.v1*copy.v1+copy.v2*copy.v2);
        if(n!=(T)0){
            copy.v1=copy.v1/n;
            copy.v2=copy.v2/n;}
        return copy;
    }
	
    void Print(){
        std::cout<<"Vector 2D = (" << v1 << " , " << v2 << ")" << std::endl;}
};
	
static VECTOR_2D<double> operator*(const double scale,const VECTOR_2D<double>& input){return VECTOR_2D<double>(scale*input.x(),scale*input.y());}
static VECTOR_2D<float> operator*(const float scale,const VECTOR_2D<float>& input){return VECTOR_2D<float>(scale*input.x(),scale*input.y());}	
	
template<class T>
class VECTOR_3D{
    T v1,v2,v3;
public:
    VECTOR_3D(const T input):v1(input),v2(input),v3(input){}
    VECTOR_3D():v1((T)0),v2((T)0),v3((T)0){}
    VECTOR_3D(T v1_input,T v2_input,T v3_input):v1(v1_input),v2(v2_input),v3(v3_input){}
    VECTOR_3D(const VECTOR_3D<T>& v_input):v1(v_input.v1),v2(v_input.v2),v3(v_input.v3){}
    VECTOR_3D(const INDEX_3D& indices):v1(indices.I()),v2(indices.J()),v3(indices.K()){}
		
    VECTOR_3D<T>& operator=(const VECTOR_3D<T>& input){v1=input.v1;v2=input.v2;v3=input.v3;return *this;}
    VECTOR_3D<T> operator-(const VECTOR_3D<T>& input){return VECTOR_3D<T>(v1-input.v1,v2-input.v2,v3-input.v3);}
    VECTOR_3D<T> operator+(const VECTOR_3D<T>& input){return VECTOR_3D<T>(v1+input.v1,v2+input.v2,v3+input.v3);}
    VECTOR_3D<T> operator*(const T scale){return VECTOR_3D<T>(scale*v1,scale*v2,scale*v3);}
		
    VECTOR_3D<T> operator-()const{return VECTOR_3D<T>(-v1,-v2,-v3);}
	
    T operator()(const int component)const {
        if(component==0) return v1;
        if (component==1) return v2;
        if (component==2) return v3;
        assert(false);
        return (T)0;
    }
	
    static VECTOR_3D<T> Standard_Basis_Vector(const int i){
        if(i==0) return VECTOR_3D<T>((T)1,0,0);
        if(i==1) return VECTOR_3D<T>(0,(T)1,0);
        if(i==2) return VECTOR_3D<T>(0,0,(T)1);
        assert(false);
        return VECTOR_3D();
    }
	
    T& operator()(const int component){
        if(component==0) return v1;
        if (component==1) return v2;
        if (component==2) return v3;
        assert(false);
        return v1;
    }
	
    VECTOR_3D<T> Sorted()//return the sorted vector without changing the original vector
    {
        VECTOR_3D<T> v=(*this);
        T k;
        for(int i = 0; i< 2; i++)
            for(int j = i+1; j < 3; j++)
                if(v(i) > v(j)){
                    k = v(i);
                    v(i) = v(j);
                    v(j)=k;}			
        return v;
    }
		
    static T Dot_Product(const VECTOR_3D<T>& u,const VECTOR_3D<T>& v){return u.Dot(v);}
    static VECTOR_3D<T> Cross_Product(const VECTOR_3D<T>& u,const VECTOR_3D<T>& v){return VECTOR_3D(u.v2*v.v3-v.v2*u.v3,u.v3*v.v1-u.v1*v.v3,u.v1*v.v2-u.v2*v.v1);}
			
    T& x() {return v1;} 
    T& y() {return v2;}
    T& z() {return v3;}
    const T& x() const {return v1;} 
    const T& y() const {return v2;}
    const T& z() const {return v3;}
		
    T Dot(const VECTOR_3D<T>& v)const {return v.v1*v1+v.v2*v2+v.v3*v3;}
    T Magnitude(){return sqrt(v1*v1+v2*v2+v3*v3);}
    T Magnitude_Squared(){return v1*v1+v2*v2+v3*v3;}
		
    void Normalize(){
        T n=sqrt(v1*v1+v2*v2+v3*v3);
        if(n!=(T)0){
            v1=v1/n;
            v2=v2/n;
            v3=v3/n;}}
	
    void Print(){
        std::cout<<"Vector 3D = (" << v1 << " , " << v2 << " , " << v3 << ")" << std::endl;}
	
	
};

static VECTOR_3D<double> operator*(const double scale,const VECTOR_3D<double>& input){return VECTOR_3D<double>(scale*input.x(),scale*input.y(),scale*input.z());}
static VECTOR_3D<float> operator*(const float scale,const VECTOR_3D<float>& input){return VECTOR_3D<float>(scale*input.x(),scale*input.y(),scale*input.z());}
	
template<class T>
class VECTOR_4D{
    T v1,v2,v3,v4;
public:
    VECTOR_4D(const T input):v1(input),v2(input),v3(input),v4(input){}
    VECTOR_4D():v1((T)0),v2((T)0),v3((T)0),v4((T)0){}
    VECTOR_4D(T v1_input,T v2_input,T v3_input,T v4_input):v1(v1_input),v2(v2_input),v3(v3_input),v4(v4_input){}
    VECTOR_4D(const VECTOR_4D<T>& v_input):v1(v_input.v1),v2(v_input.v2),v3(v_input.v3),v4(v_input.v4){}
	
    VECTOR_4D<T>& operator=(const VECTOR_4D<T>& input){v1=input.v1;v2=input.v2;v3=input.v3;v4=input.v4;return *this;}
	
    T operator()(const int component)const {
        if(component==0) return v1;
        if (component==1) return v2;
        if (component==2) return v3;
        if (component==3) return v4;
        assert(false);
        return (T)0;
    }
	
    T& operator()(const int component){
        if(component==0) return v1;
        if (component==1) return v2;
        if (component==2) return v3;
        if (component==3) return v4;
        assert(false);
        return v1;
    }
};
	
template<class T>
class MATRIX_2X2{
public:
    T a11,a21,a12,a22;//column major
    MATRIX_2X2(const T input):a11(input),a21(input),a12(input),a22(input){}
    MATRIX_2X2(const MATRIX_2X2<T>& A_input):a11(A_input.a11),a21(A_input.a21),a12(A_input.a12),a22(A_input.a22){}
    MATRIX_2X2():a11((T)0),a21((T)0),a12((T)0),a22((T)0){}
    MATRIX_2X2(const T a11_input,const T a21_input,const T a12_input,const T a22_input):a11(a11_input),a21(a21_input),a12(a12_input),a22(a22_input){}
    MATRIX_2X2(const VECTOR_2D<T>& c1,const VECTOR_2D<T>& c2):a11(c1.x_copy()),a21(c1.y_copy()),a12(c2.x_copy()),a22(c2.y_copy()){}
    // MATRIX_2X2(const VECTOR_2D<T> c1,const VECTOR_2D<T> c2):a11(c1.x_copy()),a21(c1.y_copy()),a12(c2.x_copy()),a22(c2.y_copy()){}
    MATRIX_2X2(const VECTOR_2D<T>& diagonal):a11(diagonal.x()),a21((T)0),a12((T)0),a22((T)diagonal.y()){}
	
    MATRIX_2X2<T>& operator=(const MATRIX_2X2<T>& A_input){a11=A_input.a11;a21=A_input.a21;a12=A_input.a12;a22=A_input.a22;return *this;}
    MATRIX_2X2<T> operator+(const MATRIX_2X2<T>& A_input){return MATRIX_2X2<T>(a11+A_input.a11,a21+A_input.a21,a12+A_input.a12,a22+A_input.a22);}
    MATRIX_2X2<T> operator-(const MATRIX_2X2<T>& A_input){return MATRIX_2X2<T>(a11-A_input.a11,a21-A_input.a21,a12-A_input.a12,a22-A_input.a22);}
    MATRIX_2X2<T> operator*(const T scale) const{return MATRIX_2X2<T>(scale*a11,scale*a21,scale*a12,scale*a22);}
    VECTOR_2D<T> operator*(const VECTOR_2D<T>& x) const{return VECTOR_2D<T>(x.x()*a11+x.y()*a12,x.x()*a21+x.y()*a22);}
    MATRIX_2X2<T> operator*(const MATRIX_2X2<T>& B) const{return MATRIX_2X2<T>(a11*B.a11+a12*B.a21,a21*B.a11+a22*B.a21,a11*B.a12+a12*B.a22,a21*B.a12+a22*B.a22);}
	
    VECTOR_2D<T> Column(const int number){
        assert(number>=0 && number<2);
        if(number==0) return VECTOR_2D<T>(a11,a21);
        else return VECTOR_2D<T>(a12,a22);		
    }

    T Determinant(){return a11*a22-a21*a12;}
	
    T Norm(){return std::sqrt(a11*a11+a12*a12+a21*a21+a22*a22);}

    MATRIX_2X2<T> Data(){return MATRIX_2X2<T>(a11,a21,a12,a22);}

    MATRIX_2X2<T> Transpose(){return MATRIX_2X2<T>(a11,a12,a21,a22);}
	
    MATRIX_2X2<T> Transposed() const{
        return MATRIX_2X2(a11,a12,a21,a22);
    }
	
    T Trace(){return a11+a22;}
	
    void SVD(MATRIX_2X2<T>& U,VECTOR_2D<T>& sigma,MATRIX_2X2<T>& V,const T tol=(T)1e-10)const{
        // Compute V using FtF
        MATRIX_2X2<T> F=(*this);
        MATRIX_2X2<T> Ft=this->Transposed();
        MATRIX_2X2<T> FtF=Ft*F;
        T c,s;
        if(std::abs(FtF(1,0))<tol){
            c=1;s=0;}
        else{
            T tau=(FtF(0,0)-FtF(1,1))/(2.0*FtF(1,0));
            T t1=tau+std::sqrt(1+tau*tau),t2=tau-std::sqrt(1+tau*tau),t;
            if(std::abs(t1)>std::abs(t2)) t=t2;
            else t=t1;
            c=1.0/std::sqrt(1.0+t*t);
            s=t*c;}
        V=MATRIX_2X2<T>(c,-s,s,c);

        // Sort V
        MATRIX_2X2<T> Sigma_sq=V.Transposed()*F.Transposed()*F*V;
        if(Sigma_sq(0,0)<Sigma_sq(1,1)){
            V=MATRIX_2X2<T>(s,c,-c,s);}

        // FV=QR=USigma
        MATRIX_2X2<T> FV=F*V;

        if(std::abs(FV(1,0))<tol){
            c=1;s=0;}
        else{
            T alpha=1.0/std::sqrt(FV(0,0)*FV(0,0)+FV(1,0)*FV(1,0));
            s=-FV(1,0)*alpha;
            c=FV(0,0)*alpha;}
        U=MATRIX_2X2<T>(c,-s,s,c);
        MATRIX_2X2<T> Sigma=U.Transposed()*F*V;
        sigma(0)=Sigma(0,0);
        sigma(1)=Sigma(1,1);

        // sign convention
        if((sigma(0)>0&&sigma(1)<0&&std::abs(sigma(1))>std::abs(sigma(0)))||(sigma(0)<0&&sigma(1)>0&&std::abs(sigma(0))>=std::abs(sigma(1)))||(sigma(0)<=0&&sigma(1)<=0)){
            sigma(0)=-sigma(0);
            sigma(1)=-sigma(1);
            U=(-1.0)*U;}
    }
	
    void Delta_Sigma(const MATRIX_2X2<T>& delta_F, VECTOR_2D<T>& delta_sigma)const{
        MATRIX_2X2<T> F=*this,U,V,UtDFV;
        VECTOR_2D<T> sigma;
        SVD(U,sigma,V);
        UtDFV=U.Transposed()*delta_F*V;
        delta_sigma=VECTOR_2D<T>(UtDFV(0,0),UtDFV(1,1));
    }
	
    void Delta_SVD(const MATRIX_2X2<T>& delta_F,VECTOR_2D<T>& delta_sigma,MATRIX_2X2<T>& delta_U,MATRIX_2X2<T>& delta_V)const{
        MATRIX_2X2<T> F=*this;
        MATRIX_2X2<T> U,V,UtDFV;
        VECTOR_2D<T> sigma;
        F.SVD(U,sigma,V);
        if(std::abs(sigma(0)-sigma(1))<1e-10 || std::abs(sigma(0)+sigma(1))<1e-10) std::cout<<"FATAL ERROR: Weird case detected -- repeated sigma.\n";

        UtDFV=U.Transposed()*delta_F*V;
        delta_sigma=VECTOR_2D<T>(UtDFV(0,0),UtDFV(1,1));

        VECTOR_2D<T> rhs(UtDFV(1,0),UtDFV(0,1));
        MATRIX_2X2<T> A(-sigma(0),sigma(1),sigma(1),-sigma(0));
        A.Invert();
        VECTOR_2D<T> xy=A*rhs;
        T x=xy(0);
        MATRIX_2X2<T> helperU(0,x,-x,0);
        delta_U=(helperU*(U.Transposed())).Transposed();
        T y=xy(1);
        MATRIX_2X2<T> helperV(0,y,-y,0);
        delta_V=(helperV*(V.Transposed())).Transposed();
    }
	
    void Delta_RS(const MATRIX_2X2<T>& delta_F,MATRIX_2X2<T>& delta_R,MATRIX_2X2<T>& delta_S)const{
        MATRIX_2X2<T> F=*this,U,V,R,S,Sigma;
        VECTOR_2D<T> sigma;
        F.SVD(U,sigma,V);
        if(std::abs(sigma(0)+sigma(1))<1e-10) std::cout<<"FATAL ERROR: Weird case detected -- sigma(0)+signma(1)=0, R not differentiable.\n";

        Sigma=MATRIX_2X2<T>(sigma(0),0,0,sigma(1));
        R=U*(V.Transposed());
        S=V*Sigma*(V.Transposed());

        MATRIX_2X2<T> RtdF=R.Transposed()*delta_F;
        T a=RtdF(0,0),b=RtdF(0,1),c=RtdF(1,0),d=RtdF(1,1),g=S(1,1),e=S(0,0),f=S(0,1);
        T x=(b-c)/(g+e);
        T y=a-x*f;
        T z=d+x*f;
        T w=b-x*g;

        delta_R=MATRIX_2X2<T>(-x*R(0,1),-x*R(1,1),x*R(0,0),x*R(1,0));
        delta_S=MATRIX_2X2<T>(y,w,w,z);
    }

    static MATRIX_2X2<T> Givens_Rotation(const VECTOR_2D<T>& x,const T tol=(T)1e-10){
        //This returns the matrix G=[c,-s;s,c] such that Gx=[X;0]
        T denominator=x.Magnitude();
        if(denominator>tol){
            T c=x.x()/denominator;T s=-x.y()/denominator;
            return MATRIX_2X2<T>(c,s,-s,c);}
        else return MATRIX_2X2<T>((T)1,(T)0,(T)0,(T)1);
    }
	
    static MATRIX_2X2<double> Outer_Product(const VECTOR_2D<T>& u,const VECTOR_2D<T>& v){return MATRIX_2X2(u.x_copy()*v.x_copy(),u.y_copy()*v.x_copy(),u.x_copy()*v.y_copy(),u.y_copy()*v.y_copy());}
    static T Contract(const MATRIX_2X2<T>& A,const MATRIX_2X2<T>& B){return A.a11*B.a11+A.a21*B.a21+A.a12*B.a12+A.a22*B.a22;}
    static MATRIX_2X2<double> Identity(){return MATRIX_2X2((T)1,(T)0,(T)0,(T)1);}

    void QR(MATRIX_2X2<T>&Q,MATRIX_2X2<T>&R){
        T r,s,c;
        if(a21==(double)0){
            c=1;s=0;r=a11;}
        else{
            if(fabs(a21)>fabs(a11)){
                T r=-a11/a21;s=(double)1/sqrt((double)1+r*r);
                c=s*r;}
            else {
                T r=-a21/a11;c=(double)1/sqrt((double)1+r*r);
                s=c*r;}}		
        Q=MATRIX_2X2<T>(c,s,-s,c);
        VECTOR_2D<T> c2=Q.Transpose()*this->Column(1);
        R=MATRIX_2X2<T>(r,0,c2.x_copy(),c2.y_copy());
    }
	
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
    }

    T& operator()(const int i,const int j){
        if(i==0 && j==0) return a11;
        if(i==1 && j==0) return a21;
        if(i==0 && j==1) return a12;
        if(i==1 && j==1) return a22;
        assert(false);
    }

    VECTOR_2D<T> operator*(VECTOR_2D<T>& v_input) const {return VECTOR_2D<T>(a11*v_input.x()+a12*v_input.y(),a21*v_input.x()+a22*v_input.y());}
	
    void Print(){
        std::cout<<"Matrix 2X2 = (" << a11 << " , " << a12 << std::endl;
        std::cout<<"\t\t\t\t" << a21 << " , " << a22 << ")" << std::endl;
    }
	
    MATRIX_2X2<T> Symmetric_Part(){
        return MATRIX_2X2<T>(a11,(T).5*(a12+a21),(T).5*(a12+a21),a22);
    }

    static MATRIX_2X2<T> Deformation_Gradient(const VECTOR_2D<T> &X1,const VECTOR_2D<T> &X2,const VECTOR_2D<T> &X3,const VECTOR_2D<T> &x1,const VECTOR_2D<T> &x2,const VECTOR_2D<T> &x3){
        VECTOR_2D<T> es2,es3,em2,em3;
        es2=x2-x1;es3=x3-x1;em2=X2-X1;em3=X3-X1;
        MATRIX_2X2<T> Ds(es2,es3);
        MATRIX_2X2<T> Dm(em2,em3);
        T eps=1e-15;
        assert(Dm.Determinant()>=eps||Dm.Determinant()<=-eps);
        Dm.Invert();
        return Ds*Dm;
    }

    static MATRIX_2X2<T> Deformation_Gradient(const MATRIX_2X2<T> &Dm_inv,const VECTOR_2D<T> &x1,const VECTOR_2D<T> &x2,const VECTOR_2D<T> &x3){
        VECTOR_2D<T> es2,es3;
        es2=x2-x1;es3=x3-x1;
        MATRIX_2X2<T> Ds(es2,es3);
        return Ds*Dm_inv;
    }
	
};

static MATRIX_2X2<double> operator*(const double scale,const MATRIX_2X2<double>& input){return MATRIX_2X2<double>(scale*input(0,0),scale*input(1,0),scale*input(0,1),scale*input(1,1));}
static MATRIX_2X2<float> operator*(const float scale,const MATRIX_2X2<float>& input){return MATRIX_2X2<float>(scale*input(0,0),scale*input(1,0),scale*input(0,1),scale*input(1,1));}

template<class T>
class MATRIX_3X3{
public:
    T x[9];//column major, TODO: make this private, the current implementations of the static multiply is the culprit
public:	
    MATRIX_3X3(const T input){
        for(int i=0;i<9;i++) x[i]=input;
    }
    MATRIX_3X3(const MATRIX_3X3<T>& A_input)
    {
        for(int i=0;i<9;i++) x[i]=A_input.x[i];
    }
    MATRIX_3X3(){
        for(int i=0;i<9;i++) x[i]=(T)0;
    }
    MATRIX_3X3(const T a11_input,const T a21_input,const T a31_input,const T a12_input,const T a22_input,const T a32_input,const T a13_input,const T a23_input,const T a33_input)
    {
        x[0]=a11_input;x[3]=a12_input;x[6]=a13_input;
        x[1]=a21_input;x[4]=a22_input;x[7]=a23_input;
        x[2]=a31_input;x[5]=a32_input;x[8]=a33_input;
    }
    MATRIX_3X3(const VECTOR_3D<T>& c1,const VECTOR_3D<T>& c2,const VECTOR_3D<T>& c3){
        x[0]=c1(0);x[1]=c1(1);x[2]=c1(2);
        x[3]=c2(0);x[4]=c2(1);x[5]=c2(2);
        x[6]=c3(0);x[7]=c3(1);x[8]=c3(2);
    }
	
    MATRIX_3X3(const VECTOR_3D<T>&v){
        x[0]=v(0);x[1]=0;x[2]=0;
        x[3]=0;x[4]=v(1);x[5]=0;
        x[6]=0;x[7]=0;x[8]=v(2);
    }
	
    VECTOR_3D<T> operator*(const VECTOR_3D<T>& x_in) const{
        return VECTOR_3D<T>(x[0]*x_in(0)+x[3]*x_in(1)+x[6]*x_in(2),
            x[1]*x_in(0)+x[4]*x_in(1)+x[7]*x_in(2),
            x[2]*x_in(0)+x[5]*x_in(1)+x[8]*x_in(2));
    }
	
    MATRIX_3X3<T> operator+(const MATRIX_3X3<T>& A_input)
    {return MATRIX_3X3<T>(x[0]+A_input.x[0],x[1]+A_input.x[1],x[2]+A_input.x[2],x[3]+A_input.x[3],x[4]+A_input.x[4],x[5]+A_input.x[5],x[6]+A_input.x[6],x[7]+A_input.x[7],x[8]+A_input.x[8]);}
    MATRIX_3X3<T> operator*(const T scale) const{return MATRIX_3X3<T>(scale*x[0],scale*x[1],scale*x[2],scale*x[3],scale*x[4],scale*x[5],scale*x[6],scale*x[7],scale*x[8]);}
	
    T operator()(const int i,const int j) const {
        int index=3*j+i;
        assert(index>=0&&index<9);
        return x[index];
    }
	
    T& operator()(const int i,const int j){
        int index=3*j+i;
        assert(index>=0&&index<9);
        return x[index];
    }
	
    static MATRIX_3X3<T> Identity(){return MATRIX_3X3<T>((T)1,(T)0,(T)0,(T)0,(T)1,(T)0,(T)0,(T)0,(T)1);}
    static T Trace(const MATRIX_3X3<T>& M){return M.Trace();}
	
    T Trace() const
    {return x[0]+x[4]+x[8];}
	
    T Determinant() const
    {return x[0]*(x[4]*x[8]-x[7]*x[5])+x[3]*(x[7]*x[2]-x[1]*x[8])+x[6]*(x[1]*x[5]-x[4]*x[2]);}
	
    void Transpose(){
        T temp1=x[1];T temp2=x[2];T temp5=x[5];
        x[1]=x[3];x[2]=x[6];x[5]=x[7];
        x[3]=temp1;x[6]=temp2;x[7]=temp5;
    }

    T Norm(){
        T n=0;
        for(int i=0;i<9;i++) n+=x[i]*x[i];
        return std::sqrt(n);
    }
	
    MATRIX_3X3<T> Transposed() const{
        return MATRIX_3X3(x[0],x[3],x[6],x[1],x[4],x[7],x[2],x[5],x[8]);
    }
	
    VECTOR_3D<T> Column(const int column)const {
        if(column==0) return VECTOR_3D<T>(x[0],x[1],x[2]);
        else if(column==1) return VECTOR_3D<T>(x[3],x[4],x[5]);
        else if(column==2) return VECTOR_3D<T>(x[6],x[7],x[8]);
        else{
            assert(false);
            return VECTOR_3D<T>(0,0,0);}
    }
	
    VECTOR_3D<T> Row(const int row)const {
        if(row==0) return VECTOR_3D<T>(x[0],x[3],x[6]);
        else if(row==1) return VECTOR_3D<T>(x[1],x[4],x[7]);
        else if(row==2) return VECTOR_3D<T>(x[2],x[5],x[8]);
        else{
            assert(false);
            return VECTOR_3D<T>(0,0,0);}
    }
	
    void SVD(MATRIX_3X3<T>& U,VECTOR_3D<T>& sigma,MATRIX_3X3<T>& V,const T tol=(T)1e-10,const int max_iterations=20)const {
        T c,s;
        MATRIX_2X2<T> A12,A23,A13;
        MATRIX_3X3<T> G12,G23,G13;
        MATRIX_3X3<T> F=*this;
        MATRIX_3X3<T> A=F.Transposed()*F;

        // Jacobi iterations for V
        V=MATRIX_3X3<T>::Identity();
        for(int iter=0;iter<max_iterations;iter++){
            A12=MATRIX_2X2<T>(A(0,0),A(1,0),A(0,1),A(1,1));
            if(fabs(A12(1,0))<tol){
                c=1;s=0;}
            else{
                T tal=(A12(0,0)-A12(1,1))/(2*A12(1,0));
                T t1=tal+sqrt(tal*tal+1);
                T t2=tal-sqrt(tal*tal+1);
                T t=fabs(t1)<fabs(t2)?t1:t2;
                c=1/sqrt(1+t*t);
                s=t*c;}
            G12=MATRIX_3X3<T>(c,-s,0,s,c,0,0,0,1);
            A=G12.Transposed()*A*G12;
            V=V*G12;
             
            A23=MATRIX_2X2<T>(A(1,1),A(2,1),A(1,2),A(2,2));
            if(fabs(A23(1,0))<tol){
                c=1;s=0;}
            else{
                T tal=(A23(0,0)-A23(1,1))/(2*A23(1,0));
                T t1=tal+sqrt(tal*tal+1);
                T t2=tal-sqrt(tal*tal+1);
                T t=fabs(t1)<fabs(t2) ? t1:t2;
                c=1/sqrt(1+t*t);
                s=t*c;}
            G23=MATRIX_3X3<T>(1,0,0,0,c,-s,0,s,c);
            A=G23.Transposed()*A*G23;
            V=V*G23;
 
            A13=MATRIX_2X2<T>(A(0,0),A(2,0),A(0,2),A(2,2));
            if(fabs(A13(1,0))< tol){
                c=1; s=0;}
            else{
                T tal=(A13(0,0)-A13(1,1))/(2*A13(1,0));
                T t1=tal+sqrt(tal*tal+1);
                T t2=tal-sqrt(tal*tal+1);
                T t=fabs(t1)<fabs(t2) ? t1:t2;
                c=1/sqrt(1+t*t);
                s=t*c;}
            G13=MATRIX_3X3<T>(c,0,-s,0,1,0,s,0,c);
            A=G13.Transposed()*A*G13;
            V=V*G13;}

        A=F.Transposed()*F;

        // Sort sigma squared and rearrange V
        MATRIX_3X3<T> Sigma_sq=V.Transposed()*A*V;
        if(Sigma_sq(2,2)>Sigma_sq(1,1)){
            VECTOR_3D<T> temp=V.Column(2);
            V(0,2)=-V(0,1);
            V(1,2)=-V(1,1);
            V(2,2)=-V(2,1);
            V(0,1)=temp.x();
            V(1,1)=temp.y();
            V(2,1)=temp.z();
            T temp_s=Sigma_sq(1,1);
            Sigma_sq(1,1)=Sigma_sq(2,2);
            Sigma_sq(2,2)=temp_s;}
        if(Sigma_sq(1,1)>Sigma_sq(0,0)){
            VECTOR_3D<T> temp=V.Column(0);
            V(0,0)=-V(0,1);
            V(1,0)=-V(1,1);
            V(2,0)=-V(2,1);
            V(0,1)=temp.x();
            V(1,1)=temp.y();
            V(2,1)=temp.z();
            T temp_s=Sigma_sq(1,1);
            Sigma_sq(1,1)=Sigma_sq(0,0);
            Sigma_sq(0,0)=temp_s;}
        if(Sigma_sq(2,2)>Sigma_sq(1,1)){
            VECTOR_3D<T> temp=V.Column(2);
            V(0,2)=-V(0,1);
            V(1,2)=-V(1,1);
            V(2,2)=-V(2,1);
            V(0,1)=temp.x();
            V(1,1)=temp.y();
            V(2,1)=temp.z();
            T temp_s=Sigma_sq(1,1);
            Sigma_sq(1,1)=Sigma_sq(2,2);
            Sigma_sq(2,2)=temp_s;}

        // QR to get U and Sigma
        MATRIX_3X3<T> FV=(*this)*V;
        MATRIX_3X3<T> Sigma;
        FV.QR(U,Sigma,tol);
        sigma=VECTOR_3D<T>(Sigma(0,0),Sigma(1,1),Sigma(2,2));

        // sign convention
        int N_negative=(sigma(0)<0)+(sigma(1)<0)+(sigma(2)<0);
        if(N_negative==0) return;
        if(N_negative==1){
            if(sigma(0)<0){
                sigma(0)=-sigma(0);
                sigma(2)=-sigma(2);
                for(int i=0;i<3;i++) V(i,0)=-V(i,0);
                for(int i=0;i<3;i++) V(i,2)=-V(i,2);}
            else if(sigma(1)<0){
                sigma(1)=-sigma(1);
                sigma(2)=-sigma(2);
                for(int i=0;i<3;i++) V(i,1)=-V(i,1);
                for(int i=0;i<3;i++) V(i,2)=-V(i,2);}}
        if(N_negative==2){
            if(sigma(0)>=0){
                sigma(1)=-sigma(1);
                sigma(2)=-sigma(2);
                for(int i=0;i<3;i++) V(i,1)=-V(i,1);
                for(int i=0;i<3;i++) V(i,2)=-V(i,2);}
            else if(sigma(1)>=0){
                sigma(0)=-sigma(0);
                sigma(2)=-sigma(2);
                for(int i=0;i<3;i++) V(i,0)=-V(i,0);
                for(int i=0;i<3;i++) V(i,2)=-V(i,2);}
            else if(sigma(2)>=0){
                sigma(0)=-sigma(0);
                sigma(1)=-sigma(1);
                for(int i=0;i<3;i++) V(i,0)=-V(i,0);
                for(int i=0;i<3;i++) V(i,1)=-V(i,1);}}
        if(N_negative==3){
                sigma(0)=-sigma(0);
                sigma(1)=-sigma(1);
                for(int i=0;i<3;i++) V(i,0)=-V(i,0);
                for(int i=0;i<3;i++) V(i,1)=-V(i,1);}
    }
     
    void QR(MATRIX_3X3<T>&Q,MATRIX_3X3<T>&R,const T tol){
        T c,s;
        MATRIX_3X3<T> G1(0),G2(0),G3(0);
        if(fabs(x[2])<tol){
            c=1;s=0;
            G1(0,0)=G1(1,1)=G1(2,2)=1;}
        else{
            T alpha=1/sqrt(x[1]*x[1]+x[2]*x[2]);
            s=-x[2]*alpha;
            c=x[1]*alpha;
            G1(0,0)=1;
            G1(1,1)=c;
            G1(2,1)=s;
            G1(1,2)=-s;
            G1(2,2)=c;}
         
        MATRIX_3X3<T> A1;
        A1=G1*(*this);
        if(fabs(A1(1,0))<tol){
            c=1;s=0;
            G2(0,0)=G2(1,1)=G2(2,2)=1;}
        else{
            T alpha=1/sqrt(A1(0,0)*A1(0,0)+A1(1,0)*A1(1,0));
            s=-A1(1,0)*alpha;
            c=A1(0,0)*alpha;
            G2(2,2)=1;
            G2(0,0)=c;
            G2(1,0)=s;
            G2(0,1)=-s;
            G2(1,1)=c;}
         
        MATRIX_3X3<T> A2;
        A2=G2*A1;
        if(fabs(A2(2,1))<tol){
            c=1;s=0;
            G3(0,0)=G3(1,1)=G3(2,2)=1;}
        else{
            T alpha=1/sqrt(A2(1,1)*A2(1,1)+A2(2,1)*A2(2,1));
            s=-A2(2,1)*alpha;
            c=A2(1,1)*alpha;
            G3(0,0)=1;
            G3(1,1)=c;
            G3(2,1)=s;
            G3(1,2)=-s;
            G3(2,2)=c;}
         
        Q=(G3*G2*G1).Transposed();
        R=Q.Transposed()*(*this);
    }

    void Delta_Sigma(const MATRIX_3X3<T>& delta_F, VECTOR_3D<T>& delta_sigma)const{
        MATRIX_3X3<T> F=*this,U,V,UtDFV;
        VECTOR_3D<T> sigma;
        SVD(U,sigma,V);
        UtDFV=U.Transposed()*delta_F*V;
        delta_sigma=VECTOR_3D<T>(UtDFV(0,0),UtDFV(1,1));
    }
	
    void Delta_SVD(const MATRIX_3X3<T>& delta_F,VECTOR_3D<T>& delta_sigma,MATRIX_3X3<T>& delta_U,MATRIX_3X3<T>& delta_V)const{
        MATRIX_3X3<T> F=*this,U,V,UtDFV;
        VECTOR_3D<T> sigma;
        F.SVD(U,sigma,V);
        UtDFV=U.Transposed()*delta_F*V;
        T sigma1=sigma(0),sigma2=sigma(1),sigma3=sigma(2);
        if(std::abs(sigma1-sigma2)<1e-10||std::abs(sigma1-sigma3)<1e-10||std::abs(sigma2-sigma3)<1e-10||std::abs(sigma1+sigma2)<1e-10||std::abs(sigma1+sigma3)<1e-10||std::abs(sigma2+sigma3)<1e-10) std::cout<<"FATAL ERROR: Repeated sigma detected.\n";

        delta_sigma=VECTOR_3D<T>(UtDFV(0,0),UtDFV(1,1),UtDFV(2,2));
        T a=UtDFV(0,1),b=UtDFV(0,2),c=UtDFV(1,2),d=UtDFV(1,0),e=UtDFV(2,0),f=UtDFV(2,1);
        MATRIX_2X2<T> xpM(sigma2,-sigma1,-sigma1,sigma2);xpM.Invert();
        MATRIX_2X2<T> yqM(sigma3,-sigma1,-sigma1,sigma3);yqM.Invert(); 
        MATRIX_2X2<T> zrM(sigma3,-sigma2,-sigma2,sigma3);zrM.Invert();
        VECTOR_2D<T> xprhs(a,d),yqrhs(b,e),zrrhs(c,f);
        VECTOR_2D<T> xp=xpM*xprhs,yq=yqM*yqrhs,zr=zrM*zrrhs;
        T x=xp(0),p=xp(1),y=yq(0),q=yq(1),z=zr(0),r=zr(1);
        MATRIX_3X3<T> dUtU(0,x,y,-x,0,z,-y,-z,0);
        MATRIX_3X3<T> dVtV(0,p,q,-p,0,r,-q,-r,0);
        delta_U=(dUtU*U.Transposed()).Transposed();
        delta_V=(dVtV*V.Transposed()).Transposed();
    }

    void Delta_RS(const MATRIX_3X3<T>& delta_F,MATRIX_3X3<T>& delta_R,MATRIX_3X3<T>& delta_S)const{
        MATRIX_3X3<T> F=*this,U,V,R,S,Sigma;
        VECTOR_3D<T> sigma;
        F.SVD(U,sigma,V);
        if(std::abs(sigma(0)+sigma(1))<1e-10||std::abs(sigma(1)+sigma(2))<1e-10||std::abs(sigma(0)+sigma(2))<1e-10) std::cout<<"FATAL ERROR: Weird case detected -- sigma(i)+signma(j)=0, R not differentiable.\n";

        Sigma=MATRIX_3X3<T>(sigma(0),0,0,0,sigma(1),0,0,0,sigma(2));
        R=U*(V.Transposed());
        S=V*Sigma*(V.Transposed());

        MATRIX_3X3<T> W=R.Transposed()*delta_F;
        T trS=S(0,0)+S(1,1)+S(2,2);
        MATRIX_3X3<T> SmtrS=S-trS*Identity();
        SmtrS.Invert();
        VECTOR_3D<T> ws(W(1,2)-W(2,1),W(2,0)-W(0,2),W(0,1)-W(1,0));
        VECTOR_3D<T> r=SmtrS*ws;
        T r1=r(0),r2=r(1),r3=r(2);
        MATRIX_3X3<T> rMat(0,r3,-r2,-r3,0,r1,r2,-r1,0);
        delta_R=R*rMat;
        delta_S=W-rMat*S;
    }
	
    void Invert(){
        T cofactor11=x[4]*x[8]-x[7]*x[5],cofactor12=x[7]*x[2]-x[1]*x[8],cofactor13=x[1]*x[5]-x[4]*x[2];
        T determinant=x[0]*cofactor11+x[3]*cofactor12+x[6]*cofactor13;assert(determinant!=0);T s=1/determinant;
        T a11=s*cofactor11;
        T a21=s*cofactor12;
        T a31=s*cofactor13;
        T a12=s*(x[6]*x[5]-x[3]*x[8]);
        T a22=s*(x[0]*x[8]-x[6]*x[2]);
        T a32=s*(x[3]*x[2]-x[0]*x[5]);
        T a13=s*(x[3]*x[7]-x[6]*x[4]);
        T a23=s*(x[6]*x[1]-x[0]*x[7]);
        T a33=s*(x[0]*x[4]-x[3]*x[1]);
        x[0]=a11;x[3]=a12;x[6]=a13;
        x[1]=a21;x[4]=a22;x[7]=a23;
        x[2]=a31;x[5]=a32;x[8]=a33;
    }
	
    MATRIX_3X3<T> Cofactor_Matrix(){
        T cofactor11=x[4]*x[8]-x[7]*x[5],cofactor12=x[7]*x[2]-x[1]*x[8],cofactor13=x[1]*x[5]-x[4]*x[2];
        T a11=cofactor11;
        T a21=cofactor12;
        T a31=cofactor13;
        T a12=(x[6]*x[5]-x[3]*x[8]);
        T a22=(x[0]*x[8]-x[6]*x[2]);
        T a32=(x[3]*x[2]-x[0]*x[5]);
        T a13=(x[3]*x[7]-x[6]*x[4]);
        T a23=(x[6]*x[1]-x[0]*x[7]);
        T a33=(x[0]*x[4]-x[3]*x[1]);
        T x0=a11;T x3=a12;T x6=a13;
        T x1=a21;T x4=a22;T x7=a23;
        T x2=a31;T x5=a32;T x8=a33;
        return MATRIX_3X3<T>(x0,x1,x2,x3,x4,x5,x6,x7,x8);
    }
	
    MATRIX_3X3<T> Cofactor_Matrix_Transposed() const{
        T cofactor11=x[4]*x[8]-x[7]*x[5],cofactor12=x[7]*x[2]-x[1]*x[8],cofactor13=x[1]*x[5]-x[4]*x[2];
        T a11=cofactor11;
        T a21=cofactor12;
        T a31=cofactor13;
        T a12=(x[6]*x[5]-x[3]*x[8]);
        T a22=(x[0]*x[8]-x[6]*x[2]);
        T a32=(x[3]*x[2]-x[0]*x[5]);
        T a13=(x[3]*x[7]-x[6]*x[4]);
        T a23=(x[6]*x[1]-x[0]*x[7]);
        T a33=(x[0]*x[4]-x[3]*x[1]);
        T x0=a11;T x1=a12;T x2=a13;
        T x3=a21;T x4=a22;T x5=a23;
        T x6=a31;T x7=a32;T x8=a33;
        return MATRIX_3X3<T>(x0,x1,x2,x3,x4,x5,x6,x7,x8);
    }
	
    void Print(){
        std::cout<<"Matrix 3X3 = (" << x[0] << " , " << x[3] << " , " << x[6] << std::endl;
        std::cout<<"\t\t\t\t" << x[1] << " , " << x[4] << " , " << x[7] << std::endl;
        std::cout<<"\t\t\t\t" << x[2] << " , " << x[5] << " , " << x[8] << ")" << std::endl;
    }

    static MATRIX_3X3<T> Deformation_Gradient(const VECTOR_3D<T> X1,const VECTOR_3D<T> X2,const VECTOR_3D<T> X3,const VECTOR_3D<T> X4,const VECTOR_3D<T> x1,const VECTOR_3D<T> x2,const VECTOR_3D<T> x3,const VECTOR_3D<T> x4){
        MATRIX_3X3<T> Ds(x2(0)-x1(0),x2(1)-x1(1),x2(2)-x1(2),x3(0)-x1(0),x3(1)-x1(1),x3(2)-x1(2),x4(0)-x1(0),x4(1)-x1(1),x4(2)-x1(2));
        MATRIX_3X3<T> Dm(X2(0)-X1(0),X2(1)-X1(1),X2(2)-X1(2),X3(0)-X1(0),X3(1)-X1(1),X3(2)-X1(2),X4(0)-X1(0),X4(1)-X1(1),X4(2)-X1(2));
        T eps=1e-15;
        assert(Dm.Determinant()>=eps||Dm.Determinant()<=-eps);
        Dm.Invert();
        return Ds*Dm;
    }

    static MATRIX_3X3<T> Deformation_Gradient(const MATRIX_3X3<T> &Dm_inv,const VECTOR_3D<T> &x1,const VECTOR_3D<T> &x2,const VECTOR_3D<T> &x3,const VECTOR_3D<T> &x4){
        MATRIX_3X3<T> Ds(x2(0)-x1(0),x2(1)-x1(1),x2(2)-x1(2),x3(0)-x1(0),x3(1)-x1(1),x3(2)-x1(2),x4(0)-x1(0),x4(1)-x1(1),x4(2)-x1(2));
        return Ds*Dm_inv;
    }
};   


	
static MATRIX_3X3<double> operator+(const MATRIX_3X3<double>& A,const MATRIX_3X3<double>& B){
    return MATRIX_3X3<double>(A.x[0]+B.x[0],A.x[1]+B.x[1],A.x[2]+B.x[2],A.x[3]+B.x[3],A.x[4]+B.x[4],A.x[5]+B.x[5],A.x[6]+B.x[6],A.x[7]+B.x[7],A.x[8]+B.x[8]);}
static MATRIX_3X3<float> operator+(const MATRIX_3X3<float>& A,const MATRIX_3X3<float>& B){
    return MATRIX_3X3<float>(A.x[0]+B.x[0],A.x[1]+B.x[1],A.x[2]+B.x[2],A.x[3]+B.x[3],A.x[4]+B.x[4],A.x[5]+B.x[5],A.x[6]+B.x[6],A.x[7]+B.x[7],A.x[8]+B.x[8]);}
	
static MATRIX_3X3<double> operator-(const MATRIX_3X3<double>& A,const MATRIX_3X3<double>& B){
    return MATRIX_3X3<double>(A.x[0]-B.x[0],A.x[1]-B.x[1],A.x[2]-B.x[2],A.x[3]-B.x[3],A.x[4]-B.x[4],A.x[5]-B.x[5],A.x[6]-B.x[6],A.x[7]-B.x[7],A.x[8]-B.x[8]);}
static MATRIX_3X3<float> operator-(const MATRIX_3X3<float>& A,const MATRIX_3X3<float>& B){
    return MATRIX_3X3<float>(A.x[0]-B.x[0],A.x[1]-B.x[1],A.x[2]-B.x[2],A.x[3]-B.x[3],A.x[4]-B.x[4],A.x[5]-B.x[5],A.x[6]-B.x[6],A.x[7]-B.x[7],A.x[8]-B.x[8]);}	
	
static MATRIX_3X3<double> operator*(const MATRIX_3X3<double>& A,const MATRIX_3X3<double>& B){
    return MATRIX_3X3<double>(A.x[0]*B.x[0]+A.x[3]*B.x[1]+A.x[6]*B.x[2],A.x[1]*B.x[0]+A.x[4]*B.x[1]+A.x[7]*B.x[2],A.x[2]*B.x[0]+A.x[5]*B.x[1]+A.x[8]*B.x[2],
        A.x[0]*B.x[3]+A.x[3]*B.x[4]+A.x[6]*B.x[5],A.x[1]*B.x[3]+A.x[4]*B.x[4]+A.x[7]*B.x[5],A.x[2]*B.x[3]+A.x[5]*B.x[4]+A.x[8]*B.x[5],
        A.x[0]*B.x[6]+A.x[3]*B.x[7]+A.x[6]*B.x[8],A.x[1]*B.x[6]+A.x[4]*B.x[7]+A.x[7]*B.x[8],A.x[2]*B.x[6]+A.x[5]*B.x[7]+A.x[8]*B.x[8]);}
static MATRIX_3X3<float> operator*(const MATRIX_3X3<float>& A,const MATRIX_3X3<float>& B){
    return MATRIX_3X3<float>(A.x[0]*B.x[0]+A.x[3]*B.x[1]+A.x[6]*B.x[2],A.x[1]*B.x[0]+A.x[4]*B.x[1]+A.x[7]*B.x[2],A.x[2]*B.x[0]+A.x[5]*B.x[1]+A.x[8]*B.x[2],
        A.x[0]*B.x[3]+A.x[3]*B.x[4]+A.x[6]*B.x[5],A.x[1]*B.x[3]+A.x[4]*B.x[4]+A.x[7]*B.x[5],A.x[2]*B.x[3]+A.x[5]*B.x[4]+A.x[8]*B.x[5],
        A.x[0]*B.x[6]+A.x[3]*B.x[7]+A.x[6]*B.x[8],A.x[1]*B.x[6]+A.x[4]*B.x[7]+A.x[7]*B.x[8],A.x[2]*B.x[6]+A.x[5]*B.x[7]+A.x[8]*B.x[8]);}
	
static MATRIX_3X3<double> operator*(const double scale,const MATRIX_3X3<double>& input){return MATRIX_3X3<double>(scale*input(0,0),scale*input(1,0),scale*input(2,0),scale*input(0,1),scale*input(1,1),scale*input(2,1),scale*input(0,2),scale*input(1,2),scale*input(2,2));}
static MATRIX_3X3<float> operator*(const float scale,const MATRIX_3X3<float>& input){return MATRIX_3X3<float>(scale*input(0,0),scale*input(1,0),scale*input(2,0),scale*input(0,1),scale*input(1,1),scale*input(2,1),scale*input(0,2),scale*input(1,2),scale*input(2,2));}	

template<class T>
class MATRIX_3X2{
	
public:	
	
    T x[6];
	
    MATRIX_3X2(const T input){
        for(int i=0;i<6;i++) x[i]=input;
    }
	
    MATRIX_3X2(){
        for(int i=0;i<6;i++) x[i]=(T)0;
    }
	
    MATRIX_3X2(const T x11,const T x21,const T x31,const T x12,const T x22,const T x32){
        x[0]=x11;x[3]=x12;
        x[1]=x21;x[4]=x22;
        x[2]=x31;x[5]=x32;
    }
	
    MATRIX_3X2(const VECTOR_3D<T>& x1,const VECTOR_3D<T>& x2){
        x[0]=x1.x();x[3]=x2.x();
        x[1]=x1.y();x[4]=x2.y();
        x[2]=x1.z();x[5]=x2.z();
    }
	
    VECTOR_3D<T> Column(const int column) const{
        if(column==0) return VECTOR_3D<T>(x[0],x[1],x[2]);
        else if(column==1) return VECTOR_3D<T>(x[3],x[4],x[5]);
        else{
            assert(false);
            return VECTOR_3D<T>(0,0,0);}
    }
	
    static T Contract(const MATRIX_3X2<T>& A,const MATRIX_3X2<T>& B){
        return A.x[0]*B.x[0]+A.x[1]*B.x[1]+A.x[2]*B.x[2]+A.x[3]*B.x[3]+A.x[4]*B.x[4]+A.x[5]*B.x[5];}
	
    T Frobenius_Norm_Squared() const{
        return x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+x[3]*x[3]+x[4]*x[4]+x[5]*x[5];
    }
	
    T Determinant_ABS() const{
        T a=x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
        T b=x[0]*x[3]+x[1]*x[4]+x[2]*x[5];
        T c=x[3]*x[3]+x[4]*x[4]+x[5]*x[5];
		
        return sqrt(MATRIX_2X2<T>(a,b,b,c).Determinant());
    }
	
    T operator()(const int i,const int j) const {
        int index=3*j+i;
        assert(index>=0);
        assert(index<6);
        return x[index];
    }
	
    T& operator()(const int i,const int j){
        int index=3*j+i;
        assert(index>=0);
        assert(index<6);
        return x[index];
    }
	
    T& operator()(const int i){
        assert(i>=0);
        assert(i<6);
        return x[i];
    }
	
    T operator()(const int i) const {
        assert(i>=0);
        assert(i<6);
        return x[i];
    }
	
    VECTOR_3D<T> operator*(const VECTOR_2D<T>& x_in) const{
        return VECTOR_3D<T>(x[0]*x_in(0)+x[3]*x_in(1),
            x[1]*x_in(0)+x[4]*x_in(1),
            x[2]*x_in(0)+x[5]*x_in(1));
    }
	
    MATRIX_3X2<T> operator+(const MATRIX_3X2<T>& A_input){
        return MATRIX_3X2(x[0]+A_input(0),x[1]+A_input(1),x[2]+A_input(2),x[3]+A_input(3),x[4]+A_input(4),x[5]+A_input(5));
    }
	
    MATRIX_3X2<T> operator-(const MATRIX_3X2<T>& A_input){
        return MATRIX_3X2(x[0]-A_input(0),x[1]-A_input(1),x[2]-A_input(2),x[3]-A_input(3),x[4]-A_input(4),x[5]-A_input(5));
    }
	
};

static MATRIX_3X2<double> operator*(const MATRIX_3X2<double>& A,const MATRIX_2X2<double>& B){
    VECTOR_3D<double> c1=A.Column(0);VECTOR_3D<double> c2=A.Column(1);
    double a11=B.a11;double a21=B.a21;double a12=B.a12;double a22=B.a22;
    return MATRIX_3X2<double>(a11*c1+a21*c2,a12*c1+a22*c2);}
static MATRIX_3X2<float> operator*(const MATRIX_3X2<float>& A,const MATRIX_2X2<float>& B){
    VECTOR_3D<float> c1=A.Column(0);VECTOR_3D<float> c2=A.Column(1);
    float a11=B.a11;float a21=B.a21;float a12=B.a12;float a22=B.a22;
    return MATRIX_3X2<float>(a11*c1+a21*c2,a12*c1+a22*c2);}	
static MATRIX_3X2<double> operator+(const MATRIX_3X2<double>& A,const MATRIX_3X2<double>& B){
    return MATRIX_3X2<double>(A.x[0]+B.x[0],A.x[1]+B.x[1],A.x[2]+B.x[2],A.x[3]+B.x[3],A.x[4]+B.x[4],A.x[5]+B.x[5]);}
static MATRIX_3X2<float> operator+(const MATRIX_3X2<float>& A,const MATRIX_3X2<float>& B){
    return MATRIX_3X2<float>(A.x[0]+B.x[0],A.x[1]+B.x[1],A.x[2]+B.x[2],A.x[3]+B.x[3],A.x[4]+B.x[4],A.x[5]+B.x[5]);}	
	
static MATRIX_3X2<double> operator*(const double scale,const MATRIX_3X2<double>& input){return MATRIX_3X2<double>(scale*input(0,0),scale*input(1,0),scale*input(2,0),scale*input(0,1),scale*input(1,1),scale*input(2,1));}
static MATRIX_3X2<float> operator*(const float scale,const MATRIX_3X2<float>& input){return MATRIX_3X2<float>(scale*input(0,0),scale*input(1,0),scale*input(2,0),scale*input(0,1),scale*input(1,1),scale*input(2,1));}	
	
template<class T>
class LIST{
    //This is a dynamically resizable VECTOR. It just keeps a buffer and then resizes if more info than available in the buffer is requested
	
    int n,n_extended;//n is the apparent size of the array and n_extended is the actual size of the array
    int buffer_size;//Whenever the array needs to be resized becasue we have passed n_extended, we add the buffer_size more entries than would be minimally necessary
    T* values;
public:
    LIST():n(0),buffer_size(10){
        n_extended=n+buffer_size;
        values=new T[n_extended];
        for(int i=0;i<n_extended;i++) values[i]=T();}
	
    LIST(const int n_input):n(n_input),buffer_size(10){
        n_extended=n+buffer_size;
        values=new T[n_extended];
        for(int i=0;i<n_extended;i++) values[i]=T();}
	
    ~LIST() {delete[] values;}
	
    int Size(){return n;}
	
    T& operator()(const int i){assert(i<n); return values[i];}
	
    void Resize(const int n_input){
        if(n_input>n_extended){
            delete values;
            n_extended=n_input+buffer_size;
            values=new T[n_extended];
            for(int i=0;i<n_extended;i++) values[i]=T();}
        else{
            n=n_input;}}
	
    void Append_Unique(const T& entry){
        for(int i=0;i<n;i++)
            if(values[i]==entry) return;
		
        if(n<n_extended){
            n=n+1;
            values[n-1]=entry;}
        else{
            //delete values;
            n=n+1;
            n_extended=n+buffer_size;
            T* values_new=new T[n_extended];
            for(int i=0;i<n_extended;i++) values_new[i]=T();
            for(int i=0;i<n-1;i++) values_new[i]=values[i];
            values_new[n-1]=entry;
            delete values;
            values=values_new;}
    }
	
    void Append_Element(const T& entry){
        if(n<n_extended){
            n=n+1;
            values[n-1]=entry;}
        else{
            //delete values;
            n=n+1;
            n_extended=n+buffer_size;
            T* values_new=new T[n_extended];
            for(int i=0;i<n_extended;i++) values_new[i]=T();
            for(int i=0;i<n-1;i++) values_new[i]=values[i];
            values_new[n-1]=entry;
            delete values;
            values=values_new;}
    }
};
	
template<class T>
class VECTOR{
    int n;
    T* values;
public:
    VECTOR():n(1) {
        values=new T[n];
        for(int i=0;i<=n-1;i++) values[i]=T();}

    VECTOR(const int n_input):n(n_input) {
        values=new T[n];
        for(int i=0;i<=n-1;i++) values[i]=T();}

    ~VECTOR() {delete[] values;}

    bool Resize(const int n_new){
        if (n_new==n) return true;
        assert(n_new>0);
        delete[] values;
        values=new T[n_new];
        n=n_new;
        for(int i=0;i<=n-1;i++) values[i]=T();
        return (values!=NULL);}

    T& operator()(INDEX_2D& index){
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
	
    T Min(){
        T min=FLT_MAX;
        for(int i=0;i<n;i++) if(values[i]<min) min=values[i];
        return min;
    }
	
    T Max(){
        T max=-FLT_MAX;
        for(int i=0;i<n;i++) if(values[i]>max) max=values[i];
        return max;
    }
	
    T& operator()(const int i) {assert(0<=i && i<=n-1);return values[i];}
	
    T operator()(const int i) const {assert(0<=i && i<=n-1);return values[i];}
	
    void Set_To_Zero(){for(int i=0;i<n;i++) values[i]=0;}
	
    T L_inf(){
        T max_norm=(T)0;
        for(int i=0;i<n;i++) if(fabs(values[i])>max_norm) max_norm=fabs(values[i]);
        return max_norm;}
	
    T Sum(){T sum=(T)0;for(int i=0;i<n;i++) sum+=values[i];return sum;}
	
    int Size() const {return n;}
	
    T Magnitude(){
        T sum=(T)0;
        for(int i=0;i<n;i++) sum+=(values[i]*values[i]);
        return sqrt(sum);}
	
    T Normalize(){
        T sum=(T)0;
        for(int i=0;i<n;i++) sum+=(values[i]*values[i]);
        T norm=sqrt(sum);
        for(int i=0;i<n;i++) values[i]=values[i]/norm;
        return norm;}
	
    void Enforce_Zero_Sum(){
        T sum=(T)0;
        for(int i=0;i<n;i++) sum+=values[i];
        for(int i=0;i<n;i++) values[i]-=sum/((T)n);}
	
    void operator+=(VECTOR<T>& x) {assert(n==x.Size());for(int i=0;i<n;i++) values[i]=values[i]+x(i);}
	
    void operator-=(VECTOR<T>& x) {assert(n==x.Size());for(int i=0;i<n;i++) values[i]=values[i]-x(i);}
	
    // void Write_DAT_File(std::string file){
    // 	FILE* fpointer;
    // 	fpointer=fopen(file.c_str(),"w");
    // 	for(int i=0;i<n;i++)
    // 		fprintf(fpointer,"%g\n",(double)values[i]);
    // 	fclose(fpointer);}
	
    void Print(){
        std::cout << "Vector =";
        for(int i=0;i<n;i++) std::cout << " " << values[i] << " , ";
        std::cout << "." << std::endl;}
};

template <class T>
std::ostream& operator<<(std::ostream & os, const VECTOR<T>& v){
    os<<"[";for(int i=0;i<v.Size();i++) os<<v(i)<<" ";os<<"]";
    return os;
}
template <class T>
std::ostream& operator<<(std::ostream & os, const VECTOR_2D<T>& v){
    os<<"["<<v.x()<<" "<<v.y()<<"]";
    return os;
}

template<class T>
class SPARSE_ROW{
    const int n;
    int size;
    int* indices;
    T* values;
public:
    SPARSE_ROW(const int n_input):n(n_input),size(0),indices(0),values(0) {}

    ~SPARSE_ROW() {delete[] indices;delete[] values;}
	
    void Permute_Columns(const VECTOR<int>& permutation){
        assert(permutation.Size()==n);
        for(int i=0;i<size;i++) indices[i]=permutation(indices[i]);
    }
	
    void Zero_Out_Without_Changing_Sparsity(){for(int i=0;i<size;i++) values[i]=(T)0;}
	
    bool Is_Non_Zero(const int index){for(int i=0;i<size;i++) if(indices[i]==index) return true;return false;}
	
    void operator+=(SPARSE_ROW<T>& input){
        for(int i=0;i<input.Number_Nonzero();i++){
            if(this->Value_Exists_At_Entry(input.Index(i)))
                (*this)(input.Index(i))+=input.Value_At_Sparse_Index(i);
            else
                this->Add_Entry(input.Index(i),input.Value_At_Sparse_Index(i));}
    }

    T& operator()(const int i){
        assert(0<=i && i<=n-1);
        for(int j=0;j<=size-1;j++) if(indices[j]==i) return values[j];
        assert(false);
        return values[0];}
	
    T Row_Sum(){
        T sum=(T)0;
        for(int i=0;i<=size-1;i++) sum+=values[i];
        return sum;}
	
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
	
    T Value_At_Sparse_Index(const int i_hat) const{assert(i_hat<size);return values[i_hat];}
    T& Value_At_Sparse_Index(const int i_hat){assert(i_hat<size);return values[i_hat];}
	
    T Dot_Product(VECTOR<T>& v){
        assert(v.Size()==n);
        T result=0;for(int i=0;i<=size-1;i++) result+=values[i]*v(indices[i]);
        return result;}

    void Add_Entry(const int index,const T value){
        bool found=false;int entry=0;
        for(int i=0;i<size;i++) if(indices[i]==index){found=true;entry=i;}
        if(found){
            values[entry]=value;
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

    void Zero_Out_Without_Changing_Sparsity(){
        for(int i=0;i<m;i++) rows[i]->Zero_Out_Without_Changing_Sparsity();}
	
    T& operator()(const int i,const int j){
        assert(0<=i && i<=m-1);return (*rows[i])(j);}
	
    void Column(const int j,VECTOR<T>& c){
        assert(j>=0 && j<n && c.Size()==m);
        for(int i=0;i<m;i++){
            SPARSE_ROW<T>& Ai=Row(i);
            if(Ai.Is_Non_Zero(j)) c(i)=Ai(j);
            else c(i)=(T)0;}}
	
    T Column_Sum(const int j){
        assert(j>=0 && j<n);
        T sum=(T)0;
        for(int i=0;i<m;i++){
            SPARSE_ROW<T>& Ai=Row(i);
            if(Ai.Is_Non_Zero(j)) sum+=Ai(j);}
        return sum;
    }
	
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
	
    void Scale_Rows(const VECTOR<T>& diagonal_row_scaling){
        for(int i=0;i<=m-1;i++) rows[i]->Scale(diagonal_row_scaling(i));}
	
    void Scale_Columns(const VECTOR<T>& diagonal_column_scaling){
        for(int i=0;i<=m-1;i++){
            SPARSE_ROW<T>& row=*(rows[i]);
            for(int j_index=0;j_index<row.Number_Nonzero();j_index++){
                int j=row.Index(j_index);
                row.Value_At_Sparse_Index(j_index)*=diagonal_column_scaling(j);}}
    }
	
    void Print_Sparsity_Information(){
        std::cout << "Sparse Matrix: " << std::endl;
        for(int i=0;i<m;i++){
            std::cout << "Number non-zero in row " << i << " = " << Row(i).Number_Nonzero()<<std::endl;
        }
    }
	
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
	
    void Multiply_With_Transpose(VECTOR<T>& x,VECTOR<T>& b){//as in b=A'x
        assert(b.Size()==n);
        assert(x.Size()==m);
        b.Set_To_Zero();
        for(int j=0;j<m;j++){
            SPARSE_ROW<T>& row=Row(j);
            for(int i=0;i<row.Number_Nonzero();i++){
                int index=row.Index(i);
                b(index)+=row.Value_At_Sparse_Index(i)*x(j);}}
    }
	
    void Multiply_With_Transpose(const VECTOR<T>& x,VECTOR<T>& b){//as in b=A'x
        assert(b.Size()==n);
        assert(x.Size()==m);
        b.Set_To_Zero();
        for(int j=0;j<m;j++){
            SPARSE_ROW<T>& row=Row(j);
            for(int i=0;i<row.Number_Nonzero();i++){
                int index=row.Index(i);
                b(index)+=row.Value_At_Sparse_Index(i)*x(j);}}
    }
	
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
	
    // void Write_DAT_File(std::string file,const bool sparse=false){
    // 	FILE* fpointer;
    // 	fpointer=fopen(file.c_str(),"w");
    // 	if(!sparse){
    // 		for(int i=0;i<m;i++){
    // 			SPARSE_ROW<T>& row=Row(i);
    // 			for(int j=0;j<n;j++){
    // 				bool found=false;
    // 				if(row.Value_Exists_At_Entry(j)){
    // 						fprintf(fpointer,"%g ",(*this)(i,j));
    // 					found=true;}
    // 				if(!found) fprintf(fpointer,"%g ",(T)0);}
    // 			fprintf(fpointer,"\n");}}
    // 	else{
    // 		for(int i=0;i<m;i++){
    // 			SPARSE_ROW<T>& row=Row(i);
    // 			for(int j=0;j<row.Number_Nonzero();j++){
    // 				int column=row.Index(j);
    // 				fprintf(fpointer,"%d %d %g \n",i+1,column+1,(*this)(i,column));}}}
    // 	fclose(fpointer);
    // }
	
};
	
static void Sparse_Multiply(SPARSE_MATRIX<double>& C,SPARSE_MATRIX<double>& A,SPARSE_MATRIX<double>& B){
    //C=A*B;
    assert(A.N()==B.M());
    assert(C.M()==A.M());
    assert(C.N()==B.N());
	
    for(int i=0;i<C.M();i++){
        SPARSE_ROW<double>& a_row=A.Row(i);
        SPARSE_ROW<double>& c_row=C.Row(i);
        for(int j=0;j<a_row.Number_Nonzero();j++){
            double aij=a_row.Value_At_Sparse_Index(j);
            int b_row_index=a_row.Index(j);
            SPARSE_ROW<double>& b_row=B.Row(b_row_index);
            for(int k=0;k<b_row.Number_Nonzero();k++){
                int column_index=b_row.Index(k);
                double product=aij*b_row.Value_At_Sparse_Index(k);
                if(!c_row.Value_Exists_At_Entry(column_index)){
                    c_row.Add_Entry(column_index,product);}
                else c_row(column_index)+=product;}}}
}
	
static void Sparse_Multiply(SPARSE_MATRIX<float>& C,SPARSE_MATRIX<float>& A,SPARSE_MATRIX<float>& B){
    //C=A*B;
    assert(A.N()==B.M());
    assert(C.M()==A.M());
    assert(C.N()==B.N());
	
    for(int i=0;i<C.M();i++){
        SPARSE_ROW<float>& a_row=A.Row(i);
        SPARSE_ROW<float>& c_row=C.Row(i);
        for(int j=0;j<a_row.Number_Nonzero();j++){
            float aij=a_row.Value_At_Sparse_Index(j);
            int b_row_index=a_row.Index(j);
            SPARSE_ROW<float>& b_row=B.Row(b_row_index);
            for(int k=0;k<b_row.Number_Nonzero();k++){
                int column_index=b_row.Index(k);
                float product=aij*b_row.Value_At_Sparse_Index(k);
                if(!c_row.Value_Exists_At_Entry(column_index)){
                    c_row.Add_Entry(column_index,product);}
                else c_row(column_index)+=product;}}}
}		
	
template<class T>
class MATRIX_MXN{
    const int m,n;
    VECTOR<T>** rows;
public:
    MATRIX_MXN(const int m_input,const int n_input):m(m_input),n(n_input){
        rows=new VECTOR<T>*[m];
        for(int i=0;i<=m-1;i++) rows[i]=new VECTOR<T>(n);}
	
    ~MATRIX_MXN(){
        for(int i=0;i<=m-1;i++) delete rows[i];
        delete[] rows;
    }
	
    int M(){return m;}
    int N(){return n;}
	
    void Set_To_Zero(){
        for(int i=0;i<m;i++)for(int j=0;j<n;j++) (*this)(i,j)=(T)0;
    }
	
    T& operator()(const int i,const int j){
        assert(0<=i && i<=m-1);return (*rows[i])(j);}
	
    void Transpose(MATRIX_MXN<T>& result){
        assert(result.M()==n && result.N()==m);
        for(int i=0;i<result.M();i++){
            for(int j=0;j<result.N();j++){
                result(i,j)=(*this)(j,i);}}
    }
	
    void Print(){
        for(int i=0;i<m;i++){
            for(int j=0;j<n;j++)
                std::cout<<(*this)(i,j)<<" ";
            std::cout<<std::endl;}
    }
    
};
	
static void Multiply(MATRIX_MXN<double>& A, MATRIX_MXN<double>& B, MATRIX_MXN<double>& result){
    assert(A.N()==B.M());
    assert(result.M()==A.M());
    assert(result.N()==B.N());
    for(int i=0;i<result.M();i++){
        for(int j=0;j<result.N();j++){
            result(i,j)=(double)0;
            for(int k=0;k<A.N();k++){
                result(i,j)+=A(i,k)*B(k,j);}}}
}	

static void Multiply(MATRIX_MXN<float>& A, MATRIX_MXN<float>& B, MATRIX_MXN<float>& result){
    assert(A.N()==B.M());
    assert(result.M()==A.M());
    assert(result.N()==B.N());
    for(int i=0;i<result.M();i++){
        for(int j=0;j<result.N();j++){
            result(i,j)=(float)0;
            for(int k=0;k<A.N();k++){
                result(i,j)+=A(i,k)*B(k,j);}}}
}
	
template<class T>
class MINRES{
    //TODO: Change this so that the initial guess can be non-zero.
	
    SPARSE_MATRIX<T>& A;
    VECTOR<T>& x;
    VECTOR<T>& b;
    int n;//The number of unknowns
    MATRIX_2X2<T> Gk,Gkm1,Gkm2;//The Q in the QR of Hk: Givens rotations
    int max_iterations;
    T gamma,delta,epsilon;//These are the newest entries in R from the QR of Hk
    T beta_kp1,alpha_k,beta_k;//This is the last column in Hk
    VECTOR<T> mk,mkm1,mkm2;//mk is the newest search direction in the memory friendly basis for the Krylov space, you need the other two to generate mk
    VECTOR<T> qkp1,qk,qkm1;//These are the two Lanczos vectors needed at each iteration 	
    T tk;//This is the step length in the direction of mk
    VECTOR_2D<T> last_two_components_of_givens_transformed_least_squares_rhs;//This will track the residual with just a constant number of flops per iteration
    int current_iteration;
    T tolerance;
    VECTOR<int>* dirichlet_dofs;
	
public:
    MINRES(SPARSE_MATRIX<T>& A_input,VECTOR<T>& x_input,VECTOR<T>& b_input,const int max_it_input):A(A_input),
                                                                                                  x(x_input),b(b_input),max_iterations(max_it_input),n(b.Size()),mk(n),mkm1(n),mkm2(n),
                                                                                                  qkp1(n),qk(n),qkm1(n),dirichlet_dofs(0){Set_Tolerance();}
	
    ~MINRES(){}
	
    void Set_Dirichlet_Dofs(VECTOR<int>& dirichlet){dirichlet_dofs=&dirichlet;}
	
    void Set_Tolerance(const T tolerance_input=(T)1e-10){tolerance=tolerance_input;}
	
    void Initialize(){
        //This is the zero initial guess version of MINRES
        if(!dirichlet_dofs) x.Set_To_Zero();
        else{
            VECTOR<T> dirichlet_values(dirichlet_dofs->Size());
            for(int i=0;i<dirichlet_dofs->Size();i++) dirichlet_values(i)=x((*dirichlet_dofs)(i));
            x.Set_To_Zero();
            for(int i=0;i<dirichlet_dofs->Size();i++) x((*dirichlet_dofs)(i))=dirichlet_values(i);}
        mk.Set_To_Zero();mkm1.Set_To_Zero();mkm2.Set_To_Zero();
        qkm1.Set_To_Zero();qk.Set_To_Zero();qkp1.Set_To_Zero();
        gamma=(T)0;delta=(T)0;epsilon=(T)0;
        beta_kp1=(T)0;alpha_k=(T)0;beta_k=(T)0;
        tk=(T)0;
		
        Gk=MATRIX_2X2<T>::Identity();Gkm1=MATRIX_2X2<T>::Identity();Gkm2=MATRIX_2X2<T>::Identity();
    }
	
    void Set_Lanczos_Dirichlet_Indices_To_Zero(){
        assert(dirichlet_dofs);
        for(int i=0;i<dirichlet_dofs->Size();i++){
            int index=(*dirichlet_dofs)(i);
            qkp1(index)=(T)0;}
    }
	
    int Solve(const bool verbose=false){
        Initialize();
		
        A.Residual(b,x,qkp1);//qkp1=b-Ax: Initial Krylov vector is the residual
        if(dirichlet_dofs) Set_Lanczos_Dirichlet_Indices_To_Zero();
        T two_norm_of_residual=qkp1.Normalize();//This is the zero initial guess version.
        last_two_components_of_givens_transformed_least_squares_rhs=VECTOR_2D<T>(two_norm_of_residual,(T)0);
		
        for(int k=0;k<max_iterations;k++){
            current_iteration=k;//Keep track of the iteration number.
            if(two_norm_of_residual<tolerance) return k;//Output the number of iterations.
			
            qkm1=qk;qk=qkp1;//Save the last two Lanczos vectors
            A.Multiply(qk,qkp1);//Get the important part of the next Lanczos vector: q_k+1
            if(dirichlet_dofs) Set_Lanczos_Dirichlet_Indices_To_Zero();
            alpha_k=qkp1.Dot(qk);
            for(int j=0;j<n;j++) qkp1(j)-=alpha_k*qk(j);
            beta_k=beta_kp1;
            for(int j=0;j<n;j++) qkp1(j)-=beta_k*qkm1(j);
            beta_kp1=qkp1.Normalize();//qkp1 is now complete
			
            two_norm_of_residual=Apply_All_Previous_Givens_Rotations_And_Determine_New_Givens();//This determines the newest Givens rotation and applies the previous two where appropriate
			
            if(verbose){
                std::cout<<"Iteration = " << current_iteration << std::endl;
                std::cout<<"\t Newest Lanczos vector: ";qkp1.Print();
                std::cout<<"\t Two norm of the residual = " << two_norm_of_residual << std::endl;}
			
            mkm2=mkm1;mkm1=mk;//Get ready to get the newest m
            for(int j=0;j<n;j++) mk(j)=((T)1/gamma)*(qk(j)-delta*mkm1(j)-epsilon*mkm2(j));//Three term recurence for the m's
            for(int j=0;j<n;j++) x(j)+=tk*mk(j);}
		
        return max_iterations;
    }
	
    T Apply_All_Previous_Givens_Rotations_And_Determine_New_Givens(){
		
        //QR the LHS: gamma, delta, epsilon
        Gkm2=Gkm1;Gkm1=Gk;		
        VECTOR_2D<T> epsilon_k_and_phi_k=Gkm2*VECTOR_2D<T>((T)0,beta_k);
        epsilon=epsilon_k_and_phi_k.x();
        VECTOR_2D<T> delta_k_and_zsi_k=Gkm1*VECTOR_2D<T>(epsilon_k_and_phi_k.y(),alpha_k);
        delta=delta_k_and_zsi_k.x();
        VECTOR_2D<T> temp(delta_k_and_zsi_k.y(),beta_kp1);
        Gk=MATRIX_2X2<T>::Givens_Rotation(temp);
        VECTOR_2D<T> temp2=Gk*temp;
        gamma=temp2.x();
		
        //Now deal with the RHS: tk and residual (two norm)
        last_two_components_of_givens_transformed_least_squares_rhs=Gk*last_two_components_of_givens_transformed_least_squares_rhs;
        tk=last_two_components_of_givens_transformed_least_squares_rhs.x();
        T residual=last_two_components_of_givens_transformed_least_squares_rhs.y();//This is the two norm of the residual.
        last_two_components_of_givens_transformed_least_squares_rhs=VECTOR_2D<T>(residual,(T)0);//Set up for the next iteration
        if(residual<(T)0) return -residual;
        else return residual;
    }
	
};
	
template<class T>
class GMRES{
    SPARSE_MATRIX<T>& A;
    VECTOR<T>& x;
    VECTOR<T>& b;
    VECTOR<T> lambda,givens_transformed_least_squares_rhs;	
    VECTOR<VECTOR<T>*> arnoldi_basis_vectors;//Orthonormal basis vectors for the Krylov space
    MATRIX_MXN<T> R;//Upper triangular part of the QR'd upper Hessenberg matrix arising from the Arnoldi conjugated A
    VECTOR<MATRIX_2X2<T> > givens_rotations;//Used in the QR of H
    int current_iteration,max_iterations;
    T tolerance;
    int number_dofs;
	
public:
    GMRES(SPARSE_MATRIX<T>& A_input,VECTOR<T>& x_input,VECTOR<T>& b_input,const int max_it_input):A(A_input),
                                                                                                 x(x_input),b(b_input),max_iterations(max_it_input),lambda(max_it_input),arnoldi_basis_vectors(max_it_input+1),R(max_it_input+1,max_it_input),
                                                                                                 givens_rotations(max_it_input),givens_transformed_least_squares_rhs(max_it_input+1),number_dofs(x_input.Size()){
        for(int i=0;i<max_iterations;i++) arnoldi_basis_vectors(i)=new VECTOR<T>(number_dofs);
        lambda.Set_To_Zero();
        givens_transformed_least_squares_rhs.Set_To_Zero();
        x.Set_To_Zero();//initial guess is zero
        Set_Tolerance();}
	
    void Set_Tolerance(const T tolerance_input=(T)1e-10){tolerance=tolerance_input;}
	
    ~GMRES(){for(int i=0;i<max_iterations;i++) delete arnoldi_basis_vectors(i);}
	
    MATRIX_2X2<T>& Givens(int i){return givens_rotations(i);}//Gi operates on the last two entries of a i+1 length column vector
    VECTOR<T>& Q(const int i){return *(arnoldi_basis_vectors(i));}
	
    int Solve(const bool verbose=false){
        Q(0)=b;T two_norm_of_residual=Q(0).Normalize();//This is the zero initial guess version.
        givens_transformed_least_squares_rhs(0)=two_norm_of_residual;//This is the RHS in the least squares minimization of the residual.
		
        for(int k=0;k<max_iterations;k++){
            current_iteration=k;//Keep track of the iteration number.
            if(two_norm_of_residual<tolerance){
                Solve_For_Lambda_And_Assemble_X();//We only do this on the last iteration.
                return k;}//Output the number of iterations.
			
            A.Multiply(Q(k),Q(k+1));//Get the important part of the next Arnoldi vector: q_k+1
            for(int i=0;i<k+1;i++){//Now do the Graham-Schmidt like orthogonalization of q_k+1
                T hik=Q(k+1).Dot(Q(i));
                Add_Entry_To_Last_Column_In_Upper_Hessenberg_H(i,hik);//Keep track of these components in H.
                for(int j=0;j<number_dofs;j++) Q(k+1)(j)-=hik*Q(i)(j);}
			
            T hkp1k=Q(k+1).Normalize();//This is the newest entry in the subdiagonal of the upper Hessenberg H.
            Add_Entry_To_Last_Column_In_Upper_Hessenberg_H(k+1,hkp1k);			
            Apply_All_Previous_Givens_Rotations_To_Last_Column_Of_H();//Update the last column to be consistent with previous Givens applicaitons.
            Determine_And_Apply_New_Givens_Rotation_For_Last_Column_Of_H();//Now get the lastest Givens to turn H into uppertriangular.
            two_norm_of_residual=Apply_Newest_Givens_Rotation_To_Least_Squares_RHS();//Keep track of the Givens effect on the RHS, also this yields the residual.
		
            if(verbose) Print();}
		
        return max_iterations;
    }
	
    T Apply_Newest_Givens_Rotation_To_Least_Squares_RHS(){
        MATRIX_2X2<T>& Gi=givens_rotations(current_iteration);
        VECTOR_2D<T> last_two_components(givens_transformed_least_squares_rhs(current_iteration),givens_transformed_least_squares_rhs(current_iteration+1));
        VECTOR_2D<T> transformed_last_two_components=Gi*last_two_components;
        givens_transformed_least_squares_rhs(current_iteration)=transformed_last_two_components.x();
        givens_transformed_least_squares_rhs(current_iteration+1)=transformed_last_two_components.y();
        if(givens_transformed_least_squares_rhs(current_iteration+1)<0) return -givens_transformed_least_squares_rhs(current_iteration+1);
        else return givens_transformed_least_squares_rhs(current_iteration+1);
    }
	
    void Determine_And_Apply_New_Givens_Rotation_For_Last_Column_Of_H(){
        VECTOR_2D<T> hk_to_kp1(R(current_iteration,current_iteration),R(current_iteration+1,current_iteration));
        givens_rotations(current_iteration)=MATRIX_2X2<T>::Givens_Rotation(hk_to_kp1);
        VECTOR_2D<T> x=givens_rotations(current_iteration)*hk_to_kp1;
        R(current_iteration,current_iteration)=x.x();
        R(current_iteration+1,current_iteration)=x.y();
    }
	
    void Print(){
        std::cout<< "Current iteration = " << current_iteration << std::endl;
        std::cout<< "The Arnoldi basis vectors:"<<std::endl;
        for(int i=0;i<current_iteration+1;i++){
            std::cout<<"Writing vector: " << i <<std::endl;
            Q(i).Print();}
        std::cout<<"Arnoldi basis orthogonality check: " << std::endl;
        for(int i=0;i<current_iteration+1;i++){
            for(int j=0;j<current_iteration+1;j++){
                std::cout<< "(qi,qj) = " << Q(i).Dot(Q(j)) << " , ";}
            std::cout<<std::endl;}
        std::cout<< "Upper triangular part of H = " << std::endl;
        R.Print();
        std::cout<< "The RHS for the upper triangular solve:"<<std::endl;
        givens_transformed_least_squares_rhs.Print();
    }
	
    void Apply_All_Previous_Givens_Rotations_To_Last_Column_Of_H(){
        for(int i=0;i<current_iteration;i++){
            MATRIX_2X2<T>& Gi=givens_rotations(i);
            VECTOR_2D<T> hi_to_ip1=Gi*VECTOR_2D<T>(R(i,current_iteration),R(i+1,current_iteration));
            R(i,current_iteration)=hi_to_ip1.x();
            R(i+1,current_iteration)=hi_to_ip1.y();}
    }
	
    void Add_Entry_To_Last_Column_In_Upper_Hessenberg_H(const int row_number,const T value){
        R(row_number,current_iteration)=value;
    }
	
    void Solve_For_Lambda_And_Assemble_X(){
        x.Set_To_Zero();
        for(int i=current_iteration-1;i>=0;i--){
            lambda(i)=givens_transformed_least_squares_rhs(i);
            for(int j=i+1;j<current_iteration;j++){
                lambda(i)-=R(i,j)*lambda(j);}
            lambda(i)/=R(i,i);
            for(int k=0;k<x.Size();k++) x(k)+=lambda(i)*Q(i)(k);}
    }
};
	
template<class T>
class PRECONDITIONED_CONJUGATE_GRADIENT{
    SPARSE_MATRIX<T>& A;
    VECTOR<T>& x;
    VECTOR<T>& b;
    VECTOR<T> z,r,p,q,d,temp;
    VECTOR<int>* dirichlet_dofs;
    int max_iterations;
    T tolerance,alpha,beta;
    bool use_jacobi_preconditioner;
		
public:
    PRECONDITIONED_CONJUGATE_GRADIENT(SPARSE_MATRIX<T>& A_input,VECTOR<T>& x_input,VECTOR<T>& b_input,const int max_it):A(A_input),
                                                                                                                       x(x_input),b(b_input),r(x_input.Size()),max_iterations(max_it),
                                                                                                                       p(x_input.Size()),q(x_input.Size()),d(x_input.Size()),temp(x_input.Size()),dirichlet_dofs(0),use_jacobi_preconditioner(true){
        Set_Tolerance((T)1e-14);
        if(use_jacobi_preconditioner){
            for(int i=0;i<d.Size();i++){
                if(A.Row(i).Value_Exists_At_Entry(i))
                    d(i)=(T)1/(T)sqrt(A(i,i));
                else
                    d(i)=(T)1;}}}
		
    ~PRECONDITIONED_CONJUGATE_GRADIENT(){}
	
    virtual void Apply_Preconditioner(){
        //This should compute z=(M^-1)r where M is the preconditioner
        z=r;
    }
	
    void Set_Tolerance(const T tol_input){tolerance=tol_input;}
	
    void Solve(const bool verbose){
        A.Residual(x,b,r);
        p.Set_To_Zero();
		
        for(int i=0;i<max_iterations;i++){
            Apply_Preconditioner();
			
        }
    }
};

template<class T>
class CONJUGATE_GRADIENT{
    SPARSE_MATRIX<T>& A;
    VECTOR<T>& x;
    VECTOR<T>& b;
    VECTOR<T> r,p,q,d,temp;
    VECTOR<int>* dirichlet_dofs;
    int max_iterations;
    T tolerance;
    bool use_jacobi_preconditioner;
	
public:
    CONJUGATE_GRADIENT(SPARSE_MATRIX<T>& A_input,VECTOR<T>& x_input,VECTOR<T>& b_input,const int max_it):A(A_input),
                                                                                                        x(x_input),b(b_input),r(x_input.Size()),max_iterations(max_it),
                                                                                                        p(x_input.Size()),q(x_input.Size()),d(x_input.Size()),temp(x_input.Size()),dirichlet_dofs(0),use_jacobi_preconditioner(true){
        Set_Tolerance((T)1e-14);
        if(use_jacobi_preconditioner){
            for(int i=0;i<d.Size();i++){
                if(A.Row(i).Value_Exists_At_Entry(i))
                    d(i)=(T)1/(T)sqrt(A(i,i));
                else
                    d(i)=(T)1;}
        }
    }
	
    ~CONJUGATE_GRADIENT(){}
	
    void Set_Tolerance(const T tol_input){tolerance=tol_input;}
	
    void Set_Dirichlet_Dofs(VECTOR<int>& dirichlet){
        dirichlet_dofs=&dirichlet;
    }
	
    void Zero_Dirichlet_Residual(){
        if(dirichlet_dofs){
            for(int i=0;i<dirichlet_dofs->Size();i++) r((*dirichlet_dofs)(i))=(T)0;}
    }
	
    void Update_Residual(){
        if(!use_jacobi_preconditioner) A.Residual(b,x,r);
        else{
            for(int i=0;i<b.Size();i++) temp(i)=d(i)*x(i);
            A.Residual(b,temp,r);
            for(int i=0;i<b.Size();i++) r(i)=d(i)*r(i);}
    }
	
    void Multiply_To_Get_q_From_p(){
        if(!use_jacobi_preconditioner) A.Multiply(p,q);
        else{
            for(int i=0;i<p.Size();i++) temp(i)=d(i)*p(i);
            A.Multiply(temp,q);
            for(int i=0;i<p.Size();i++) q(i)=d(i)*q(i);}
    }
	
    int Exit_CG(const int iteration){
        if(!use_jacobi_preconditioner) return iteration;
        else{
            for(int i=0;i<p.Size();i++) x(i)=x(i)*d(i);
            return iteration;}
    }
	
    int Solve(const bool verbose=false){
		
        if(dirichlet_dofs && use_jacobi_preconditioner){
            for(int i=0;i<dirichlet_dofs->Size();i++){
                int index=(*dirichlet_dofs)(i);
                x(index)=x(index)/d(index);}}
		
        Update_Residual();
        Zero_Dirichlet_Residual();
        if(r.L_inf()<tolerance) return 0;
        p=r;
        Multiply_To_Get_q_From_p();
        T r_dot_r=r.Dot(r);
        T alpha=r_dot_r/p.Dot(q);
		
        for(int it=0;it<max_iterations;it++){
			
            for(int i=0;i<r.Size();i++){
                x(i)+=alpha*p(i);
                r(i)-=alpha*q(i);}
			
            Zero_Dirichlet_Residual();
			
            if(verbose) std::cout << "Residual at iteration " << it+1 << " = " << r.L_inf() << std::endl;
            if(r.L_inf()<tolerance) return Exit_CG(it+1);
			
            T r_dot_r_new=r.Dot(r);
            T beta=r_dot_r_new/r_dot_r;
            r_dot_r=r_dot_r_new;
			
            for(int i=0;i<p.Size();i++) p(i)=beta*p(i)+r(i);
			
            Multiply_To_Get_q_From_p();
            alpha=r_dot_r/p.Dot(q);}
		
        return Exit_CG(max_iterations);
    }
	
};
}
#endif
