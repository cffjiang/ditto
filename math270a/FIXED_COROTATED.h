/*
 *  FIXED_COROTATED.h
 *
 *  Created by Chenfanfu Jiang on 10/31/12.
 *
 */
#ifndef _fixed_corotated_
#define _fixed_corotated_

#include "ALGEBRA.h"
#include "GEOMETRY.h"
#include "DEFORMABLE_OBJECTS.h"

using namespace GEOMETRY;
using namespace ALGEBRA;

template <class T>
class FIXED_COROTATED_3D:public HYPERELASTICITY_CONSTITUTIVE_MODEL_3D<T>{
    T lambda;
    T mu;
	
  public:
    FIXED_COROTATED_3D(const T youngs_modulus,const T poisson_ratio){
        lambda=youngs_modulus*poisson_ratio/(((T)1+poisson_ratio)*((T)1-(T)2*poisson_ratio));
        mu=youngs_modulus/((T)2*((T)1+poisson_ratio));
    }
	
    MATRIX_3X3<T> P(const MATRIX_3X3<T>& F){
        MATRIX_3X3<T> U,V,Sigma,R,S,F_inv,F_inv_t;
        VECTOR_3D<T> sigma; 
        T J;
        F.SVD(U,sigma,V);
        Sigma=MATRIX_3X3<T>(sigma(0),0,0,0,sigma(1),0,0,0,sigma(2));
        R=U*(V.Transposed());
        S=V*Sigma*(V.Transposed());
        F_inv=F;
        F_inv.Invert();
        F_inv_t=F_inv.Transposed();
        J=F.Determinant();
        return 2*mu*(F-R)+lambda*(J-1)*J*F_inv_t;
    }
	
    MATRIX_3X3<T> dP(const MATRIX_3X3<T>& F, const MATRIX_3X3<T>& dF){
        MATRIX_3X3<T> dR,dS,F_inv,F_inv_t;
        F.Delta_RS(dF,dR,dS);
        T J=F.Determinant();
        F_inv=F;
        F_inv.Invert();
        F_inv_t=F_inv.Transposed();
        T a=F(0,0),b=F(1,0),c=F(2,0),d=F(0,1),e=F(1,1),f=F(2,1),g=F(0,2),h=F(1,2),i=F(2,2);
        T da=dF(0,0),db=dF(1,0),dc=dF(2,0),dd=dF(0,1),de=dF(1,1),df=dF(2,1),dg=dF(0,2),dh=dF(1,2),di=dF(2,2);
        T dJ=da*e*i+a*de*i+a*e*di-da*h*f-a*dh*f-a*h*df+dd*h*c+d*dh*c+d*h*dc-dd*b*i-d*db*i-d*b*di+dg*b*f+g*db*f+g*b*df-dg*e*c-g*de*c-g*e*dc;
        MATRIX_3X3<T> dJFiT(de*i+e*di-dh*f-h*df,dg*f+g*df-dd*i-d*di,dd*h+d*dh-dg*e-g*de,dh*c+h*dc-db*i-b*di,da*i+a*di-dg*c-g*dc,dg*b+g*db-da*h-a*dh,db*f+b*df-de*c-e*dc,dd*c+d*dc-da*f-a*df,da*e+a*de-dd*b-d*db);
        return 2*mu*(dF-dR)+lambda*(dJ*J*F_inv_t+(J-1)*dJFiT);
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


#endif
