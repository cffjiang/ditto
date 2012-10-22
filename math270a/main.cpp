#include <ALGEBRA.h>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace MATH270A;

int main()
{
    typedef double T;
    srand ( time(NULL) );

    //##############################################################################################
    // 2D TEST(hmwk) Deformation Gradient
    //##############################################################################################
    if(1){
        VECTOR_2D<T> X0(0,0),X1(1,0),X2(0,1),x0(0,1),x1(0,1),x2(0,0);
        MATRIX_2X2<T> F=MATRIX_2X2<T>::Deformation_Gradient(X0,X1,X2,x0,x1,x2);
        F.Print();
    }

    //##############################################################################################
    // 3D TEST(hmwk) Deformation Gradient
    //##############################################################################################
    if(1){
        VECTOR_3D<T> X0(0,0,0),X1(2,0,0),X2(0,1,0),X3(0,0,3),x0(0,0,0),x1(1,0,0),x2(0,1,0),x3(0,0,1);
        MATRIX_3X3<T> F=MATRIX_3X3<T>::Deformation_Gradient(X0,X1,X2,X3,x0,x1,x2,x3);
        F.Print();
    }

    //##############################################################################################
    // 2D TEST for SVD and differentiating SVD
    //##############################################################################################
    if(0){
        MATRIX_2X2<T> F,U,V,UtU,VtV,UUt,VVt,DeltaF,DeltaU,DeltaV,DeltaSigma;
        VECTOR_2D<T> sigma,delta_sigma;
        int random_integer;

        for(int test=0;test<1000000;test++){
            std::cout<<"test 2d svd "<<test<< " \n";

            for(int i=0;i<2;i++){
                for(int j=0;j<2;j++){
                    random_integer =-2+rand()%6;
                    F(i,j)=random_integer;
                    random_integer =-4+rand()%8;
                    DeltaF(i,j)=random_integer;}}

            // F=MATRIX_2X2<T>(1,2,3,4);
            // DeltaF=MATRIX_2X2<T>(1,0,0,1);

            F.SVD(U,sigma,V);
            F.Delta_SVD(DeltaF,delta_sigma,DeltaU,DeltaV);

            UtU=U.Transposed()*U;
            VtV=V.Transposed()*V;
            UUt=U*U.Transposed();
            VVt=V*V.Transposed();
            DeltaSigma=MATRIX_2X2<T>(delta_sigma(0),0,0,delta_sigma(1));
            MATRIX_2X2<T> Sigma(sigma(0),0,0,sigma(1));
            MATRIX_2X2<T> LeftDifferentiate=U.Transposed()*DeltaF*V;
            MATRIX_2X2<T> RightDifferentiate=U.Transposed()*DeltaU*Sigma+DeltaSigma+Sigma*(DeltaV.Transposed())*V;
                
            MATRIX_2X2<T> Iden=MATRIX_2X2<T>::Identity();
            MATRIX_2X2<T> utu_error=UtU-Iden;
            MATRIX_2X2<T> vtv_error=VtV-Iden;
            MATRIX_2X2<T> uut_error=UUt-Iden;
            MATRIX_2X2<T> vvt_error=VVt-Iden;
            MATRIX_2X2<T> restore=U*Sigma*V.Transposed();;
            MATRIX_2X2<T> error=restore-F;
            MATRIX_2X2<T> differentiate_error=RightDifferentiate-LeftDifferentiate;
            MATRIX_2X2<T> FtF=F.Transposed()*F;
            MATRIX_2X2<T> FV=F*V;

            // std::cout<<"F "; F.Print();
            // std::cout<<"DelataF "; DeltaF.Print();
            // std::cout<<"FV "; FV.Print();
            // std::cout<<"U "; U.Print();
            // std::cout<<"UtU "; UtU.Print();
            // std::cout<<"V ";    V.Print();
            // std::cout<<"VtV "; VtV.Print();
            // std::cout<<"sigma ";    sigma.Print();
            // std::cout<<"restore ";    restore.Print();
            // std::cout<<"DeltaSigma "; DeltaSigma.Print();
            // std::cout<<"DeltaU "; DeltaU.Print();
            // std::cout<<"DeltaV "; DeltaV.Print();
            // std::cout<<"error ";    error.Print();
            // std::cout<<"differentiate_error ";    differentiate_error.Print();
            // std::cout<<"u determinant "<< U.Determinant()<<std::endl;
            // std::cout<<"V determinant "<< V.Determinant()<<std::endl;
            // throw(0);


            if(std::abs(U.Determinant()-1)>1e-10 || 
                abs(V.Determinant()-1)>1e-10 ||
                utu_error.Norm()>1e-10 || 
                vtv_error.Norm()>1e-10 || 
                uut_error.Norm()>1e-10 || 
                vvt_error.Norm()>1e-10 ||
                error.Norm() >1e-10 ||
                (differentiate_error.Norm() > 1e-8 && std::abs(sigma(0)-sigma(1))>1e-10 && std::abs(sigma(0)+sigma(1))>1e-10) ||
                (sigma(0)<0&&sigma(1)<0) || 
                (sigma(0)<0 && std::abs(sigma(0))>std::abs(sigma(1))) || 
                (sigma(1)<0 && std::abs(sigma(1))>std::abs(sigma(0))) ||
                !(sigma(0)>=sigma(1)) ||
                !(std::abs(sigma(0))>=std::abs(sigma(1))) ){

                std::cout<<"F "; F.Print();
                std::cout<<"DelataF "; DeltaF.Print();
                std::cout<<"FV "; FV.Print();
                std::cout<<"U "; U.Print();
                std::cout<<"UtU "; UtU.Print();
                std::cout<<"V ";    V.Print();
                std::cout<<"VtV "; VtV.Print();
                std::cout<<"sigma ";    sigma.Print();
                std::cout<<"DeltaSigma "; DeltaSigma.Print();
                std::cout<<"DeltaU "; DeltaU.Print();
                std::cout<<"DeltaV "; DeltaV.Print();
                std::cout<<"restore ";    restore.Print();
                std::cout<<"error ";    error.Print();
                std::cout<<"differentiate_error ";    differentiate_error.Print();
                std::cout<<"u determinant "<< U.Determinant()<<std::endl;
                std::cout<<"V determinant "<< V.Determinant()<<std::endl;
                throw(0);}

        }
    }
    //##############################################################################################
    // 3D TEST for SVD and differentiating SVD
    //##############################################################################################
    if(0){
        MATRIX_3X3<T> F,U,V,UtU,VtV,UUt,VVt,DeltaF,DeltaSigma,DeltaU,DeltaV;
        VECTOR_3D<T> sigma,delta_sigma;
        int random_integer;

        for(int test=0;test<1000000;test++){
            std::cout<<"test 3d svd "<<test<< " \n";

            for(int i=0;i<3;i++){
                for(int j=0;j<3;j++){
                    random_integer =-2+rand()%6;
                    F(i,j)=random_integer;
                    random_integer =-4+rand()%8;
                    DeltaF(i,j)=random_integer;}}

            // F=MATRIX_3X3<T>(1,2,3,4,5,6,7,8,9);
            // DeltaF=MATRIX_3X3<T>(1,0,0,0,1,0,0,0,1);

            F.SVD(U,sigma,V);
            F.Delta_SVD(DeltaF,delta_sigma,DeltaU,DeltaV);
            MATRIX_3X3<T> Sigma(sigma(0),0,0,0,sigma(1),0,0,0,sigma(2));
            MATRIX_3X3<T> DeltaSigma(delta_sigma(0),0,0,0,delta_sigma(1),0,0,0,delta_sigma(2));

            MATRIX_3X3<T> Iden=MATRIX_3X3<T>::Identity();
            UtU=U.Transposed()*U;
            VtV=V.Transposed()*V;
            UUt=U*U.Transposed();
            VVt=V*V.Transposed();
            MATRIX_3X3<T> utu_error=UtU-Iden;
            MATRIX_3X3<T> vtv_error=VtV-Iden;
            MATRIX_3X3<T> uut_error=UUt-Iden;
            MATRIX_3X3<T> vvt_error=VVt-Iden;
            MATRIX_3X3<T> restore=U*Sigma*V.Transposed();;
            MATRIX_3X3<T> error=restore-F;
            MATRIX_3X3<T> LeftDifferentiate=U.Transposed()*DeltaF*V;
            MATRIX_3X3<T> RightDifferentiate=U.Transposed()*DeltaU*Sigma+DeltaSigma+Sigma*(DeltaV.Transposed())*V;
            MATRIX_3X3<T> differentiate_error=RightDifferentiate-LeftDifferentiate;
            MATRIX_3X3<T> FtF=F.Transposed()*F;
            MATRIX_3X3<T> FV=F*V;

            // std::cout<<"F "; F.Print();
            // std::cout<<"DeltaF "; DeltaF.Print();
            // std::cout<<"FV "; FV.Print();
            // std::cout<<"U "; U.Print();
            // std::cout<<"UtU "; UtU.Print();
            // std::cout<<"UUt "; UUt.Print();
            // std::cout<<"V ";    V.Print();
            // std::cout<<"VtV "; VtV.Print();
            // std::cout<<"VVt "; VVt.Print();
            // std::cout<<"sigma : "; sigma.Print();
            // std::cout<<"restore ";    restore.Print();
            // std::cout<<"DeltaSigma "; DeltaSigma.Print();
            // std::cout<<"DeltaU "; DeltaU.Print();
            // std::cout<<"DeltaV "; DeltaV.Print();
            // std::cout<<"error ";    error.Print();
            // std::cout<<"differentiate_error ";    differentiate_error.Print();
            // std::cout<<"u determinant "<< U.Determinant()<<std::endl;
            // std::cout<<"V determinant "<< V.Determinant()<<std::endl;
            // throw(0);

            if(
                std::abs(U.Determinant()-1)>1e-7 || 
                abs(V.Determinant()-1)>1e-7 ||
                utu_error.Norm()>1e-7 || 
                vtv_error.Norm()>1e-7 || 
                uut_error.Norm()>1e-7 || 
                vvt_error.Norm()>1e-7 ||
                error.Norm() >1e-7 ||
                (differentiate_error.Norm() > 1e-5 && std::abs(sigma(0)-sigma(1))>1e-10 && std::abs(sigma(0)-sigma(2))>1e-10 && std::abs(sigma(1)-sigma(2))>1e-10 && std::abs(sigma(0)+sigma(1))>1e-10 && std::abs(sigma(0)+sigma(2))>1e-10 && std::abs(sigma(1)+sigma(2))>1e-10) ||
                !(sigma(0)>=sigma(1)-1e-10&&sigma(1)>=sigma(2)-1e-10) ||
                !(std::abs(sigma(0))>=std::abs(sigma(1))-1e-10&&std::abs(sigma(1))>=std::abs(sigma(2))-1e-10) ||
                !(sigma(0)>=0&&sigma(1)>=0)
            ){
                std::cout<<"F "; F.Print();
                std::cout<<"DelataF "; DeltaF.Print();
                std::cout<<"FV "; FV.Print();
                std::cout<<"U "; U.Print();
                std::cout<<"UtU "; UtU.Print();
                std::cout<<"UUt "; UUt.Print();
                std::cout<<"V ";    V.Print();
                std::cout<<"VtV "; VtV.Print();
                std::cout<<"VVt "; VVt.Print();
                std::cout<<"sigma : "; sigma.Print();
                std::cout<<"restore ";    restore.Print();
                std::cout<<"DeltaSigma "; DeltaSigma.Print();
                std::cout<<"DeltaU "; DeltaU.Print();
                std::cout<<"DeltaV "; DeltaV.Print();

                std::cout<<"error ";    error.Print();
                std::cout<<"differentiate_error ";    differentiate_error.Print();
                std::cout<<"u determinant "<< U.Determinant()<<std::endl;
                std::cout<<"V determinant "<< V.Determinant()<<std::endl;
                throw(0);}

        }
    }


    return 0;
}
