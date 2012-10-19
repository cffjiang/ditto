#include <ALGEBRA.h>
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;
using namespace MATH270A;

int main()
{
    typedef double T;
    srand ( time(NULL) );

    MATRIX_2X2<T> F,U,V,UtU,VtV;
    VECTOR_2D<T> sigma;

    for(int test=0;test<100000;test++){
        std::cout<<"test "<<test<<std::endl;
        // Randomly build F from -2, 0, 1, 2
        random_integer = -3+rand()%5+1;
        F(0,0)=random_integer;
        random_integer = -3+rand()%5+1;
        F(0,1)=random_integer;
        random_integer = -3+rand()%5+1;
        F(1,0)=random_integer;
        random_integer = -3+rand()%5+1;
        F(1,1)=random_integer;

        F.SVD(U,sigma,V);
        MATRIX_2X2<T> Iden=MATRIX_2X2<T>::Identity();
        UtU=U.Transposed()*U;
        VtV=V.Transposed()*V;
        MATRIX_2X2<T> utu_error=UtU-Iden;
        MATRIX_2X2<T> vtv_error=VtV-Iden;
        MATRIX_2X2<T> Sigma(sigma(0),0,0,sigma(1));
        MATRIX_2X2<T> restore=U*Sigma*V.Transposed();;
        MATRIX_2X2<T> error=restore-F;
        if(utu_error.Norm()>1e-10 || vtv_error.Norm()>1e-10 || error.Norm() >1e-10 || (sigma(0)<0&&sigma(1)<0) || (sigma(0)<0 && std::abs(sigma(0))>std::abs(sigma(1))) || (sigma(1)<0 && std::abs(sigma(1))>std::abs(sigma(0)))){
            std::cout<<"F "; F.Print();
            std::cout<<"U "; U.Print();
            std::cout<<"UtU "; UtU.Print();
            std::cout<<"V ";    V.Print();
            std::cout<<"VtV "; VtV.Print();
            std::cout<<"sigma ";    sigma.Print();
            std::cout<<"restore ";    restore.Print();
            std::cout<<"error ";    error.Print();
            throw(0);
        }
    }

    return 0;
}
