#include <iostream>
#include "equality_least_squares_with_pairwise_equality_constraints_and_lower_bounds.hpp"



int main()
{
    using namespace sake::optimization::equality_least_squares_with_pairwise_equality_constraints_and_lower_bounds;
    int num_constraints = 3;
    int i[3]={0, 1, 2};
    int j[3]={3,3,3};
    double a[3]={1, 1, 4};
    double b[3]={3, 1, 4};
    int s=4;
    double phi_dist[4] = {3*1.414/4, 0.5, 0.5, -1.414/4};
    double phi[4];
    double phi_lb[4] = {5,5,5,-70};
    int result =solve(i,j,a,b,num_constraints,phi_dist,phi_lb,s,1e-15, 1e-12,phi);

    if (result<0){
        std::cout<< "optimization failed!\n";

        throw(0);
    }
    else{

        std::cout<<"result: ";
        for(int k=0;k<s;k++)
            std::cout<<std::setprecision(15)<<phi[k]<<"\n";
        std::cout<<std::endl;


    }


    return 0;
}
