//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
// equality_least_squares_with_pairwise_equality_constraints_and_lower_bounds.hpp
// Created by Chenfanfu Jiang in 2012
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
#ifndef EQUALITY_LEAST_SQUARES_WITH_PAIRWISE_EQUALITY_CONSTRAINTS_AND_LOWER_BOUNDS_HPP
#define EQUALITY_LEAST_SQUARES_WITH_PAIRWISE_EQUALITY_CONSTRAINTS_AND_LOWER_BOUNDS_HPP

// 
//  Solves the following optimization problem using nlopt library:
//
//    solve for the x that minimizes sum||x[i]-y[i]||_2^2
//    while satisfying the equality constraints: a[i]x[i]+b[j]x[j]=0 for some i,j
//    and also satisfying the inequality constraints: x[i] >= lb[i]
//
//  It guarantees to give a solution when the following conditions are satisfied:
//    1. the number of constraints is not more than the number of x
//    2. equality constraints are linear independent (constraint matrix is full rank)
//    3. the lower bound input by user makes sense -- this means the user is
//       responsible for making sure the solution to the problem exists.
//
//  usage:  
//    int result = sake::optimiation::equality_least_squares_with_pairwise_equality_constraints_and_lower_bounds::solve(
//                         const int *i,                             // see description above
//                         const int *j,                             // see description above
//                         const double *a,                          // see description above
//                         const double* b,                          // see description above
//                         const int num_constraints,                // number of equality constraints
//                         double *y,                                // see description above
//                         const double *lb,                         // see description above
//                         const int number_of_DOFs,                 // number of entries in x (and y and lb)
//                         const double constraint_eps,              // the criteria for satisfying equality constraint 
//                         const double solution_convergence_eps,    // the criteria for saying x converges 
//
//                         double *x                                 // the output x
//                       ) 
//
//    return value:   positive means success, negative means failure.
//                    for detailed meanings for debugging see 
//                    http://ab-initio.mit.edu/wiki/index.php/NLopt_Reference#Local.2Fsubsidiary_optimization_algorithm
//

#include <iostream>
#include <iomanip>
#include <vector>

#include <nlopt.hpp>

namespace sake{
namespace optimization{
namespace equality_least_squares_with_pairwise_equality_constraints_and_lower_bounds{


typedef struct {
    double *phi_dist_approximate;
    int s;
} my_function_data;


typedef struct {
    double a, b;
    int i,j;
    int s;
} my_constraint_data;


double myfunc(unsigned n, const double* x, double* grad, void* my_func_data_in)
{
    my_function_data *my_func_data = (my_function_data *)my_func_data_in;

    if(grad){
        for(int i=0;i<my_func_data->s;i++){
            grad[i] = 2*(x[i]-my_func_data->phi_dist_approximate[i]);
        }
    }

    double value=0;
    for(int i=0;i<my_func_data->s;i++){
        double diff = x[i]-my_func_data->phi_dist_approximate[i];
        value += diff*diff;
    }
        
    return value;
}

double myconstraint(unsigned n, const double* x, double *grad, void *data)
{
    my_constraint_data *d = (my_constraint_data *) data;
    double a=d->a, b=d->b;
    if (grad) {
        for(int i=0;i<d->s;i++){
            if(i == d->i) grad[i] = a;
            else if(i == d->j) grad[i] = b;
            else grad[i] = double(0);
        }
    }

    return a*x[d->i]+b*x[d->j];
}

int solve(const int *i, const int *j, const double *a, const double* b,const int num_constraints,   double *phi_dist_approximate, const double *phi_lb, const int s , const double constraint_eps, const double solution_convergence_eps, double * phi_result ) 
{
    if(num_constraints > s){
        std::cout<< "SHIT! too many constraints!";
        throw(0);
    }

    nlopt_opt opt;

    // use local optimization
    opt = nlopt_create(NLOPT_LD_SLSQP, s);
    //opt = nlopt_create(NLOPT_GN_ISRES, s);

    // set lower bound to not less than phi before damge evolution
    double *lb = new double[s];
    for(int k=0;k<s;k++) lb[k]=phi_lb[k];

    nlopt_set_lower_bounds(opt, lb);

    my_function_data func_data={phi_dist_approximate, s};
    nlopt_set_min_objective(opt, myfunc, &func_data);

    std::vector<my_constraint_data> data(num_constraints);
    for(int k=0;k<num_constraints;k++) {
        data[k].a=a[k];
        data[k].b=b[k];
        data[k].i=i[k];
        data[k].j=j[k];
        data[k].s=s;
        nlopt_add_equality_constraint(opt,myconstraint, &data[k],constraint_eps);
    }
    nlopt_set_xtol_rel(opt, solution_convergence_eps);

    double *x = new double[s];

    for(int k=0;k<s;k++) {
        x[k]= (phi_dist_approximate[k]<lb[k]) ? lb[k]:phi_dist_approximate[k];  // need to let the initialial guess in bound!!
    }
    double minf;
    
    int result= nlopt_optimize(opt,x,&minf);

    if(result > 0){
        for(int k=0; k<s;k++){
            phi_result[k]=x[k];
        }
    }
    else{
        std::cout<<"FATAL ERROR! Can't solve optimization.\n";
        std::cout<<result<<std::endl;
        std::cout<<minf<<std::endl;
    }

    nlopt_destroy(opt);


    delete []lb;
    delete []x;
    return result;

}

}}}

#endif
