
#include <iostream>
#include <iomanip>
#include <vector>

#include <nlopt.hpp>

using namespace std;


int s=4;
double phi_dist_approximate[4] = {3.0*1.414/4.0, 0.5, 0.5, -1.414/4.0};
double phi_result[4];

typedef struct {
    double a, b;
    int i,j;
} my_constraint_data;


double myfunc(unsigned n, const double* x, double* grad, void* my_func_data)
{
    if(grad){
        for(int i=0;i<s;i++){
            grad[i] = 2*(x[i]-phi_dist_approximate[i]);
        }
    }

    double value=0;
    for(int i=0;i<s;i++){
        double diff = x[i]-phi_dist_approximate[i];
        value += diff*diff;
    }
        
    return value;
}

double myconstraint(unsigned n, const double* x, double *grad, void *data)
{
    my_constraint_data *d = (my_constraint_data *) data;
    double a=d->a, b=d->b;
    if (grad) {
        for(int i=0;i<s;i++){
            if(i == d->i) grad[i] = a;
            else if(i == d->j) grad[i] = b;
            else grad[i] = double(0);
        }
    }

    return a*x[d->i]+b*x[d->j];
}


// solve: build phi_result
// constraints are in the form of i*phi[a]+j*phi[b]=0
bool solve_problem_with_constraints(const int *i, const int *j, const double *a, const double* b,const int num_constraints) 
{
    nlopt_opt opt;

    // use local optimization
    opt = nlopt_create(NLOPT_LD_SLSQP, s);

    // set lower bound to not less than phi before damge evolution
    double lb[4] = { 1.0, 0.32, 0.32, -0.5};
    nlopt_set_lower_bounds(opt, lb);

    nlopt_set_min_objective(opt, myfunc, NULL);

    std::vector<my_constraint_data> data(num_constraints);
    for(int k=0;k<num_constraints;k++) {
        data[k].a=a[k];
        data[k].b=b[k];
        data[k].i=i[k];
        data[k].j=j[k];
        nlopt_add_equality_constraint(opt,myconstraint, &data[k],1e-15);
    }
    nlopt_set_xtol_rel(opt, 1e-12);

    double x[4];
    for(int k=0;k<s;k++) x[k]=phi_dist_approximate[k];

    double minf;
    

    if (nlopt_optimize(opt, x, &minf)<0) {
        std::cout<< "optimization failed!\n";
        nlopt_destroy(opt);
        return false;
    }
    else{
        for(int k=0; k<s; k++)
            phi_result[k] = x[k];
        nlopt_destroy(opt);

        return true;
    }
}

int main()
{
    int num_constraints = 3;
    int i[3]={0, 1, 2};
    int j[3]={3,3,3};
    double a[3]={1, 1, 4};
    double b[3]={3, 1, 4};

    solve_problem_with_constraints(i,j,a,b,num_constraints);

    std::cout<<"result: ";
    for(int k=0;k<s;k++)
        std::cout<<std::setprecision(15)<<phi_result[k]<<"\n";
    std::cout<<std::endl;


    return 0;
}
