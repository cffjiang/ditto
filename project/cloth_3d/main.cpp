//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// ditto/project/cloth_3d/main.cpp
// Copyright 2012, Chenfanfu Jiang
// Description: Cloth simulation in 3d space.
// Supporting discretization: 
//         FVM with explicity Euler
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

#include <string>
#include <GL/glut.h>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#include <ctime>
#include <ditto/public_library/geometry/triangle_mesh_3d.h>
#include <ditto/public_library/cloth_3d/neo_hookean_cloth_3d_fvm_explicit.h>

int main(int argc, char ** argv)
{
#ifdef _OPENMP
    omp_set_num_threads(4);
#endif

    typedef double T;
    int test = 1;
    T dt;
    T E;
    T rho;
    bool use_gravity;
    T ballv;
    T ball_spring_constant;
    int frame;
    int last_frame = 5000;
   
    if (test == 1) { // ball pass several hang clothes
        dt = 1.0/1000.0;
        E = 1000;
        rho = 10;
        use_gravity = true;
        ballv = 2;
        ball_spring_constant = 1e4;
        int num_of_clothes = 5;
        ditto::geometry::Triangle_Mesh_3d<T> tm;

        tm.initialize_parellel_clothes(num_of_clothes ,11,11, -0.5,0.5, -0.5, 0.5,-0.3,0.3, rho);
        // tm.initialize_regular_mesh(11,11, -0.5,0.5, -0.5, 0.5, rho);

        ditto::cloth_3d::Neo_Hookean_Cloth_3d_Fvm_Explicit<T, ditto::geometry::Triangle_Mesh_3d<T> > cloth(tm, dt, E, 0.3, 0.001, use_gravity);
        cloth.aware_of_num_of_clothes(num_of_clothes);
        cloth.set_dirichlet_with_a_bounding_box(-100,100, 0.4, 0.6, -100, 100, 0.0, 0.0, 0.0);
        frame = 1;
        double very_begin= omp_get_wtime();
        while (1) {
            double begin= omp_get_wtime();
    
            std::cout << "--- Frame " << frame <<"  -------------------------------------" << std::endl;
            std::cout << "  simulating time: " << (frame-1)*dt <<  "\n";
            cloth.compute_elasticity();
            cloth.add_gravity();
            cloth.add_ball_collision(0, 0, 0.7-(frame-1)*dt*ballv, 0.3, ball_spring_constant);
            cloth.switch_self_collision(true, 1, 0.01);
            cloth.update_one_step();

            cloth.write_output(frame++);
            // cloth.write_vtk(frame++);

            if (frame > last_frame) break;

            double end = omp_get_wtime();
            double cost = end - begin;
            double time_passed = end - very_begin;
            double average_cost = time_passed/frame;
            std::cout << "  Time cost for this frame: " << cost << " s"<< std::endl;
            std::cout << "  Time passed in total: " << time_passed << " s" << std::endl;
            std::cout << "  Average time cost for each frame: " << average_cost << " s" << std::endl;
            std::cout << "  To finish " << last_frame << " frames, still need " << average_cost*(last_frame-frame)/60.0 << " minutes.\n" << std::endl;

        } 
    }
    else if (test == 2) {
        dt = 1.0/1000.0;
        E = 1000;
        rho = 10;
        use_gravity = true;
        ballv = 2;
        ball_spring_constant = 1e6;
        int num_of_clothes = 1;
        ditto::geometry::Triangle_Mesh_3d<T> tm;

        tm.initialize_parellel_clothes(num_of_clothes ,11,21, -0.5,0.5, 1, 3,-0.1,0.1, rho);

        ditto::cloth_3d::Neo_Hookean_Cloth_3d_Fvm_Explicit<T, ditto::geometry::Triangle_Mesh_3d<T> > cloth(tm, dt, E, 0.3, 0.001, use_gravity);
        cloth.aware_of_num_of_clothes(num_of_clothes);
        // cloth.set_dirichlet_with_a_bounding_box(-100,100, 0.4, 0.6, -100, 100, 0.0, 0.0, 0.0);
        frame = 1;
        while (1) {
            std::cout << "--- Frame " << frame <<"  -------------------------------------" << std::endl;
            std::cout << "  simulating time: " << (frame-1)*dt <<  "\n\n";
            cloth.compute_elasticity();
            cloth.add_gravity();
            cloth.add_ball_collision(0, 0.5, 0, 0.5, ball_spring_constant);
            cloth.add_ground_collision(-0.6, 1e4, 10);
            cloth.switch_self_collision(true, 2, 0.01);
            cloth.update_one_step();

            cloth.write_output(frame++);
            // cloth.write_vtk(frame++); 
        } 
    }
  
    return 0;
}


