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
#include <ditto/public_library/algebra/linear_algebra.h>
#include <ditto/public_library/geometry/triangle_mesh_3d.h>
#include <ditto/public_library/cloth_3d/neo_hookean_cloth_3d_fvm_explicit.h>

int main(int argc, char ** argv)
{
#ifdef _OPENMP
    omp_set_num_threads(4);
#endif
    typedef double T;
    typedef ditto::algebra::VECTOR_3D<T> Vec3;
    int test = 5;

    T dt;
    T E;
    T rho;
    bool use_gravity;
    Vec3 ballV;
    Vec3 ballP;
    T ball_spring_constant;
    int frame;

    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*   
    // Test 1: ball through 1 hanged cloth
    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    if (test == 1) { 
        dt = 1.0/1000.0;
        int last_frame = 8000;
        E = 5000;
        rho = 10;
        use_gravity = true;
        ball_spring_constant = 1e4;
        ditto::geometry::Triangle_Mesh_3d<T> tm;
        tm.initialize_regular_mesh(21, 21, -0.5, 0.5, -0.4, 0.6, rho);
        ditto::cloth_3d::Neo_Hookean_Cloth_3d_Fvm_Explicit<T, ditto::geometry::Triangle_Mesh_3d<T> > cloth(tm, dt, E, 0.3, 0.0001, use_gravity);
        cloth.set_dirichlet_with_a_bounding_box(-100,100, 0.58, 0.62, -100, 100, 0.0, 0.0, 0.0);
        cloth.switch_self_collision(true, 2, 0.01);
        frame = 1;
        ballP.Set_Value(0, 0, 0.4);
        ballV.Set_Value(0, 0, -1.5);
        T ballR = 0.3;
        double very_begin= omp_get_wtime();
        while (1) {
            double begin= omp_get_wtime();
            std::cout << "--- Frame " << frame <<"  -------------------------------------" << std::endl;
            std::cout << "  simulating time: " << (frame-1)*dt <<  "\n";
            cloth.compute_elasticity();
            cloth.add_gravity(0, -9.8,0);
            cloth.add_ball_collision(ballP(0), ballP(1), ballP(2), ballR, ball_spring_constant);
            ballP = ballP + ballV*dt;
            cloth.update_one_step();
            // cloth.write_output(frame++);
            cloth.write_vtk(frame++);
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

    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*   
    // Test 2: ball through 5 hanged clothes
    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    else if (test == 2) {
        dt = 1.0/1000.0;
        int last_frame = 5000;
        E = 1000;
        rho = 10;
        use_gravity = true;
        ball_spring_constant = 1e4;
        int num_of_clothes = 5;
        ditto::geometry::Triangle_Mesh_3d<T> tm;
        tm.initialize_parellel_clothes(num_of_clothes ,21,21, -0.5,0.5, -0.5, 0.5,-0.3,0.3, rho);
        ditto::cloth_3d::Neo_Hookean_Cloth_3d_Fvm_Explicit<T, ditto::geometry::Triangle_Mesh_3d<T> > cloth(tm, dt, E, 0.3, 0.001, use_gravity);
        cloth.aware_of_num_of_clothes(num_of_clothes);
        cloth.set_dirichlet_with_a_bounding_box(-100,100, 0.48, 0.52, -100, 100, 0.0, 0.0, 0.0);
        cloth.switch_self_collision(true, 2, 0.01);
        frame = 1;
        ballP.Set_Value(0, 0, 0.7);
        ballV.Set_Value(0, 0, -2);
        T ballR = 0.3;
        double very_begin= omp_get_wtime();
        while (1) {
            double begin= omp_get_wtime();
            std::cout << "--- Frame " << frame <<"  -------------------------------------" << std::endl;
            std::cout << "  simulating time: " << (frame-1)*dt <<  "\n";
            cloth.compute_elasticity();
            cloth.add_gravity(0, -9.8,0);
            cloth.add_ball_collision(ballP(0), ballP(1), ballP(2), ballR, ball_spring_constant);
            ballP = ballP + ballV*dt;
            cloth.update_one_step();
            // cloth.write_output(frame++);
            cloth.write_vtk(frame++);
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
 
    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*   
    // Test 3: cloth(es) slide over ball, drop to ground
    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    else if (test == 3) { 
        dt = 1.0/5000.0;
        int last_frame = 80000;
        E = 1000;
        rho = 10;
        use_gravity = true;
        ball_spring_constant = 1e4;
        int num_of_clothes = 2;
        ditto::geometry::Triangle_Mesh_3d<T> tm;
        // tm.initialize_regular_mesh_horizental(96, 96, -0.5, 0.5, -0.5, 0.5, rho);
        tm.initialize_parellel_clothes_horizental(num_of_clothes, 96, 96, -0.5, 0.5, 0, 2.0, -0.5, 0.5, rho);
        ditto::cloth_3d::Neo_Hookean_Cloth_3d_Fvm_Explicit<T, ditto::geometry::Triangle_Mesh_3d<T> > cloth(tm, dt, E, 0.3, 0.001, use_gravity);
        cloth.aware_of_num_of_clothes(num_of_clothes);
        cloth.switch_self_collision(true, 2, 0.005);
        frame = 1;
        ballP.Set_Value(0, -0.5, 0.2);
        ballV.Set_Value(0, 0, 0);
        T ballR = 0.5;
        double very_begin= omp_get_wtime();
        while (1) {
            double begin= omp_get_wtime();
            std::cout << "--- Frame " << frame <<"  -------------------------------------" << std::endl;
            std::cout << "  simulating time: " << (frame-1)*dt <<  "\n";
            cloth.compute_elasticity();
            cloth.add_gravity(0, -9.8,0);
            cloth.add_ball_collision(ballP(0), ballP(1), ballP(2), ballR, ball_spring_constant);
            cloth.add_ground_collision(-1, 1e3, 0.5);
            cloth.update_one_step();
            // cloth.write_output(frame++);
            cloth.write_vtk(frame++);
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

    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*   
    // Test 4: clothes drop on small ball
    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    else if (test == 4) { 
        dt = 1.0/1000.0;
        int last_frame = 5000;
        E = 1000;
        rho = 10;
        use_gravity = true;
        ball_spring_constant = 1e4;
        int num_of_clothes = 5;
        ditto::geometry::Triangle_Mesh_3d<T> tm;
        tm.initialize_regular_mesh_horizental(41, 41, -0.5, 0.5, -0.5, 0.5, rho);
        // tm.initialize_parellel_clothes_horizental(num_of_clothes, 41, 41, -0.5, 0.5, 0, 0.4, -0.5, 0.5, rho);
        ditto::cloth_3d::Neo_Hookean_Cloth_3d_Fvm_Explicit<T, ditto::geometry::Triangle_Mesh_3d<T> > cloth(tm, dt, E, 0.3, 0.001, use_gravity);
        // cloth.aware_of_num_of_clothes(num_of_clothes);
        cloth.switch_self_collision(true, 1, 0.005);
        frame = 1;
        ballP.Set_Value(0, -0.2, 0);
        ballV.Set_Value(0, 0, 0);
        T ballR = 0.1;
        double very_begin = omp_get_wtime();
        while (1) {
            double begin= omp_get_wtime();
            std::cout << "--- Frame " << frame <<"  -------------------------------------" << std::endl;
            std::cout << "  simulating time: " << (frame-1)*dt <<  "\n";
            cloth.compute_elasticity();
            cloth.add_gravity(0, -9.8,0);
            cloth.add_ball_collision(ballP(0), ballP(1), ballP(2), ballR, ball_spring_constant);
            cloth.update_one_step();
            // cloth.write_output(frame++);
            cloth.write_vtk(frame++);
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

    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*   
    // Test 5 clothes piling
    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    else if (test == 5) { 
        dt = 1.0/1000.0;
        int last_frame = 20000;
        E = 1000;
        rho = 10;
        use_gravity = true;
        ball_spring_constant = 1e4;
        int num_of_clothes = 75;
        ditto::geometry::Triangle_Mesh_3d<T> tm;
        tm.initialize_parellel_clothes_horizental(num_of_clothes, 11, 11, -0.5, 0.5, 0, 6, -0.5, 0.5, rho);
        ditto::cloth_3d::Neo_Hookean_Cloth_3d_Fvm_Explicit<T, ditto::geometry::Triangle_Mesh_3d<T> > cloth(tm, dt, E, 0.3, 0.001, use_gravity);
        cloth.aware_of_num_of_clothes(num_of_clothes);
        cloth.switch_self_collision(true, 2, 0.015);
        frame = 1;
        double very_begin = omp_get_wtime();
        while (1) {
            double begin= omp_get_wtime();
            std::cout << "--- Frame " << frame <<"  -------------------------------------" << std::endl;
            std::cout << "  simulating time: " << (frame-1)*dt <<  "\n";
            cloth.compute_elasticity();
            cloth.add_gravity(0, -9.8,0);
            cloth.add_ground_collision(-1, 1e3, 0.5);
            cloth.update_one_step();
            // cloth.write_output(frame++);
            cloth.write_vtk(frame++);
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

    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*   
    // Test 6: folding cloth
    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    else if (test == 6) { 
        dt = 1.0/1000.0;
        int last_frame = 8000;
        E = 1000;
        rho = 10;
        use_gravity = true;
        ditto::geometry::Triangle_Mesh_3d<T> tm;
        tm.initialize_regular_mesh_horizental(21, 21, -0.5, 0.5, -0.5, 0.5, rho);
        ditto::cloth_3d::Neo_Hookean_Cloth_3d_Fvm_Explicit<T, ditto::geometry::Triangle_Mesh_3d<T> > cloth(tm, dt, E, 0.3, 0.001, use_gravity);
        T dirichlet_v = 1;
        T dirichlet_new = 0;
        cloth.switch_self_collision(true, 2, 0.01);
        frame = 1;
        double very_begin= omp_get_wtime();
        while (1) {
            double begin= omp_get_wtime();
            std::cout << "--- Frame " << frame <<"  -------------------------------------" << std::endl;
            std::cout << "  simulating time: " << (frame-1)*dt <<  "\n";

            if (frame % 500 == 0) dirichlet_v *= -1;
            Vec3 dirichlet_displacement( 0, 0,  dirichlet_new);
            dirichlet_new += dt*dirichlet_v;

            cloth.clear_previous_dirichlet_conditions();
            
            if (frame < 1100) {
                cloth.set_dirichlet_with_a_bounding_box(-0.52, -0.48,-100,100, 0.48, 0.52, dirichlet_displacement(0), dirichlet_displacement(1), dirichlet_displacement(2));
                cloth.set_dirichlet_with_a_bounding_box(0.48, 0.52,-100,100, 0.48, 0.52, dirichlet_displacement(0), dirichlet_displacement(1), dirichlet_displacement(2)); }

            cloth.compute_elasticity();
            cloth.add_gravity(0, -9.8,0);
            cloth.add_ground_collision(-0.3, 1e4, 0.5);
            cloth.update_one_step();
            // cloth.write_output(frame++);
            cloth.write_vtk(frame++);
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
  
  
    return 0;
}


