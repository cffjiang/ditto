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

#include <ditto/public_library/geometry/triangle_mesh_3d.h>
#include <ditto/public_library/cloth_3d/neo_hookean_cloth_3d_fvm_explicit.h>

int main(int argc, char ** argv)
{
    typedef double T;
    int test = 2;
    T dt;
    T E;
    T rho;
    bool use_gravity;
    T ballv;
    T ball_spring_constant;
    int frame;
    
    if (test == 1) { // ball pass hang cloth
        dt = 1.0/1000.0;
        E = 1000;
        rho = 10;
        use_gravity = true;
        ballv = 1;
        ball_spring_constant = 1e3;
        ditto::geometry::Triangle_Mesh_3d<T> tm(21,21, -0.5,0.5, -0.5, 0.5, rho);
        ditto::cloth_3d::Neo_Hookean_Cloth_3d_Fvm_Explicit<T, ditto::geometry::Triangle_Mesh_3d<T> > cloth(tm, dt, E, 0.3, 0.001, use_gravity);
        cloth.set_dirichlet_with_a_bounding_box(-100,100, 0.4, 0.6, -100, 100, 0.0, 0.0, 0.0);
        cloth.write_output(1);
        frame = 2;
        while (1) {
            std::cout << "--- Frame " << frame <<"  -------------------------------------" << std::endl;
            std::cout << "  simulating time: " << (frame-1)*dt <<  "\n\n";
            cloth.compute_elasticity();
            cloth.add_gravity();
            cloth.add_ball_collision(0, 0.0, 0.5-(frame-1)*dt*ballv, 0.3, ball_spring_constant);
            cloth.update_one_step();
            cloth.write_output(frame++); } }
    else if (test == 2) {  // ball through corner fixed cloth (fracture)
        dt = 1.0/1000.0;
        E = 1000;
        rho = 10;
        use_gravity = true;
        ballv = 1;
        ball_spring_constant = 1e4;
        ditto::geometry::Triangle_Mesh_3d<T> tm(21,21, -0.5,0.5, -0.5, 0.5, rho);
        ditto::cloth_3d::Neo_Hookean_Cloth_3d_Fvm_Explicit<T, ditto::geometry::Triangle_Mesh_3d<T> > cloth(tm, dt, E, 0.3, 0.001, use_gravity);
        cloth.set_dirichlet_with_a_bounding_box(0.45, 0.55, 0.45, 0.55, -100, 100, 0.0, 0.0, 0.0);
        cloth.set_dirichlet_with_a_bounding_box(0.45, 0.55, -0.55, -0.45, -100, 100, 0.0, 0.0, 0.0);
        cloth.set_dirichlet_with_a_bounding_box(-0.55, -0.45, 0.45, 0.55, -100, 100, 0.0, 0.0, 0.0);
        cloth.set_dirichlet_with_a_bounding_box(-0.55 ,-0.45, -0.55, -0.45, -100, 100, 0.0, 0.0, 0.0);
        cloth.write_output(1);
        frame = 2;
        while (1) {
            std::cout << "--- Frame " << frame <<"  -------------------------------------" << std::endl;
            std::cout << "  simulating time: " << (frame-1)*dt <<  "\n\n";
            cloth.compute_elasticity();
            cloth.add_gravity();
            cloth.add_ball_collision(0, 0.0, 0.5-(frame-1)*dt*ballv, 0.3, ball_spring_constant);
            cloth.update_one_step();
            cloth.write_output(frame++); } }
    else if (test == 3) { // ball pass several hang clothes
        dt = 1.0/1000.0;
        E = 1000;
        rho = 10;
        use_gravity = true;
        ballv = 1;
        ball_spring_constant = 1e3;
        int num_of_clothes = 5;

        ditto::geometry::Triangle_Mesh_3d<T> tm;
        tm.initialize_parellel_clothes(num_of_clothes ,21,21, -0.5,0.5, -0.5, 0.5,-0.3,0.3, rho);
        ditto::cloth_3d::Neo_Hookean_Cloth_3d_Fvm_Explicit<T, ditto::geometry::Triangle_Mesh_3d<T> > cloth(tm, dt, E, 0.3, 0.001, use_gravity);
        cloth.aware_of_num_of_clothes(num_of_clothes);
        cloth.set_dirichlet_with_a_bounding_box(-100,100, 0.4, 0.6, -100, 100, 0.0, 0.0, 0.0);
        cloth.write_output(1);
        frame = 2;
        while (1) {
            std::cout << "--- Frame " << frame <<"  -------------------------------------" << std::endl;
            std::cout << "  simulating time: " << (frame-1)*dt <<  "\n\n";
            cloth.compute_elasticity();
            cloth.add_gravity();
            cloth.add_ball_collision(0, 0.0, 0.7-(frame-1)*dt*ballv, 0.3, ball_spring_constant);
            cloth.update_one_step();
            cloth.write_output(frame++); } }

    return 0;
}


