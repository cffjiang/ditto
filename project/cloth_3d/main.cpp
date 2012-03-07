//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// ditto/project/cloth_3d/main.cpp
// Copyright 2012, Chenfanfu Jiang
//
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

    T dt = 1.0/1000.0;
    T E = 1000;
    T rho = 10;
    bool use_gravity = true;

    T ballv = 1;
    
    ditto::geometry::Triangle_Mesh_3d<T> tm(21,21, -0.5,0.5, -0.5, 0.5, rho);
    ditto::cloth_3d::Neo_Hookean_Cloth_3d_Fvm_Explicit<T, ditto::geometry::Triangle_Mesh_3d<T> > cloth(tm, dt, E, 0.3, 0.001, use_gravity);

    cloth.set_dirichlet_with_a_bounding_box(-100,100, 0.4, 0.6, -100, 100, 0.0, 0.0, 0.0);
    cloth.write_output(1);
    int frame = 2;
    while (1) {
        std::cout << "--- Frame " << frame <<"  -------------------------------------" << std::endl;
        std::cout << "  simulating time: " << (frame-1)*dt <<  "\n\n";
        cloth.compute_elasticity();
        cloth.add_gravity();
        cloth.add_ball_collision(0, 0.0, 0.5-(frame-1)*dt*ballv, 0.3, 1e3);
        cloth.update_one_step();
        cloth.write_output(frame++); 
    }

    return 0;
}


