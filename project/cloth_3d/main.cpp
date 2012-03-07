#include <string>
#include <GL/glut.h>
#include <cstdlib>
#include <cmath>

#include <ditto/public_library/geometry/triangle_mesh_3d.h>
#include <ditto/public_library/cloth_3d/neo_hookean_cloth_3d_fvm_explicit.h>

int main(int argc, char ** argv)
{
    typedef double T;

    T dt = 0.0001;
    T E = 1000;
    T rho = 1;
    bool use_gravity = true;
    
    ditto::geometry::Triangle_Mesh_3d<T> tm(12,12, -0.5,0.5, -0.5, 0.5, rho);
    ditto::cloth_3d::Neo_Hookean_Cloth_3d_Fvm_Explicit<T, ditto::geometry::Triangle_Mesh_3d<T> > cloth(tm, dt, E, 0.3, use_gravity);

    cloth.set_dirichlet_with_a_bounding_box(0.48,0.62, 0.48, 0.62, -100, 100, 0.0, 0.0, 0.0);
    cloth.set_dirichlet_with_a_bounding_box(0.45,0.65,  -0.62,-0.48,  -100, 100, 0.0, 0.0, 0.0);
    cloth.set_dirichlet_with_a_bounding_box(-0.65, -0.45, 0.48, 0.62, -100, 100, 0.0, 0.0, 0.0);
    cloth.set_dirichlet_with_a_bounding_box(-0.65,-0.45, -0.62,-0.45, -100, 100, 0.0, 0.0, 0.0);

        
    cloth.write_output(1);
    int frame = 1;
    while (1) {
        cloth.update();
        cloth.write_output(++frame); 
    }

    return 0;
}


