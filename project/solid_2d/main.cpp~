#include <string>
#include <GL/glut.h>
#include <cstdlib>
#include <cmath>

#include <ditto/public_library/geometry/triangle_mesh_2d.h>
#include <ditto/public_library/solid_3d/invertible_finite_element_cloth.h>

int main(int argc, char ** argv)
{
    typedef double T;

    T dt = 0.001;
    T E = 10;
    T rho = 50;

    ditto::geometry::Triangle_Mesh_2d<T> tm(12,12, -0.5,0.5, -0.5, 0.5, rho);
    ditto::solid_3d::Invertible_Finite_Element_Cloth<T, ditto::geometry::Triangle_Mesh_2d<T> > cloth("neo-hookean", "backward euler", tm, dt, E, 0.3);

    cloth.write_output(1);
    int frame = 1;
    while (1) {
        cloth.update(1e-3, 100, 1e-8, 2000);
        cloth.write_output(++frame); }

    return 0;
}


