#include <string>
#include <GL/glut.h>
#include <cstdlib>
#include <cmath>

#include <ditto/public_library/geometry/triangle_mesh_2d.h>
#include <ditto/public_library/solid_2d/neo_hookean_2d_fem_implicit.h>

int main(int argc, char ** argv)
{
    typedef double T;

    T dt = 0.001;
    T E = 10;
    T rho = 50;

    ditto::geometry::Triangle_Mesh_2d<T> tm(12,12, -0.5,0.5, -0.5, 0.5, rho);
    ditto::solid_2d::Neo_Hookean_2d_Fem_Implicit<T, ditto::geometry::Triangle_Mesh_2d<T> > cloth(tm, dt, E, 0.3);

    cloth.write_output(1);
    int frame = 1;
    while (1) {
        cloth.update(1e-3, 100, 1e-8, 2000);
        cloth.write_output(++frame); }

    return 0;
}


