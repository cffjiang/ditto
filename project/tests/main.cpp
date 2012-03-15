
#include <string>
#include <GL/glut.h>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#include <ctime>
#include <ditto/public_library/geometry/triangle_mesh_3d.h>
#include <ditto/public_library/geometry/box_hierarchy.h>

int main(int argc, char ** argv)
{
    typedef double T;

    ditto::geometry::Box_Hierarchy<T> tree;

    for (int leafs = 4; leafs <=8 ; leafs++) {
        std::cout << " --------------------- test with " << leafs << " leafs\n";

        int boxes = 0;
        for (int num_current_level = leafs; num_current_level >= 1; num_current_level/=2) {
            boxes += num_current_level; }

        tree.build_hierarchy_structure_based_on_mesh_connectivity(leafs, boxes);

        for (unsigned int i=0; i<tree.childrens.size(); i++) {
            std::cout << i << " : ";
            for (unsigned int j=0; j<tree.childrens[i].size(); j++) {
                std::cout << tree.childrens[i][j] << " "; }
            std::cout << std::endl; }

    }  


    return 0;
}


