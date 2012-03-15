#include <vector>
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

    ditto::geometry::Triangle_Mesh_3d<T> mesh;
    ditto::geometry::Box_Hierarchy<T> hierarchy;

    mesh.initialize_parellel_clothes(5 ,11,11, -0.5,0.5, -0.5, 0.5,-0.3,0.3, 1);

    hierarchy.build_tree(mesh.elements, mesh.nodes, 0.01);

    std::cout << "------- childrens -----------------\n";
    for (unsigned int i=0; i<hierarchy.childrens.size(); i++) {
        std::cout << i << " : ";
        for (unsigned int j=0; j<hierarchy.childrens[i].size(); j++) {
            std::cout << hierarchy.childrens[i][j] << " "; }
        std::cout << std::endl; }

    std::cout << "\n----------- boxes ---------------\n";
    for (unsigned int i=0; i<hierarchy.boxes.size(); i++) {
        std::cout << i << " : (" << hierarchy.boxes[i].Pmin(0) <<","<<hierarchy.boxes[i].Pmin(1)<<","<<hierarchy.boxes[i].Pmin(2)<<") ("<<hierarchy.boxes[i].Pmax(0) <<","<<hierarchy.boxes[i].Pmax(1)<<","<<hierarchy.boxes[i].Pmax(2) << ")";
        std::cout << std::endl;}

    std::cout << "\n------------ point query test ----------------\n";

    ditto::algebra::VECTOR_3D<T> P;

    while(1) {
        hierarchy.update_box_positions(mesh.elements, mesh.nodes);
        std::vector<int> intersection_list;
        std::cin >> P(0) >> P(1) >> P(2);
        hierarchy.query_point(P, intersection_list);

        std::cout << "The potential intersection list : \n";
        for (unsigned int i=0; i<intersection_list.size(); i++) {
            std::cout << intersection_list[i] << " "; }
        std::cout << std::endl;
    }

    return 0;
}


