//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-**-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Trying to use PhysBAM Tetrahedral mesh hierarchy
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-**-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

#include <iostream>

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>

int main(int argc, char ** argv)
{
    typedef float T;

    PhysBAM::TETRAHEDRALIZED_VOLUME<T>* tv = PhysBAM::TETRAHEDRALIZED_VOLUME<T>::Create();
    int N_elements = 3;
    int N_particles = 9;

    // initialize tetrahedralized volume (tv): particles and mesh.elements
    tv->Clean_Memory();
    tv->particles.Delete_All_Elements();
    tv->particles.Add_Elements(N_particles);
    tv->mesh.elements.Exact_Resize(N_elements);

    // manually assign elements
    tv->mesh.elements(0).Set(0, 1, 2, 3);
    tv->mesh.elements(1).Set(0, 2, 1, 4);
    tv->mesh.elements(2).Set(5, 6, 7, 8);

    // manually assign particles
    tv->particles.X(0) = PhysBAM::VECTOR<T,3>(1, 0, 0);
    tv->particles.X(1) = PhysBAM::VECTOR<T,3>(0, 1, 0);
    tv->particles.X(2) = PhysBAM::VECTOR<T,3>(0, 0, 0);
    tv->particles.X(3) = PhysBAM::VECTOR<T,3>(0, 0, 1);
    tv->particles.X(4) = PhysBAM::VECTOR<T,3>(0, 0, -1);
    tv->particles.X(5) = PhysBAM::VECTOR<T,3>(.5, 0, 0);
    tv->particles.X(6) = PhysBAM::VECTOR<T,3>(0, .5, 0);
    tv->particles.X(7) = PhysBAM::VECTOR<T,3>(0, 0, .5);
    tv->particles.X(8) = PhysBAM::VECTOR<T,3>(0, 0, 1.5);

    tv->Update_Number_Nodes(); // for safety

    // initialize box hierarchy
    tv->Initialize_Hierarchy();

    // test points
    while(true)
    {
        PhysBAM::ARRAY<int> intersection_list;
        PhysBAM::VECTOR<T,3> P;

        std::cout << "Please input test point position: ";
        std::cin >> P(0) >> P(1) >> P(2);

        tv->hierarchy->Intersection_List(P, intersection_list);

        std::cout << "The potential intersecting tetrahedrals are: ";
    }
    delete tv;
    return 0;
}


