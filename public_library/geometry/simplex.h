//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// ditto/public_library/geometry/simplex.h
// Copyright 2012, Chenfanfu Jiang
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

#ifndef DITTO_PUBLIC_LIBRARY_GEOMETRY_SIMPLEX_H
#define DITTO_PUBLIC_LIBRARY_GEOMETRY_SIMPLEX_H

#include <cmath>

#include <ditto/public_library/algebra/linear_algebra.h>

namespace ditto { namespace geometry {

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Class: Triangle_2d
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T>
class Triangle_2d {
public:
    T x[3][2];
    
    Triangle_2d(T input[]) {
        int k = 0;
        for (int i=0; i<3; i++) {
            for (int j=0; j<2; j++) {
                x[i][j] = input[k++];
            }
        }
    }

    T get_area() {
        return 0.5 * ( x[0][0]*x[1][1] + x[1][0]*x[2][1] + x[2][0]*x[0][1] - x[0][0]*x[2][1] - x[1][0]*x[0][1] - x[2][0]*x[1][1] );
    }
};



} } // end namespaces

#endif
