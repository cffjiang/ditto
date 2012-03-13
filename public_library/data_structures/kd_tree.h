//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// ditto/public_library/data_structures/kd_tree.h
// Copyright 2012, Chenfanfu Jiang
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

#ifndef KD_TREE_H
#define KD_TREE_H

#include <vector>
#include <algorithm>

namespace ditto {

template<int dimensions=3, class VectorT> //VectorT needs to have a [] access member, a - operator, and a lengthSquared() method (distance metric)
class KD_Tree {
public:
    // constructor
    KD_Tree() {
        kdtree = 0;
        kdtreesize = 0; }

    // destructor
    ~KD_Tree() {
        if (kdtree) delete [] kdtree; }

    KD_Tree &operator=(const KD_Tree

};










}
#endif
