//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// ditto/public_library/visualization/vtk_writer.h
// Copyright 2012, Chenfanfu Jiang
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

#ifndef DITTO_PUBLIC_LIBRARY_VISUALIZATION_VTK_WRITER_H
#define DITTO_PUBLIC_LIBRARY_VISUALIZATION_VTK_WRITER_H

#include <cstdlib>
#include <cmath>
#include <fstream>
#include <vector>
#include <set>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <string>

namespace ditto{ namespace visualization {

template<typename T>
void write_vtk_sphere(T x, T y, T z, T r, std::string filename) 
{
    // input file
    std::string input_file = "../../public_data/vtk_data/sphere.vtk";
    std::ifstream in(input_file.c_str());
    std::istringstream ss;
    std::string s;
    
    // output file
    std::fstream out(filename.c_str(), std::ios_base::out); 
    
    int currentline = 1;
    while (getline(in, s)) {
        ss.clear();
        ss.str(s);
        
        if (currentline >= 6 && currentline <= 125) {
            T pts[9];
            for (int i=0; i<9; i++) {
                ss >> pts[i];
                pts[i] *= r; }
            for (int i=0; i<=6; i+=3) {
                pts[i] += x; }
            for (int i=1; i<=7; i+=3) {
                pts[i] += y; }
            for (int i=2; i<=8; i+=3) {
                pts[i] += z; }
            
            for (int i=0; i<9; i++) {
                out << pts[i] << " "; }
            out << std::endl;
        } else if (currentline == 126) {
            T pts[6];
            for (int i=0; i<6; i++) {
                ss >> pts[i];
                pts[i] *= r; }
            for (int i=0; i<=3; i+=3) {
                pts[i] += x; }
            for (int i=1; i<=4; i+=3) {
                pts[i] += y; }
            for (int i=2; i<=5; i+=3) {
                pts[i] += z; }
            
            for (int i=0; i<6; i++) {
                out << pts[i] << " "; }
            out << std::endl;
        } else {
            out << s << std::endl;
        }
        currentline++;
    }
}


} }

#endif
