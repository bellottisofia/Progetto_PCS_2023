#ifndef __FINDDELUNAY_H
#define __FINDDELUNAY_H

#include "iostream"
#include "list"
#include "Eigen/Eigen"
#include "map"
#include "delaunay.hpp"
using namespace std;
using namespace Eigen;

namespace ProjectLibrary
{ Triangulation triangulation;
  Triangle triangle;
  Triangle triangle_max;
  triangle_max = triangle.findMaximumTriangle(points);
  Triangulation triangulation;

}



#endif // __FINDDELUNAY_H
