#ifndef __DELAUNAY_H
#define __DELAUNAY_H

#include <iostream>
#include "list"
#include "Eigen/Eigen"
#include "map"
#include <cmath>
#include "sorting.hpp"


using namespace std;
using namespace Eigen;
namespace ProjectLibrary
{


  constexpr double max_tolerance(const double& x, const double& y)
  {
    return x > y ? x : y;
  }


struct Point {
    double x;
    double y;
    unsigned int id;
    Point* succ = nullptr;
    Point* prec = nullptr;

    static constexpr double geometricTol = 1.0e-12;
    static constexpr double geometricTol_Squared = max_tolerance(Point::geometricTol * Point::geometricTol,
                                                                 numeric_limits<double>::epsilon());
    Point(const double& x, const double& y, const unsigned int& id)
        : x(x), y(y), id(id)
    {
    }

    Point(const Point& p)
        : x(p.x), y(p.y), id(p.id)
    {
    }
};

double normSquared(const double& x, const double& y)
{
    return x * x + y * y;
}

bool operator==(const Point& p1, const Point& p2)
{
    return (normSquared(p1.x - p2.x, p1.y - p2.y) <=
            Point::geometricTol * Point::geometricTol *
                std::max(normSquared(p1.x, p1.y), normSquared(p2.x, p2.y)));
}

bool operator!=(const Point& p1, const Point& p2)
{
    return !(p1 == p2);
}

std::ostream& operator<<(std::ostream& os, const Point& p2)
{
    os << p2.id;
    return os;
}

bool operator>(const Point& p1, const Point& p2)
{
    return p1.x > p2.x + Point::geometricTol * std::max(p1.x, p2.x);
}

bool operator<=(const Point& p1, const Point& p2)
{
    return !(p1 > p2);
}

Point operator-(const Point& p1, const Point& p2)
{
    return Point(p1.x - p2.x, p1.y - p2.y, 0);
}

Point operator+(const Point& p1, const Point& p2)
{
    return Point(p1.x + p2.x, p1.y + p2.y, 0);
}

  class Triangle {
      Point p1, p2, p3;
      findAdjacentTriangles(const std::vector<Triangle>& triangles);

  // Function to calculate the area of a triangle
  double calculateArea(const Point& p1, const Point& p2, const Point& p3) {
      return 0.5 * ((p1.x - p3.x) * (p2.y - p3.y) - (p1.y - p3.y) * (p2.x - p3.x));
  }

  // Function to find the triangle with the maximum area usando l'approccio divide et impera
  Triangle findMaxAreaTriangle(const std::vector<Point>& points) ;
// ordina i vertici di un triangolo in senso antiorario
  bool SortVertices(const Point& p1,
                    const Point& p2,
                    const Point& p3);


  // Function to check if a point is inside a triangle
  Point getCircleCenter();
  bool isPointInsideSquare(const Point& Q, const Triangle& T);
  bool isInsideCircumcircle(const Point& point, const Triangle& triangle);
  bool isPointInsideTriangle(const Point& Q, const Triangle& T);
  double calculateAngle(const Point& A, const Point& B, const Point& C);
  bool verifyDelaunayHypothesis(const Point& pointA, const Point& pointB, const Point& pointC, const Point& pointD);


bool isPointInPolygon(const Point& Q, const vector<Triangle>& triangulation);



// Triangulation data structure
class Triangulation {
public:
    // Constructor
    Triangulation(const std::vector<Triangle>& initialTriangles);
    buildAdjacency();
    updateAdjacency();




private:
    std::vector<Triangle> triangles;
    std::map<Triangle, std::vector<Triangle>> adjacency;


};




#endif // __DELAUNAY_H
