#ifndef __DELAUNAY_H
#define __DELAUNAY_H

#include <iostream>
#include "list"
#include "Eigen/Eigen"
#include "map"
#include <unordered_map>

#include<vector>
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


class Point {
public:

        double x;
        double y;
        unsigned int id;

        Point *succ = nullptr;
        Point *prec = nullptr;

        static constexpr double geometricTol = 1.0e-12;
        static constexpr double geometricTol_Squared = max_tolerance(Point::geometricTol * Point::geometricTol,
                                                                     numeric_limits<double>::epsilon());

        Point(const double& x,
              const double& y,
              const unsigned int& id);

        Point(const Point& p);
        Point& operator=(const Point& other);
        double calculateAngle( const Point& B, const Point& C);//funge
        bool isPointOnSegment( const Point& q, const Point& r);
        int calculateOrientation(const Point& p1, const Point& p2);


    };


inline double normSquared(const double& x, const double& y)
{
  return x * x + y * y;
}

inline bool operator==(const Point& p1, const Point& p2)
{
  return (normSquared(p1.x - p2.x, p1.y - p2.y) <=
          Point::geometricTol * Point::geometricTol *
          max(normSquared(p1.x, p1.y), normSquared(p2.x, p2.y)));
}

inline bool operator!=(const Point& p1, const Point& p2)
{
  return !(p1 == p2);
}

inline ostream& operator<<(ostream& os, const Point& p2)
{
  os << p2.id;
  return os;
}

inline bool operator>(const Point& p1, const Point& p2)
{
  return p1.x > p2.x + Point::geometricTol * max(p1.x, p2.x);
}

inline bool operator<=(const Point& p1, const Point& p2)
{
  return !(p1 > p2);
}



inline Point operator-(const Point& p1, const Point& p2)
{
    return Point(p1.x - p2.x, p1.y - p2.y, 0);
}

inline Point operator+(const Point& p1, const Point& p2)
{
    return Point(p1.x + p2.x, p1.y + p2.y, 0);
}


  class Triangle {
      public:
      Point p1, p2, p3;
      std::vector<Triangle*> neighbors;



      Triangle(const Point& point1, const Point& point2, const Point& point3); //testato
      bool isInsideCircumcircle(const Point& point); //testato

      bool isPointInsideTriangle(const Point& Q); //testato

      bool doSegmentsIntersect(const Point& Q);

      bool isVertexShared(const Point& vertex); //testato
      bool isAdjacent(const Triangle& other);  //testato


      double calculateArea(); //testato
      // ordina i vertici di un triangolo in senso antiorario
       bool IsVerticesSort();  //testato
       void SortVertices();







  // Function to find the triangle with the maximum area usando l'approccio divide et impera
  Triangle findMaximumTriangleArea(const std::vector<Point>& points, int start, int end);
  Triangle findMaximumTriangle(const std::vector<Point>& points);
};


class Triangulation {
public:

    std::vector<Triangle> DelunayTriangles;
    std::unordered_map<const Triangle*, std::vector<const Triangle*>> adjacencyList;
    Triangulation(const std::vector<Triangle>& DelunayTriangles);

        void addTriangle(const Triangle& DelunayTriangles);
        void addPointToTriangulation(const Point& Q);

        const std::vector<Triangle>& getTriangles() const ;

        const std::vector<const Triangle*>& getAdjacentTriangles(const Triangle& DelunayTriangles) const;

        void updateAdjacency(const Triangle& DelunayTriangles);



    bool verifyDelaunayHypothesis(const Point& pointA, const Point& pointB, const Point& pointC, const Point& pointD);





};





#endif // __DELAUNAY_H
