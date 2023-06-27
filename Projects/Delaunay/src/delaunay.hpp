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
        Point(): x(0.0),y(0.0),id(0){}

        Point(const double& x,
              const double& y,
              const unsigned int& id);
        Point(const Point& p);


        Point& operator=(const Point& other);
        double calculateAngle( const Point& B, const Point& C);//funge
        bool isPointOnSegment( const Point& q, const Point& r)const;
        //int calculateOrientation(const Point& p1, const Point& p2);
        //int calculateOrientation(const Point& p1, const Point& p2, const Point& p3);
        bool arePointsCollinear(const Point& p2,const Point& p3);
        bool doSegmentsIntersect( const Point& p2, const Point& p3, const Point& p4);
        bool isPointInVector( const std::vector<Point>& points)const;




    };

int calculateOrientation(const Point& p1, const Point& p2, const Point& p3);
void SortVertices( Point& p1, Point& p2, Point& p3);
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
  os << "id: "<< p2.id << " x: "<< p2.x<< " y: "<< p2.y;
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

struct MinMax {
    double minX;
    double maxX;
    double minY;
    double maxY;

};
MinMax findMinMax(const std::vector<Point>& points);
  class Triangle {
      public:
      Point p1, p2, p3;
      unsigned int id;

      Triangle() = default;
      Triangle(const Point& point1, const Point& point2, const Point& point3, const unsigned int& id); //testato
      bool isInsideCircumcircle(const Point& point)const; //testato
      bool isPointInsideTriangle(const Point& Q)const; //testato

      std::vector<Point> getVertices() const ;
      const Point& getVertexA() const ;
      const Point& getVertexB() const ;
      const Point& getVertexC() const ;
      bool isPointOnEdge(const Point& Q) const;


      bool doSegmentsIntersect(const Point& Q);

      bool isVertexShared(const Point& vertex); //testato
      bool isAdjacent(const Triangle& other);  //testato


      double calculateArea()const; //testato
      // ordina i vertici di un triangolo in senso antiorario
       bool IsVerticesSort();  //testato
       void SortVertices();

  // Function to find the triangle with the maximum area usando l'approccio divide et impera
  Triangle findMaximumTriangleArea(const std::vector<Point>& points, int start, int end);
  Triangle findMaximumTriangle(const std::vector<Point>& points);
  bool checkTriangleOverlap(Triangle& other);
  bool hasSharedPartialEdge(const Triangle& other);
  std::vector<Point> isPointOnEdge(const Point& Q);
  int areTrianglesDelaunay(Triangle& triangle1);
  std::vector<Point> findIntersection(const Point& q);
  void flip(Triangle& triangle1);
  std::vector<Point> PointOnEdge(const Point& Q);

  bool operator!=(const Triangle& other) const {
      // Confronta i vertici dei triangoli
      return !(p1 == other.p1 && p2 == other.p2 && p3 == other.p3);
  }

  bool operator==(const Triangle& other) const {
      // Confronta i vertici dei triangoli
      return (p1 == other.p1 && p2 == other.p2 && p3 == other.p3);
  }
  void print() const {
          std::cout << "Triangle " << id << ": ";
          std::cout << p1 << ", ";
          std::cout << p2 << ", ";
          std::cout << p3 << std::endl;
      }
};


        class Triangulation {
        private:
            unsigned int maxTriangleId = 0;


        public:
            std::vector<Triangle> DelunayTriangles;
            std::unordered_map<unsigned int, std::vector<unsigned int>> adjacencyList;
            Triangulation() = default;
            Triangulation(const std::vector<Triangle>& triangles);


            void addTriangle(const Triangle& triangle);

            void addAdjacentTriangle(const unsigned int& triangleId,const unsigned int& adjacentTriangleId) ;

            const std::vector<unsigned int>& getAdjacentTriangles(const unsigned int& triangleId) ;

            int PointInsideTriangulation(const Point& Q);
            void createSubtriangulation(const Point& Q, const unsigned int& triangleId);
            void connectPointToVertices(const Point& Q);
            vector<int> connectPointOnEdge(const Triangle& t, const Point& Q, const vector<Point>& edge);
            void connectPointOnEdgeToTriangulation(const Triangle& triangle, const Point& Q, const vector<Point>& edge, const MinMax& minMax);
            void flipTrianglesIfNotDelaunay();
            void flipAndUpdateAdjacency(Triangle& triangle, Triangle& adjacentTriangle);
            void updateAdjacency(Triangle& triangle, Triangle& adjacentTriangle);
            bool areAllTrianglesDelaunay();
            bool isBoundaryTriangle(const Triangle& triangle);
            void addPointToTriangulation(const Point& Q, const MinMax& minMax);
            Triangle findTriangleById(const unsigned int& triangleId);
            void removeAdjacentTriangle(unsigned int triangleId1, unsigned int triangleId2);
            bool isIdInAdjacency(unsigned int triangleId,vector<unsigned int> adjacency);
            bool isBoundaryEdge(const Triangle& triangle,const vector<Point>& edge, const MinMax& minMax);

            bool isBoundaryPoint(const Point& point, const MinMax& minMax);
            inline unsigned int getMaxTriangleId() const {
                    return maxTriangleId;
                }

                inline void incrementMaxTriangleId(unsigned int increment = 1) {
                    maxTriangleId += increment;
                }


            bool operator==(const std::vector<Triangle>& other) const {
                    // Verifica le dimensioni dei vettori
                    if (DelunayTriangles.size() != other.size()) {
                        return false;
                    }

                    // Confronta gli elementi dei vettori
                    for (size_t i = 0; i < DelunayTriangles.size(); i++) {
                        if (DelunayTriangles[i] != other[i]) {
                            return false;
                        }
                    }

                    return true;
                }

            void print() const {
                    for (const auto& triangle : DelunayTriangles) {
                        triangle.print();
                    }
                }


};

Triangulation DelunayTriangulation(const std::vector<Point>& points, const MinMax& minMax);
void findMinMax(const std::vector<Point>& sortedPoints, double& minX, double& maxX, double& minY, double& maxY);
}

#endif // __DELAUNAY_H
