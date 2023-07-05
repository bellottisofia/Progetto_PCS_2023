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


struct Point {


        double x;
        double y;
        unsigned int id;

        static constexpr double geometricTol = 1.0e-12;
        static constexpr double geometricTol_Squared = max_tolerance(Point::geometricTol * Point::geometricTol,
                                                                     numeric_limits<double>::epsilon());
        Point(): x(0.0),y(0.0),id(0){}

        Point(const double& x,
              const double& y,
              const unsigned int& id);
        Point(const Point& p);


        Point& operator=(const Point& other);
        double calculateAngle( const Point& B, const Point& C);
        bool isPointOnSegment( const Point& q, const Point& r)const;
        bool arePointsCollinear(const Point& p2,const Point& p3);
        bool doSegmentsIntersect( const Point& p2, const Point& p3, const Point& p4)const;
        bool isPointInVector( const std::vector<Point>& points)const;




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


  class Triangle {
      mutable Point p1, p2, p3;
      mutable unsigned int id;
  public:
      Triangle() = default;
      Triangle(const Point& point1, const Point& point2, const Point& point3, const unsigned int& id);
      bool isInsideCircumcircle(const Point& point)const;
      bool isPointInsideTriangle(const Point& Q)const;

      //std::vector<Point> getVertices() const ;
      const Point getVertex1() const  ;
      void setVertex1(const Point& newVertex1);
      const Point getVertex2() const ;
      void setVertex2(const Point& newVertex2) ;
      const Point getVertex3() const ;
      void setVertex3(const Point& newVertex3) ;
       unsigned int getId() const ;
      void setId(unsigned int newId) const;
      bool isPointOnEdge(const Point& Q)const;
      friend class Triangulation;



      bool isVertexShared(const Point& vertex) const;
      bool isAdjacent(const Triangle& other);


      double calculateArea()const;
      // ordina i vertici di un triangolo in senso antiorario
       bool IsVerticesSort();
       void SortVertices();

  // Function to find the triangle with the maximum area usando l'approccio divide et impera
  Triangle findMaximumTriangleArea(const std::vector<Point>& points, int start, int end);
  Triangle findMaximumTriangle(const std::vector<Point>& points);
  bool checkTriangleOverlap(Triangle& other);
  bool hasSharedPartialEdge(const Triangle& other);

  int areTrianglesDelaunay(Triangle& triangle1);
  std::vector<Point> findIntersection(const Point& q);
  void flip(Triangle& triangle1);
  std::vector<Point> PointOnEdge(const Point& Q);

  bool operator!=(const Triangle& other) const {
    // Confronta i vertici dei triangoli
    unsigned int shared = 0;
    if (isVertexShared(other.p1)) shared++;
    if (isVertexShared(other.p2)) shared++;
    if (isVertexShared(other.p3)) shared++;
    return shared < 3;
  }

  bool operator==(const Triangle& other) const {
    // Confronta i vertici dei triangoli
    unsigned int shared = 0;
    if (isVertexShared(other.p1)) shared++;
    if (isVertexShared(other.p2)) shared++;
    if (isVertexShared(other.p3)) shared++;
    return shared > 2 ;
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
            std::vector<Triangle> DelaunayTriangles;
            std::unordered_map<unsigned int, std::vector<unsigned int>> adjacencyList;
            Triangulation() = default;
            Triangulation(const std::vector<Triangle>& triangles);
            void addTriangle(const Triangle& triangle);

            void addAdjacentTriangle(const unsigned int& triangleId,const unsigned int& adjacentTriangleId) ;

            const std::vector<unsigned int>& getAdjacentTriangles(const unsigned int& triangleId) ;

            int PointInsideTriangulation(const Point& Q);
            void createSubtriangulation(const Point& Q, const unsigned int& triangleId);
            void connectPointToVertices(const Point& Q);
            void connectPointOnEdgeInside(Triangle& t, const Point& Q, const vector<Point>& edge);
            void connectPointOnEdgeToTriangulation(Triangle& triangle, const Point& Q, const vector<Point>& edge);
            void flipTrianglesIfNotDelaunay();
            void flipAndUpdateAdjacency(Triangle& triangle, Triangle& adjacentTriangle);
            void updateAdjacency(const Triangle& triangle, const Triangle& adjacentTriangle);
            bool areAllTrianglesDelaunay();
            bool isBoundaryTriangle(const Triangle& triangle);
            void addPointToTriangulation(const Point& Q);
            void removeAdjacentTriangle(unsigned int triangleId1, unsigned int triangleId2);
            bool isIdInAdjacency(unsigned int triangleId,vector<unsigned int> adjacency);
            bool isBoundaryEdge(const Triangle& triangle, const std::vector<Point>& edge);
            inline unsigned int getMaxTriangleId() const {
                    return maxTriangleId;
                }

            inline void incrementMaxTriangleId(unsigned int increment = 1) {
                maxTriangleId += increment;
            }


            bool operator==(const std::vector<Triangle>& other) const {
                    // Verifica le dimensioni dei vettori
                    if (DelaunayTriangles.size() != other.size()) {
                        return false;
                    }

                    // Confronta gli elementi dei vettori
                    for (size_t i = 0; i < DelaunayTriangles.size(); i++) {
                        if (DelaunayTriangles[i] != other[i]) {
                            return false;
                        }
                    }

                    return true;
                }

            void print() const {
                    for (const auto& triangle : DelaunayTriangles) {
                        triangle.print();
                    }
                }


};

Triangulation DelaunayTriangulation(const std::vector<Point>& points);
}

#endif // __DELAUNAY_H
