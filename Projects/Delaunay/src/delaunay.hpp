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
        bool isPointOnSegment( const Point& q, const Point& r);
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


  class Triangle {
      public:
      Point p1, p2, p3;
      unsigned int id;
      std::vector<unsigned int> adjacentTriangles;


      Triangle() = default;
      Triangle(const Point& point1, const Point& point2, const Point& point3, const unsigned int& id); //testato
      bool isInsideCircumcircle(const Point& point)const; //testato
      bool isPointInsideTriangle(const Point& Q)const; //testato

      std::vector<Point> getVertices() const ;
      const Point& getVertexA() const ;
      const Point& getVertexB() const ;
      const Point& getVertexC() const ;



      bool doSegmentsIntersect(const Point& Q);

      bool isVertexShared(const Point& vertex); //testato
      bool isAdjacent(const Triangle& other);  //testato


      double calculateArea()const; //testato
      // ordina i vertici di un triangolo in senso antiorario
       bool IsVerticesSort();  //testato
       void SortVertices();

  // Function to find the triangle with the maximum area usando l'approccio divide et impera
  //Triangle findMaximumTriangleArea(const std::vector<Point>& points, int start, int end);
  //Triangle findMaximumTriangle(const std::vector<Point>& points);
  Triangle findMaximumTriangleArea(const std::vector<Point>& points);
  bool areTrianglesDelaunay(Triangle& triangle1);
  std::vector<Point> findIntersection(const Point& q);
  void flip(Triangle& triangle1);
};

/*
class Triangulation {
public:

    std::vector<Triangle> DelunayTriangles;
    std::unordered_map<int, std::vector<int>> adjacencyList;
    Triangulation(const std::vector<Triangle>& DelunayTriangles);



       void addTriangle(const Triangle& triangle) {
           triangles.push_back(triangle);
           adjacencyList[triangle.id] = {}; // Inizializza la lista di adiacenza per il triangolo
       }

       void addAdjacentTriangle(int triangleId, int adjacentTriangleId) {
           adjacencyList[triangleId].push_back(adjacentTriangleId);
           adjacencyList[adjacentTriangleId].push_back(triangleId);
       }

       const std::vector<int>& getAdjacentTriangles(int triangleId) {
           return adjacencyList[triangleId];
       }
        void addPointToTriangulation(const Point&){}



        void updateAdjacency(const Triangle& DelunayTriangles){}
*/
        class Triangulation {
        public:
            std::vector<Triangle> DelunayTriangles;
            std::unordered_map<unsigned int, std::vector<unsigned int>> adjacencyList;

            Triangulation(const std::vector<Triangle>& DelunayTriangles);

            void addTriangle(const Triangle& triangle);

            void addAdjacentTriangle(int triangleId, int adjacentTriangleId) ;

            const std::vector<unsigned int>& getAdjacentTriangles(int triangleId) ;

            int PointInsideTriangulation(const Point& Q);
            void createSubtriangulation(const Point& Q, int triangleId);
            void connectPointToVertices(const Point& Q);
            void addPointToTriangulation(const Point&);

            void updateAdjacency(const std::vector<Triangle>& DelunayTriangles);/* {
                this->DelunayTriangles = DelunayTriangles;
                buildAdjacencyList();
            }
*/
        private:
            void buildAdjacencyList();


};



}

#endif // __DELAUNAY_H
