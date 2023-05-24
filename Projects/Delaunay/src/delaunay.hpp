#ifndef __DELAUNAY_H
#define __DELAUNAY_H

#include <iostream>
#include "list"
#include "Eigen/Eigen"
#include "map"
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
      bool isVertexShared(const Point& vertex); //testato
      bool isAdjacent(const Triangle& other);  //testato
      void addAdjacentTriangle(const Triangle* triangle);
      void addAdjacentTriangles(std::vector<Triangle>& triangles);
      std::vector<Triangle*>& getAdjacentTriangles();
      void createTriangularMesh(std::vector<Triangle>& triangles) {
          for (size_t i = 0; i < triangles.size(); ++i) {
              Triangle& currentTriangle = triangles[i];

              // Controlla l'adiacenza con i triangoli precedenti
              for (size_t j = 0; j < i; ++j) {
                  Triangle& previousTriangle = triangles[j];
                  if (currentTriangle.isAdjacent(previousTriangle)) {
                      currentTriangle.addAdjacentTriangle(&previousTriangle);
                      previousTriangle.addAdjacentTriangle(&currentTriangle);
                  }
              }
          }
      }

     double calculateAngle(const Point& A, const Point& B, const Point& C);

      void printAdjacentTriangles() ;
      /*
      // Funzione per confrontare i puntatori ai triangoli
      bool operator==(const Triangle& other) const {
          return (p1 == other.p1 && p2 == other.p2 && p3 == other.p3) ||
                 (p1 == other.p2 && p2 == other.p3 && p3 == other.p1) ||
                 (p1 == other.p3 && p2 == other.p1 && p3 == other.p2);
      }
      */

      // Function to calculate the area of a triangle
     // double calculateArea(const Point& p1, const Point& p2, const Point& p3);
      double calculateArea(); //testato
      // ordina i vertici di un triangolo in senso antiorario
       bool IsVerticesSort();  //testato
       void SortVertices();

/*
    std::unordered_map<const Triangle*, std::vector<const Triangle*>> adjacencyMap;
};

// Durante la creazione della mesh triangolare, quando crei un nuovo triangolo, aggiorna le informazioni di adiacenza
void createTriangularMesh(std::vector<Triangle>& triangles) {
    for (size_t i = 0; i < triangles.size(); ++i) {
        Triangle& currentTriangle = triangles[i];

        // Controlla l'adiacenza con i triangoli precedenti
        for (size_t j = 0; j < i; ++j) {
            Triangle& previousTriangle = triangles[j];
            if (currentTriangle.isAdjacent(previousTriangle)) {
                currentTriangle.addAdjacentTriangle(&previousTriangle);
                previousTriangle.addAdjacentTriangle(&currentTriangle);
            }
        }
    }
}
Nell'esempio sopra, Triangle contiene un membro adjacencyMap che è un std::unordered_map. La chiave della mappa è un puntatore al triangolo corrente, e il valore è un vettore di puntatori ai triangoli adiacenti.

Durante la creazione della mesh triangolare, la funzione createTriangularMesh itera attraverso tutti i triangoli e controlla l'adiacenza con i triangoli precedenti. Se due triangoli sono adiacenti, vengono aggiunti come adiacenti l'uno all'altro utilizzando la funzione addAdjacentTriangle.

Assicurati di implementare correttamente la classe Triangle con le funzioni isAdjacent e addAdjacentTriangle, in modo da gestire correttamente l'aggiornamento delle informazioni di adiacenza.






*/






  // Function to find the triangle with the maximum area usando l'approccio divide et impera
  Triangle findMaximumTriangleArea(const std::vector<Point>& points, int start, int end);
  Triangle findMaximumTriangle(const std::vector<Point>& points);
  friend std::ostream& operator<<(std::ostream& os, const Triangle& triangle);

};
  inline std::ostream& operator<<(std::ostream& os, const Triangle& triangle) {
      os << "Triangolo adiacente: " << triangle.p1 << ", " << triangle.p2 << ", " << triangle.p3;
      return os;
  }
}
/*
// Triangulation data structure
class Triangulation {
public:
    // Constructor
    Triangulation(const std::vector<Triangle>& initialTriangles);
    buildAdjacency();
    updateAdjacency();
    double calculateAngle(const Point& A, const Point& B, const Point& C);
    bool verifyDelaunayHypothesis(const Point& pointA, const Point& pointB, const Point& pointC, const Point& pointD);





private:
    std::vector<Triangle> triangles;


};


*/


#endif // __DELAUNAY_H
