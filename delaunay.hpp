#ifndef __DELAUNAY_H
#define __DELAUNAY_H

#include <iostream>
#include <list>
#include <Eigen/Eigen>
#include <map>
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

        Point(const double& x, const double& y, const unsigned int& id) : x(x), y(y), id(id) {}

        Point(const Point& p) : x(p.x), y(p.y), id(p.id) {}
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

    inline std::ostream& operator<<(std::ostream& os, const Point& p2)
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
        Triangle(const Point& point1, const Point& point2, const Point& point3)
            : p1(point1), p2(point2), p3(point3) {}

        // TROVO IL TRIANGOLO DI AREA MASSIMA
        double calculateArea(const Point& p1, const Point& p2, const Point& p3);
        Triangle findMaximumTriangleArea(const std::vector<Point>& points, int start, int end);
        Triangle findMaximumTriangle(const std::vector<Point>& points);

        // CAPISCO SE IL PUNTO INSERITO è INTERNO ALLA CIRCONFERENZA CIRCOSCRITTA
        bool isInsideCircumcircle(const Point& Q);

        // CAPISCO SE UN PUNTO è INTERNO AL TRIANGOLO
        double calculateAngle(const Point& p1, const Point& p2, const Point& p3);
        bool isPointInsideTriangle(const Point& Q);

        std::vector<Triangle> findAdjacentTriangles(const std::vector<Triangle>& triangles);
        bool sortVertices(const Point& p1, const Point& p2, const Point& p3);
        Point getCircleCenter();


        bool verifyDelaunayHypothesis(const Point& pointA, const Point& pointB, const Point& pointC,
                                      const Point& pointD);
        bool isPointInPolygon(const Point& Q, const std::vector<Triangle>& triangulation);
    };

    class Triangulation {
    public:
        Triangulation(const std::vector<Triangle>& initialTriangles);
        void buildAdjacency();
        void updateAdjacency();

    private:
        std::vector<Triangle> triangles;
        std::map<Triangle, std::vector<Triangle>> adjacency;
    };
}

#endif // __DELAUNAY_H
