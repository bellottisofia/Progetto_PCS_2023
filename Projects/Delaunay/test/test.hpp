#ifndef __TEST_H
#define __TEST_H

#include <gtest/gtest.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include "delaunay.hpp"

using namespace testing;
using namespace std;
using namespace ProjectLibrary;

// Test case per la funzione calculateArea
TEST(DelaunayTest, calculateArea) {
    // Crea un vettore di punti
    vector<Point> points = {
        Point(0, 0, 0),
        Point(1, 1, 1),
        Point(1, 0, 2),
        Point(3, 3, 3),
        Point(4, 4, 4)
    };

    // Crea un oggetto Triangle utilizzando i primi tre punti del vettore
    Triangle triangle(points[0], points[1], points[2],1);

    // Verifica se l'area calcolata corrisponde al valore atteso di 0.5
    EXPECT_EQ(triangle.calculateArea(), 0.5);
}
/*
// Test case per la funzione IsVerticesSort
TEST(DelaunayTest, IsVerticesSort) {
    // Crea un vettore di punti
    vector<Point> points = {
        Point(0, 0, 0),
        Point(1, 1, 1),
        Point(1,0,2)

    };

    // Crea un oggetto Triangle utilizzando i primi tre punti del vettore
    Triangle triangle(points[0], points[1], points[2]);


    // Verifica se i vertici sono ordinati in senso antiorario
    EXPECT_EQ(IsVerticesSort(points[0], points[1], points[2]), 0);

}
/*
// Test case per la funzione SortVertices
TEST(DelaunayTest, SortVertices) {
    // Crea un vettore di punti
    vector<Point> points = {
        Point(0, 0, 0),
        Point(1, 1, 1),
        Point(1,0,2)

    };

    // Crea un oggetto Triangle utilizzando i primi tre punti del vettore
    Triangle triangle(points[0], points[1], points[2]);
    triangle.SortVertices();


    // Verifica se i vertici sono ordinati in senso antiorario
    EXPECT_EQ(triangle.IsVerticesSort(), 1);

}
*/

// Test case per la funzione SortVertices
TEST(DelaunayTest, SortVertices) {
    // Crea un vettore di punti
    vector<Point> points = {
        Point(0, 0, 0),
        Point(1, 1, 1),
        Point(1,0,2)

    };
    Triangle t(points[0],points[1],points[2],1);


    EXPECT_TRUE(t.IsVerticesSort());
    // Crea un oggetto Triangle utilizzando i primi tre punti del vettore
    t.SortVertices();
    // Verifica se i vertici sono ordinati in senso antiorario
    EXPECT_TRUE(t.IsVerticesSort());

}

// Test case per la funzione isVertexShared
TEST(DelaunayTest, isVertexShared) {
    // Crea un vettore di punti
    vector<Point> points = {
        Point(0, 0, 0),
        Point(1, 1, 1),
        Point(1,0,2)

    };

    // Crea un oggetto Triangle utilizzando i primi tre punti del vettore
    Triangle triangle(points[0], points[1], points[2],1);
    vector<Point> vertex={Point(1,1,1),Point(3,3,3),Point(1,0,2)};


    // Verifica se il vertcice fornito è uno dei vertici del triangolo
    EXPECT_EQ(triangle.isVertexShared(vertex[0]), 1);
    EXPECT_EQ(triangle.isVertexShared(vertex[1]), 0);


}
// Test case per la funzione isAdjacent(const Triangle& other)
TEST(DelaunayTest, isAdjacent) {
    // Crea un vettore di punti
    vector<Point> points = {
        Point(0, 0, 0),
        Point(1, 1, 1),
        Point(1,0,2)

    };

    // Crea un oggetto Triangle utilizzando i primi tre punti del vettore
    Triangle triangle(points[0], points[1], points[2],1);
    // Creazione di un altro triangolo adiacente
        vector<Point> adjacentVertex = {Point(1, 1, 1), Point(2, 2, 2), Point(1, 0, 2)};
        Triangle adjacentTriangle(adjacentVertex[0], adjacentVertex[1], adjacentVertex[2],2);

        // Verifica se i triangoli sono adiacenti
        EXPECT_EQ(triangle.isAdjacent(adjacentTriangle), true);
}
/*
// Test case per la funzione addAdjacentTriangles(std::vector<Triangle>& triangles)
TEST(DelaunayTest, addAdjacentTriangles) {
    // Crea un vettore di punti
    vector<Point> points = {
        Point(0, 0, 0),
        Point(1, 1, 1),
        Point(1,0,2),
         Point(2,2,3),
         Point(2,0,4),
         Point(0,3,5),
         Point(0,4,6),
         Point(4,0,7),
         Point(5,1,8)

    };

    // Crea un oggetto Triangle utilizzando i primi tre punti del vettore
    Triangle triangle(points[0], points[1], points[2]);
    std::vector<Triangle*> neighbors;
    Triangle x(points[0], points[4], points[7]);
    Triangle y(points[2], points[3], points[8]);
    Triangle z(points[0], points[1], points[4]);
    vector<Triangle> triangles;
    triangles.push_back(x);
    triangles.push_back(y);
    triangles.push_back(z);

    triangle.addAdjacentTriangles(triangles);
    triangle.addAdjacentTriangle(z);

    triangle.getAdjacentTriangles();
    //triangle.printAdjacentTriangles();




        // Verifica la dimensione dei triangoli adiacenti
        ASSERT_EQ(triangle.getAdjacentTriangles().size(), 1);

        // Verifica che i triangoli adiacenti siano stati correttamente aggiunti
        ASSERT_EQ(triangle.neighbors[0], &z);

    // Verifica che i triangoli adiacenti siano stati correttamente aggiunti
       ASSERT_EQ(triangle.neighbors.size(), 1);

      // ASSERT_NEQ(triangle.getAdjacentTriangles()[1], &x);  // non dovrebbe avere x come triangolo adiacente

    }

/*
TEST(CalculateAngleTest, CorrectAngle) {
    Point p1(0.0, 0.0, 1);
    Point p2(1.0, 0.0, 2);
    Point p3(0.0, 1.0, 3);

    Triangle triangle(p1, p2, p3);

    double expectedAngle = 45;

    double angledegrees = triangle.calculateAngle(p1, p2, p3);

    EXPECT_DOUBLE_EQ(angledegrees, expectedAngle);
}
*/
/*


TEST(DelaunayTest, addAdjacentTriangles2) {
    // Crea un vettore di punti
    std::vector<Point> points = {
        Point(0, 0, 0),
        Point(1, 1, 1),
        Point(1, 0, 2),
        Point(2, 2, 3),
        Point(2, 0, 4),
        Point(0, 3, 5),
        Point(0, 4, 6),
        Point(4, 0, 7),
        Point(5, 1, 8)
    };

    // Crea un oggetto Triangle utilizzando i primi tre punti del vettore
    Triangle triangle(points[0], points[1], points[2]);

    // Crea altri triangoli
    Triangle x(points[0], points[1], points[7]);
    Triangle y(points[2], points[3], points[8]);
    Triangle z(points[0], points[1], points[4]);

    // Crea un vettore di puntatori ai triangoli
    std::vector<Triangle> triangles;
    triangles.push_back(x);
    triangles.push_back(y);
    triangles.push_back(z);

    // Aggiungi i triangoli adiacenti a triangle
    triangle.addAdjacentTriangles(triangles);

    // Verifica la lunghezza dei triangoli adiacenti
    int expectedSize = 3;
    int actualSize = triangle.getAdjacentTriangles().size();
    EXPECT_EQ(actualSize, expectedSize);

    // Stampa i triangoli adiacenti
    triangle.printAdjacentTriangles();
}




TEST(DelaunayTest, createTriangularMesh) {
    // Crea un vettore di punti
    std::vector<Point> points = {
        Point(0, 0, 0),
        Point(1, 1, 1),
        Point(1, 0, 2),
        Point(2, 2, 3),
        Point(2, 0, 4),
        Point(0, 3, 5),
        Point(0, 4, 6),
        Point(4, 0, 7),
        Point(5, 1, 8)
    };

    // Crea la mesh triangolare
    std::vector<Triangle> mesh ;
    Triangle::createTriangularMesh(points);

    // Verifica la dimensione della mesh
    int expectedSize = 12;  // Numero di triangoli per una mesh completa di 9 punti
    int actualSize = mesh.size();
    ASSERT_EQ(actualSize, expectedSize);

    // Verifica l'adiacenza dei triangoli
    for (const Triangle& triangle : mesh) {
        for (const Triangle* adjacent : triangle.getAdjacentTriangles()) {
            // Verifica che il triangolo corrente sia adiacente a tutti i suoi triangoli adiacenti
            ASSERT_TRUE(adjacent->isAdjacent(triangle));
        }
    }
}
*/
// Test case per la funzione findMaximumTriangleArea
TEST(DelaunayTest, FindMaximumTriangleArea) {
    std::vector<Point> points = {
                Point(0, 0, 0),
                Point(1, 1, 1),
                Point(2, 2, 2),
                Point(3, 3, 3),
                Point(4, 4, 4)
    };

    Triangle triangle(points[0], points[1], points[2],1);
    Triangle maxTriangle = triangle.findMaximumTriangle(points);

    // Assert per verificare il risultato atteso
    EXPECT_EQ(maxTriangle.p1, points[0]);
    EXPECT_EQ(maxTriangle.p2, points[1]);
    EXPECT_EQ(maxTriangle.p3, points[2]);
}


// Test case per verificare se un punto è all'interno della circonferenza circoscritta al triangolo
TEST(TriangleTest, PointInsideCircumcircle) {
    Point p1(0.0, 0.0, 1);
    Point p2(2.0, 0.0, 2);
    Point p3(1.0, 2.0, 3);
    Triangle triangle(p1, p2, p3,1);

    Point pointInside(1.0, 1.0, 4);

    bool isInside = triangle.isInsideCircumcircle(pointInside);

    EXPECT_TRUE(isInside);
}

// Test case per verificare se un punto è all'esterno della circonferenza circoscritta al triangolo
TEST(TriangleTest, PointOutsideCircumcircle) {
    Point p1(0.0, 0.0, 1);
    Point p2(2.0, 0.0, 2);
    Point p3(1.0, 2.0, 3);
    Triangle triangle(p1, p2, p3,1);

    Point pointOutside(0.0, 2.0, 4);

    bool isInside = triangle.isInsideCircumcircle(pointOutside);

    EXPECT_FALSE(isInside);
}


// Test per la funzione isPointInsideTriangle
TEST(IsPointInsideTriangleTest, PointInside) {
    Point p1(0.0, 0.0, 1);
    Point p2(1.0, 0.0, 2);
    Point p3(0.0, 1.0, 3);

    Triangle triangle(p1, p2, p3,1);


    Point Q(0.5, 0.5, 4);

    bool result = triangle.isPointInsideTriangle(Q);

        EXPECT_TRUE(result);

}

TEST(IsPointInsideTriangleTest, PointOutside) {
    Point p1(0.0, 0.0, 1);
    Point p2(0.0, 1.0, 2);
    Point p3(1.0, 0.0, 3);

    Triangle triangle(p1, p2, p3,1);

    Point Q(0.0, 2.0, 4);

    bool result = triangle.isPointInsideTriangle(Q);

    EXPECT_FALSE(result);
}





 // Verifica se il punto q si trova sul segmento p-r

TEST(isPointOnSegment, PointOUT) {
    Point p1(2.0, 3.0, 1);
    Point p2(3.0, 0.0, 2);
    Point p3(1.0, 0.0, 3);

    EXPECT_FALSE(p1.isPointOnSegment(p2, p3));

}

TEST(isPointOnSegment, PointON) {
    Point p1(2.0, 0.0, 1);
    Point p2(3.0, 0.0, 2);
    Point p3(1.0, 0.0, 3);

    EXPECT_TRUE(p1.isPointOnSegment(p2, p3));

}

TEST(doSegmentIntersect, no) {

    Point p1(0.0, 0.0, 1);
    Point p2(1.0, 0.0, 2);
    Point p3(0.0, 1.0, 3);
    Point p4(1.0, 1.0, 6);
    EXPECT_FALSE(p1.doSegmentsIntersect( p2,  p3,  p4));
}
TEST(doSegmentIntersect, yes) {

    Point p1(0.0, 0.0, 1);
    Point p2(1.0, 2.0, 2);
    Point p3(0.0, 1.0, 3);
    Point p4(1.0, 0.0, 6);
    EXPECT_TRUE(p1.doSegmentsIntersect( p2,  p3,  p4));
}

   // bool Triangle::findIntersection(const Point& q)

TEST(findIntersection, yes) {

    Point p1(0.0, 0.0, 1);
    Point p2(1.0, 2.0, 2);
    Point p3(0.0, 1.0, 3);
    Point p4(1.0, 0.0, 6);
    Triangle triangle(p1,p2,p3,1);

  EXPECT_TRUE(p3.isPointInVector(triangle.findIntersection( p4)));
    // devo scrivere p4 in
}


TEST(AreTrianglesDelaunayTest, TrianglesSatisfyCriterion) {
  // Create triangles that satisfy the Delaunay criterion
    Point p1(0.0, 0.0, 1);
    Point p2(1.0, 0.0, 2);
    Point p3(0.0, 1.0, 3);
    Point p4(1.0, 1.0, 6);



    Triangle triangle1(p1, p2, p3,1);
    Triangle triangle2(p1, p2, p4,2);

  // Check if the triangles satisfy the Delaunay criterion
    EXPECT_TRUE(triangle1.areTrianglesDelaunay(triangle2));
    EXPECT_TRUE(triangle2.areTrianglesDelaunay(triangle1));
}

// Test case for when a triangle does not satisfy the Delaunay criterion
TEST(AreTrianglesDelaunayTest, TriangleDoesNotSatisfyCriterion) {
  // Create triangles where one does not satisfy the Delaunay criterion
    Point p1(0.0, 0.0, 1);
    Point p2(1.0, 0.0, 2);
    Point p3(0.0, 1.0, 3);
    Point p4(2.0, 0.0, 6);



    Triangle triangle1(p1, p2, p3,1);
    Triangle triangle3(p1, p2, p4,2);
  // Check if a triangle doesn't satisfy the Delaunay criterion with another triangle
  EXPECT_FALSE(triangle1.areTrianglesDelaunay(triangle3));
}
TEST(AreTrianglesDelaunayTest, Triangles_not_adjacent) {
  // Create triangles where one does not satisfy the Delaunay criterion
    Point p1(0.0, 0.0, 1);
    Point p2(1.0, 0.0, 2);
    Point p3(0.0, 1.0, 3);
    Point p4(15.0,7.0, 6);
    Point p5(9.0, 7.0, 5);
    Point p6(20.0, 14.0, 5);



    Triangle triangle1(p1, p2, p3,1);
    Triangle triangle3(p3, p4, p5,2);
  // Check if a triangle doesn't satisfy the Delaunay criterion with another triangle
  EXPECT_EQ(triangle1.areTrianglesDelaunay(triangle3),-1);
}
// Test case per la funzione flip: vertici scambiati correttamente
TEST(TriangleTest, FlipVerticesSwapped)
{
    // Crea i triangoli di prova
    Point p1(0, 0, 1);
    Point p2(1, 0, 2);
    Point p3(0, 1, 3);
    Triangle triangle1(p1, p2, p3,1);

    Point p4(0.7, 0.7, 4);
    Triangle triangle2(p2, p4, p3,2);

    triangle1.flip(triangle2);

    // Verifica che i vertici siano stati scambiati correttamente
    EXPECT_EQ(triangle1.p1, p1);
    EXPECT_EQ(triangle1.p2, p2);
    EXPECT_EQ(triangle1.p3, p4);

    EXPECT_EQ(triangle2.p1, p1);
    EXPECT_EQ(triangle2.p2, p4);
    EXPECT_EQ(triangle2.p3, p3);
}
//punto esterno
TEST(connettiPuntoEsterno, interseca) {
  // Create triangles where one does not satisfy the Delaunay criterion
    Point p1(0.0, 0.0, 1);
    Point p2(1.0, 0.0, 2);
    Point p3(0.0, 1.0, 3);
    Point Q(1.0, 1.0, 6);



    Triangle triangle1(p1, p2, p3,1);
    std::vector<Triangle> triangles;
    triangles.push_back(triangle1);
    Triangulation Delunay;
    Delunay.DelunayTriangles=triangles;


    Delunay.addPointToTriangulation(Q);

    Triangle triangle(p2, Q, p3,2);
    EXPECT_TRUE(triangle.isAdjacent(triangle1));


std::vector<Triangle> ciao;
ciao.push_back(triangle1);
ciao.push_back(triangle);
Delunay.print();
  EXPECT_EQ(Delunay,ciao);
}
#endif // __TEST_H
