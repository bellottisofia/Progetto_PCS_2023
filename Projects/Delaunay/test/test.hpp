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
    Triangle triangle(points[0], points[1], points[2]);

    // Verifica se l'area calcolata corrisponde al valore atteso di 0.5
    EXPECT_EQ(triangle.calculateArea(), 0.5);
}
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
    EXPECT_EQ(triangle.IsVerticesSort(), 0);

}
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
// Test case per la funzione isVertexShared
TEST(DelaunayTest, isVertexShared) {
    // Crea un vettore di punti
    vector<Point> points = {
        Point(0, 0, 0),
        Point(1, 1, 1),
        Point(1,0,2)

    };

    // Crea un oggetto Triangle utilizzando i primi tre punti del vettore
    Triangle triangle(points[0], points[1], points[2]);
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
    Triangle triangle(points[0], points[1], points[2]);
    // Creazione di un altro triangolo adiacente
        vector<Point> adjacentVertex = {Point(1, 1, 1), Point(2, 2, 2), Point(1, 0, 2)};
        Triangle adjacentTriangle(adjacentVertex[0], adjacentVertex[1], adjacentVertex[2]);

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




*/
/*
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

    Triangle triangle(points[0], points[1], points[2]);
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
    Triangle triangle(p1, p2, p3);

    Point pointInside(1.0, 1.0, 4);

    bool isInside = triangle.isInsideCircumcircle(pointInside);

    EXPECT_TRUE(isInside);
}

// Test case per verificare se un punto è all'esterno della circonferenza circoscritta al triangolo
TEST(TriangleTest, PointOutsideCircumcircle) {
    Point p1(0.0, 0.0, 1);
    Point p2(2.0, 0.0, 2);
    Point p3(1.0, 2.0, 3);
    Triangle triangle(p1, p2, p3);

    Point pointOutside(0.0, 2.0, 4);

    bool isInside = triangle.isInsideCircumcircle(pointOutside);

    EXPECT_FALSE(isInside);
}

TEST(CalculateAngleTest, CorrectAngle) {
    Point p1(0.0, 0.0, 1);
    Point p2(0.0, 1.0, 2);
    Point p3(1.0, 0.0, 3);

    Triangle triangle(p1, p2, p3);

    double expectedAngleRadians = M_PI * 0.25;

    double angleRadians = abs(p1.calculateAngle( p2, p3));

    EXPECT_DOUBLE_EQ(angleRadians, expectedAngleRadians);
}
// Test per la funzione isPointInsideTriangle
TEST(IsPointInsideTriangleTest, PointInside) {
    Point p1(0.0, 0.0, 1);
    Point p2(0.0, 1.0, 2);
    Point p3(1.0, 0.0, 3);

    Triangle triangle(p1, p2, p3);

    Point Q(0.5, 0.5, 4);

    bool result = triangle.isPointInsideTriangle(Q);

    EXPECT_TRUE(result);
}

TEST(IsPointInsideTriangleTest, PointOutside) {
    Point p1(0.0, 0.0, 1);
    Point p2(0.0, 1.0, 2);
    Point p3(1.0, 0.0, 3);

    Triangle triangle(p1, p2, p3);

    Point Q(0.0, 2.0, 4);

    bool result = triangle.isPointInsideTriangle(Q);

    EXPECT_FALSE(result);
}

#endif // __TEST_H
