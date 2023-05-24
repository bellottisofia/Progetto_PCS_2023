#ifndef __TEST_H
#define __TEST_H

#include <gtest/gtest.h>
#include "delaunay.hpp"



using namespace testing;
using namespace std;
using namespace ProjectLibrary;

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

    double angleRadians = abs(triangle.calculateAngle(p1, p2, p3));

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
