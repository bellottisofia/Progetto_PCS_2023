#include "delaunay.hpp"
#include <math.h>
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;

namespace ProjectLibrary
{


    double Triangle::calculateArea(const Point& p1, const Point& p2, const Point& p3) {
        return std::abs((p1.x * (p2.y - p3.y) + p2.x * (p3.y - p1.y) + p3.x * (p1.y - p2.y)) / 2);
    }

    Triangle Triangle::findMaximumTriangleArea(const std::vector<Point>& points, int start, int end) {
        if (start == end) {
            // Non ci sono abbastanza punti per formare un triangolo
            return Triangle(Point(0, 0, 0), Point(0, 0, 0), Point(0, 0, 0));
        }
        if (end - start == 1) {
            // Nessun triangolo in due punti
            return Triangle(Point(0, 0, 0), Point(0, 0, 0), Point(0, 0, 0));
        }
        if (end - start == 2) {
            // Calcolo l'area del triangolo con i tre punti
            return Triangle(points[start], points[start + 1], points[start + 2]);
        }

        // Calcolo il punto mediano
        int mid = (start + end) / 2;

        // Calcolo l'area massima del triangolo a sinistra
        Triangle leftMaxTriangle = findMaximumTriangleArea(points, start, mid);

        // Calcolo l'area massima del triangolo a destra
        Triangle rightMaxTriangle = findMaximumTriangleArea(points, mid + 1, end);

        // Calcolo l'area massima del triangolo che attraversa la divisione
        double crossMaxArea = 0.0;
        Triangle crossMaxTriangle(Point(0, 0, 0), Point(0, 0, 0), Point(0, 0, 0));
        for (int i = mid; i >= start; i--) {
            for (int j = mid + 1; j <= end; j++) {
                double area = calculateArea(points[i], points[mid], points[j]);
                if (area > crossMaxArea) {
                    crossMaxArea = area;
                    crossMaxTriangle = Triangle(points[i], points[mid], points[j]);
                }
            }
        }

        // Restituisco il triangolo con l'area massima tra le tre opzioni
        if (calculateArea(leftMaxTriangle.p1,leftMaxTriangle.p2,leftMaxTriangle.p3) >= calculateArea(rightMaxTriangle.p1,rightMaxTriangle.p2,rightMaxTriangle.p3) && calculateArea(leftMaxTriangle.p1,leftMaxTriangle.p2,leftMaxTriangle.p3) >= crossMaxArea) {
            return leftMaxTriangle;
        } else if (calculateArea(rightMaxTriangle.p1,rightMaxTriangle.p2,rightMaxTriangle.p3) >= calculateArea(leftMaxTriangle.p1,leftMaxTriangle.p2,leftMaxTriangle.p3) &&calculateArea(rightMaxTriangle.p1,rightMaxTriangle.p2,rightMaxTriangle.p3) >= crossMaxArea) {
            return rightMaxTriangle;
        } else {
            return crossMaxTriangle;
        }
    }


    Triangle Triangle::findMaximumTriangle(const std::vector<Point>& points) {
        return findMaximumTriangleArea(points, 0, points.size() - 1);
    }

    bool Triangle::isInsideCircumcircle(const Point& point)  {
        // Calcola il determinante della matrice per verificare se il punto è dentro la circonferenza circoscritta
            double a = p1.x - point.x;
            double b = p1.y - point.y;
            double c = a*a + b*b;
            double d = p2.x - point.x;
            double e = p2.y - point.y;
            double f = d*d + e*e;
            double g = p3.x - point.x;
            double h = p3.y - point.y;
            double i = g*g + h*h;

            double det = a * (e * i - h * f) - b * (d * i - f * g) + c * (d * h - g * e);

            return det > 0.0;
    }


    double Triangle::calculateAngle(const Point& A, const Point& B, const Point& C) {
        double ABx = B.x - A.x;
        double ABy = B.y - A.y;
        double BCx = C.x - B.x;
        double BCy = C.y - B.y;

        // Calculate the dot product of vectors AB and BC
        double dotProduct = ABx * BCx + ABy * BCy;
        // Calculate the magnitudes of vectors AB and BC
        double magnitudeAB = std::sqrt(ABx * ABx + ABy * ABy);
        double magnitudeBC = std::sqrt(BCx * BCx + BCy * BCy);

        // Calculate the angle in radians using the dot product and magnitudes
        double cosAngle = dotProduct / (magnitudeAB * magnitudeBC);
        double angleRadians = std::acos(cosAngle);
        return M_PI - angleRadians;
    }

    //usa il determinante per controllare se Q è dentro il cerchio circoscritto
    bool Triangle::isPointInsideTriangle(const Point& Q) {
        double angleSum = 0.0;

        // Calcola gli angoli tra il punto Q e ciascun vertice del triangolo ABC
        double angleA = calculateAngle(p3, Q, p1);
        double angleB = calculateAngle(p1, Q, p2);
        double angleC = calculateAngle(p2, Q, p3);

        // Calcola la somma degli angoli
        angleSum = angleA + angleB + angleC;

        // Verifica se la somma degli angoli è uguale a 360 gradi
        return (abs(angleSum - 2*M_PI) < 1e-12);
    }
/*
//opzione se ho punti ordinati
double findMaximumTriangleArea(const std::vector<Point>& points) {
    int n = points.size();
    double maxArea = 0.0;

    for (int i = 0; i < n - 2; i++) {
        int j = i + 1;
        int k = j + 1;

        while (k < n) {
            double area = calculateTriangleArea(points[i], points[j], points[k]);
            maxArea = std::max(maxArea, area);
            k++;
        }
    }

    return maxArea;
}

// Find adjacent triangles
    std::vector<Triangle> findAdjacentTriangles(const std::vector<Triangle>& triangles) const {
        std::vector<Triangle> adjacentTriangles;

        for (const auto& triangle : triangles) {
            if (triangle != *this && areAdjacent(triangle)) {
                adjacentTriangles.push_back(triangle);
            }
        }

        return adjacentTriangles;
    }

/*
 * Per verificare se un punto Q è esterno al quadrato circoscritto o al cerchio circoscritto di un triangolo T,
 * possiamo utilizzare le seguenti considerazioni:

Quadrato circoscritto:

Troviamo il punto più a sinistra (minX), il punto più a destra (maxX), il punto più in alto (maxY) e
il punto più in basso (minY) tra i vertici del triangolo T.
Se il punto Q ha una coordinata x minore di minX o maggiore di maxX,
 o una coordinata y minore di minY o maggiore di maxY, allora il punto Q è esterno al quadrato circoscritto.
Cerchio circoscritto:

Calcoliamo il centro del cerchio circoscritto (xc, yc) utilizzando le coordinate dei vertici del triangolo T.
Calcoliamo il raggio del cerchio circoscritto (r) utilizzando le coordinate dei vertici del triangolo T.
Calcoliamo la distanza tra il centro del cerchio e il punto Q utilizzando la formula della distanza euclidea.
Se la distanza tra il centro del cerchio e il punto Q è maggiore del raggio r, allora il punto Q è esterno al cerchio circoscritto.


bool isPointInsideSquare(const Point& Q, const Triangle& T) {
    double minX = std::min({ T.a.x, T.b.x, T.c.x });
    double maxX = std::max({ T.a.x, T.b.x, T.c.x });
    double minY = std::min({ T.a.y, T.b.y, T.c.y });
    double maxY = std::max({ T.a.y, T.b.y, T.c.y });

    return (Q.x >= minX && Q.x <= maxX && Q.y >= minY && Q.y <= maxY);
}





bool SortVertices(const Point& p1,
               const Point& p2,
               const Point& p3)
{
  double dx1 = p2.x - p1.x;
  double dy1 = p2.y - p1.y;
  double dx2 = p3.x - p1.x;
  double dy2 = p3.y - p1.y;

  if (dx1 * dy2 > dy1 * dx2)
    return true;
  else if (dx1 * dy2 < dy1 * dx2)
    return false;
  else
    return (dx1 * dx2 + dy1 * dy2 >= 0);
}


bool isPointInPolygon(const Point& Q, const vector<Triangle>& triangulation) {
    int count = 0;

    for (const Triangle& T : triangulation) {
        if (Q == T.p1 || Q == T.p2 || Q == T.p3) {
            return true;  // Il punto Q è un vertice della triangolazione (di bordo)
        }

        if ((Q.y > T.p1.y && Q.y <= T.p2.y) || (Q.y > T.p2.y && Q.y <= T.p1.y)) {
            double xIntersection = (Q.y - T.p1.y) * (T.p2.x - T.p1.x) / (T.p2.y - T.p1.y) + T.p1.x;
            if (Q.x < xIntersection) {
                count++;
            }
        }
    }

    return (count % 2 == 1);  // true se il punto Q è interno, false se è di bordo
}






    bool verifyDelaunayHypothesis(const Point& pointA, const Point& pointB, const Point& pointC, const Point& pointD) {
        // Calcola gli angoli opposti al lato di adiacenza BC
        double angleABC = calculateAngle(pointA, pointB, pointC);
        double angleBDC = calculateAngle(pointB, pointD, pointC);

        // Verifica se la somma degli angoli è minore di °
        double sumAngles = angleABC + angleBDC;
        const double maxSumAngles = 180.0 - 0.0001;  // ° - tolleranza
        return sumAngles < maxSumAngles;
    }

    bool doSegmentsIntersect(const Point& p1, const Point& q1, const Point& p2, const Point& q2) {
        // Calcola le direzioni dei segmenti
        double d1 = direction(p2, q2, p1);
        double d2 = direction(p2, q2, q1);
        double d3 = direction(p1, q1, p2);
        double d4 = direction(p1, q1, q2);

        // Controlla se i segmenti si intersecano
        if ((d1 > 0.0 && d2 < 0.0 || d1 < 0.0 && d2 > 0.0) && (d3 > 0.0 && d4 < 0.0 || d3 < 0.0 && d4 > 0.0))

         {
            return true;
        } else if (d1 == 0 && isPointOnSegment(p2, q2, p1)) {
            return true;
        } else if (d2 == 0 && isPointOnSegment(p2, q2, q1)) {
            return true;
        } else if (d3 == 0 && isPointOnSegment(p1, q1, p2)) {
            return true;
        } else if (d4 == 0 && isPointOnSegment(p1, q1, q2)) {
            return true;
        }

        return false;
    }


    // Calcola la direzione del punto r rispetto alla linea formata da p e q
    double direction(const Point& p, const Point& q, const Point& r) {
        return (r.x - p.x) * (q.y - p.y) - (r.y - p.y) * (q.x - p.x);
    }

    // Verifica se il punto q si trova sul segmento p-r
    bool isPointOnSegment(const Point& p, const Point& q, const Point& r) {
        return q.x >= min(p.x, r.x) && q.x <= max(p.x, r.x) &&
               q.y >= min(p.y, r.y) && q.y <= max(p.y, r.y);
    }

    void addPointToTriangulation(const Point& Q, std::vector<Triangle>& triangulation) {
        Triangle containingTriangle;

        // Verifica se il punto Q è interno alla triangolazione
        bool isInside = false;
        for (const Triangle& triangle : triangulation) {
            if (isPointInsideTriangle(Q, triangle)) {
                isInside = true;
                containingTriangle = triangle;
                break;
            }
        }

        if (isInside) {
            // Se il punto Q è interno alla triangolazione, crea una sottotriangolazione unendo Q con i vertici del triangolo T
            Triangle newTriangle1(Q, containingTriangle.A, containingTriangle.B);
            Triangle newTriangle2(Q, containingTriangle.B, containingTriangle.C);
            Triangle newTriangle3(Q, containingTriangle.C, containingTriangle.A);

            // Aggiungi i nuovi triangoli alla triangolazione
            triangulation.push_back(newTriangle1);
            triangulation.push_back(newTriangle2);
            triangulation.push_back(newTriangle3);

            // Rimuovi il triangolo T dalla triangolazione
            triangulation.erase(std::remove(triangulation.begin(), triangulation.end(), containingTriangle), triangulation.end());
        } else {
            // Il punto Q è esterno alla triangolazione, unisce Q con tutti i vertici escludendo le intersezioni
            for (size_t i = 0; i < triangulation.size(); ++i) {
                Triangle& triangle = triangulation[i];
                Triangle QTriangle(Q, triangle.A, triangle.B);
                if (!doSegmentsIntersect(triangle.A, triangle.B, triangle.C, Q)) {
                    // Aggiungi il nuovo triangolo QTriangle solo se non ha intersezioni con il segmento AB
                    triangulation.push_back(QTriangle);
                }

                QTriangle = Triangle(Q, triangle.B, triangle.C);
                if (!doSegmentsIntersect(triangle.A, triangle.B, triangle.C, Q)) {
                    // Aggiungi il nuovo triangolo QTriangle solo se non ha intersezioni con il segmento BC
                    triangulation.push_back(QTriangle);
                }

                QTriangle = Triangle(Q, triangle.C, triangle.A);
                if (!doSegmentsIntersect(triangle.A, triangle.B, triangle.C, Q)) {
                    // Aggiungi il nuovo triangolo QTriangle solo se non ha intersezioni con il segmento CA
                    triangulation.push_back(QTriangle);
                }

                // Rimuovi il triangolo corrente dalla triangolazione
                triangulation.erase(triangulation.begin() + i);
                --i;
            }
        }
    }

    Triangulation::Triangulation(const std::vector<Triangle>& initialTriangles)
            : triangles(initialTriangles) {}
/*
 * // Tipo per l'identificatore dei triangoli
using TriangleId = size_t;

// Struttura per l'adiacenza dei triangoli
struct TriangleAdjacency {
    std::vector<TriangleId> adjacentTriangles;
};

// Funzione per creare una nuova entrata di adiacenza per un triangolo
TriangleAdjacency createTriangleAdjacency() {
    TriangleAdjacency adjacency;
    adjacency.adjacentTriangles.clear();
    return adjacency;
}

// Funzione per aggiungere l'adiacenza tra due triangoli
void addTriangleAdjacency(std::vector<TriangleAdjacency>& adjacencyList, TriangleId triangleId1, TriangleId triangleId2) {
    // Aggiungi l'adiacenza tra triangleId1 e triangleId2
    adjacencyList[triangleId1].adjacentTriangles.push_back(triangleId2);
    adjacencyList[triangleId2].adjacentTriangles.push_back(triangleId1);
}

// Funzione per calcolare il prodotto vettoriale tra due vettori
double crossProduct(const Point& A, const Point& B, const Point& C) {
    return (B.x - A.x) * (C.y - A.y) - (C.x - A.x) * (B.y - A.y);
}

// Funzione per ordinare i vertici del triangolo in senso antiorario
vector<Point> orderVerticesCounterclockwise(const Point& A, const Point& B, const Point& C) {
    // Calcola il baricentro del triangolo
    Point centroid;
    centroid.x = (A.x + B.x + C.x) / 3;
    centroid.y = (A.y + B.y + C.y) / 3;

    // Calcola il prodotto vettoriale tra il vettore del vertice e il vettore del baricentro per ogni vertice
    double crossAB = crossProduct(A, B, centroid);
    double crossBC = crossProduct(B, C, centroid);
    double crossCA = crossProduct(C, A, centroid);

    // Crea un vettore per i vertici ordinati
    vector<Point> orderedVertices;

    // Ordina i vertici in base al segno del prodotto vettoriale
    if ((crossAB < 0 && crossBC < 0 && crossCA < 0) || (crossAB > 0 && crossBC > 0 && crossCA > 0)) {
        // Caso in cui i vertici sono già in senso antiorario
        orderedVertices.push_back(A);
        orderedVertices.push_back(B);
        orderedVertices.push_back(C);
    } else {
        // Caso in cui i vertici non sono in senso antiorario, quindi si fa lo scambio per ottenere l'ordinamento corretto
        orderedVertices.push_back(A);
        orderedVertices.push_back(C);
        orderedVertices.push_back(B);
    }

    return orderedVertices;
}

 */

}
