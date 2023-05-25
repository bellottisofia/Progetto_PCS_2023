#include "delaunay.hpp"

namespace ProjectLibrary
{

Point::Point(const double& x, const double& y, const unsigned int& id)
    : x(x), y(y), id(id)
{
}
Point::Point(const Point& p)
    : x(p.x), y(p.y), id(p.id)
{
}
Point& Point::operator=(const Point& other) {
  if (this != &other) {
    // Effettua una copia profonda dei membri
    x = other.x;
    y = other.y;
    id = other.id;
  }
  return *this;
}
double Point::calculateAngle( const Point& B, const Point& C) {
        double ABx = B.x - x;
        double ABy = B.y - y;
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

/*
double cross(const Point& a, const Point& b) {
    return a.x * b.y - a.y * b.x;
}

double area(const Point& a, const Point& b, const Point& c) {
    return abs(cross(b - a, c - a)) / 2;
}
*/
Triangle::Triangle(const Point& point1, const Point& point2, const Point& point3)
            : p1(point1), p2(point2), p3(point3) {}
/*
double calculateArea(){
    return 0.5 * ((p1.x - p3.x) * (p2.y - p3.y) - (p1.y - p3.y) * (p2.x - p3.x));
}
*/
double Triangle::calculateArea(){
    return abs(0.5 * ((p1.x - p3.x) * (p2.y - p3.y) - (p1.y - p3.y) * (p2.x - p3.x)));
}
bool Triangle::IsVerticesSort()
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
void Triangle::SortVertices(){
    if (IsVerticesSort()==true)
    return;
    else
    {
        Point temp = p2;
            p2 = p3;
            p3 = temp;
    }
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


bool Triangle::isPointInsideTriangle(const Point& Q) {
    // Calcola le aree dei triangoli formati dal punto Q e ciascun lato del triangolo ABC
    Triangle t1(p1, p2, p3);
    Triangle t2(Q, p2, p3);
    Triangle t3(p1, Q, p3);
    Triangle t4(p1, p2, Q);

    double areaABC = t1.calculateArea();
    double areaPBC = t2.calculateArea();
    double areaPCA = t3.calculateArea();
    double areaPAB = t4.calculateArea();

    // Verifica se la somma delle aree dei triangoli interni è uguale all'area totale del triangolo ABC
    return abs(areaABC - (areaPBC + areaPCA + areaPAB)) < 1e-12;
}


bool Triangle::isVertexShared(const Point& vertex)  {
    return (vertex == p1) || (vertex == p2) || (vertex == p3);
}

bool Triangle::isAdjacent(const Triangle& other) {
    int sharedVertices = 0;

    if (isVertexShared(other.p1))
        sharedVertices++;
    if (isVertexShared(other.p2))
        sharedVertices++;
    if (isVertexShared(other.p3))
        sharedVertices++;

    if( sharedVertices >= 2){
        return true;
    }
    else
    {return false;
    }

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
                double area = crossMaxTriangle.calculateArea();
                if (area > crossMaxArea) {
                    crossMaxArea = area;
                    crossMaxTriangle = Triangle(points[i], points[mid], points[j]);
                }
            }
        }

        // Restituisco il triangolo con l'area massima tra le tre opzioni
        if (leftMaxTriangle.calculateArea() >= rightMaxTriangle.calculateArea() && leftMaxTriangle.calculateArea() >= crossMaxArea) {
            return leftMaxTriangle;
        } else if (rightMaxTriangle.calculateArea() >= leftMaxTriangle.calculateArea() && rightMaxTriangle.calculateArea() >= crossMaxArea) {
            return rightMaxTriangle;
        } else {
            return crossMaxTriangle;
        }
    }
    Triangle Triangle::findMaximumTriangle(const std::vector<Point>& points) {
        return findMaximumTriangleArea(points, 0, points.size() - 1);
    }





    Triangulation::Triangulation(const std::vector<Triangle>& DelunayTriangles);
    DelunayTriangles(DelunayTriangles){};

    void Triangulation::addPointToTriangulation(const Point& Q);

    const Triangulation::std::vector<Triangle>& getTriangles() const {
        return DelunayTriangles;
    }

    const Triangulation::std::vector<const Triangle*>& getAdjacentTriangles(const Triangle& DelunayTriangles) const {
        return adjacencyList.at(&DelunayTriangles);
    }

    void Triangulation::updateAdjacency(const Triangle& DelunayTriangles) {
        std::vector<const Triangle*> adjacentTriangles;

        for (const Triangle& other : DelunayTriangles) {
            if (DelunayTriangles.isAdjacent(other)) {
                adjacentTriangles.push_back(&other);
                adjacencyList[&other].push_back(&DelunayTriangles);
            }
        }

        adjacencyList[&DelunayTriangles] = adjacentTriangles;
    }






    int Point::calculateOrientation(const Point& p1, const Point& p2) {
        double value = (p2.y - p1.y) * (x - p2.x) - (p2.x - p1.x) * (y - p2.y);

        if (value == 0) {
            return 0;  // I punti sono allineati
        } else if (value > 0) {
            return 1;  // Orientamento orario
        } else {
            return -1; // Orientamento antiorario
        }
    }


    // Verifica se il punto q si trova sul segmento p-r
    bool Point::isPointOnSegment( const Point& q, const Point& r) {
        return x >= min(p.x, r.x) && q.x <= max(x, r.x) &&
               y >= min(p.y, r.y) && q.y <= max(y, r.y);
    }

    bool Triangle::doSegmentsIntersect( const Point& Q) {
        // Calcola l'orientamento dei punti rispetto ai segmenti
        int orientation1 = -1;
        int orientation2 = calculateOrientation(p1, p2, Q);
        int orientation3 = calculateOrientation(p3, Q, p1);
        int orientation4 = calculateOrientation(p3, Q, p2);

        // Controlla se i segmenti si intersecano
        if (orientation1 != orientation2 && orientation3 != orientation4) {
            return true;
        }

        // Controlla i casi speciali in cui i segmenti sono sovrapposti
        if (orientation1 == 0 && isPointOnSegment(p1, p2, p3)) {
            return true;
        }
        if (orientation2 == 0 && isPointOnSegment(p1, p2, p4)) {
            return true;
        }
        if (orientation3 == 0 && isPointOnSegment(p3, p4, p1)) {
            return true;
        }
        if (orientation4 == 0 && isPointOnSegment(p3, p4, p2)) {
            return true;
        }

        return false;
    }



    void Triangulation::addPointToTriangulation(const Point& Q) {
        Triangle containingTriangle;

        // Verifica se il punto Q è interno alla triangolazione
        bool isInside = false;
        for (const Triangle& triangle : DelunayTriangles) {
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
            DelunayTriangles.push_back(newTriangle1);
            DelunayTriangles.push_back(newTriangle2);
            DelunayTriangles.push_back(newTriangle3);

            // Rimuovi il triangolo T dalla triangolazione
            DelunayTriangles.erase(std::remove(DelunayTriangles.begin(), DelunayTriangles.end(), containingTriangle), DelunayTriangles.end());
        } else {
            // Il punto Q è esterno alla triangolazione, unisce Q con tutti i vertici escludendo le intersezioni
            for (size_t i = 0; i < DelunayTriangles.size(); ++i) {
                Triangle& triangle = DelunayTriangles[i];
                Triangle QTriangle(Q, triangle.A, triangle.B);
                if (!doSegmentsIntersect(triangle.A, triangle.B, triangle.C, Q)) {
                    // Aggiungi il nuovo triangolo QTriangle solo se non ha intersezioni con il segmento AB
                    DelunayTriangles.push_back(QTriangle);
                }

                QTriangle = Triangle(Q, triangle.B, triangle.C);
                if (!doSegmentsIntersect(triangle.A, triangle.B, triangle.C, Q)) {
                    // Aggiungi il nuovo triangolo QTriangle solo se non ha intersezioni con il segmento BC
                    DelunayTriangles.push_back(QTriangle);
                }

                QTriangle = Triangle(Q, triangle.C, triangle.A);
                if (!doSegmentsIntersect(triangle.A, triangle.B, triangle.C, Q)) {
                    // Aggiungi il nuovo triangolo QTriangle solo se non ha intersezioni con il segmento CA
                    DelunayTriangles.push_back(QTriangle);
                }

                // Rimuovi il triangolo corrente dalla triangolazione
                DelunayTriangles.erase(DelunayTriangles.begin() + i);
                --i;
            }
        }
    }



