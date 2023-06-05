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
/*
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


int calculateOrientation(const Point& p1, const Point& p2, const Point& p3) {
    double value = (p2.y - p1.y) * (p3.x - p2.x) - (p2.x - p1.x) * (p3.y - p2.y);

    if (value == 0) {
        return 0;  // I punti sono allineati
    } else if (value > 0) {
        return 1;  // Orientamento orario
    } else {
        return -1; // Orientamento antiorario
    }
}
*/
bool Point::isPointInVector(const std::vector<Point>& points) const {
    return std::find(points.begin(), points.end(), *this) != points.end();
}

/*
bool Point::arePointsCollinear(const Point& p2,const Point& p3){


    if ((p2.y - y) * (p3.x - p2.x) - (p2.x - x) * (p3.y - p2.y)){
        return true;}
    else
        {
            return false;
        }

}
*/


Triangle::Triangle(const Point& point1, const Point& point2, const Point& point3, const unsigned int& id):
    p1(point1),
    p2(point2),
    p3(point3),
    id(id)
{}

double Triangle::calculateArea()const{
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
    if (IsVerticesSort())
    return;
    else
    {
        Point temp = p2;
            p2 = p3;
            p3 = temp;
    }
}

bool Triangle::isInsideCircumcircle(const Point& point) const {
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

           double det = a * (e * i - g * f) - b * (d * i - f * g) + c * (d * h - g * e);

           return det > 0.0;
   }


bool Triangle::isPointInsideTriangle(const Point& Q)const {
    if (!isInsideCircumcircle(Q))
    {return false;}
    else
    {
    // Calcola le aree dei triangoli formati dal punto Q e ciascun lato del triangolo ABC
    //Triangle t1(p1, p2, p3);
    Triangle t2(Q, p2, p3,id+1);
    Triangle t3(p1, Q, p3,id+2);
    Triangle t4(p1, p2, Q,id+3);

    //double areaABC = t1.calculateArea();
    double areaABC = calculateArea();
    double areaPBC = t2.calculateArea();
    double areaPCA = t3.calculateArea();
    double areaPAB = t4.calculateArea();

    // Verifica se la somma delle aree dei triangoli interni è uguale all'area totale del triangolo ABC
    return abs(areaABC - (areaPBC + areaPCA + areaPAB)) < 1e-12;
    }
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
/*
Triangle Triangle::findMaximumTriangleArea(const std::vector<Point>& points, int start, int end) {
        if (start == end) {
            // Non ci sono abbastanza punti per formare un triangolo
            return Triangle(Point(0, 0, 0), Point(0, 0, 0), Point(0, 0, 0),0);
        }
        if (end - start == 1) {
            // Nessun triangolo in due punti
            return Triangle(Point(0, 0, 0), Point(0, 0, 0), Point(0, 0, 0),0);
        }
        if (end - start == 2) {
            // Calcolo l'area del triangolo con i tre punti
            return Triangle(points[start], points[start + 1], points[start + 2],0);
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
                    crossMaxTriangle = Triangle(points[i], points[mid], points[j],0);
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
*/
Triangle Triangle::findMaximumTriangleArea(const std::vector<Point>& points) {
    if (points.size() < 3) {
        // Non ci sono abbastanza punti per formare un triangolo
        return Triangle(Point(0, 0,0), Point(0, 0,0), Point(0,0,0),0);
    }

    double maxArea = 0.0;
    Triangle maxTriangle(Point(0, 0,0), Point(0, 0,0), Point(0, 0,0),0);

    for (size_t i = 0; i < points.size() - 2; ++i) {
        const Point& p1 = points[i];
        const Point& p2 = points[i + 1];
        const Point& p3 = points[i + 2];

        double area = Triangle(p1, p2, p3,1).calculateArea();

        if (area > maxArea) {
            maxArea = area;
            maxTriangle = Triangle(p1, p2, p3,1);
        }
    }

    return maxTriangle;
}

    std::vector<Point> Triangle::getVertices() const {
        return {p1, p2, p3};
    }

    const Point& Triangle::getVertexA() const {
        return p1;
    }

    const Point& Triangle::getVertexB() const {
        return p2;
    }

    const Point& Triangle::getVertexC() const {
        return p3;
    }


    Triangulation::Triangulation(const std::vector<Triangle>& DelunayTriangles) : DelunayTriangles(DelunayTriangles) {
        buildAdjacencyList();
    }

    void Triangulation::addTriangle(const Triangle& triangle) {
        DelunayTriangles.push_back(triangle);
        adjacencyList[triangle.id] = {}; // Inizializza la lista di adiacenza per il triangolo
    }
    void Triangulation::addAdjacentTriangle(int triangleId, int adjacentTriangleId) {
        adjacencyList[triangleId].push_back(adjacentTriangleId);
        adjacencyList[adjacentTriangleId].push_back(triangleId);
    }
    const std::vector<unsigned int>& Triangulation::getAdjacentTriangles(int triangleId) {
        return adjacencyList[triangleId];
    }


/*
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





*/

    // Verifica se il punto q si trova sul segmento p-r



    bool Point::isPointOnSegment(const Point& p,  const Point& r) {
        // Verifica se il punto q si trova sul segmento p-r
        return x >= std::min(p.x, r.x) && x <= std::max(p.x, r.x) &&
               y >= std::min(p.y, r.y) && y <= std::max(p.y, r.y);
    }




    bool Point::doSegmentsIntersect( const Point& p2, const Point& p3, const Point& p4) {
        // Calcola i coefficienti delle rette
        double m1 = (p2.y -y) / (p2.x - x);
        double m2 = (p4.y - p3.y) / (p4.x - p3.x);

        // Calcola l'intercetta delle rette
        double q1 = y - m1 * x;
        double q2 = p3.y - m2 * p3.x;

        // Calcola le coordinate dell'intersezione
        double x_coord = (q2 - q1) / (m1 - m2);
        double y_coord = m1 * x_coord + q1;

        // Restituisci il punto di intersezione
        Point intersection;
        intersection.x = x_coord;
        intersection.y = y_coord;
        intersection.id = -1;
        if(intersection.isPointOnSegment(p3,p4))
        return true;
        else
            return false;
    }

    std::vector<Point> Triangle::findIntersection(const Point& q) {
        std::vector<Point> intersections;

        // Controlla l'intersezione con i lati del triangolo
        if (p3.doSegmentsIntersect(p1, p2, q)) {
            intersections.push_back(p2);
            intersections.push_back(q);
        }
        else if (p1.doSegmentsIntersect(p2, p3, q)) {
            intersections.push_back(p3);
            intersections.push_back(q);
        }
        else if (p2.doSegmentsIntersect(p3, p1, q)) {
            intersections.push_back(p1);
            intersections.push_back(q);
        }

        return intersections;
    }



bool Triangle::areTrianglesDelaunay(Triangle& triangle1) {
            // Verifica se i punti del primo triangolo sono all'interno della circonferenza circoscritta del secondo triangolo
            bool point1InsideCircle = isInsideCircumcircle(triangle1.p1);
            bool point2InsideCircle = isInsideCircumcircle(triangle1.p2);
            bool point3InsideCircle = isInsideCircumcircle(triangle1.p3);

            // Verifica se i punti del secondo triangolo sono all'interno della circonferenza circoscritta del primo triangolo
            bool point4InsideCircle = triangle1.isInsideCircumcircle(p1);
            bool point5InsideCircle = triangle1.isInsideCircumcircle(p2);
            bool point6InsideCircle = triangle1.isInsideCircumcircle(p3);

            // L'ipotesi di Delaunay è verificata solo se tutti i punti sono all'esterno delle rispettive circonferenze circoscritte
            return !point1InsideCircle && !point2InsideCircle && !point3InsideCircle &&
                   !point4InsideCircle && !point5InsideCircle && !point6InsideCircle;
        }
void Triangle::flip(Triangle& triangle1){
            // Verifica se i due triangoli sono adiacenti
            int commonVerticesCount = 0;
            Point commonVertices[2];

            // Trova i vertici comuni tra this e triangle1

            if (isVertexShared(triangle1.p1)) {
                commonVertices[commonVerticesCount++] = triangle1.p1;
            }
            if (isVertexShared(triangle1.p2)) {
                commonVertices[commonVerticesCount++] = triangle1.p2;
            }
            if (isVertexShared(triangle1.p3)) {
                commonVertices[commonVerticesCount++] = triangle1.p3;
            }

            //int count = 0;
            Point nonCommonVertex[2];
            if (commonVerticesCount == 2) {
                // Scambia i vertici dei due triangoli


                if (p1 != commonVertices[0] && p1 != commonVertices[1]) {
                    nonCommonVertex[0] = p1;
                }
                else if (p2 != commonVertices[0] && p2 != commonVertices[1]) {
                    nonCommonVertex[0] = p2;
                }
                else {
                    nonCommonVertex[0] = p3;
                }
                if (triangle1.p1 != commonVertices[0] && triangle1.p1 != commonVertices[1]) {
                    nonCommonVertex[1] = triangle1.p1;
                }
                else if (triangle1.p2 != commonVertices[0] && triangle1.p2 != commonVertices[1]) {
                    nonCommonVertex[1] = triangle1.p2;
                }
                else {
                    nonCommonVertex[1] = triangle1.p3;
                }

            }

            p1 = nonCommonVertex[0];
            p2 = commonVertices[0];
            p3 = nonCommonVertex[1];

            triangle1.p1= nonCommonVertex[0];
            triangle1.p2= nonCommonVertex[1];
            triangle1.p3= commonVertices[1];

    }

void Triangulation::buildAdjacencyList() {
    adjacencyList.clear();

    for (const Triangle& triangle : DelunayTriangles) {
        adjacencyList[triangle.id] = {}; // Inizializza la lista di adiacenza per il triangolo
    }

    // Aggiorna la lista di adiacenza basata sulle relazioni tra i triangoli
    // Implementazione specifica in base alle regole della triangolazione
}
/*
void Triangulation::addAdjacentTriangle(int triangleId, int adjacentTriangleId) {
    triangles[triangleId].adjacentTriangles.push_back(adjacentTriangleId);
    triangles[adjacentTriangleId].adjacentTriangles.push_back(triangleId);
}
*/
int Triangulation::PointInsideTriangulation(const Point& Q) {
    // Itera attraverso tutti i triangoli della triangolazione
    for (const Triangle& triangle : DelunayTriangles) {
        if (!triangle.isInsideCircumcircle(Q)) {
            return -1;
        }
        // Verifica se il punto Q è interno al triangolo corrente
        else if (triangle.isPointInsideTriangle(Q)) {
            return triangle.id;  //  restituisco l'id del triangolo a cui è interno
        }
    }

    return -1;  // Il punto non è interno alla triangolazione
}



void Triangulation::createSubtriangulation(const Point& Q, int triangleId) {
    // Ottieni il triangolo T dalla lista dei triangoli
    Triangle& T = DelunayTriangles[triangleId];

    // Crea il nuovo triangolo collegando il punto Q ai vertici di T
    Triangle newTriangle(Q, T.p1, T.p2,triangleId );

    // Aggiungi il nuovo triangolo alla lista dei triangoli
    DelunayTriangles[triangleId] = newTriangle;

    // Aggiorna la lista di adiacenza per il nuovo triangolo
    addAdjacentTriangle(triangleId, triangleId);

    // Aggiungi gli altri due triangoli della sottotriangolazione alla lista dei triangoli
    Triangle t1(Q, T.p2, T.p3,triangleId +1);
    Triangle t2(Q, T.p3, T.p1,triangleId +2);
    DelunayTriangles.push_back(t1);
    DelunayTriangles.push_back(t2);

    // Aggiorna la lista di adiacenza per i nuovi triangoli
    addAdjacentTriangle(triangleId, t1.id);
    addAdjacentTriangle(triangleId, t2.id);
}

 void Triangulation::connectPointToVertices(const Point& Q){
// Il punto Q è esterno alla triangolazione, unisce Q con tutti i vertici escludendo le intersezioni
for (size_t i = 0; i < DelunayTriangles.size(); ++i) {
    Triangle& triangle = DelunayTriangles[i];
    std::vector<Point> intersection = triangle.findIntersection(Q);

    if (intersection.empty()) {
    int id=triangle.id;
    Triangle t1(Q,triangle.p1,triangle.p2,id);
    t1.SortVertices();
    Triangle t2(Q,triangle.p1,triangle.p3,id );
    t2.SortVertices();
    Triangle t3(Q,triangle.p2,triangle.p3,id);
    t3.SortVertices();
    if(!t1.isPointInsideTriangle(triangle.p3))
    {
    id+=1;
    DelunayTriangles.push_back(t1);}
    if(!t2.isPointInsideTriangle(triangle.p2))
    {id+=1;
    DelunayTriangles.push_back(t2);
    }
     if(!t3.isPointInsideTriangle(triangle.p1))
     {id+=1;
    DelunayTriangles.push_back(t3);
     }
    }
    else{
        if (triangle.p1.isPointInVector(intersection)){
           DelunayTriangles.push_back(Triangle(Q,triangle.p2,triangle.p3,triangle.id +1));
        }
        else if (triangle.p2.isPointInVector(intersection)){
            DelunayTriangles.push_back(Triangle(Q,triangle.p1,triangle.p3,triangle.id +1));
         }
        else {
              DelunayTriangles.push_back(Triangle(Q,triangle.p1,triangle.p2,triangle.id +1));
                }



}}}

void Triangulation::addPointToTriangulation(const Point& Q) {

    if (PointInsideTriangulation(Q)>0) {
        int triangleId = PointInsideTriangulation(Q);
        createSubtriangulation(Q, triangleId);
    } else {
        connectPointToVertices(Q);

    }
}

/*
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

/*
 *
    void Triangulation::addPoint(const Point& point) {
        // Check if the point is already in the triangulation
        if (std::find_if(points.begin(), points.end(), [&](const Point& p) { return p.id == point.id; }) != points.end())
            return;

        points.push_back(point);

        // Check if the point is inside the circumcircle of any existing triangle
        for (auto it = DelaunayTriangles.begin(); it != DelaunayTriangles.end(); ++it) {
            if (it->isInsideCircumcircle(point)) {
                // Create three new triangles by connecting the new point with each edge of the existing triangle
                Triangle t1(point, it->p1, it->p2);
                Triangle t2(point, it->p2, it->p3);
                Triangle t3(point, it->p3, it->p1);

                // Remove the existing triangle from the triangulation
                DelaunayTriangles.erase(it);

                // Add the three new triangles to the triangulation
                DelaunayTriangles.push_back(t1);
                DelaunayTriangles.push_back(t2);
                DelaunayTriangles.push_back(t3);

                // Update the adjacency information of the triangles
                updateTriangleAdjacency(t1);
                updateTriangleAdjacency(t2);
                updateTriangleAdjacency(t3);

                // Check for legalizing edges and flip adjacent triangles if necessary
                legalizeEdges(t1, point);
                legalizeEdges(t2, point);
                legalizeEdges(t3, point);

                break;
            }
        }
    }

    void Triangulation::updateTriangleAdjacency(Triangle& triangle) {
        // Find adjacent triangles by checking if they share any vertices
        for (auto& t : DelaunayTriangles) {
            if (triangle.isVertexShared(t)) {
                // Update the adjacency information of both triangles
                triangle.adjacentTriangles.push_back(&t);
                t.adjacentTriangles.push_back(&triangle);
            }
        }
    }

    void Triangulation::legalizeEdges(Triangle& triangle, const Point& point) {
        for (auto& adjTriangle : triangle.adjacentTriangles) {
            // Check if the point is inside the circumcircle of the adjacent triangle
            if (adjTriangle->isInsideCircumcircle(point)) {
                // Flip the common edge if necessary
                flipEdge(triangle, *adjTriangle, point);

                // Recursively legalize edges for the adjacent triangle
                legalizeEdges(triangle, point);
                legalizeEdges(*adjTriangle, point);
            }
        }
    }

    void Triangulation::flipEdge(Triangle& triangle1, Triangle& triangle2, const Point& point) {
        // Find the common edge between the two triangles
        std::pair<int, int> commonEdge = findCommonEdge(triangle1, triangle2);

        // Get the indices of the vertices of the common edge
        int v1 = commonEdge.first;
        int v2 = commonEdge.second;

        // Get the indices of the remaining vertices
        int v3 = getThirdVertexIndex(triangle1, v1, v2);
        int v4 = getThirdVertexIndex(triangle2, v1, v2);

        // Get the remaining vertices
        Point& p3 = points[v3];
        Point& p4 = points[v4];

        // Remove the two triangles from the triangulation
        removeTriangle(triangle1);
        removeTriangle(triangle2);

        // Create two new triangles by connecting the remaining vertices with the point
        Triangle newTriangle1(p3, p4, point);
        Triangle newTriangle2(p4, p3, point);

        // Add the two new triangles to the triangulation
        DelaunayTriangles.push_back(newTriangle1);
        DelaunayTriangles.push_back(newTriangle2);

        // Update the adjacency information of the new triangles
        updateTriangleAdjacency(newTriangle1);
        updateTriangleAdjacency(newTriangle2);

        // Check for legalizing edges and flip adjacent triangles if necessary
        legalizeEdges(newTriangle1, point);
        legalizeEdges(newTriangle2, point);
    }

    std::pair<int, int> Triangulation::findCommonEdge(const Triangle& triangle1, const Triangle& triangle2) {
        std::vector<int> vertices1 = { triangle1.p1.id, triangle1.p2.id, triangle1.p3.id };
        std::vector<int> vertices2 = { triangle2.p1.id, triangle2.p2.id, triangle2.p3.id };

        // Find the common vertices between the two triangles
        std::vector<int> commonVertices;
        for (int v1 : vertices1) {
            if (std::find(vertices2.begin(), vertices2.end(), v1) != vertices2.end())
                commonVertices.push_back(v1);
        }

        // The common vertices form the common edge
        return std::make_pair(commonVertices[0], commonVertices[1]);
    }

    int Triangulation::getThirdVertexIndex(const Triangle& triangle, int v1, int v2) {
        std::vector<int> vertices = { triangle.p1.id, triangle.p2.id, triangle.p3.id };

        // Find the index of the third vertex that is not v1 or v2
        for (int v : vertices) {
            if (v != v1 && v != v2)
                return v;
        }

        return -1; // Invalid index
    }

    void Triangulation::removeTriangle(Triangle& triangle) {
        // Remove the triangle from the triangulation
        auto it = std::find_if(DelaunayTriangles.begin(), DelaunayTriangles.end(),
            [&](const Triangle& t) { return t.p1.id == triangle.p1.id && t.p2.id == triangle.p2.id && t.p3.id == triangle.p3.id; });
        if (it != DelaunayTriangles.end())
            DelaunayTriangles.erase(it);

        // Remove the triangle from the adjacency information of other triangles
        for (auto& t : DelaunayTriangles) {
            auto adjIt = std::find_if(t.adjacentTriangles.begin(), t.adjacentTriangles.end(),
                [&](const Triangle* adjTriangle) { return adjTriangle->p1.id == triangle.p1.id && adjTriangle->p2.id == triangle.p2.id && adjTriangle->p3.id == triangle.p3.id; });
            if (adjIt != t.adjacentTriangles.end())
                t.adjacentTriangles.erase(adjIt);
        }
    }

    std::vector<Triangle> Triangulation::getDelaunayTriangles() const {
        return DelaunayTriangles;
    }

    std::vector<Point> Triangulation::getPoints() const {
        return points;
    }

*/

}
