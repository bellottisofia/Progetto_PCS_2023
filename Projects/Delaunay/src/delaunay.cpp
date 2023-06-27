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
    double ax = p1.x - point.x;
        double ay = p1.y - point.y;
        double bx = p2.x - point.x;
        double by = p2.y - point.y;
        double cx = p3.x - point.x;
        double cy = p3.y - point.y;

        double ab = ax * by - ay * bx;
        double bc = bx * cy - by * cx;
        double ca = cx * ay - cy * ax;

        double dot = ax * ax + ay * ay;
        double det = bx * bx + by * by;
        double det2 = cx * cx + cy * cy;

        double radiusSquared = dot * bc + det * ca + det2 * ab;

        return radiusSquared > 0.0;
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

bool Triangle::hasSharedPartialEdge(const Triangle& other) {
    // Controlla ogni possibile coppia di segmenti tra i lati dei due triangoli
    if (p1.doSegmentsIntersect(p2, other.p1, other.p2) ||
        p1.doSegmentsIntersect(p2, other.p1, other.p3) ||
        p1.doSegmentsIntersect(p2, other.p2, other.p3))
    {
        return true;
    }

    if (p2.doSegmentsIntersect(p3, other.p1, other.p2) ||
        p2.doSegmentsIntersect(p3, other.p1, other.p3) ||
        p2.doSegmentsIntersect(p3, other.p2, other.p3))
    {
        return true;
    }

    if (p3.doSegmentsIntersect(p1, other.p1, other.p2) ||
        p3.doSegmentsIntersect(p1, other.p1, other.p3) ||
        p3.doSegmentsIntersect(p1, other.p2, other.p3))
    {
        return true;
    }

    // Nessuna intersezione tra segmenti, quindi i due triangoli non condividono solo una parte del lato
    return false;
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
                Triangle crossMaxTriangle(Point(0, 0, 0), Point(0, 0, 0), Point(0, 0, 0),0);
                for (int i = mid; i >= start; i--) {
                    for (int j = mid + 1; j <= end; j++) {
                        Triangle triangle = Triangle(points[i], points[mid], points[j],0);
                        double area = triangle.calculateArea();
                        if (area > crossMaxArea) {
                            crossMaxArea = area;
                            crossMaxTriangle = triangle;
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

            std::vector<Point> sortedPoints = points;

            // Ordina il vettore di punti in base all'asse x utilizzando il MergeSort
            SortLibrary::MergeSortx(sortedPoints, 0, sortedPoints.size() -1 );

            return findMaximumTriangleArea(sortedPoints, 0, sortedPoints.size() -1 );

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


    Triangulation::Triangulation(const std::vector<Triangle>& triangles)
    : DelunayTriangles(triangles) {}



    void Triangulation::addTriangle(const Triangle& triangleAgg) {

        DelunayTriangles.push_back(triangleAgg);
        adjacencyList.emplace(triangleAgg.id, std::vector<unsigned int>());// Inizializza la lista di adiacenza per il triangolo
    }

    void Triangulation::addAdjacentTriangle(const unsigned int& triangleId,const unsigned int& adjacentTriangleId) {
        adjacencyList[triangleId].push_back(adjacentTriangleId);
        adjacencyList[adjacentTriangleId].push_back(triangleId);
    }
    const std::vector<unsigned int>& Triangulation::getAdjacentTriangles(const unsigned int& triangleId) {
        return adjacencyList[triangleId];
    }




    // Verifica se il punto q si trova sul segmento p-r



    bool Point::isPointOnSegment(const Point& p,  const Point& r) const{
        // Controlla se il punto Q ha la stessa pendenza rispetto ai punti p e r
        double crossProduct = (y - p.y) * (r.x - p.x) - (x - p.x) * (r.y - p.y);
        if (std::abs(crossProduct) > std::numeric_limits<double>::epsilon()) {
            return false;  // I punti non hanno la stessa pendenza, quindi Q non giace sul segmento
        }

        // Controlla se il punto Q rientra nell'intervallo tra i punti p e r
        if (x < std::min(p.x, r.x) || x > std::max(p.x, r.x) || y < std::min(p.y, r.y) || y > std::max(p.y, r.y)) {
            return false;  // Il punto Q non rientra nell'intervallo tra p e r
        }

        return true;  // Il punto Q giace sul segmento p-r
    }




    bool Point::doSegmentsIntersect( const Point& p2, const Point& p3, const Point& p4) {
            // Verifica se uno dei segmenti è verticale
                if (x == p2.x) {
                    // Calcola l'intersezione quando il primo segmento è verticale
                    double x_coord = x;
                    double m2 = (p4.y - p3.y) / (p4.x - p3.x);
                    double q2 = p3.y - m2 * p3.x;
                    double y_coord = m2 * x_coord + q2;

                    // Crea il punto di intersezione
                    Point intersection;
                    intersection.x = x_coord;
                    intersection.y = y_coord;
                    intersection.id = -1;

                    // Verifica se il punto di intersezione rientra all'interno dei segmenti
                    if (intersection.isPointOnSegment(p3, p4)) {
                        return true;  // C'è intersezione
                    }
                } else if (p3.x == p4.x) {
                    // Calcola l'intersezione quando il secondo segmento è verticale
                    double x_coord = p3.x;
                    double m1 = (p2.y - y) / (p2.x - x);
                    double q1 = y - m1 * x;
                    double y_coord = m1 * x_coord + q1;

                    // Crea il punto di intersezione
                    Point intersection;
                    intersection.x = x_coord;
                    intersection.y = y_coord;
                    intersection.id = -1;

                    // Verifica se il punto di intersezione rientra all'interno dei segmenti
                    if (intersection.isPointOnSegment(p3, p4)) {
                        return true;  // C'è intersezione
                    }
                } else {

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


                    if(intersection.isPointOnSegment(p3,p4) && intersection.isPointOnSegment(*this, p2))
                        return true;
                    else
                        return false;
                }
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



    int Triangle::areTrianglesDelaunay(Triangle& triangle1) {

        if (isAdjacent(triangle1)) {
            // Verifica se i punti del primo triangolo sono all'interno della circonferenza circoscritta del secondo triangolo
            bool point1InsideCircle = isInsideCircumcircle(triangle1.p1);
            bool point2InsideCircle = isInsideCircumcircle(triangle1.p2);
            bool point3InsideCircle = isInsideCircumcircle(triangle1.p3);

            // Verifica se i punti del secondo triangolo sono all'interno della circonferenza circoscritta del primo triangolo
            bool point4InsideCircle = triangle1.isInsideCircumcircle(p1);
            bool point5InsideCircle = triangle1.isInsideCircumcircle(p2);
            bool point6InsideCircle = triangle1.isInsideCircumcircle(p3);

            // L'ipotesi di Delaunay è verificata solo se tutti i punti sono all'esterno delle rispettive circonferenze circoscritte
            if ( !(point1InsideCircle || point2InsideCircle || point3InsideCircle ||
                     point4InsideCircle || point5InsideCircle || point6InsideCircle) ) {
                return 1;    // 1 = i triangoli sono di delaunay
            }
            else{
                return 0;    // 0 = i non triangoli sono di delaunay
            }
        } else {
            cerr << "Triangoli non adiacenti, non ha senso controllare l'ipotesi di Delaunay" << endl;
            return -1;
        }
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
            SortVertices();

            triangle1.p1= nonCommonVertex[0];
            triangle1.p2= nonCommonVertex[1];
            triangle1.p3= commonVertices[1];
            triangle1.SortVertices();

    }



    // Aggiorna la lista di adiacenza basata sulle relazioni tra i triangoli
    // Implementazione specifica in base alle regole della triangolazione

/*
void Triangulation::addAdjacentTriangle(int triangleId, int adjacentTriangleId) {
    triangles[triangleId].adjacentTriangles.push_back(adjacentTriangleId);
    triangles[adjacentTriangleId].adjacentTriangles.push_back(triangleId);
}
*/
bool Triangle::isPointOnEdge(const Point& Q) const {
    return Q.isPointOnSegment(p1, p2) || Q.isPointOnSegment(p1, p3) || Q.isPointOnSegment(p3, p2);
}

std::vector<Point> Triangle::PointOnEdge(const Point& Q){
    std::vector<Point> edge;

    if (Q.isPointOnSegment(p1,p2)){
        edge.push_back(p1);
        edge.push_back(p2);
    }
    else if (Q.isPointOnSegment(p1,p3)){
        edge.push_back(p1);
        edge.push_back(p3);
    }
    else if(Q.isPointOnSegment(p3,p2)){
        edge.push_back(p3);
        edge.push_back(p2);
    }
    else{
        Point point;
        edge.push_back(point);
        edge.push_back(point);

    }
    return edge;
}

int Triangulation::PointInsideTriangulation(const Point& Q) {
    // Itera attraverso tutti i triangoli della triangolazione
    for (const Triangle& triangle : DelunayTriangles) {
        if (!triangle.isInsideCircumcircle(Q)) {
            //return -1;
        }
        // Verifica se il punto Q è interno al triangolo corrente
        else if (triangle.isPointInsideTriangle(Q)|| triangle.isPointOnEdge(Q)) {

        return triangle.id;  //  restituisco l'id del triangolo a cui è interno
        }
    }

    return -1;  // Il punto non è interno alla triangolazione
}

Triangle Triangulation::findTriangleById(const unsigned int& triangleId) {
    for (const Triangle& triangle : DelunayTriangles) {
        if (triangle.id == triangleId) {
            return triangle;
        }
    }
    // Restituisce un triangolo vuoto se l'ID non viene trovato
    return Triangle();
}

void Triangulation::removeAdjacentTriangle(unsigned int triangleId1, unsigned int triangleId2) {
    // Rimuovi il triangoloId2 dalla lista di adiacenza del triangoloId1
    auto it = std::find(adjacencyList[triangleId1].begin(), adjacencyList[triangleId1].end(), triangleId2);
    if (it != adjacencyList[triangleId1].end()) {
        adjacencyList[triangleId1].erase(it);
    }

    // Rimuovi il triangoloId1 dalla lista di adiacenza del triangoloId2
    it = std::find(adjacencyList[triangleId2].begin(), adjacencyList[triangleId2].end(), triangleId1);
    if (it != adjacencyList[triangleId2].end()) {
        adjacencyList[triangleId2].erase(it);
    }
}


void Triangulation::createSubtriangulation(const Point& Q, const unsigned int& triangleId) {
    // Ottieni il triangolo T dalla lista dei triangoli
    Triangle& T = DelunayTriangles[triangleId];
    Point a = T.p1;
    Point b = T.p2;
    Point c = T.p3;
    unsigned int Tid = T.id;

    // Crea il nuovo triangolo collegando il punto Q ai vertici di T
    Triangle newTriangle(Q, a, b, triangleId);
    newTriangle.SortVertices();

    // Aggiungi il nuovo triangolo alla lista dei triangoli
    DelunayTriangles[triangleId] = newTriangle;

    // Ottieni la lista di adiacenza originale di T
    std::vector<unsigned int> originalAdjacencyList = adjacencyList[Tid];


    // Rimuovi i triangoli adiacenti che non sono più adiacenti a newTriangle
    for (unsigned int adjTriangleId : adjacencyList[triangleId]) {
        Triangle& adjTriangle = DelunayTriangles[adjTriangleId];
        if (!adjTriangle.isAdjacent(newTriangle)) {
            removeAdjacentTriangle(triangleId, adjTriangleId);
        }
    }


/*
    // Aggiorna la lista di adiacenza per il nuovo triangolo
    for (unsigned int adjTriangleId : adjacencyList[triangleId]) {
        Triangle& adjTriangle = DelunayTriangles[adjTriangleId];
        if (adjTriangle.isAdjacent(newTriangle) && !newTriangle.isAdjacent(adjTriangle)) {
            addAdjacentTriangle(triangleId, adjTriangleId);
        }
    }
*/


    // Aggiungi gli altri due triangoli della sottotriangolazione alla lista dei triangoli
    unsigned int newTriangleId = getMaxTriangleId();
    unsigned int t1Id = newTriangleId + 1;
    unsigned int t2Id = newTriangleId + 2;
    Triangle t1(Q, b, c, t1Id);
    Triangle t2(Q, c, a, t2Id);
    t1.SortVertices();
    t2.SortVertices();
    addTriangle(t1);
    addTriangle(t2);
    incrementMaxTriangleId(2);

    // Controlla l'adiacenza di t1 con i triangoli originariamente adiacenti a T
        for (unsigned int adjTriangleId : originalAdjacencyList) {
            Triangle& adjTriangle = DelunayTriangles[adjTriangleId];
            if (adjTriangle.isAdjacent(t1)) {
                addAdjacentTriangle(adjTriangleId, t1Id);

            }
        }

        // Controlla l'adiacenza di t2 con i triangoli originariamente adiacenti a T
        for (unsigned int adjTriangleId : originalAdjacencyList) {
            Triangle& adjTriangle = DelunayTriangles[adjTriangleId];
            if (adjTriangle.isAdjacent(t2) ) {
                addAdjacentTriangle(adjTriangleId, t2Id);

            }
        }

    addAdjacentTriangle(triangleId, t1Id);
    addAdjacentTriangle(triangleId, t2Id);
    addAdjacentTriangle(t1Id, t2Id);
}


void Triangulation::connectPointToVertices(const Point& Q) {
    unsigned int n = DelunayTriangles.size();
    std::vector<Triangle> newTriangles;  // Triangoli da aggiungere

    // Genera tutti i triangoli possibili
    for (size_t i = 0; i < n; i++) {
        Triangle& triangle = DelunayTriangles[i];
        std::vector<Point> intersection = triangle.findIntersection(Q);

        if (intersection.empty()) {
            // Genera i triangoli combinando il nuovo punto con i vertici del triangolo esistente
            Triangle t1(Q, triangle.p1, triangle.p2, getMaxTriangleId());
            t1.SortVertices();
            newTriangles.push_back(t1);

            Triangle t2(Q, triangle.p1, triangle.p3, getMaxTriangleId());
            t2.SortVertices();
            newTriangles.push_back(t2);

            Triangle t3(Q, triangle.p2, triangle.p3, getMaxTriangleId());
            t3.SortVertices();
            newTriangles.push_back(t3);
        }
        else {
                  if (triangle.p1.isPointInVector(intersection)) {
                      unsigned int id = getMaxTriangleId();
                      Triangle t1(Q, triangle.p2, triangle.p3, id);
                      ++id;
                      t1.SortVertices();
                      if (t1.calculateArea() > 0.0) {
                          newTriangles.push_back(t1);
                      }
                  }
                  else if (triangle.p2.isPointInVector(intersection)) {
                      unsigned int id = getMaxTriangleId();
                      Triangle t2(Q, triangle.p1, triangle.p3, id);
                      ++id;
                      t2.SortVertices();
                      if (t2.calculateArea() > 0.0) {
                          newTriangles.push_back(t2);
                      }
                  }
                  else {
                      unsigned int id = getMaxTriangleId();
                      Triangle t3(Q, triangle.p1, triangle.p2, id);
                      ++id;
                      t3.SortVertices();
                      if (t3.calculateArea() > 0.0) {
                          newTriangles.push_back(t3);
                      }
                  }
        }
    }

    // Verifica e aggiunge solo i triangoli che non si sovrappongono ai triangoli esistenti
    for (Triangle& newTriangle : newTriangles) {
        bool isOverlapping = false;

        for (Triangle& existingTriangle : DelunayTriangles) {
            if (newTriangle.checkTriangleOverlap(existingTriangle)) {
                isOverlapping = true;
                break;
            }
        }

        if (!isOverlapping) {
            newTriangle.id = getMaxTriangleId() +1;
            addTriangle(newTriangle);
            incrementMaxTriangleId();

            // Aggiorna l'adiacenza dei triangoli adiacenti
            unsigned int n = DelunayTriangles.size();
            for (size_t i = 0; i < n; i++) {
                Triangle adjacentTriangle = DelunayTriangles[i];
                if (adjacentTriangle.isAdjacent(newTriangle) && adjacentTriangle.id != newTriangle.id){
                    addAdjacentTriangle(adjacentTriangle.id, newTriangle.id);
                }
            }
        }
    }
}

bool Triangle::checkTriangleOverlap(Triangle& other) {
    if(isAdjacent(other)){
    // Verifica l'intersezione dei segmenti tra i triangoli
        int intersectCount = 0;

        // Verifica l'intersezione dei segmenti tra i lati del triangolo corrente e gli altri lati
        if ( p1.doSegmentsIntersect(p2, other.p2, other.p3) && ((p1 != other.p2) && (p1 != other.p3) && (p2 != other.p2) && (p2 != other.p3))) intersectCount++;
        if ( p1.doSegmentsIntersect(p2, other.p1, other.p2) && ((p1 != other.p1) && (p1 != other.p2) && (p2 != other.p1) && (p2 != other.p2))) intersectCount++;
        if ( p1.doSegmentsIntersect(p2, other.p1, other.p3) && ((p1 != other.p1) && (p1 != other.p3) && (p2 != other.p1) && (p2 != other.p3))) intersectCount++;
        if ( p2.doSegmentsIntersect(p3, other.p2, other.p3) && ((p2 != other.p2) && (p2 != other.p3) && (p3 != other.p2) && (p3 != other.p3))) intersectCount++;
        if ( p2.doSegmentsIntersect(p3, other.p1, other.p2) && ((p2 != other.p1) && (p2 != other.p2) && (p3 != other.p1) && (p3 != other.p2))) intersectCount++;
        if ( p2.doSegmentsIntersect(p3, other.p1, other.p3) && ((p2 != other.p1) && (p2 != other.p3) && (p3 != other.p1) && (p3 != other.p3))) intersectCount++;
        if ( p3.doSegmentsIntersect(p1, other.p2, other.p3) && ((p3 != other.p2) && (p3 != other.p3) && (p1 != other.p2) && (p1 != other.p3))) intersectCount++;
        if ( p3.doSegmentsIntersect(p1, other.p1, other.p2) && ((p3 != other.p1) && (p3 != other.p2) && (p1 != other.p1) && (p1 != other.p2))) intersectCount++;
        if ( p3.doSegmentsIntersect(p1, other.p1, other.p3) && ((p3 != other.p1) && (p3 != other.p3) && (p1 != other.p1) && (p1 != other.p3))) intersectCount++;
        // Restituisci true se ci sono almeno due intersezioni, altrimenti false
        return intersectCount >= 1;}
    else{return 0;}
}

bool Triangulation::isBoundaryTriangle(const Triangle& triangle) {
    unsigned int numAdjacentTriangles = adjacencyList[triangle.id].size();

    // If the number of adjacent triangles is less than three, the triangle is a boundary triangle
    if (numAdjacentTriangles < 3) {
        return true; // The triangle is a boundary triangle
    }

    return false; // The triangle is not a boundary triangle
}

vector<int> Triangulation::connectPointOnEdge(const Triangle& t, const Point& Q, const vector<Point>& edge) {
    vector<unsigned int> oldadj = adjacencyList[t.id];
    Triangle t1;
    Triangle t2;
    if (!t.p1.isPointInVector(edge)) {
        t1 = Triangle(Q, t.p1, t.p2, t.id);
        t1.SortVertices();
        t2 = Triangle(Q, t.p1, t.p3, getMaxTriangleId() + 1);
        t2.SortVertices();

    } else if (!t.p2.isPointInVector(edge)) {
        t1 = Triangle(Q, t.p2, t.p1, t.id);
        t1.SortVertices();
        t2 = Triangle(Q, t.p2, t.p3, getMaxTriangleId() + 1);
        t2.SortVertices();

    } else {
        t1 = Triangle(Q, t.p3, t.p1, t.id);
        t1.SortVertices();
        t2 = Triangle(Q, t.p3, t.p2, getMaxTriangleId() + 1);
        t2.SortVertices();

    }

    // Replace the previous triangle with t1 in DelunayTriangles
    DelunayTriangles[t.id] = t1;

    // Add the new triangle t2 to DelunayTriangles
    addTriangle(t2);
    incrementMaxTriangleId();

    // Update adjacency
    addAdjacentTriangle(t1.id, t2.id);


    // Update adjacency for t2
    for (unsigned int adjId : oldadj) {

            if (DelunayTriangles[adjId].isAdjacent(t2)) {

                    addAdjacentTriangle(adjId, t2.id);

        }
    }
    vector<int> ids;
    int id1=t1.id;
    int id2=t2.id;
    ids.push_back(id1);
    ids.push_back(id2);
    return ids;
}


void Triangulation::connectPointOnEdgeToTriangulation(const Triangle& triangle, const Point& Q, const vector<Point>& edge, const MinMax& minMax) {


    // Connect the point to the vertex opposite the edge
    vector<int>ids;
    ids=connectPointOnEdge(triangle, Q, edge); // Create two subtriangles from one triangle

    // If the triangle is not a boundary triangle, check the other adjacent triangle on the same edge
    vector<int>ids2;
    if (!isBoundaryEdge(triangle, edge, minMax)) {


            unsigned int n = adjacencyList[triangle.id].size();
            for (size_t i = 0; i < n; i++) {
                Triangle t = DelunayTriangles[ adjacencyList[triangle.id][i]];

                vector<Point> edge2 = t.PointOnEdge(Q); //lato in cui giace il punto Q

                if (edge2==edge){

                    ids2=connectPointOnEdge(t, Q,edge); // Create two subtriangles from one triangle


                        }
                    }

                }
    std::vector<int> mergedIds;
    mergedIds.reserve(ids.size() + ids2.size());
    std::merge(ids.begin(), ids.end(), ids2.begin(), ids2.end(), std::back_inserter(mergedIds));
    for (size_t i = 0; i < mergedIds.size(); i++) {
        unsigned int triangleId = mergedIds[i];
        int n=adjacencyList[triangleId].size();
        vector<unsigned int> oldaj=adjacencyList[triangleId];
        for (size_t k = 0; k < n; k++) {
            unsigned int adjId = oldaj[k];
            Triangle& adjTriangle = DelunayTriangles[adjId];

            if (!DelunayTriangles[triangleId].isAdjacent(adjTriangle)) {
                removeAdjacentTriangle(triangleId, adjId);
            }
        }
    }
if (mergedIds.size()==4){
if(DelunayTriangles[mergedIds[0]].isAdjacent(DelunayTriangles[mergedIds[2]])){
    if(!isIdInAdjacency(mergedIds[0],adjacencyList[mergedIds[2]])){
    addAdjacentTriangle(mergedIds[0],mergedIds[2]);}}
if(DelunayTriangles[mergedIds[0]].isAdjacent(DelunayTriangles[mergedIds[3]])){
    if(!isIdInAdjacency(mergedIds[0],adjacencyList[mergedIds[3]])){
    addAdjacentTriangle(mergedIds[0],mergedIds[3]);}}
if(DelunayTriangles[mergedIds[1]].isAdjacent(DelunayTriangles[mergedIds[2]])){
    if(!isIdInAdjacency(mergedIds[1],adjacencyList[mergedIds[2]])){
    addAdjacentTriangle(mergedIds[1],mergedIds[2]);}}
if(DelunayTriangles[mergedIds[1]].isAdjacent(DelunayTriangles[mergedIds[3]])){
    if(!isIdInAdjacency(mergedIds[1],adjacencyList[mergedIds[3]])){
    addAdjacentTriangle(mergedIds[1],mergedIds[3]);}}
}
}





void Triangulation::addPointToTriangulation(const Point& Q, const MinMax& minMax) {

    if (PointInsideTriangulation(Q) >= 0) {
        int triangleId = PointInsideTriangulation(Q);
        vector<Point> edge;
        Triangle triangle(findTriangleById(triangleId));
        edge = triangle.PointOnEdge(Q); //lato in cui giace il punto Q
        if (!(edge[0].x==0.0 && edge[0].y==0.0 && edge[1].x==0.0 && edge[1].y==0.0)){  //se è interno è il vettore (0,0)
            connectPointOnEdgeToTriangulation(triangle,Q,edge,minMax);
        }
        else
        createSubtriangulation(Q, triangleId);
    } else {
        connectPointToVertices(Q);
    }
}

void Triangulation::flipTrianglesIfNotDelaunay() {
    bool hasNonDelaunayTriangles = true;

    while (hasNonDelaunayTriangles) {
        hasNonDelaunayTriangles = false;

        for (size_t i = 0; i < DelunayTriangles.size(); ++i) {
            Triangle& triangle = DelunayTriangles[i];

            // Crea una copia della lista di adiacenza del triangolo
            unsigned int id = triangle.id;
            vector<unsigned int> adjacency(adjacencyList[id].begin(), adjacencyList[id].end());
            unsigned int n = adjacency.size();

            bool triangleFlipped = false;

            for (size_t j = 0; j < n; ++j) {
                unsigned int adjTriangleId = adjacency[j];
                Triangle& adjacentTriangle = DelunayTriangles[adjTriangleId];

                if (adjacentTriangle.areTrianglesDelaunay(triangle) == 0) {
                    flipAndUpdateAdjacency(triangle, adjacentTriangle);
                    triangleFlipped = true;
                    break;
                }
            }

            if (triangleFlipped) {
                hasNonDelaunayTriangles = true;
                flipTrianglesIfNotDelaunay();  // Recursive call to restart the flipping process
                break;
            }
        }
    }
}


void Triangulation::updateAdjacency(Triangle& triangle, Triangle& adjacentTriangle) {
    // Concatena i due vettori
    vector<unsigned int> oldAdjacency;
    vector<unsigned int> adjacency(adjacencyList[triangle.id].begin(), adjacencyList[triangle.id].end());
    vector<unsigned int> adjacency2(adjacencyList[adjacentTriangle.id].begin(), adjacencyList[adjacentTriangle.id].end());
    oldAdjacency.reserve(adjacency.size() + adjacencyList[adjacentTriangle.id].size());
    std::merge(adjacency.begin(), adjacency.end(), adjacencyList[adjacentTriangle.id].begin(), adjacencyList[adjacentTriangle.id].end(), std::back_inserter(oldAdjacency));

                // Rimuovi i duplicati
                auto last = std::unique(oldAdjacency.begin(), oldAdjacency.end());
                oldAdjacency.erase(last, oldAdjacency.end());

                // Aggiorna le liste di adiacenza dei triangoli dopo il flip

                // AGGIORNO LE MODIFICHE PER TRIANGLE
                for (size_t i = 0; i < oldAdjacency.size(); i++) {
                    unsigned int adjId = oldAdjacency[i];
                    Triangle& adjTriangle = DelunayTriangles[adjId];
                    if (adjTriangle.isAdjacent(triangle) && (adjId != triangle.id) ){
                        if( !(isIdInAdjacency(adjId, adjacency)))
                            addAdjacentTriangle(adjId,triangle.id);
                        }
                    else{
                        if( isIdInAdjacency(adjId, adjacency) && (adjId != triangle.id))
                            removeAdjacentTriangle(adjId,triangle.id);
                    }
                }

                // Aggiornamento per adjacentTriangle
                for (size_t i = 0; i < oldAdjacency.size(); i++) {
                    unsigned int adjId = oldAdjacency[i];
                    Triangle& adjTriangle = DelunayTriangles[adjId];
                    if (adjTriangle.isAdjacent(adjacentTriangle) && (adjId != adjacentTriangle.id)) {
                        if (!isIdInAdjacency(adjId, adjacency2))
                            addAdjacentTriangle(adjId, adjacentTriangle.id);
                    } else {
                        if (isIdInAdjacency(adjId, adjacency2) && (adjId != adjacentTriangle.id))
                            removeAdjacentTriangle(adjId, adjacentTriangle.id);
                    }
                }

}

bool Triangulation::isIdInAdjacency(unsigned int triangleId, vector<unsigned int> adjacency) {
    for (unsigned int adjId : adjacency) {
        if (adjId == triangleId) {
            return true;
        }
    }
return false;
}
/*
// QUESTA SOTTO NON L'HO UTILIZZATA PER ORA
bool Triangulation::areAllTrianglesDelaunay() {
    for (Triangle& triangle : DelunayTriangles) {
        for (unsigned int adjTriangleId : triangle.adjacentTriangles) {
            Triangle& adjacentTriangle = DelunayTriangles[adjTriangleId];

            if (triangle.areTrianglesDelaunay(adjacentTriangle) == 0) {
                return false;
            }
        }
    }
    return true;
}
*/

void Triangulation::flipAndUpdateAdjacency(Triangle& triangle, Triangle& adjacentTriangle) {


        if (triangle.areTrianglesDelaunay(adjacentTriangle) == 0) {

            // Esegui il flip dei triangoli
            triangle.flip(adjacentTriangle);

            // Aggiorna i triangoli nella lista DelunayTriangles
            DelunayTriangles[triangle.id] = triangle;
            DelunayTriangles[adjacentTriangle.id] = adjacentTriangle;

            // Aggiorna le liste di adiacenza dei triangoli
            updateAdjacency(triangle, adjacentTriangle);
        }

}


bool Triangulation::isBoundaryEdge(const Triangle& triangle, const std::vector<Point>& edge ,const MinMax& minMax) {
    // Verifica se il vettore edge rappresenta un lato di bordo della triangolazione
    // Verifica se entrambi i punti del lato sono punti di bordo
    // Restituisce true se il lato è di bordo, altrimenti false

    for (const Point& point : edge) {
        // Controlla se il punto non è un punto di bordo
        if (!isBoundaryPoint(point,minMax)) {
            return false;
        }
    }

    return true;
}





bool Triangulation::isBoundaryPoint(const Point& point, const MinMax& minMax) {
    // Verifica se il punto point è un punto di bordo della triangolazione
    // Restituisci true se il punto è di bordo, altrimenti false

    if (point.x == minMax.minX || point.x == minMax.maxX || point.y == minMax.minY || point.y == minMax.maxY) {
        return true;
    }

    return false;
}




MinMax findMinMax(const std::vector<Point>& points) {
    MinMax minMax;

    if (points.empty()) {
        // Restituisci valori di default se il vettore di punti è vuoto
        minMax.minX = 0.0;
        minMax.maxX = 0.0;
        minMax.minY = 0.0;
        minMax.maxY = 0.0;
        return minMax;
    }

    // Crea una copia dei punti da ordinare
    std::vector<Point> sortedPoints = points;

    // Ordina i punti per coordinata x utilizzando MergeSort
    SortLibrary::MergeSortx(sortedPoints, 0, sortedPoints.size() - 1);

    // Trova il valore minimo e massimo di x
    minMax.minX = sortedPoints[0].x;
    minMax.maxX = sortedPoints[sortedPoints.size() - 1].x;

    // Ordina i punti per coordinata y utilizzando MergeSort
    SortLibrary::MergeSorty(sortedPoints, 0, sortedPoints.size() - 1);

    // Trova il valore minimo e massimo di y
    minMax.minY = sortedPoints[0].y;
    minMax.maxY = sortedPoints[sortedPoints.size() - 1].y;

    return minMax;
}

Triangulation DelunayTriangulation(const std::vector<Point>& points, const MinMax& minMax) {
    Triangle triangle_default;
    Triangle triangle_max;



    triangle_max = triangle_default.findMaximumTriangle(points);


    std::vector<Triangle> DelunayTriangles;


    Triangulation triangulation(DelunayTriangles);
    triangulation.addTriangle(triangle_max);

    std::vector<Point> pointsCopy;
    pointsCopy = points;

    pointsCopy.erase(std::remove(pointsCopy.begin(), pointsCopy.end(), triangle_max.p1), pointsCopy.end());
    pointsCopy.erase(std::remove(pointsCopy.begin(), pointsCopy.end(), triangle_max.p2), pointsCopy.end());
    pointsCopy.erase(std::remove(pointsCopy.begin(), pointsCopy.end(), triangle_max.p3), pointsCopy.end());
    for (const auto& point : pointsCopy) {
        cout << point << endl;
    }

    for (size_t i = 0; i <7; i++) {
        triangulation.addPointToTriangulation(pointsCopy[i],minMax);

        unsigned int n = triangulation.DelunayTriangles.size();
        for (size_t j = 0; j < n ; j++) {
            //Triangle& triangle = triangulation.DelunayTriangles[j];
            // const std::vector<unsigned int>& adjacency = triangulation.adjacencyList[triangle.id];
            triangulation.flipTrianglesIfNotDelaunay();
        }

        // Rimuovi i triangoli duplicati
        auto last = std::unique(triangulation.DelunayTriangles.begin(), triangulation.DelunayTriangles.end());
        triangulation.DelunayTriangles.erase(last, triangulation.DelunayTriangles.end());
    }

    return triangulation;
}}
