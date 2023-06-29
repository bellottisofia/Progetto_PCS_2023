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

bool Point::isPointInVector(const std::vector<Point>& points) const {
    return std::find(points.begin(), points.end(), *this) != points.end();
}




Triangle::Triangle(const Point& point1, const Point& point2, const Point& point3, const unsigned int& id):
    p1(point1),
    p2(point2),
    p3(point3),
    id(id)
{}

double Triangle::calculateArea()const{
    return abs(0.5 * ((getVertex1().x - getVertex3().x) * (getVertex2().y - getVertex3().y) - (getVertex1().y - getVertex3().y) * (getVertex2().x - getVertex3().x)));
}


bool Triangle::IsVerticesSort()
{
  double dx1 = getVertex2().x - getVertex1().x;
  double dy1 = getVertex2().y - getVertex1().y;
  double dx2 = getVertex3().x - getVertex1().x;
  double dy2 = getVertex3().y - getVertex1().y;

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
        getVertex1();
        getVertex3();
        getVertex1();
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
    //Triangle t1(getVertex1(), getVertex2(), getVertex3());
    Triangle t2(Q, getVertex2(), getVertex3(),getId()+1);
    Triangle t3(getVertex1(), Q, getVertex3(),getId()+2);
    Triangle t4(getVertex1(), getVertex2(), Q,getId()+3);

    //double areaABC = t1.calculateArea();
    double areaABC = calculateArea();
    double areaPBC = t2.calculateArea();
    double areaPCA = t3.calculateArea();
    double areaPAB = t4.calculateArea();

    // Verifica se la somma delle aree dei triangoli interni è uguale all'area totale del triangolo ABC
    return abs(areaABC - (areaPBC + areaPCA + areaPAB)) < 1e-12;

    }
}


bool Triangle::isVertexShared(const Point& vertex) const {
    return (vertex == getVertex1()) || (vertex == getVertex2()) || (vertex == getVertex3());
}




bool Triangle::isAdjacent(const Triangle& other) {
    int sharedVertices = 0;

    if (isVertexShared(other.getVertex1()))
        sharedVertices++;
    if (isVertexShared(other.getVertex2()))
        sharedVertices++;
    if (isVertexShared(other.getVertex3()))
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


    Triangulation::Triangulation(const std::vector<Triangle>& triangles)
    : DelaunayTriangles(triangles) {}



    void Triangulation::addTriangle(const Triangle& triangleAgg) {

        DelaunayTriangles.push_back(triangleAgg);
        adjacencyList.emplace(triangleAgg.getId(), std::vector<unsigned int>());// Inizializza la lista di adiacenza per il triangolo
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
        if (std::abs(crossProduct) > 1e-12) {
            return false;  // I punti non hanno la stessa pendenza, quindi Q non giace sul segmento
        }

        // Controlla se il punto Q rientra nell'intervallo tra i punti p e r
        if (x < std::min(p.x, r.x) || x > std::max(p.x, r.x) || y < std::min(p.y, r.y) || y > std::max(p.y, r.y)) {
            return false;  // Il punto Q non rientra nell'intervallo tra p e r
        }

        return true;  // Il punto Q giace sul segmento p-r
    }




    bool Point::doSegmentsIntersect( const Point& p2, const Point& p3, const Point& p4) const{
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

    const Point Triangle::getVertex1() const {
            return p1;
        }
    void Triangle::setVertex1(const Point& newVertex1) { p1 = newVertex1; }
    const Point Triangle::getVertex2() const {
        return p2;
    }
    void Triangle::setVertex2(const Point& newVertex2) { p2 = newVertex2; }
    const Point Triangle::getVertex3() const {
        return p3;
    }
    void Triangle::setVertex3(const Point& newVertex3) { p3 = newVertex3; }
    unsigned int Triangle::getId() const {
        return id;
    }
    void Triangle::setId(unsigned int newId) const { id = newId; }

    std::vector<Point> Triangle::findIntersection(const Point& q) {
        std::vector<Point> intersections;

        // Controlla l'intersezione con i lati del triangolo
        if (getVertex3().doSegmentsIntersect(getVertex1(), getVertex2(), q)) {
            intersections.push_back(getVertex2());
            intersections.push_back(q);
        }
        else if (getVertex1().doSegmentsIntersect(getVertex2(), getVertex3(), q)) {
            intersections.push_back(getVertex3());
            intersections.push_back(q);
        }
        else if (getVertex2().doSegmentsIntersect(getVertex3(), getVertex1(), q)) {
            intersections.push_back(getVertex1());
            intersections.push_back(q);
        }

        return intersections;
    }




    int Triangle::areTrianglesDelaunay(Triangle& triangle1) {

        if (isAdjacent(triangle1)) {
            // Verifica se i punti del primo triangolo sono all'interno della circonferenza circoscritta del secondo triangolo
            bool point1InsideCircle = isInsideCircumcircle(triangle1.getVertex1());
            bool point2InsideCircle = isInsideCircumcircle(triangle1.getVertex2());
            bool point3InsideCircle = isInsideCircumcircle(triangle1.getVertex3());

            // Verifica se i punti del secondo triangolo sono all'interno della circonferenza circoscritta del primo triangolo
            bool point4InsideCircle = triangle1.isInsideCircumcircle(getVertex1());
            bool point5InsideCircle = triangle1.isInsideCircumcircle(getVertex2());
            bool point6InsideCircle = triangle1.isInsideCircumcircle(getVertex3());

            // L'ipotesi di Delaunay è verificata solo se tutti i punti sono all'esterno delle rispettive circonferenze circoscritte
            if ( !(point1InsideCircle || point2InsideCircle || point3InsideCircle ||
                     point4InsideCircle || point5InsideCircle || point6InsideCircle) ) {
                return 1;    // 1 = i triangoli sono di delaunay
            }
            else{
                return 0;    // 0 = i non triangoli sono di delaunay
            }
        } else {
                return -1;
        }
    }

void Triangle::flip(Triangle& triangle1){
            // Verifica se i due triangoli sono adiacenti
            int commonVerticesCount = 0;
            Point commonVertices[2];

            // Trova i vertici comuni tra this e triangle1

            if (isVertexShared(triangle1.getVertex1())) {
                commonVertices[commonVerticesCount++] = triangle1.getVertex1();
            }
            if (isVertexShared(triangle1.getVertex2())) {
                commonVertices[commonVerticesCount++] = triangle1.getVertex2();
            }
            if (isVertexShared(triangle1.getVertex3())) {
                commonVertices[commonVerticesCount++] = triangle1.getVertex3();
            }

            //int count = 0;
            Point nonCommonVertex[2];
            if (commonVerticesCount == 2) {
                // Scambia i vertici dei due triangoli


                if (getVertex1() != commonVertices[0] && getVertex1() != commonVertices[1]) {
                    nonCommonVertex[0] = getVertex1();
                }
                else if (getVertex2() != commonVertices[0] && getVertex2() != commonVertices[1]) {
                    nonCommonVertex[0] = getVertex2();
                }
                else {
                    nonCommonVertex[0] = getVertex3();
                }
                if (triangle1.getVertex1() != commonVertices[0] && triangle1.getVertex1() != commonVertices[1]) {
                    nonCommonVertex[1] = triangle1.getVertex1();
                }
                else if (triangle1.getVertex2() != commonVertices[0] && triangle1.getVertex2() != commonVertices[1]) {
                    nonCommonVertex[1] = triangle1.getVertex2();
                }
                else {
                    nonCommonVertex[1] = triangle1.getVertex3();
                }

            }


            p1=nonCommonVertex[0];
            p2= commonVertices[0];
            p3=nonCommonVertex[1];
            SortVertices();

            triangle1.p1=nonCommonVertex[0];
            triangle1.p2=nonCommonVertex[1];
            triangle1.p3=commonVertices[1];
            triangle1.SortVertices();

    }




bool Triangle::isPointOnEdge(const Point& Q)const {
    return Q.isPointOnSegment(getVertex1(), getVertex2()) || Q.isPointOnSegment(getVertex1(), getVertex3()) || Q.isPointOnSegment(getVertex3(), getVertex2());
}

std::vector<Point> Triangle::PointOnEdge(const Point& Q){
    std::vector<Point> edge;

    if (Q.isPointOnSegment(getVertex1(),getVertex2())){
        edge.push_back(getVertex1());
        edge.push_back(getVertex2());
    }
    else if (Q.isPointOnSegment(getVertex1(),getVertex3())){
        edge.push_back(getVertex1());
        edge.push_back(getVertex3());
    }
    else if(Q.isPointOnSegment(getVertex3(),getVertex2())){
        edge.push_back(getVertex3());
        edge.push_back(getVertex2());
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
    for (const Triangle& triangle : DelaunayTriangles) {
        if (!triangle.isInsideCircumcircle(Q)) {
            //return -1;
        }
        // Verifica se il punto Q è interno al triangolo corrente
        else if (triangle.isPointInsideTriangle(Q)|| triangle.isPointOnEdge(Q)) {

        return triangle.getId();  //  restituisco l'id del triangolo a cui è interno
        }
    }

    return -1;  // Il punto non è interno alla triangolazione
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
    Triangle& T = DelaunayTriangles[triangleId];
    Point a = T.getVertex1();
    Point b = T.getVertex2();
    Point c = T.getVertex3();
    unsigned int Tid = T.getId();

    // Crea il nuovo triangolo collegando il punto Q ai vertici di T
    Triangle newTriangle(Q,a,b, triangleId);
    newTriangle.SortVertices();

    // Aggiungi il nuovo triangolo alla lista dei triangoli
    DelaunayTriangles[triangleId] = newTriangle;

    // Ottieni la lista di adiacenza originale di T
    std::vector<unsigned int> originalAdjacencyList = adjacencyList[Tid];


    // Rimuovi i triangoli adiacenti che non sono più adiacenti a newTriangle
    for (unsigned int adjTriangleId : adjacencyList[triangleId]) {
        Triangle& adjTriangle = DelaunayTriangles[adjTriangleId];
        if (!adjTriangle.isAdjacent(newTriangle)) {
            removeAdjacentTriangle(triangleId, adjTriangleId);
        }
    }


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
            Triangle& adjTriangle = DelaunayTriangles[adjTriangleId];
            if (adjTriangle.isAdjacent(t1)) {
                addAdjacentTriangle(adjTriangleId, t1Id);

            }
        }

        // Controlla l'adiacenza di t2 con i triangoli originariamente adiacenti a T
        for (unsigned int adjTriangleId : originalAdjacencyList) {
            Triangle& adjTriangle = DelaunayTriangles[adjTriangleId];
            if (adjTriangle.isAdjacent(t2) ) {
                addAdjacentTriangle(adjTriangleId, t2Id);

            }
        }

    addAdjacentTriangle(triangleId, t1Id);
    addAdjacentTriangle(triangleId, t2Id);
    addAdjacentTriangle(t1Id, t2Id);
}


void Triangulation::connectPointToVertices(const Point& Q) {
    unsigned int n = DelaunayTriangles.size();
    std::vector<Triangle> newTriangles;  // Triangoli da aggiungere

    // Genera tutti i triangoli possibili
    for (size_t i = 0; i < n; i++) {
        Triangle& triangle = DelaunayTriangles[i];
        std::vector<Point> intersection = triangle.findIntersection(Q);

        if (intersection.empty()) {
            // Genera i triangoli combinando il nuovo punto con i vertici del triangolo esistente
            Triangle t1(Q, triangle.getVertex1(), triangle.getVertex2(), getMaxTriangleId());
            t1.SortVertices();
            newTriangles.push_back(t1);

            Triangle t2(Q, triangle.getVertex1(), triangle.getVertex3(), getMaxTriangleId());
            t2.SortVertices();
            newTriangles.push_back(t2);

            Triangle t3(Q, triangle.getVertex2(), triangle.getVertex3(), getMaxTriangleId());
            t3.SortVertices();
            newTriangles.push_back(t3);
        }
        else {
                  if (triangle.getVertex1().isPointInVector(intersection)) {
                      unsigned int id = getMaxTriangleId();
                      Triangle t1(Q, triangle.getVertex2(), triangle.getVertex3(), id);
                      ++id;
                      t1.SortVertices();
                      if (t1.calculateArea() > 0.0) {
                          newTriangles.push_back(t1);
                      }
                  }
                  else if (triangle.getVertex2().isPointInVector(intersection)) {
                      unsigned int id = getMaxTriangleId();
                      Triangle t2(Q, triangle.getVertex1(), triangle.getVertex3(), id);
                      ++id;
                      t2.SortVertices();
                      if (t2.calculateArea() > 0.0) {
                          newTriangles.push_back(t2);
                      }
                  }
                  else {
                      unsigned int id = getMaxTriangleId();
                      Triangle t3(Q, triangle.getVertex1(), triangle.getVertex2(), id);
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

        for (Triangle& existingTriangle : DelaunayTriangles) {
            if (newTriangle.checkTriangleOverlap(existingTriangle)) {
                isOverlapping = true;
                break;
            }
        }

        if (!isOverlapping) {

            newTriangle.setId(getMaxTriangleId() +1);
                   addTriangle(newTriangle);
                   incrementMaxTriangleId();
            // Aggiorna l'adiacenza dei triangoli adiacenti
            unsigned int n = DelaunayTriangles.size();
            for (size_t i = 0; i < n; i++) {
                Triangle adjacentTriangle = DelaunayTriangles[i];
                if (adjacentTriangle.isAdjacent(newTriangle) && adjacentTriangle.getId() != newTriangle.getId()){
                    addAdjacentTriangle(adjacentTriangle.getId(), newTriangle.getId());
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
            if ( getVertex1().doSegmentsIntersect(getVertex2(), other.getVertex2(), other.getVertex3()) && ((getVertex1() != other.getVertex2()) && (getVertex1() != other.getVertex3()) && (getVertex2() != other.getVertex2()) && (getVertex2() != other.getVertex3()))) intersectCount++;
            if ( getVertex1().doSegmentsIntersect(getVertex2(), other.getVertex1(), other.getVertex2()) && ((getVertex1() != other.getVertex1()) && (getVertex1() != other.getVertex2()) && (getVertex2() != other.getVertex1()) && (getVertex2() != other.getVertex2()))) intersectCount++;
            if ( getVertex1().doSegmentsIntersect(getVertex2(), other.getVertex1(), other.getVertex3()) && ((getVertex1() != other.getVertex1()) && (getVertex1() != other.getVertex3()) && (getVertex2() != other.getVertex1()) && (getVertex2() != other.getVertex3()))) intersectCount++;
            if ( getVertex2().doSegmentsIntersect(getVertex3(), other.getVertex2(), other.getVertex3()) && ((getVertex2() != other.getVertex2()) && (getVertex2() != other.getVertex3()) && (getVertex3() != other.getVertex2()) && (getVertex3() != other.getVertex3()))) intersectCount++;
            if ( getVertex2().doSegmentsIntersect(getVertex3(), other.getVertex1(), other.getVertex2()) && ((getVertex2() != other.getVertex1()) && (getVertex2() != other.getVertex2()) && (getVertex3() != other.getVertex1()) && (getVertex3() != other.getVertex2()))) intersectCount++;
            if ( getVertex2().doSegmentsIntersect(getVertex3(), other.getVertex1(), other.getVertex3()) && ((getVertex2() != other.getVertex1()) && (getVertex2() != other.getVertex3()) && (getVertex3() != other.getVertex1()) && (getVertex3() != other.getVertex3()))) intersectCount++;
            if ( getVertex3().doSegmentsIntersect(getVertex1(), other.getVertex2(), other.getVertex3()) && ((getVertex3() != other.getVertex2()) && (getVertex3() != other.getVertex3()) && (getVertex1() != other.getVertex2()) && (getVertex1() != other.getVertex3()))) intersectCount++;
            if ( getVertex3().doSegmentsIntersect(getVertex1(), other.getVertex1(), other.getVertex2()) && ((getVertex3() != other.getVertex1()) && (getVertex3() != other.getVertex2()) && (getVertex1() != other.getVertex1()) && (getVertex1() != other.getVertex2()))) intersectCount++;
            if ( getVertex3().doSegmentsIntersect(getVertex1(), other.getVertex1(), other.getVertex3()) && ((getVertex3() != other.getVertex1()) && (getVertex3() != other.getVertex3()) && (getVertex1() != other.getVertex1()) && (getVertex1() != other.getVertex3()))) intersectCount++;

            if ((other.isPointInsideTriangle(getVertex1()) || other.isPointInsideTriangle(getVertex2()) ||  other.isPointInsideTriangle(getVertex3()))) intersectCount++;
            if ((isPointInsideTriangle(other.getVertex1()) || isPointInsideTriangle(other.getVertex2()) ||  isPointInsideTriangle(other.getVertex3()))) intersectCount++;
            // Restituisci true se ci sono almeno due intersezioni, altrimenti false
            return intersectCount >= 1;}
        else{return 0;}
}


bool Triangulation::isBoundaryTriangle(const Triangle& triangle) {
    unsigned int numAdjacentTriangles = adjacencyList[triangle.getId()].size();

    // If the number of adjacent triangles is less than three, the triangle is a boundary triangle
    if (numAdjacentTriangles < 3) {
        return true; // The triangle is a boundary triangle
    }

    return false; // The triangle is not a boundary triangle
}

void Triangulation::connectPointOnEdgeInside(Triangle& t, const Point& Q, const vector<Point>& edge) {
    vector<unsigned int> oldadj = adjacencyList[t.getId()];
    Triangle t1;
    Triangle t2;
    for (unsigned int i : oldadj){

        Triangle adjT = DelaunayTriangles[i];
        Point q1;
        Point q2;

        if (((adjT.getVertex1() == edge[0]) || (adjT.getVertex2() == edge[0]) || (adjT.getVertex3() == edge[0])) && ((adjT.getVertex1() == edge[1]) || (adjT.getVertex2() == edge[1]) || (adjT.getVertex3() == edge[1]))) {
           if(t.getVertex1() != edge[0] && t.getVertex1() != edge[1]) q1 = t.getVertex1();
           else if (t.getVertex2() != edge[0] && t.getVertex2() != edge[1]) q1 = t.getVertex2();
           else q1 = t.getVertex3();

           if(adjT.getVertex1() != edge[0] && adjT.getVertex1() != edge[1]) q2 = adjT.getVertex1();
           else if (adjT.getVertex2() != edge[0] && adjT.getVertex2() != edge[1]) q2 = adjT.getVertex2();
           else q2 = adjT.getVertex3();
/*
           t.setVertex1(edge[0]); t.setVertex2(Q); t.setVertex3(q1);
           t1.setVertex1(edge[1]); t1.setVertex2(q1); t1.setVertex3(Q); t1.setId(getMaxTriangleId() + 1);
           t2.setVertex1(edge[0]); t2.setVertex2(Q); t2.setVertex3(q2); t2.setId(getMaxTriangleId() + 2);
           adjT.setVertex1(edge[1]); adjT.setVertex2( Q); adjT.setVertex3(q2);*/
           t.p1=edge[0]; t.p2=Q; t.p3=q1;
           t1.p1=edge[1]; t1.p2=q1; t1.p3=Q; t1.setId(getMaxTriangleId() + 1);
           t2.p1=edge[0]; t2.p2=Q; t2.p3=q2; t2.setId(getMaxTriangleId() + 2);
           adjT.p1=edge[1]; adjT.p2= Q; adjT.setVertex3(q2);
           t.SortVertices();
           t1.SortVertices();
           t2.SortVertices();
           adjT.SortVertices();
           // Replace the previous triangle with t1 in DelaunayTriangles
           DelaunayTriangles[t.getId()] = t;
           DelaunayTriangles[adjT.getId()] = adjT;

           addTriangle(t1);
           addTriangle(t2);
           incrementMaxTriangleId(2);
           // Update adjacency
           addAdjacentTriangle(t.getId(), t2.getId());
           addAdjacentTriangle(t.getId(), t1.getId());
           addAdjacentTriangle(adjT.getId(), t1.getId());
           addAdjacentTriangle(adjT.getId(), t2.getId());
           removeAdjacentTriangle(adjT.getId(), t.getId());

           // Update adjacency for t2
           for (unsigned int adjId : oldadj) {

               if (DelaunayTriangles[adjId].isAdjacent(t2)){
                   if (!(isIdInAdjacency(adjId, adjacencyList[t2.getId()])))
                      addAdjacentTriangle(adjId, t2.getId());
               }
               else removeAdjacentTriangle(adjId, t2.getId());
           }

           // Update adjacency for t1
           for (unsigned int adjId : oldadj) {

                   if (DelaunayTriangles[adjId].isAdjacent(t1)) {
                       if (!(isIdInAdjacency(adjId, adjacencyList[t1.getId()])))

                           addAdjacentTriangle(adjId, t1.getId());

               }
                   else removeAdjacentTriangle(adjId, t1.getId());
           }

           // Update adjacency for t1
           for (unsigned int adjId : oldadj) {

                   if (DelaunayTriangles[adjId].isAdjacent(t) && adjId != t.getId() ) {
                       if (!(isIdInAdjacency(adjId, adjacencyList[t.getId()])))

                           addAdjacentTriangle(adjId, t.getId());

               }
                   else removeAdjacentTriangle(adjId, t.getId());
           }

           // Update adjacency for t1
           for (unsigned int adjId : oldadj) {

                   if (DelaunayTriangles[adjId].isAdjacent(adjT) && adjId != adjT.getId()) {
                       if (!(isIdInAdjacency(adjId, adjacencyList[adjT.getId()])))

                           addAdjacentTriangle(adjId, adjT.getId());

               }
                   else removeAdjacentTriangle(adjId, adjT.getId());
           }




        }
    }



}


void Triangulation::connectPointOnEdgeToTriangulation(Triangle& t, const Point& Q, const vector<Point>& edge) {


    if (!isBoundaryEdge(t, edge)) {
        connectPointOnEdgeInside(t, Q,edge); // Create two subtriangles from one triangle
    }

    else{
        Point q1;
        vector<unsigned int> oldadj = adjacencyList[t.getId()];
        if(t.getVertex1() != edge[0] && t.getVertex1() != edge[1]) q1 = t.getVertex1();
        else if (t.getVertex2() != edge[0] && t.getVertex2() != edge[1]) q1 = t.getVertex2();
        else q1 = t.getVertex3();
        Triangle t1;
        Point p1;
        t.setVertex1(edge[0]); t.setVertex2(Q); t.setVertex3(q1);


        t1.setVertex1(edge[1]); t1.setVertex2(Q); t1.setVertex3(q1); t1.setId(getMaxTriangleId() + 1);
        t.SortVertices();
        t1.SortVertices();
        DelaunayTriangles[t.getId()] = t;
        addTriangle(t1);
        incrementMaxTriangleId();
        addAdjacentTriangle(t.getId(),t1.getId());

        for (unsigned int adjId : oldadj) {

                if (DelaunayTriangles[adjId].isAdjacent(t)) {
                    if (!(isIdInAdjacency(adjId, adjacencyList[t.getId()])))
                        addAdjacentTriangle(adjId, t.getId());

            }
                else removeAdjacentTriangle(adjId, t.getId());
        }

        for (unsigned int adjId : oldadj) {

                if (DelaunayTriangles[adjId].isAdjacent(t1)) {
                    if (!(isIdInAdjacency(adjId, adjacencyList[t1.getId()])))

                        addAdjacentTriangle(adjId, t1.getId());

            }
                else removeAdjacentTriangle(adjId, t1.getId());
        }


    }

}





void Triangulation::addPointToTriangulation(const Point& Q) {

    if (PointInsideTriangulation(Q) >= 0) {
        int triangleId = PointInsideTriangulation(Q);
        vector<Point> edge;
        Triangle triangle(DelaunayTriangles[triangleId]);
        edge = triangle.PointOnEdge(Q); //lato in cui giace il punto Q
        if (!(edge[0].x==0.0 && edge[0].y==0.0 && edge[1].x==0.0 && edge[1].y==0.0)){  //se è interno è il vettore (0,0)
            connectPointOnEdgeToTriangulation(triangle,Q,edge);
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

        for (size_t i = 0; i < DelaunayTriangles.size(); ++i) {
            Triangle& triangle = DelaunayTriangles[i];

            // Crea una copia della lista di adiacenza del triangolo
            unsigned int id = triangle.getId();
            vector<unsigned int> adjacency(adjacencyList[id].begin(), adjacencyList[id].end());
            unsigned int n = adjacency.size();

            bool triangleFlipped = false;

            for (size_t j = 0; j < n; ++j) {
                unsigned int adjTriangleId = adjacency[j];
                Triangle& adjacentTriangle = DelaunayTriangles[adjTriangleId];

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
    vector<unsigned int> adjacency(adjacencyList[triangle.getId()].begin(), adjacencyList[triangle.getId()].end());
    vector<unsigned int> adjacency2(adjacencyList[adjacentTriangle.getId()].begin(), adjacencyList[adjacentTriangle.getId()].end());
    oldAdjacency.reserve(adjacency.size() + adjacencyList[adjacentTriangle.getId()].size());
    std::merge(adjacency.begin(), adjacency.end(), adjacencyList[adjacentTriangle.getId()].begin(), adjacencyList[adjacentTriangle.getId()].end(), std::back_inserter(oldAdjacency));

                // Rimuovi i duplicati
                auto last = std::unique(oldAdjacency.begin(), oldAdjacency.end());
                oldAdjacency.erase(last, oldAdjacency.end());

                // Aggiorna le liste di adiacenza dei triangoli dopo il flip

                // AGGIORNO LE MODIFICHE PER TRIANGLE
                for (size_t i = 0; i < oldAdjacency.size(); i++) {
                    unsigned int adjId = oldAdjacency[i];
                    Triangle& adjTriangle = DelaunayTriangles[adjId];
                    if (adjTriangle.isAdjacent(triangle) && (adjId != triangle.getId()) ){
                        if( !(isIdInAdjacency(adjId, adjacency)))
                            addAdjacentTriangle(adjId,triangle.getId());
                        }
                    else{
                        if( isIdInAdjacency(adjId, adjacency) && (adjId != triangle.getId()))
                            removeAdjacentTriangle(adjId,triangle.getId());
                    }
                }

                // Aggiornamento per adjacentTriangle
                for (size_t i = 0; i < oldAdjacency.size(); i++) {
                    unsigned int adjId = oldAdjacency[i];
                    Triangle& adjTriangle = DelaunayTriangles[adjId];
                    if (adjTriangle.isAdjacent(adjacentTriangle) && (adjId != adjacentTriangle.getId())) {
                        if (!isIdInAdjacency(adjId, adjacency2))
                            addAdjacentTriangle(adjId, adjacentTriangle.getId());
                    } else {
                        if (isIdInAdjacency(adjId, adjacency2) && (adjId != adjacentTriangle.getId()))
                            removeAdjacentTriangle(adjId, adjacentTriangle.getId());
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


void Triangulation::flipAndUpdateAdjacency(Triangle& triangle, Triangle& adjacentTriangle) {


        if (triangle.areTrianglesDelaunay(adjacentTriangle) == 0) {

            // Esegui il flip dei triangoli
            triangle.flip(adjacentTriangle);

            // Aggiorna i triangoli nella lista DelaunayTriangles
            DelaunayTriangles[triangle.getId()] = triangle;
            DelaunayTriangles[adjacentTriangle.getId()] = adjacentTriangle;

            // Aggiorna le liste di adiacenza dei triangoli
            updateAdjacency(triangle, adjacentTriangle);
        }

}






bool Triangulation::isBoundaryEdge(const Triangle& triangle, const std::vector<Point>& edge) {
    const std::vector<unsigned int>& adjList = adjacencyList[triangle.getId()];

    // Verifica che entrambi gli estremi del lato siano condivisi con un triangolo nella lista di adiacenza
    bool sharedEndpoints = false;
    for (unsigned int adjId : adjList) {
        const Triangle& adjTriangle = DelaunayTriangles[adjId];

        // Verifica se uno dei vertici del lato è condiviso con il triangolo adiacente
        if ((adjTriangle.isVertexShared(edge[0]) && adjTriangle.isVertexShared(edge[1]))) {
            sharedEndpoints = true;
            break;
        }
    }

    // Se entrambi gli estremi sono condivisi, il lato non è un lato di bordo
    return !sharedEndpoints;
}



Triangulation DelaunayTriangulation(const std::vector<Point>& points) {
    Triangle triangle_default;
    Triangle triangle_max;



    triangle_max = triangle_default.findMaximumTriangle(points);


    std::vector<Triangle> DelaunayTriangles;


    Triangulation triangulation(DelaunayTriangles);
    triangulation.addTriangle(triangle_max);

    std::vector<Point> pointsCopy;
    pointsCopy = points;

    pointsCopy.erase(std::remove(pointsCopy.begin(), pointsCopy.end(), triangle_max.getVertex1()), pointsCopy.end());
    pointsCopy.erase(std::remove(pointsCopy.begin(), pointsCopy.end(), triangle_max.getVertex2()), pointsCopy.end());
    pointsCopy.erase(std::remove(pointsCopy.begin(), pointsCopy.end(), triangle_max.getVertex3()), pointsCopy.end());
    for (const auto& point : pointsCopy) {
        cout << point << endl;
    }

    for (size_t i = 0; i < pointsCopy.size(); i++) {
        triangulation.addPointToTriangulation(pointsCopy[i]);

        unsigned int n = triangulation.DelaunayTriangles.size();
        for (size_t j = 0; j < n ; j++) {
            //Triangle& triangle = triangulation.DelaunayTriangles[j];
            // const std::vector<unsigned int>& adjacency = triangulation.adjacencyList[triangle.getId()];
            triangulation.flipTrianglesIfNotDelaunay();
        }

        // Rimuovi i triangoli duplicati
        auto last = std::unique(triangulation.DelaunayTriangles.begin(), triangulation.DelaunayTriangles.end());
        triangulation.DelaunayTriangles.erase(last, triangulation.DelaunayTriangles.end());
    }

    return triangulation;
}}
