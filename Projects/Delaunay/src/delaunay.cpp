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
            SortLibrary::MergeSort(sortedPoints, 0, sortedPoints.size() -1 );

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



    void Triangulation::addTriangle(const Triangle& triangle) {

        DelunayTriangles.push_back(triangle);
        adjacencyList.clear();
        adjacencyList[triangle.id] = {}; // Inizializza la lista di adiacenza per il triangolo
    }

    void Triangulation::addAdjacentTriangle(int triangleId, int adjacentTriangleId) {
        adjacencyList[triangleId].push_back(adjacentTriangleId);
        adjacencyList[adjacentTriangleId].push_back(triangleId);
    }
    const std::vector<unsigned int>& Triangulation::getAdjacentTriangles(int triangleId) {
        return adjacencyList[triangleId];
    }




    // Verifica se il punto q si trova sul segmento p-r



    bool Point::isPointOnSegment(const Point& p,  const Point& r) const{
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
            return !(point1InsideCircle || point2InsideCircle || point3InsideCircle ||
                     point4InsideCircle || point5InsideCircle || point6InsideCircle);
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

            triangle1.p1= nonCommonVertex[0];
            triangle1.p2= nonCommonVertex[1];
            triangle1.p3= commonVertices[1];

    }



    // Aggiorna la lista di adiacenza basata sulle relazioni tra i triangoli
    // Implementazione specifica in base alle regole della triangolazione

/*
void Triangulation::addAdjacentTriangle(int triangleId, int adjacentTriangleId) {
    triangles[triangleId].adjacentTriangles.push_back(adjacentTriangleId);
    triangles[adjacentTriangleId].adjacentTriangles.push_back(triangleId);
}
*/
std::vector<Point> Triangle::PointOnEdge(const Point& Q){
std::vector<Point> edge;

if (Q.isPointOnSegment(p1,p2))
{   edge.push_back(p1);
    edge.push_back(p2);


}
else if (Q.isPointOnSegment(p1,p3))
{
    edge.push_back(p1);
    edge.push_back(p3);

}
else if(Q.isPointOnSegment(p3,p2))
{    edge.push_back(p3);
    edge.push_back(p2);


}
else
{Point point;
    edge.push_back(point);
    edge.push_back(point);

}
return edge;}

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

Triangle Triangulation::findTriangleById(const unsigned int triangleId) {
    for (const Triangle& triangle : DelunayTriangles) {
        if (triangle.id == triangleId) {
            return triangle;
        }
    }
    // Restituisce un triangolo vuoto se l'ID non viene trovato
    return Triangle();
}


void Triangulation::createSubtriangulation(const Point& Q, int triangleId) {
    // Ottieni il triangolo T dalla lista dei triangoli
    Triangle& T = DelunayTriangles[triangleId];

    // Crea il nuovo triangolo collegando il punto Q ai vertici di T
    Triangle newTriangle(Q, T.p1, T.p2,triangleId );
    newTriangle.SortVertices();

    // Aggiungi il nuovo triangolo alla lista dei triangoli
    DelunayTriangles[triangleId] = newTriangle;
  /*
    //cancelliamo T
    DelunayTriangles.erase(T);
    //aggiungiamo newTriangle al posto di T
    addTriangle(newTriangle);
    */

    // Aggiorna la lista di adiacenza per il nuovo triangolo
    //ddAdjacentTriangle(triangleId, triangleId);

    // Aggiungi gli altri due triangoli della sottotriangolazione alla lista dei triangoli
    Triangle t1(Q, T.p2, T.p3,triangleId +1);
    Triangle t2(Q, T.p3, T.p1,triangleId +2);
    t1.SortVertices();
    t2.SortVertices();
    addTriangle(t1);
    addTriangle(t2);


    // Aggiorna la lista di adiacenza per i nuovi triangoli
    unsigned int n = adjacencyList[T.id].size();
    for(size_t i = 0; i < n; i++)
    {Triangle triangle;
        triangle=findTriangleById(adjacencyList[T.id][i]);
        if(triangle.isAdjacent(newTriangle)){
            addAdjacentTriangle(triangle.id,newTriangle.id);}
        else if (triangle.isAdjacent(t1)){
                addAdjacentTriangle(triangle.id,t1.id);}
        else if (triangle.isAdjacent(t2)){
                addAdjacentTriangle(triangle.id,t2.id);

        }
        }

    addAdjacentTriangle(triangleId, t1.id);
    addAdjacentTriangle(triangleId, t2.id);
    addAdjacentTriangle(t1.id, t2.id);

    }

 void Triangulation::connectPointToVertices(const Point& Q){
// Il punto Q è esterno alla triangolazione, unisce Q con tutti i vertici escludendo le intersezioni
unsigned int n = DelunayTriangles.size();
for (size_t i = 0; i < n ; i++) {
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
    addTriangle(t1);
    addAdjacentTriangle(id, triangle.id);
    }
    if(!t2.isPointInsideTriangle(triangle.p2))
    {id+=1;
    addTriangle(t2);
    //DelunayTriangles.push_back(t2);
    addAdjacentTriangle(id, triangle.id);
    }
    if(!t3.isPointInsideTriangle(triangle.p1))
    {id+=1;
    addTriangle(t3);
    addAdjacentTriangle(id, triangle.id);
    }

    }
    else{
        if (triangle.p1.isPointInVector(intersection)){
            Triangle t1(Q,triangle.p2,triangle.p3,triangle.id +1);
            t1.SortVertices();
           addTriangle(t1);
           addAdjacentTriangle(t1.id, triangle.id);
        }
        else if (triangle.p2.isPointInVector(intersection)){
            Triangle t2(Q,triangle.p1,triangle.p3,triangle.id +1);
            t2.SortVertices();
            addTriangle(t2);
            addAdjacentTriangle(t2.id, triangle.id);
         }
        else {
              Triangle t3(Q,triangle.p1,triangle.p2,triangle.id +1);
              t3.SortVertices();
              addTriangle(t3);
              addAdjacentTriangle(t3.id, triangle.id);
                }



}}}
void Triangulation::connectPointOnEdge( const Triangle& t,const Point& Q, const vector<Point> edge){
    if (!t.p1.isPointInVector(edge)){
        Triangle t1(Q,t.p1,t.p2,t.id );
        t1.SortVertices();
       DelunayTriangles[t.id] =t1;
       Triangle t2(Q,t.p1,t.p3,t.id +1);
       t2.SortVertices();
       addTriangle(t2);

    }
    else if (!t.p2.isPointInVector(edge)){
        Triangle t1(Q,t.p2,t.p1,t.id );
        t1.SortVertices();
        DelunayTriangles[t.id] =t1;//sostituisco
        Triangle t2(Q,t.p2,t.p3,t.id+1 );
        t2.SortVertices();
        addTriangle(t2);//aggiungo
     }
    else {
        Triangle t1(Q,t.p3,t.p1,t.id );
        t1.SortVertices();
          DelunayTriangles[t.id] =t1;
          Triangle t2(Q,t.p3,t.p2,t.id+1);
          t2.SortVertices();
          addTriangle(t2);

            }
    //aggiorno adiacenza
    addAdjacentTriangle(t.id,t.id+1);
    unsigned int n = adjacencyList[t.id].size();
    for(size_t i = 0; i < n; i++)
    {Triangle triangle;
        triangle=findTriangleById(adjacencyList[t.id][i]);
        if(triangle.isAdjacent(DelunayTriangles[t.id])){
            addAdjacentTriangle(triangle.id,DelunayTriangles[t.id].id);}
        if(triangle.isAdjacent(DelunayTriangles[t.id+1])){
            addAdjacentTriangle(triangle.id,DelunayTriangles[t.id+1].id);}

}
}

 void Triangulation::connectPointOnEdgeToTriangulation(const Triangle& triangle, const Point& Q, vector<Point>& edge){
     //connetto il punto al vertice opposto al lato
     connectPointOnEdge(triangle,Q, edge); //da un triangolo creo due sottotriangoli

     //se è su un lato è "interno" a due triangoli adiacenti su quel lato,
     //cerco l'altro triangolo
     unsigned int n = adjacencyList[triangle.id].size();
     for(size_t i = 0; i < n ; i++)
     {
         Triangle t(findTriangleById(adjacencyList[triangle.id][i]));
         if(t.isPointInsideTriangle(Q)){
             connectPointOnEdge(t,Q, edge); //da un triangolo creo due sottotriangoli

         }
     }

 }


void Triangulation::addPointToTriangulation(const Point& Q) {

    if (PointInsideTriangulation(Q)>0) {
        int triangleId = PointInsideTriangulation(Q);
        vector<Point> edge;
        Triangle triangle(findTriangleById(triangleId));
         edge=triangle.PointOnEdge(Q); //lato in cui giace il punto Q
        if (edge[0]!=Point() && edge[1]!=Point()){//se è interno è il vettore (0,0)
            connectPointOnEdgeToTriangulation(triangle,Q,edge);
        }
        else
        createSubtriangulation(Q, triangleId);
    } else {
        connectPointToVertices(Q);

    }
}

Triangulation DelunayTriangulation(const std::vector<Point>& points) {
    Triangle triangle_default;
    Triangle triangle_max;

    triangle_max = triangle_default.findMaximumTriangle(points);

    std::vector<Triangle> DelunayTriangles;
    DelunayTriangles.push_back(triangle_max);
    Triangulation triangulation(DelunayTriangles);

    std::vector<Point> pointsCopy;
    pointsCopy = points;

    pointsCopy.erase(std::remove(pointsCopy.begin(), pointsCopy.end(), triangle_max.p1), pointsCopy.end());
    pointsCopy.erase(std::remove(pointsCopy.begin(), pointsCopy.end(), triangle_max.p2), pointsCopy.end());
    pointsCopy.erase(std::remove(pointsCopy.begin(), pointsCopy.end(), triangle_max.p3), pointsCopy.end());

    for (size_t i = 0; i < pointsCopy.size(); i++) {
        triangulation.addPointToTriangulation(pointsCopy[i]);

        for (size_t j = 0; j < triangulation.DelunayTriangles.size(); j++) {
            Triangle& triangle = triangulation.DelunayTriangles[j];

            if (triangulation.adjacencyList.find(triangle.id) != triangulation.adjacencyList.end()) {
                const std::vector<unsigned int>& adjacency = triangulation.adjacencyList[triangle.id];
                for (size_t k = 0; k < adjacency.size(); k++) {
                    Triangle t;
                    t = triangulation.findTriangleById(adjacency[k]);

                    if (triangle.areTrianglesDelaunay(t) != 1) {
                        triangle.flip(t);
                    }
                }
            }
        }
    }

    return triangulation;
}}

