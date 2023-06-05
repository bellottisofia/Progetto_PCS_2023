#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>


#include "delaunay.hpp"

using namespace std;
using namespace ProjectLibrary;

std::ostream& operator<<(std::ostream& os, const Triangle& triangle)
{
    os << "Triangle " << triangle.id << ":" << std::endl;
    os << "Point 1 - x: " << triangle.p1.x << ", y: " << triangle.p1.y << ", id: " << triangle.p1.id << std::endl;
    os << "Point 2 - x: " << triangle.p2.x << ", y: " << triangle.p2.y << ", id: " << triangle.p2.id << std::endl;
    os << "Point 3 - x: " << triangle.p3.x << ", y: " << triangle.p3.y << ", id: " << triangle.p3.id << std::endl;
    // Add any additional information you want to output about the triangle
    return os;
}



bool ImportData(const string& inputFilePath, vector<Point>& points);

int main()
{
    string inputFileName = "C:/Users/giorg/OneDrive/Desktop/Progetto_PCS_2023/Projects/Delaunay/Dataset/Test1.csv"; // Inserisci il percorso del tuo file
    vector<Point> points;

    if (!ImportData(inputFileName, points))
    {
      cerr<< "Something goes wrong with import"<< endl;
      return -1;
    }

    Triangle triangle;
    Triangle triangle_max;
    triangle_max = triangle.findMaximumTriangle(points);

    for (const auto& point : points) {
            std::cout << point << std::endl;
        }
    cout << triangle_max;

    return 0;
}


bool ImportData(const string& inputFilePath, vector<Point>& points)
{
    ifstream file;
    file.open(inputFilePath);

    if (!file.is_open())
    {
        cerr << "File open failed" << endl;
        return false;
    }

    string line;
    getline(file, line);
    while (getline(file, line))
    {
        istringstream iss(line);
        unsigned int currentId;
        double currentX, currentY;
        if (iss >> currentId >> currentX >> currentY)
        {
            points.push_back(Point( currentX, currentY, currentId));
        }
        else
        {
            cerr << "Errore durante la lettura dei dati dalla riga: " << line << endl;
            file.close();
            return false;
        }
    }

    file.close();
    return true;
}
