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
    os << "Triangle " << triangle.getId() << ":" << std::endl;
    os << "Point 1 - x: " << triangle.getVertex1().x << ", y: " << triangle.getVertex1().y << ", id: " << triangle.getVertex1().id << std::endl;
    os << "Point 2 - x: " << triangle.getVertex2().x << ", y: " << triangle.getVertex2().y << ", id: " << triangle.getVertex2().id << std::endl;
    os << "Point 3 - x: " << triangle.getVertex3().x << ", y: " << triangle.getVertex3().y << ", id: " << triangle.getVertex3().id << std::endl;
    // Add any additional information you want to output about the triangle
    return os;
}



bool ImportData(const string& inputFilePath, vector<Point>& points);
bool ExportResult(const std::string& outputFilePath, ProjectLibrary::Triangulation& triangulation);

int main()
{
    string inputFileName = "../Delaunay/Dataset/Test1.csv"; // Inserisci il percorso del tuo file
    vector<Point> points;


    if (!ImportData(inputFileName, points))
    {
      cerr<< "Something goes wrong with import"<< endl;
      return -1;
    }

    Triangulation Delaunaytriangulation;
    Delaunaytriangulation=DelaunayTriangulation( points);
    string outputFileName = "../Delaunay.txt";
        if (!ExportResult(outputFileName, Delaunaytriangulation))
        {
            cerr << "Something went wrong with export" << endl;
            return -1;
        }
        else
        {
            cout << "Export successful" << endl;
        }

        return 0;
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
bool ExportResult(const std::string& outputFilePath, ProjectLibrary::Triangulation& triangulation)

{
    ofstream file(outputFilePath);
    if (!file)
    {
        cerr << "File open failed" << endl;
        return false;
    }

    for (unsigned int i = 0; i < triangulation.DelaunayTriangles.size(); i++)
    {
        // Get the current Triangle object
        Triangle& triangle = triangulation.DelaunayTriangles[i];

        // Write the Triangle object's data to the file
        file << triangle<< endl;

    }

    file.close();
    return true;
}
