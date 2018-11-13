#include <QTime>
#include "MTriangulation.hpp"


MTriangulation::MTriangulation(InsertStyle is) : Delaunay(), iStyle(is){}

int MTriangulation::moveBrownian(float rMax){
    QTime timer;
    timer.start();

    std::vector<Point_2> newPoints;
    std::vector<Face_handle> hints;
    for(auto it = points_begin(); it != points_end(); ++it){
        Point_2 p = *it;
        float r = 2*(static_cast <float> (rand()) / static_cast <float> (RAND_MAX))-1; // random [-1, 1]   
        float theta = 2*3.1415*(static_cast <float> (rand()) / static_cast <float> (RAND_MAX)); // random [0, 2pi] 
        Vector_2 v(cos(theta), sin(theta));
        v *= r*rMax;
        newPoints.push_back(p+v);
    }
    clear();

    insert_move(newPoints, hints);

    return timer.elapsed();
}

void MTriangulation::insert_move(const std::vector<Point_2>& newPoints, const std::vector<Face_handle>& hints){
    switch(iStyle){
    case NAIVE:
        insert_naive(newPoints);
        break;
    case HINT:
        insert_hint(newPoints, hints);
        break;
    }
}

void MTriangulation::setInsertStyle(InsertStyle is){
    iStyle = is;
}

void MTriangulation::insert_naive(const std::vector<Point_2>& newPoints){
    insert(newPoints.begin(), newPoints.end());
}

void MTriangulation::insert_hint(const std::vector<Point_2>& newPoints, const std::vector<Face_handle>& hints){
    
}