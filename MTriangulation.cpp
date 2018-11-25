#include <QTime>
#include <CGAL/nearest_neighbor_delaunay_2.h>
#include "MTriangulation.hpp"


MTriangulation::MTriangulation(InsertStyle is) : Delaunay(), iStyle(is){}

int MTriangulation::moveBrownian(float rMax){
    QTime timer;
    timer.start();

    std::vector<Point_2> newPoints;
    std::vector<Vertex_handle> hints;
    for(auto it = all_vertices_begin(); it != all_vertices_end(); ++it){
        if(not is_infinite(it)){
            if(iStyle == HINT){
                hints.push_back(CGAL::nearest_neighbor(*this, it));
            }
            Point_2 p = it->point();
            float r = 2*(static_cast <float> (rand()) / static_cast <float> (RAND_MAX))-1; // random [-1, 1]   
            float theta = 2*3.1415*(static_cast <float> (rand()) / static_cast <float> (RAND_MAX)); // random [0, 2pi] 
            Vector_2 v(cos(theta), sin(theta));
            v *= r*rMax;
            if(iStyle == MOVE_CGAL){
                move(it, p+v);
            }else{
                newPoints.push_back(p+v);
            }
        }
    }
    if(iStyle != MOVE_CGAL){
        clear();
        insert_move(newPoints, hints);
    }

    return timer.elapsed();
}

void MTriangulation::insert_move(const std::vector<Point_2>& newPoints, const std::vector<Vertex_handle>& hints){
    switch(iStyle){
    case NAIVE:
        insert_naive(newPoints);
        break;
    case HINT:{

        break;
    }
    case MOVE_CGAL:
    break;
    }
}

void MTriangulation::setInsertStyle(InsertStyle is){
    iStyle = is;
}

void MTriangulation::insert_naive(const std::vector<Point_2>& newPoints){
    insert(newPoints.begin(), newPoints.end());
}

void MTriangulation::insert_hint(const std::vector<Point_2>& newPoints, const std::vector<Vertex_handle>& hints){
    
}