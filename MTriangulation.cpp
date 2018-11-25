#include <QTime>
#include <CGAL/nearest_neighbor_delaunay_2.h>
#include "MTriangulation.hpp"

std::size_t Hash_point::operator()(const Point_2& p) const noexcept{
    return std::hash<double>()(p.x()) ^ (std::hash<double>()(p.y()));
}


MTriangulation::MTriangulation(InsertStyle is) : Delaunay(), iStyle(is){}

MTriangulation::Vertex_handle MTriangulation::insert(const Point_2& p, const Face_handle& f){
    Vertex_handle vh = Delaunay::insert(p, f);
    nearest_neight[vh] = CGAL::nearest_neighbor(*this, vh);
    auto circulator = vh->incident_vertices();
    if(circulator != nullptr){
        auto begin = circulator;
        do{
            nearest_neight[circulator] = CGAL::nearest_neighbor(*this, circulator);
            circulator++;
        }while(circulator != begin);
    }
    return vh;
}

void MTriangulation::clear(){
    nearest_neight.clear();
    Delaunay::clear();
}

int MTriangulation::moveBrownian(float rMax){
    QTime timer;
    timer.start();

    std::vector<Point_2> newPoints;
    Hint_insertion is_nn_of;
    for(auto it = all_vertices_begin(); it != all_vertices_end(); ++it){
        if(not is_infinite(it)){
            Point_2 p = it->point();
             if(iStyle == HINT){
                 Vertex_handle nn_handle = nearest_neight[it];
                 assert(p != nn_handle->point());
                 is_nn_of.insert({nn_handle->point(), p});
             }
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
        insert_move(newPoints, is_nn_of);
    }

    return timer.elapsed();
}

void MTriangulation::insert_move(std::vector<Point_2>& newPoints, const Hint_insertion& is_nn_of){
    switch(iStyle){
    case NAIVE:
        insert_naive(newPoints);
        break;
    case HINT:
        insert_hint(newPoints, is_nn_of);
        break;
    case MOVE_CGAL:
    break;
    }
}

void MTriangulation::setInsertStyle(InsertStyle is){
    iStyle = is;
}

void MTriangulation::insert_naive(std::vector<Point_2>& points){
    CGAL::spatial_sort(points.begin(), points.end());
    Face_handle face_hint = Face_handle();
    for(Point_2 p : points){
        Vertex_handle v_handle = insert(p, face_hint);
        face_hint = v_handle->face();
    }
}

void MTriangulation::insert_hint(std::vector<Point_2>& newPoints, const Hint_insertion& is_nn_of){
    CGAL::spatial_sort(newPoints.begin(), newPoints.end());
    std::unordered_map<Point_2, Vertex_handle, Hash_point> hints;
    Face_handle face_hint = Face_handle();
    for(Point_2 p : newPoints){
         if(hints.find(p) != hints.end()){
             face_hint = hints[p]->face();
         }
        Vertex_handle v_handle = insert(p, face_hint);
        face_hint = v_handle->face();
         auto range = is_nn_of.equal_range(p);
         for(auto it = range.first; it != range.second; ++it){
             Point_2 will_be_hinted = it->second;
             hints[will_be_hinted] = v_handle;
         }
    }
}