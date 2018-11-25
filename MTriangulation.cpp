#include <QTime>
#include <CGAL/nearest_neighbor_delaunay_2.h>
#include "MTriangulation.hpp"

std::size_t Hash_point::operator()(const Point_2& p) const noexcept{
    return std::hash<double>()(p.x()) ^ (std::hash<double>()(p.y()));
}

double sqdist(const Point_2& p1, const Point_2& p2){
    return (p1.x()-p2.x())*(p1.x()-p2.x())+(p1.y()-p2.y())*(p1.y()-p2.y());
}


MTriangulation::MTriangulation(InsertStyle is) : Delaunay(), iStyle(is){}

MTriangulation::Vertex_handle MTriangulation::insert(const Point_2& p, const Face_handle& f){
    Vertex_handle vh = Delaunay::insert(p, f);

    nearest_neight[vh] = Vertex_handle();
    nearest_neight_sqdistance[vh] = -1;

    Vertex_circulator circ = vh->incident_vertices();
    if(circ == nullptr){
        return vh;
    }
    Vertex_circulator begin = circ;

    double min_dist = -1;

    do{
        double dist_neight = sqdist(vh->point(), circ->point());
        if(min_dist == -1 or dist_neight < min_dist){
            min_dist = dist_neight;
            nearest_neight[vh] = circ;
        }
        if(nearest_neight_sqdistance[circ] == -1 or dist_neight < nearest_neight_sqdistance[circ]){
            nearest_neight_sqdistance[circ] = dist_neight;
            nearest_neight[circ] = vh;
        }
        circ++;
    }while(circ != begin);
    nearest_neight_sqdistance[vh] = min_dist;
    return vh;
}

void MTriangulation::clear(){
    nearest_neight.clear();
    nearest_neight_sqdistance.clear();
    Delaunay::clear();
}

int MTriangulation::moveBrownian(float rMax){
    QTime timer;
    timer.start();

    std::vector<Point_2> newPoints;
    Hint_insertion is_nn_of;
    std::unordered_map<Point_2, Vertex_handle, Hash_point> point_to_vertex;
    std::unordered_map<Vertex_handle, Point_2> vertex_to_point;

    for(auto it = all_vertices_begin(); it != all_vertices_end(); ++it){
        if(not is_infinite(it)){
            Point_2 p = it->point();
            float r = 2*(static_cast <float> (rand()) / static_cast <float> (RAND_MAX))-1; // random [-1, 1]   
            float theta = 2*3.1415*(static_cast <float> (rand()) / static_cast <float> (RAND_MAX)); // random [0, 2pi] 
            Vector_2 v(cos(theta), sin(theta));
            v *= r*rMax;

            Point_2 nPoint = p+v;

            if(iStyle == MOVE_CGAL){
                move(it, p+v);
                continue;
            }else{
                newPoints.push_back(nPoint);
            }

             if(iStyle == HINT){
                 Vertex_handle nn_handle = nearest_neight[it];
                 assert(p != nn_handle->point());
                 is_nn_of.insert({nn_handle, it});
                 point_to_vertex[nPoint] = it;
                 vertex_to_point[it] = nPoint;
             }

        }
    }
    if(iStyle != MOVE_CGAL){
        clear();
        insert_move(newPoints, point_to_vertex, vertex_to_point, is_nn_of);
    }

    return timer.elapsed();
}

void MTriangulation::insert_move(std::vector<Point_2>& newPoints, 
    std::unordered_map<Point_2, Vertex_handle, Hash_point>& point_to_vertex,
    std::unordered_map<Vertex_handle, Point_2>& vertex_to_point,
    const Hint_insertion& is_nn_of){
    switch(iStyle){
    case NAIVE:
        insert_naive(newPoints);
        break;
    case HINT:
        insert_hint(newPoints, point_to_vertex, vertex_to_point, is_nn_of);
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

void MTriangulation::insert_hint(std::vector<Point_2>& newPoints, 
    std::unordered_map<Point_2, Vertex_handle, Hash_point>& point_to_vertex,
    std::unordered_map<Vertex_handle, Point_2>& vertex_to_point,
    const Hint_insertion& is_nn_of){

    unsigned long long number_improvement = 0;

    CGAL::spatial_sort(newPoints.begin(), newPoints.end());

    std::unordered_map<Vertex_handle, Vertex_handle> hints;

    Vertex_handle vertex_hint = Vertex_handle();
    for(const Point_2& p : newPoints){
        Vertex_handle old_vertex = point_to_vertex[p];

        if(hints.find(old_vertex) != hints.end()){ //if we better chose the hint, we do so
            if(vertex_hint == Vertex_handle() or sqdist(hints[old_vertex]->point(), p) < sqdist(p, vertex_hint->point())){
                vertex_hint = hints[old_vertex];
                number_improvement++;
            }
        }

        Vertex_handle v_handle = insert(p, vertex_hint == nullptr? Face_handle() : vertex_hint->face() );
        vertex_hint = v_handle;

        auto range = is_nn_of.equal_range(old_vertex);
        for(auto it = range.first; it != range.second; ++it){
            Vertex_handle will_be_hinted = it->second;
            hints[will_be_hinted] = v_handle;
        }
    }
    std::cout << "Fraction of improvement : " << (double)number_improvement/newPoints.size() << std::endl;
}