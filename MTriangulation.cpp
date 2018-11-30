#include <QTime>
#include <CGAL/nearest_neighbor_delaunay_2.h>
#include "MTriangulation.hpp"

std::size_t Hash_point::operator()(const Point_2& p) const noexcept{
    return std::hash<double>()(p.x()) ^ (std::hash<double>()(p.y()));
}

double sqdist(const Point_2& p1, const Point_2& p2){
    return (p1.x()-p2.x())*(p1.x()-p2.x())+(p1.y()-p2.y())*(p1.y()-p2.y());
}

struct CompXVertexMoveHint{
    bool operator()(const VertexMoveHint& v1, const VertexMoveHint& v2){
        return v1.new_position.x() < v2.new_position.x();
    }
};

struct CompYVertexMoveHint{
    bool operator()(const VertexMoveHint& v1, const VertexMoveHint& v2){
        return v1.new_position.y() < v2.new_position.y();
    }
};

struct VertexMoveHintCompTrait{
    typedef VertexMoveHint Point_2;
    typedef CompXVertexMoveHint Less_x_2;
    typedef CompYVertexMoveHint Less_y_2;
    Less_x_2 less_x_2_object() const{
        return Less_x_2();
    }
    Less_y_2 less_y_2_object() const{
        return Less_y_2();
    }
};

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
    std::vector<VertexMoveHint> newPointsHint;

    for(auto it = all_vertices_begin(); it != all_vertices_end(); ++it){
        if(not is_infinite(it)){
            Point_2 p = it->point();
            float r = 2*(static_cast <float> (rand()) / static_cast <float> (RAND_MAX))-1; // random [-1, 1]   
            float theta = 2*3.1415*(static_cast <float> (rand()) / static_cast <float> (RAND_MAX)); // random [0, 2pi] 
            Vector_2 v(cos(theta), sin(theta));
            v *= r*rMax;

            Point_2 nPoint = p+v;

            switch(iStyle){
                case MOVE_CGAL:
                    move(it, nPoint);
                    continue;
                break;
                case NAIVE:
                    newPoints.push_back(nPoint);
                    break;
                case HINT:
                    newPointsHint.emplace_back(it, nPoint);
                break;
            }
        }
    }

    switch(iStyle){
        case MOVE_CGAL:
        break;
        case NAIVE:
            clear();
            insert_naive(newPoints);
        break;
        case HINT:
            std::unordered_map<Vertex_handle, Vertex_handle> old_nn(nearest_neight);
            clear();
            insert_hint(newPointsHint, old_nn);
        break;
    }

    return timer.elapsed();
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

void MTriangulation::insert_hint(std::vector<VertexMoveHint>& newPointsHint, std::unordered_map<Vertex_handle, Vertex_handle>& old_nn){
    using namespace std;
    unsigned long long number_improvement = 0;
    VertexMoveHintCompTrait trait;
    CGAL::spatial_sort(newPointsHint.begin(), newPointsHint.end(), trait);

    std::unordered_map<Vertex_handle, Vertex_handle> newVertexs;

    Vertex_handle vertex_hint = Vertex_handle();

    for(const VertexMoveHint& h : newPointsHint){
        Vertex_handle old_vertex = h.handle;
        Point_2 nPos = h.new_position;

        Vertex_handle prev_nn = old_nn[old_vertex];

        auto hint_it = newVertexs.find(prev_nn);

        if(hint_it != newVertexs.end()){ //if there is an hint
            Vertex_handle hint_nn = hint_it->second;
            if(vertex_hint == Vertex_handle() or sqdist(hint_nn->point(), nPos) < sqdist(vertex_hint->point(), nPos)){ //if we better chose the hint, we do so
                vertex_hint = hint_nn;
                number_improvement++;
            }
        }

        Vertex_handle v_handle = insert(nPos, vertex_hint == Vertex_handle()? Face_handle() : vertex_hint->face() );
        vertex_hint = v_handle;
        newVertexs[old_vertex] = v_handle;
    }
    std::cout << "Fraction of improvement : " << (double)number_improvement/newPointsHint.size() << std::endl;
}