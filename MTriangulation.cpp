#define Dn(x) std::cout << #x << " : " << (x) << std::endl;
#define D(x) std::cout << (x) << std::endl;
#include <QTime>
#include <cmath>
#include <CGAL/nearest_neighbor_delaunay_2.h>
#include <CGAL/centroid.h>
#include "MTriangulation.hpp"

/**
 * @brief Hash functor in order to use unordered_map of Point_2
 * 
 */
struct Hash_point{
    std::size_t operator()(const Point_2& p) const noexcept{
        return std::hash<double>()(p.x()) ^ (std::hash<double>()(p.y()));
    }
};

double sqdist(const Point_2& p1, const Point_2& p2){
    return (p1.x()-p2.x())*(p1.x()-p2.x())+(p1.y()-p2.y())*(p1.y()-p2.y());
}

/**
 * @brief Compairison functors in order to be able to do a geometric sort of VertexMoveHint
 * 
 */
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

MTriangulation::MTriangulation(InsertStyle is, MovingStyle ms) : Delaunay(), iStyle(is), mStyle(ms), current_insert(){}

MTriangulation::Vertex_handle MTriangulation::insert(const Point_2& p, const Face_handle& f, Vector_2 ball){
    // use the proper structure to insert in
    Nn_map *nn = &(nearest_neight[current_insert]);
    Nn_dist_map *nn_dist = &(nearest_neight_sqdistance[current_insert]);
    Ball_map *vBall = &(vertex_ball[current_insert]);

    // Vector_2(0, 0) is a "sentinel value" in order to say "sample this ball"
    if(ball == Vector_2(0, 0)){
        float r = 2*(static_cast <float> (rand()) / static_cast <float> (RAND_MAX))-1; // random [-1, 1]   
        float theta = 2*3.1415*(static_cast <float> (rand()) / static_cast <float> (RAND_MAX)); // random [0, 2pi] 
        ball = r*Vector_2(cos(theta), sin(theta));
    }

    Vertex_handle vh = Delaunay::insert(p, f);

    // maintains the list of nn, Vertex_handle() and -1 are sentinel values
    vBall->insert({vh, ball});
    nn->insert({vh, Vertex_handle()});
    nn_dist->insert({vh, -1});

    update_nn(vh);

    return vh;
}

void MTriangulation::update_nn(Vertex_handle vh){
    Nn_map *nn = &(nearest_neight[current_insert]);
    Nn_dist_map *nn_dist = &(nearest_neight_sqdistance[current_insert]);

    Vertex_circulator circ = vh->incident_vertices();
    if(circ == nullptr){
        return;
    }

    Vertex_circulator begin = circ;

    double min_dist = -1;

    // for each neighboord, either vh is their new NN or their previous NN is still connected to them, since the NN of a point is still connected to it in a Delaunay
    do{
        // if the NN has not been inserted into nn (implementation artifact)
        if(nn->find(circ) == nn->end()){
            nn->insert({circ, Vertex_handle()});
            nn_dist->insert({circ, -1});
        }

        double dist_neight = sqdist(vh->point(), circ->point());
        
        //update of the nn of the current vertex
        if(min_dist == -1 or dist_neight < min_dist){
            min_dist = dist_neight;
            nn->at(vh) = circ;
        }

        //update of the nn of the current neighboord
        if(nn_dist->at(circ) == -1 or dist_neight < nn_dist->at(circ)){
            nn_dist->at(circ) = dist_neight;
            nn->at(circ) = vh;
        }
        circ++;
    }while(circ != begin);
    nn_dist->at(vh) = min_dist;
}

void MTriangulation::clear(){
    nearest_neight[current_insert].clear();
    nearest_neight_sqdistance[current_insert].clear();
    Delaunay::clear();
}

Point_2 MTriangulation::brownianStep(Point_2 start, float rMax){
    // just sample a vector and add it
    float r = 2*(static_cast <float> (rand()) / static_cast <float> (RAND_MAX))-1; // random [-1, 1]   
    float theta = 2*3.1415*(static_cast <float> (rand()) / static_cast <float> (RAND_MAX)); // random [0, 2pi] 
    Vector_2 v(cos(theta), sin(theta));
    v *= r*rMax;
    return start + v;
}

Point_2 MTriangulation::bouncingBallStep(Vertex_handle start, float speed, QRectF rect){
    Vector_2 vector = vertex_ball[current_insert][start];
    // get the new point by adding the vector
    Point_2 nPoint = start->point()+speed*vector;

    // Now check the position of the new point and reflect the vector if needed
    if(nPoint.x() < rect.x()){ //too far to the left
        nPoint = Point_2(2*rect.x()-nPoint.x(), nPoint.y());
        vector = Vector_2(-vector.x(), vector.y());
    }
    if(nPoint.x() > rect.right()){ //too far to the right
        nPoint = Point_2(2*rect.right()-nPoint.x(), nPoint.y());
        vector = Vector_2(-vector.x(), vector.y());
    }
    if(nPoint.y() > rect.bottom()){ // too far to the bottom
        nPoint = Point_2(nPoint.x(), 2*rect.bottom()-nPoint.y());
        vector = Vector_2(vector.x(), -vector.y());
    }
    if(nPoint.y() < rect.y()){ //top far up
        nPoint = Point_2(nPoint.x(), 2*rect.y()-nPoint.y());
        vector = Vector_2(vector.x(), -vector.y());
    }
    vertex_ball[current_insert][start] = vector;

    return nPoint;
}


Point_2 MTriangulation::lloydStep(Vertex_handle vh, QRectF rect,  std::vector<Point_2>& points_convex_hull){
    // routine used everywhere else in this function
    auto contained_in_window= [rect](Point_2 p){
        return rect.contains(p.x(), p.y());
    };

    if(not contained_in_window(vh->point())){ // we do not move point outside the window
        return vh->point();
    }

    // segments representing the sides of the view
    Segment_2 topWindow(Point_2(rect.left(), rect.top()), Point_2(rect.right(), rect.top()));
    Segment_2 bottomWindow(Point_2(rect.left(), rect.bottom()), Point_2(rect.right(), rect.bottom()));
    Segment_2 leftWindow(Point_2(rect.left(), rect.top()), Point_2(rect.left(), rect.bottom()));
    Segment_2 rightWindow(Point_2(rect.right(), rect.top()), Point_2(rect.right(), rect.bottom()));

    std::unordered_map<Point_2, bool, Hash_point> points_already_added;

    // we only want to add the points ONCE
    auto add_if_not_adlready_done= [&points_already_added, &points_convex_hull, vh](Point_2 s){
        if(points_already_added.find(s) == points_already_added.end()){
            points_already_added[s] = true;
            points_convex_hull.push_back(s);
        }
    };

    // some intersection routine to prevent code repetition
    auto insert_intersection_segseg = [add_if_not_adlready_done](Segment_2 seg, Segment_2 target){
        auto intersect = CGAL::intersection(seg, target);
        if(intersect){
            Point_2 *intersect_point = boost::get<Point_2>(&*intersect);
            add_if_not_adlready_done(*intersect_point);
        }
    };
    auto insert_intersection_rayseg = [add_if_not_adlready_done](Ray_2 r, Segment_2 target){
        auto intersect = CGAL::intersection(r, target);
        if(intersect){
            Point_2 *intersect_point = boost::get<Point_2>(&*intersect);
            add_if_not_adlready_done(*intersect_point);
        }
    };
    auto insert_intersection_lineseg = [add_if_not_adlready_done](Line_2 l, Segment_2 target){
        auto intersect = CGAL::intersection(l, target);
        if(intersect){
            Point_2 *intersect_point = boost::get<Point_2>(&*intersect);
            add_if_not_adlready_done(*intersect_point);
        }
    };

    // We circulate around the edges of the current vertex
    Edge_circulator begin = vh->incident_edges();
    Edge_circulator current = begin;

    do{
        if(is_infinite(current)){ // we do not want to deal with infinite edges
            current++;
            continue;
        }
        auto d = dual(*current);
        Segment_2 seg;
        Line_2 lin;
        Ray_2 ray;
        if(assign(seg, d)){ // if the dual is a segment
            if(contained_in_window(seg.target())){ // we only add the target and the source if they are inside the view
                add_if_not_adlready_done(seg.target());
            }else{
                insert_intersection_segseg(seg, topWindow);
                insert_intersection_segseg(seg, bottomWindow);
                insert_intersection_segseg(seg, rightWindow);
                insert_intersection_segseg(seg, leftWindow);
            }
            if(contained_in_window(seg.source())){
                add_if_not_adlready_done(seg.source());
            }else{
                insert_intersection_segseg(seg, topWindow);
                insert_intersection_segseg(seg, bottomWindow);
                insert_intersection_segseg(seg, rightWindow);
                insert_intersection_segseg(seg, leftWindow);
            }

        }else if(assign(lin, d)){ // a line
            insert_intersection_lineseg(lin, topWindow);
            insert_intersection_lineseg(lin, bottomWindow);
            insert_intersection_lineseg(lin, leftWindow);
            insert_intersection_lineseg(lin, rightWindow);
        }else if(assign(ray, d)){ // a ray
            if(contained_in_window(ray.source())){
                add_if_not_adlready_done(ray.source());
            }
            insert_intersection_rayseg(ray, topWindow);
            insert_intersection_rayseg(ray, bottomWindow);
            insert_intersection_rayseg(ray, leftWindow);
            insert_intersection_rayseg(ray, rightWindow);
        }else{ // We should never go here
            D("Error") 
        }
        
        current++;
    }while(current != begin);
    
    Point_2 c = CGAL::centroid(points_convex_hull.begin(), points_convex_hull.end());
    return c;
}


int MTriangulation::move_step(QRectF rect){
    float rMax = (rect.height() + rect.width())/(2*100);
    QTime timer;
    timer.start();

    std::vector<VertexMoveHint> newPointsHint;

    Vertex_handle closest_tl, closest_tr, closest_bl, closest_br;

    // compute the NN of the corner of the viewpoint : it will have to be taken in account for Lloyd step
    /// @TODO Does not work in some case : if we zoomed and the NN in outside the viewpoint
    if(mStyle == LLOYD){
        closest_tl = nearest_vertex(Point_2(rect.left(), rect.top()));
        closest_bl = nearest_vertex(Point_2(rect.left(), rect.bottom()));
        closest_tr = nearest_vertex(Point_2(rect.right(), rect.top()));
        closest_br = nearest_vertex(Point_2(rect.right(), rect.bottom()));
    }

    //main iteration
    for(auto it = all_vertices_begin(); it != all_vertices_end(); ++it){
        if(not is_infinite(it)){
            Point_2 p = it->point();
            Point_2 nPoint;

            switch(mStyle){
            case BROWNIAN:
                nPoint = brownianStep(p, rMax);
                break;
            case BOUNCING_BALL:
                nPoint = bouncingBallStep(it, rMax, rect);
                break;
            case LLOYD:
                // manage the case of the corners
                std::vector<Point_2> corner_points;
                if(it == closest_bl){
                    corner_points.push_back(Point_2(rect.left(), rect.bottom()));
                }
                if(it == closest_tl){
                    corner_points.push_back(Point_2(rect.left(), rect.top()));
                }
                if(it == closest_br){
                    corner_points.push_back(Point_2(rect.right(), rect.bottom()));
                }
                if(it == closest_tr){
                    corner_points.push_back(Point_2(rect.right(), rect.top()));
                }
                nPoint = lloydStep(it, rect, corner_points);
                break;
            }

            switch(iStyle){
                case MOVE_CGAL:
                    move(it, nPoint);
                    update_nn(it);
                break;
                case CLASSIC:
                case HINT:
                    newPointsHint.emplace_back(it, nPoint);
                break;
            }
        }
    }

    double improv = 0;


    switch(iStyle){
        case MOVE_CGAL:
        break;
        // If we are using the double buffering, we need to clear it
        case CLASSIC:
            clear();
            current_insert = 1-current_insert;

            insert_classic(newPointsHint);

            vertex_ball[1-current_insert].clear();
            nearest_neight[1-current_insert].clear();
            nearest_neight_sqdistance[1-current_insert].clear();

        break;
        case HINT:
            Delaunay::clear();
            current_insert = 1-current_insert;

            improv = insert_hint(newPointsHint);

            nearest_neight[1-current_insert].clear();
            nearest_neight_sqdistance[1-current_insert].clear();
            vertex_ball[1-current_insert].clear();
        break;
    }

    // Print the improvement, if we use the hint method
    int elapsed = timer.elapsed();
    if(iStyle == HINT){
        std::cout << "Fraction of improvement : " << 100*improv << "%" << std::endl;
    }

    return elapsed;
}

void MTriangulation::setInsertStyle(InsertStyle is){
    iStyle = is;
}

void MTriangulation::setMovingStyle(MovingStyle ms){
    mStyle = ms;
}

void MTriangulation::insert_classic(std::vector<Point_2>& points){
    CGAL::spatial_sort(points.begin(), points.end());
    Face_handle face_hint = Face_handle();
    for(Point_2 p : points){
        Vertex_handle v_handle = insert(p, face_hint);
        face_hint = v_handle->face();
    }
}

void MTriangulation::insert_classic(std::vector<VertexMoveHint>& points){
    CGAL::spatial_sort(points.begin(), points.end(), VertexMoveHintCompTrait());
    Face_handle face_hint = Face_handle();
    for(VertexMoveHint vh : points){
        Point_2 p = vh.new_position;
        Vertex_handle v_handle = insert(p, face_hint, vertex_ball[1-current_insert][vh.handle]);
        face_hint = v_handle->face();
    }
}


double MTriangulation::insert_hint(std::vector<VertexMoveHint>& newPointsHint){
    using namespace std;
    unsigned long long number_improvement = 0;

    // First, we spatialy sort the points
    CGAL::spatial_sort(newPointsHint.begin(), newPointsHint.end(), VertexMoveHintCompTrait());

    std::unordered_map<Vertex_handle, Vertex_handle> newVertexs;

    Vertex_handle vertex_hint = Vertex_handle();

    for(const VertexMoveHint& h : newPointsHint){
        Vertex_handle old_vertex = h.handle;
        Point_2 nPos = h.new_position;

        Vertex_handle prev_nn = nearest_neight[1-current_insert][old_vertex];

        auto hint_it = newVertexs.find(prev_nn);

        // if our previous NN has already entered the triangulation
        if(hint_it != newVertexs.end()){
            Vertex_handle hint_nn = hint_it->second;
            //if we better chose it as hint, we do so
            if(vertex_hint == Vertex_handle() or sqdist(hint_nn->point(), nPos) < sqdist(vertex_hint->point(), nPos)){
                vertex_hint = hint_nn;
                number_improvement++;
            }
        }

        // Here we insert using the hint (if it is a nullptr, we do not use any hint)
        Vertex_handle v_handle = insert(nPos, 
                                    vertex_hint == Vertex_handle()? Face_handle() : vertex_hint->face(), 
                                    vertex_ball[1-current_insert][old_vertex]);
        //the default hint for the next vertex is the current one (because of spacial sorting)
        vertex_hint = v_handle;
        // This vertex is now inserted, we will be able to use it as hint
        newVertexs[old_vertex] = v_handle;
    }
    return (double)number_improvement/newPointsHint.size();
}