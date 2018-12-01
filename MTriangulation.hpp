#include <unordered_map>
#include <tuple>
#include <array>
#include <functional>
#include <QtCore>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>


enum InsertStyle{
    NAIVE, HINT, MOVE_CGAL
};

enum MovingStyle{
    BROWNIAN, JUMPING_BALL, LLOYD
};

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Vector_2 Vector_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;

typedef CGAL::Delaunay_triangulation_2<K> Delaunay;

typedef CGAL::Delaunay_triangulation_adaptation_traits_2<Delaunay> AdaptTrait;
typedef CGAL::Voronoi_diagram_2<Delaunay, AdaptTrait> Voronoi;

struct Hash_point{
    std::size_t operator()(const Point_2& p) const noexcept;
};

struct VertexMoveHint;

class MTriangulation : public Delaunay{

    typedef std::map<Vertex_handle, Vertex_handle> Nn_map;
    typedef std::map<Vertex_handle, double> Nn_dist_map;

    typedef std::map<Vertex_handle, Vector_2> Ball_map;
public:
    MTriangulation(InsertStyle, MovingStyle);

    int move_step(QRectF);

    Point_2 brownianStep(Point_2, float);
    Point_2 jumpBallStep(Vertex_handle, float, QRectF);
    Point_2 lloydStep(Vertex_handle, QRectF, Voronoi&);

    void setInsertStyle(InsertStyle);
    void setMovingStyle(MovingStyle);

    void insert_naive(std::vector<Point_2>&);
    void insert_naive(std::vector<VertexMoveHint>&);

    Vertex_handle insert(const Point_2&, const Face_handle& f = Face_handle(), Vector_2 ball = Vector_2(0, 0));

    void clear();

private:

    void update_nn(Vertex_handle);

    double insert_hint(std::vector<VertexMoveHint>&);

    InsertStyle iStyle;
    MovingStyle mStyle;

    int current_insert;

    std::array<Nn_map, 2> nearest_neight;
    std::array<Nn_dist_map, 2> nearest_neight_sqdistance;

    std::array<Ball_map, 2> vertex_ball;
};

struct VertexMoveHint{
    VertexMoveHint(const MTriangulation::Vertex_handle& h, const Point_2& p) : handle(h), new_position(p){
    }
    MTriangulation::Vertex_handle handle;
    Point_2 new_position;
};
