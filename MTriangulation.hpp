#include <unordered_map>
#include <tuple>
#include <functional>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>


enum InsertStyle{
    NAIVE, HINT, MOVE_CGAL
};

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Vector_2 Vector_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;

typedef CGAL::Delaunay_triangulation_2<K> Delaunay;

struct Hash_point{
    std::size_t operator()(const Point_2& p) const noexcept;
};

struct VertexMoveHint;

class MTriangulation : public Delaunay{
typedef std::unordered_multimap<Vertex_handle, Vertex_handle> Hint_insertion;

public:
    MTriangulation(InsertStyle);

    int moveBrownian(float rMax);

    void setInsertStyle(InsertStyle);

    void insert_naive(std::vector<Point_2>&);

    Vertex_handle insert(const Point_2&, const Face_handle& f = Face_handle());

    void clear();

private:

    void insert_hint(std::vector<VertexMoveHint>&, std::unordered_map<Vertex_handle, Vertex_handle>& old_nn);

    InsertStyle iStyle;

    std::unordered_map<Vertex_handle, Vertex_handle> nearest_neight;
    std::unordered_map<Vertex_handle, double> nearest_neight_sqdistance;
};

struct VertexMoveHint{
    VertexMoveHint(const MTriangulation::Vertex_handle& h, const Point_2& p) : handle(h), new_position(p){
    }
    MTriangulation::Vertex_handle handle;
    Point_2 new_position;
};
