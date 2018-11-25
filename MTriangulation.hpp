#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Vector_2 Vector_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;

typedef CGAL::Delaunay_triangulation_2<K> Delaunay;

enum InsertStyle{
    NAIVE, HINT, MOVE_CGAL
};

class MTriangulation : public Delaunay{
public:
    MTriangulation(InsertStyle);

    int moveBrownian(float rMax);

    void setInsertStyle(InsertStyle);

private:

    void insert_move(const std::vector<Point_2>&, const std::vector<Vertex_handle>&);

    void insert_naive(const std::vector<Point_2>&);
    void insert_hint(const std::vector<Point_2>&, const std::vector<Vertex_handle>&);

    InsertStyle iStyle;
};