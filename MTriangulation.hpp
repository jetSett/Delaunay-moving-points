#include <unordered_map>
#include <tuple>
#include <functional>
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

struct Hash_point{
    std::size_t operator()(const Point_2& p) const noexcept;
};

using Hint_insertion = std::unordered_multimap<Point_2, Point_2, Hash_point>;

class MTriangulation : public Delaunay{
public:
    MTriangulation(InsertStyle);

    int moveBrownian(float rMax);

    void setInsertStyle(InsertStyle);

private:

    void insert_move(std::vector<Point_2>&, const Hint_insertion&);

    void insert_naive(std::vector<Point_2>&);
    void insert_hint(std::vector<Point_2>&, const Hint_insertion&);

    InsertStyle iStyle;
};