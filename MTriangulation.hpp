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
#include <CGAL/intersection_2.h>
#include <CGAL/Object.h>

enum InsertStyle{
    CLASSIC, HINT, MOVE_CGAL
};

enum MovingStyle{
    BROWNIAN, BOUNCING_BALL, LLOYD
};

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Vector_2 Vector_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;
typedef K::Segment_2 Segment_2;
typedef K::Ray_2 Ray_2;
typedef K::Line_2 Line_2;

typedef CGAL::Delaunay_triangulation_2<K> Delaunay;

typedef CGAL::Delaunay_triangulation_adaptation_traits_2<Delaunay> AdaptTrait;
typedef CGAL::Voronoi_diagram_2<Delaunay, AdaptTrait> Voronoi;

struct VertexMoveHint;

/**
 * @brief Implémentation of the Delaunay triangulation in order to make point move
 * 
 */

class MTriangulation : public Delaunay{
    /**
     * @brief Types used to store the nearest neightboard and their distance to each vertex
     * @note Vertex_handle is a pointer in disguise, so very easily hashable
     */
    typedef std::unordered_map<Vertex_handle, Vertex_handle> Nn_map;
    typedef std::unordered_map<Vertex_handle, double> Nn_dist_map;

    typedef std::unordered_map<Vertex_handle, Vector_2> Ball_map;
public:
    MTriangulation(InsertStyle, MovingStyle);

    /**
     * @brief Execute a move step
     * @arg sceneBox the current scenebox
     * @return int time elapsed in ms
     */
    int move_step(QRectF sceneBox);

    void setInsertStyle(InsertStyle);
    void setMovingStyle(MovingStyle);

    /**
     * @brief Insert points in the "classical" way : spacial sorting + insertion with previous point as hint
     * 
     */
    void insert_classic(std::vector<Point_2>&);

    /**
     * @brief Insert a new point in the triangulation
     * 
     * @param hint hint for the insertion
     * @param ball the ball attributed to the point for bouncing ball (if 0, will be chosen randomly)
     * @return Vertex_handle a handle to the newly created point
     */
    Vertex_handle insert(const Point_2&, const Face_handle& hint = Face_handle(), Vector_2 ball = Vector_2(0, 0));

    void clear();

private:
    /**
     * @brief does a jumping ball movement step
     * 
     * @param rMax the maximum step of each movement
     * @param sceneBox the current scenebox
     * @return Point_2 the new position of the vertex
     */
    Point_2 brownianStep(Point_2, float rMax);

    /**
     * @brief does a bouncing ball movement step
     * 
     * @param speed the maximum step of each movement
     * @param sceneBox the current scenebox (used for bouncing)
     * @return Point_2 the new position of the vertex
     */
    Point_2 bouncingBallStep(Vertex_handle, float speed, QRectF sceneBox);

    /**
     * @brief does a lloyd movement step
     * 
     * @param sceneBox the current scenebox
     * @param points_convex_hull offset for the points into the convex hull (used to handle corner of the viewport)
     * @return Point_2 the new position of the vertex
     */
    Point_2 lloydStep(Vertex_handle, QRectF sceneBox,  std::vector<Point_2>& points_convex_hull);

    /**
     * @brief Insert points in a "classic" way : spacial sorting + insertion with previous hint (with Hint for implementation purpose)
     * 
     */
    void insert_classic(std::vector<VertexMoveHint>&);

    /**
     * @brief Update the NN of a vertex and the vertex itself in order to maintain nearest_neight and nearest_neight_sqdistance
     * 
     */
    void update_nn(Vertex_handle);


    /**
     * @brief Insert the point using the "hint method" : first spacial sorting of the point, then when it is useful, use the previous NN's new position as hint
     * 
     * @return double the fraction of "New hint" used
     */
    double insert_hint(std::vector<VertexMoveHint>&);

    InsertStyle iStyle;
    MovingStyle mStyle;

    int current_insert;
    /**
     * @brief Those arrays are here to work as double buffering : we swap them at each new moving (when not using CGAL::move)
     * 
     */
    std::array<Nn_map, 2> nearest_neight;
    std::array<Nn_dist_map, 2> nearest_neight_sqdistance;
    /**
     * @brief The direction and intensity of movement of the points for MOVING_POINT moving style
     * 
     */
    std::array<Ball_map, 2> vertex_ball;
};

/**
 * @brief Structure agregating a Vertex_handle (the previous vertex pointer) and the new position of this vertex
 * 
 */
struct VertexMoveHint{
    VertexMoveHint(const MTriangulation::Vertex_handle& h, const Point_2& p) : handle(h), new_position(p){
    }
    MTriangulation::Vertex_handle handle;
    Point_2 new_position;
};
