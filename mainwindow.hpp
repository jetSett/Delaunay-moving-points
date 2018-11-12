#ifndef MAINWINDOW_HPP
#define MAINWINDOW_HPP

#include <QMainWindow>
#include <QTimer>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/Qt/TriangulationGraphicsItem.h>
#include <CGAL/Qt/VoronoiGraphicsItem.h>

#include <CGAL/Qt/GraphicsViewNavigation.h>
#include "pointinput.hpp"

namespace Ui {
class MainWindow;
}

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Vector_2 Vector_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;

typedef CGAL::Delaunay_triangulation_2<K> Delaunay;


class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();


public slots:
    void processInput(CGAL::Object);

    void on_addRandomPointButton_clicked(bool checked);
    void on_clearPushButton_clicked(bool checked);

    void on_startButton_clicked(bool checked);
    void on_stopButton_clicked(bool checked);

    void move();

signals:
    void changed();
private:    

    void addNavigation(QGraphicsView* graphicsView);
    
    void log(const QString&);

    void moveBrownian();    
    void addRandomPoints(unsigned int number);

    QTimer timer;

    TriangulationPointInputAndConflictZone<Delaunay> * pi;

    CGAL::Qt::GraphicsViewNavigation* navigation;
    QGraphicsScene scene;

    Ui::MainWindow *ui;

    Delaunay triangulation;

    std::vector<Point_2> points_triangulation;

    CGAL::Qt::TriangulationGraphicsItem<Delaunay>* tri_graphics;

};

#endif // MAINWINDOW_HPP
