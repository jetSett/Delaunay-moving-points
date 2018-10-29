#ifndef MAINWINDOW_HPP
#define MAINWINDOW_HPP

#include <QMainWindow>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/Qt/TriangulationGraphicsItem.h>
#include <CGAL/Qt/VoronoiGraphicsItem.h>

#include "pointinput.hpp"

namespace Ui {
class MainWindow;
}

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
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

signals:
    void changed();
private:    

    TriangulationPointInputAndConflictZone<Delaunay> * pi;


    QGraphicsScene scene;

    Ui::MainWindow *ui;

    Delaunay triangulation;

    CGAL::Qt::TriangulationGraphicsItem<Delaunay>* tri_graphics;

};

#endif // MAINWINDOW_HPP
