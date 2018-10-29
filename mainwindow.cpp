#include <iostream>
#include "mainwindow.hpp"
#include "ui_mainwindow.h"

using namespace std;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    triangulation.insert(Point_2(100, 100));
    triangulation.insert(Point_2(200, 100));
    triangulation.insert(Point_2(100, 200));
    triangulation.insert(Point_2(300, 200));

    tri_graphics = new CGAL::Qt::TriangulationGraphicsItem<Delaunay>(&triangulation);

    scene.addItem(tri_graphics);

    ui->graphicsView->setScene(&scene);

    QObject::connect(this, SIGNAL(changed()),
             tri_graphics, SLOT(modelChanged()));


    pi = new TriangulationPointInputAndConflictZone<Delaunay>(&scene, &triangulation, this );

    scene.installEventFilter(pi);

    QObject::connect(pi, SIGNAL(generate(CGAL::Object)),
             this, SLOT(processInput(CGAL::Object)));

}

MainWindow::~MainWindow()
{
    delete tri_graphics;
    delete ui;
}

void MainWindow::processInput(CGAL::Object obj){
    cout << "Input !" << endl;
    Point_2 p;
    if(CGAL::assign(p, obj)){
      cout << p.x() << " : " << p.y() << endl;
      triangulation.insert(p);
    }
    Q_EMIT( changed());
}
