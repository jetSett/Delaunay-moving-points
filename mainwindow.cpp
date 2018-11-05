#include <iostream>
#include "mainwindow.hpp"
#include "ui_mainwindow.h"

using namespace std;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    tri_graphics = new CGAL::Qt::TriangulationGraphicsItem<Delaunay>(&triangulation);

    scene.addItem(tri_graphics);

    ui->graphicsView->setScene(&scene);

    QObject::connect(this, SIGNAL(changed()),
             tri_graphics, SLOT(modelChanged()));


    pi = new TriangulationPointInputAndConflictZone<Delaunay>(&scene, &triangulation, this );

    scene.installEventFilter(pi);

    scene.setSceneRect(-100, -100, 100, 100);


    QObject::connect(pi, SIGNAL(generate(CGAL::Object)),
             this, SLOT(processInput(CGAL::Object)));

    // Navigation
      this->addNavigation(ui->graphicsView);


}

MainWindow::~MainWindow()
{
    delete tri_graphics;
    delete ui;
}

void MainWindow::processInput(CGAL::Object obj){
    Point_2 p;
    if(CGAL::assign(p, obj)){
      triangulation.insert(p);
    }
    Q_EMIT( changed());
}

void MainWindow::addNavigation(QGraphicsView* graphicsView){
  navigation = new CGAL::Qt::GraphicsViewNavigation();
  graphicsView->viewport()->installEventFilter(navigation);
  graphicsView->installEventFilter(navigation);
}