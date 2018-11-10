#include <iostream>
#include <QTime>

#include <QMessageBox>

#include <CGAL/point_generators_2.h>

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

void MainWindow::on_addRandomPointButton_clicked(bool){
    addPoints(ui->addRandomPointSpinBox->value());
}

void MainWindow::addPoints(unsigned int number){
    QRectF rect = CGAL::Qt::viewportsBbox(&scene);
    CGAL::Qt::Converter<K> convert;  
    CGAL::Iso_rectangle_2<K> isor = convert(rect);
    CGAL::Random_points_in_iso_rectangle_2<Point_2> pg((isor.min)(), (isor.max)());

    QApplication::setOverrideCursor(Qt::WaitCursor);
    std::vector<Point_2> points;
    points.reserve(number);
    for(unsigned i = 0; i < number; ++i){
        points.push_back(*pg++);
    }
    QTime timer;
    timer.start();
    triangulation.insert(points.begin(), points.end());
    ui->outputText->insertPlainText(QString("Inserted ") + QString::number(number) + 
        QString(" points in ") + QString::number(timer.elapsed()) + QString("ms\n"));
    QApplication::restoreOverrideCursor();
    Q_EMIT( changed());

}

void MainWindow::on_clearPushButton_clicked(bool){
    triangulation.clear();
    ui->outputText->insertPlainText("Triangulation cleared\n");
    Q_EMIT(changed());
}

void MainWindow::on_moveButton_clicked(bool){
    QMessageBox::critical(this, "Not implemented", "This functionality was not implemented yet");
}