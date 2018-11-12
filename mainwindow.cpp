#include <iostream>
#include <cstdlib>
#include <cmath>
#include <QTime>

#include <QMessageBox>
#include <QScrollBar>

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

      QObject::connect(&timer, SIGNAL(timeout()), this, SLOT(move()));
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
      points_triangulation.push_back(p);
    }
    Q_EMIT( changed());
}

void MainWindow::addNavigation(QGraphicsView* graphicsView){
  navigation = new CGAL::Qt::GraphicsViewNavigation();
  graphicsView->viewport()->installEventFilter(navigation);
  graphicsView->installEventFilter(navigation);
}

void MainWindow::on_addRandomPointButton_clicked(bool){
    addRandomPoints(ui->addRandomPointSpinBox->value());
}

void MainWindow::addRandomPoints(unsigned int number){
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
    for(auto &p : points){
        points_triangulation.push_back(p);
    }

    triangulation.insert(points.begin(), points.end());
    log(QString("Inserted ") + QString::number(number) + 
        QString(" points in ") + QString::number(timer.elapsed()) + QString("ms\n"));
    QApplication::restoreOverrideCursor();
    Q_EMIT( changed());

}

void MainWindow::on_clearPushButton_clicked(bool){
    triangulation.clear();
    points_triangulation.clear();
    log("Triangulation cleared\n");
    Q_EMIT(changed());
}

void MainWindow::move(){
    if(ui->brownianRadioButton->isChecked()){
        for(int i = 0; i<ui->numberIterationsSpinBox->value(); ++i){
            moveBrownian();
        }
    }else{
        QMessageBox::critical(this, "Not implemented", "This functionality was not implemented yet");
        on_stopButton_clicked(true);
    }

}

void MainWindow::on_startButton_clicked(bool){
    ui->startButton->setEnabled(false);
    ui->stopButton->setEnabled(true);
    timer.start(200);
}

void MainWindow::on_stopButton_clicked(bool){
    ui->startButton->setEnabled(true);
    ui->stopButton->setEnabled(false);
    timer.stop();
}

void MainWindow::moveBrownian(){
    QRectF rect = CGAL::Qt::viewportsBbox(&scene);;
    float maxStep = (rect.height() + rect.width())/(2*50); // We are not going to go too far
    std::vector<Point_2> newPoints;
    //std::vectr<Point_2_> nn_before;
    for(const Point_2 &p : points_triangulation){
        //nn_before.push_back(triangulation.nearest_vertex(p));
        float r = 2*(static_cast <float> (rand()) / static_cast <float> (RAND_MAX))-1; // random [-1, 1]   
        float theta = 2*3.1415*(static_cast <float> (rand()) / static_cast <float> (RAND_MAX)); // random [0, 2pi] 
        Vector_2 v(cos(theta), sin(theta));
        v *= r*maxStep;
        newPoints.push_back(p+v);
    }
    triangulation.clear();
    points_triangulation.clear();

    QTime timer;
    timer.start();
    triangulation.insert(newPoints.begin(), newPoints.end());
    log(QString("Movement in ") + QString::number(timer.elapsed()) + QString(" ms\n"));
    for(auto &p : newPoints){
        points_triangulation.push_back(p);
    }
    Q_EMIT(changed());
}

void MainWindow::log(const QString& text){
    ui->outputText->insertPlainText(text);
    ui->outputText->verticalScrollBar()->setValue(ui->outputText->verticalScrollBar()->maximum());
}