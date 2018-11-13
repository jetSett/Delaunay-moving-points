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
    ui(new Ui::MainWindow), triangulation(NAIVE)
{
    ui->setupUi(this);

    tri_graphics = new CGAL::Qt::TriangulationGraphicsItem<MTriangulation>(&triangulation);

    scene.addItem(tri_graphics);

    ui->graphicsView->setScene(&scene);

    QObject::connect(this, SIGNAL(changed()),
             tri_graphics, SLOT(modelChanged()));


    pi = new TriangulationPointInputAndConflictZone<MTriangulation>(&scene, &triangulation, this );

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

    triangulation.insert(points.begin(), points.end());
    log(QString("Inserted ") + QString::number(number) + 
        QString(" points in ") + QString::number(timer.elapsed()) + QString("ms\n"));
    QApplication::restoreOverrideCursor();
    Q_EMIT( changed());

}

void MainWindow::on_clearPushButton_clicked(bool){
    triangulation.clear();
    log("Triangulation cleared\n");
    Q_EMIT(changed());
}

void MainWindow::move(){
    int elapsed = 0;
    if(ui->brownianRadioButton->isChecked()){
        QRectF rect = CGAL::Qt::viewportsBbox(&scene);;
        float maxStep = (rect.height() + rect.width())/(2*50); // We are not going to go too far
        elapsed = triangulation.moveBrownian(maxStep);
    }else{
        QMessageBox::critical(this, "Not implemented", "This functionality was not implemented yet");
        on_stopButton_clicked(true);
    }
    log(QString("Movement in ") + QString::number(elapsed) + QString(" ms\n"));
    Q_EMIT(changed());

}

void MainWindow::on_delaySpinBox_valueChanged(int ms){
    timer.setInterval(ms);
}

void MainWindow::on_startButton_clicked(bool){
    ui->startButton->setEnabled(false);
    ui->stopButton->setEnabled(true);
    timer.start(ui->delaySpinBox->value());
}

void MainWindow::on_stopButton_clicked(bool){
    ui->startButton->setEnabled(true);
    ui->stopButton->setEnabled(false);
    timer.stop();
}

void MainWindow::log(const QString& text){
    ui->outputText->verticalScrollBar()->setValue(ui->outputText->verticalScrollBar()->maximum());
    ui->outputText->insertPlainText(text);
}

void MainWindow::on_naiveRadioButton_clicked(){
    triangulation.setInsertStyle(NAIVE);
}

void MainWindow::on_hintRadioButton_clicked(){
    triangulation.setInsertStyle(HINT);
}