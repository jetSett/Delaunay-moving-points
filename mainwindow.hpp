#ifndef MAINWINDOW_HPP
#define MAINWINDOW_HPP

#include <QMainWindow>
#include <QTimer>

#include <CGAL/Qt/TriangulationGraphicsItem.h>
#include <CGAL/Qt/VoronoiGraphicsItem.h>

#include <CGAL/Qt/GraphicsViewNavigation.h>

#include "MTriangulation.hpp"

#include "pointinput.hpp"

namespace Ui {
class MainWindow;
}

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

    void on_delaySpinBox_valueChanged(int ms);

    void on_naiveRadioButton_clicked();
    void on_hintRadioButton_clicked();
    void on_cgalMoveButton_clicked();

    void on_brownianRadioButton_clicked();
    void on_jumpingBallRadioButton_clicked();
    void on_lloydRadioButton_clicked();

    void on_showVoronoiCheckBox_toggled(bool);

    void move();

signals:
    void changed();
private:    

    void addNavigation(QGraphicsView* graphicsView);
    
    void log(const QString&);

    void addRandomPoints(unsigned int number);

    QTimer timer;

    TriangulationPointInputAndConflictZone<MTriangulation> * pi;

    CGAL::Qt::GraphicsViewNavigation* navigation;
    QGraphicsScene scene;

    Ui::MainWindow *ui;

    MTriangulation triangulation;

    CGAL::Qt::TriangulationGraphicsItem<MTriangulation>* tri_graphics;
    CGAL::Qt::VoronoiGraphicsItem<Delaunay>* vor_graphics;

};

#endif // MAINWINDOW_HPP
