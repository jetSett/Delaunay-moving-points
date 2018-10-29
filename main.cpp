#include "mainwindow.hpp"
#include <QApplication>
#include <CGAL/Qt/resources.h>

int main(int argc, char *argv[])
{

    CGAL_QT_INIT_RESOURCES;

    QApplication a(argc, argv);

    MainWindow w;
    w.show();

    return a.exec();
}
