#include "mainwindow.hpp"
#include <QApplication>
#include <cstdlib>
#include <ctime>
#include <CGAL/Qt/resources.h>

int main(int argc, char *argv[])
{
    srand(time(nullptr));
    CGAL_QT_INIT_RESOURCES;

    QApplication a(argc, argv);

    MainWindow w;
    w.show();

    return a.exec();
}
