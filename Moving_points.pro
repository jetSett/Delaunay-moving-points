#-------------------------------------------------
#
# Project created by QtCreator 2018-10-27T17:35:52
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = Moving_points
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp\
        MTriangulation.cpp

HEADERS  += mainwindow.hpp \
    pointinput.hpp

FORMS    += mainwindow.ui

LIBS        += -lCGAL
LIBS        += -lCGAL_Qt5
LIBS        += -lgmp
LIBS        += -lmpfr
QMAKE_CXXFLAGS += -frounding-math -O3 -std=c++14
