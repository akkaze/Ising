TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    model_ising2dsqrmet.cpp \
    simrun.cpp \
    isingspin.cpp \
    sched.cpp \
    qstatepotts.cpp

HEADERS += \
    model_ising2dsqrmet.hpp \
    simrun.hpp \
    sim_datastruct.hpp \
    msgtag.hpp \
    isingspin.hpp \
    sched.hpp \
    qstatepotts.h

INCLUDEPATH += /usr/include/openmpi
LIBS += -L/usr/lib/ -lgsl -lgslcblas -lboost_mpi -lmpi -lboost_serialization -lmpi_cxx -lboost_system -lboost_timer
