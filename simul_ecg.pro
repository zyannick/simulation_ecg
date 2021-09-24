TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS += -O3
QMAKE_CXXFLAGS += -DNDEBUG

SOURCES += \
        main.cc \
        opt.c

#INCLUDEPATH += /home/neil/Workspace/CPP/simul_ecg/csv_include

LIBS += -L/usr/local/lib
LIBS += -lstdc++fs
LIBS += -lpthread

DISTFILES += \
    00README \
    index.html

HEADERS += \
    dfour1.h \
    dir_utils.hh \
    opt.h \
    options.hh \
    rand1.h \
    single_include/csv.hpp
