######################################################################
# Automatically generated by qmake (1.07a) Wed Jan 13 16:14:23 2010
######################################################################

TEMPLATE = app
CONFIG -= moc qt
DEPENDPATH += . ../../include
INCLUDEPATH += . ../../include

unix{
   INCLUDEPATH += /home/pompeu/local64/include
   LIBS += -L/home/pompeu/local64/lib -lfftw3 -lsndfile
}
macx{
   CONFIG -= app_bundle
   INCLUDEPATH += /Users/mtcheou/local/include
   LIBS += -L/Users/mtcheou/local/lib -lfftw3 -lsndfile
}

LIBS += -L../../lib/ -ladsp

QMAKE_CFLAGS_RELEASE   = -pipe -D_FORTIFY_SOURCE=2 \
                         -fomit-frame-pointer
QMAKE_CXXFLAGS_RELEASE = -pipe -D_FORTIFY_SOURCE=2 \
                         -fomit-frame-pointer


# Input
HEADERS += ./src/arithmetic_codec.h\
           ./src/encode.h\
           ./src/adspdefines.h\
           ./src/adspsort.h\
           ./src/qstep.h
           
SOURCES += ./src/main.cpp \
           ./src/arithmetic_codec.cpp\
           ./src/encode.cpp\
           ./src/adspsort.cpp\
           ./src/qstep.cpp

OBJECTS_DIR  = obj
TARGET       = encode

target.path = ../../bin
INSTALLS += target
 
