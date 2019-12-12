TEMPLATE = app
GRAPHICS = x11
CONFIG += console
CONFIG += release   #default +
CONFIG -= debug     #default -
CONFIG -= app_bundle


#CONFIG += optimize_full # trying - not doing anything
#CONFIG += warn_off  # this is to suppress all warnings
# QMAKE_CFLAGS -= -O2
# QMAKE_CFLAGS -= -O1
# QMAKE_CFLAGS -= -O
# QMAKE_CFLAGS *= -m64 -O3
#
# QMAKE_CXXFLAGS -= -O2
# QMAKE_CXXFLAGS -= -O1
# QMAKE_CXXFLAGS -= -O
# QMAKE_CXXFLAGS *= -m64 -O3
#
# QMAKE_LFLAGS *= -m64 -O3

# QMAKE_CXXFLAGS_RELEASE -= -O
# QMAKE_CXXFLAGS_RELEASE -= -O1
# QMAKE_CXXFLAGS_RELEASE -= -O2
# QMAKE_CXXFLAGS_RELEASE *= -O3

# the following removes specific warnings
QMAKE_CXXFLAGS_WARN_ON += -Wno-unused-parameter
QMAKE_CXXFLAGS_WARN_ON += -Wno-unused-variable
QMAKE_CXXFLAGS_WARN_ON += -Wno-ignored-qualifiers
QMAKE_CXXFLAGS_WARN_ON += -Wno-unused-result
QMAKE_CXXFLAGS_WARN_ON += -Wno-write-strings
QMAKE_CXXFLAGS_WARN_ON += -Wno-misleading-indentation

QMAKE_CXXFLAGS += -std=c++11

TARGET = cell_evolution
MAINFILE = $$join(TARGET, " ", , ".cpp" )

message( $$MAINFILE )
message( $$TARGET )
# Input
HEADERS += ca.h \
	   hull.h \
		       gene.h \
					 genome.h \
           cell.h \
           conrec.h \
           dish.h \
           graph.h \
           info.h \
           misc.h \
           output.h \
           parameter.h \
           parse.h \
           pde.h \
           intplane.h \
           random.h \
           sqr.h \
           sticky.h \
       	   crash.h \
       	   warning.h \

SOURCES += ca.cpp \
	   hull.cpp \
		       gene.cpp \
					 genome.cpp \
           cell.cpp \
           conrec.cpp \
           dish.cpp \
           info.cpp \
           misc.cpp \
           output.cpp \
           parameter.cpp \
           parse.cpp \
           pde.cpp \
           intplane.cpp \
           random.cpp \
           crash.cpp \
           warning.cpp \

SOURCES += $$MAINFILE

#QMAKE_CXXFLAGS_RELEASE += -fexceptions
#QMAKE_CXXFLAGS_DEBUG += -fexceptions
#QMAKE_LFLAGS_RELEASE += -O4
#QMAKE_CXXFLAGS_RELEASE += -O4

contains( GRAPHICS, qt ) {
   message( "Building Qt executable" )
   SOURCES += qtgraph.cpp
   HEADERS += qtgraph.h
   QMAKE_CXXFLAGS_RELEASE += -DQTGRAPHICS
   QMAKE_CXXFLAGS_DEBUG += -DQTGRAPHICS
   QT += qt3support
   #unix {
   #   system(rm $$TARGET.o)
   #}
   win32 {
     QMAKE_LFLAGS += -L "\"C:\Program Files\GnuWin32\lib\"" -lpng -lz
     QMAKE_CXXFLAGS += -I "\"C:\Program Files\GnuWin32\include\""
   }
   #LIBS += -lpng
}

contains( GRAPHICS, qt3 ) {
   message( "Building Qt executable" )
   SOURCES += qt3graph.cpp
   HEADERS += qt3graph.h
   QMAKE_CXXFLAGS_RELEASE += -DQTGRAPHICS
   QMAKE_CXXFLAGS_DEBUG += -DQTGRAPHICS
   #unix {
   #   system(rm $$TARGET.o)
   #}
   win32 {
     QMAKE_LFLAGS += -L "C:\Program Files\GnuWin32\lib" -lpng -lz
     QMAKE_CXXFLAGS += -I "C:\Program Files\GnuWin32\include"
   }
   LIBS += -lpng
}
contains( GRAPHICS, x11 ) {
   !unix {
     error("X11 graphics only available on Unix systems.")
   }
   message("Building X11 executable")
   SOURCES += x11graph.cpp
   HEADERS += x11graph.h
   QMAKE_CXXFLAGS_RELEASE += -DX11GRAPHICS
   QMAKE_CXXFLAGS_DEBUG += -DX11GRAPHICS
   #unix {
   #   system(rm $$TARGET.o)
   #}
   CONFIG -= qt
   CONFIG += x11
   unix:LIBS += -lpng
}


#The following line was inserted by qt3to4
QT +=
