CXX = g++

ROOTFLAGS = `root-config --cflags`
ROOTLIBS  = `root-config --libs`

all: $(shell ls *.cpp | sed "s/\.cpp/\.o/g")

.cpp.o:
	$(CXX) $(ROOTFLAGS) -c $< -o $@
Style.o       : Style.h      Style.cpp
Canvas.o      : Canvas.h     Canvas.cpp
MChain.o      : MChain.h     MChain.cpp     MCut_array.h MCut_array.cpp MCut.h MCut.cpp
MCut.o        : MCut.h       MCut.cpp
MCut_array.o  : MCut_array.h MCut_array.cpp MCut.h     MCut.cpp
FitFunc.o     : FitFunc.h    FitFunc.cpp
Stat.o        : Stat.h       Stat.cpp
Manip.o       : Manip.h      Manip.cpp

.PHONY: clean
clean:
	rm $(shell ls *.cpp | sed "s/\.cpp/\.o/g")
