CXX = g++

ROOTFLAGS = `root-config --cflags`
ROOTLIBS  = `root-config --libs`

all: $(shell ls *.cpp | sed "s/\.cpp/\.o/g")
other = $(shell ls ..//Util/*.o)

.cpp.o:
	$(CXX) $(ROOTFLAGS) -c $< -o $@
Branch.o                : Branch.h                 Branch.cpp                $(other)
Nominal_cut_selection.o : Nominal_cut_selection.h  Nominal_cut_selection.cpp $(other)
makeCut.o               : makeCut.h                makeCut.cpp               $(other)

.PHONY: clean
clean:
	rm $(shell ls *.cpp | sed "s/\.cpp/\.o/g")
