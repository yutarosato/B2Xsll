CXX = g++

ROOTFLAGS = `root-config --cflags`
ROOTLIBS  = `root-config --libs`
ROOTLIBS_ROOFIT = `root-config --libs`
ROOTLIBS_ROOFIT += -lRooFit -lRooFitCore -lMinuit 

all: util set $(shell ls *.cpp | sed "s/\.cpp//g")
other = $(shell ls ../Util/*.cpp | sed "s/\.cpp/.o/g") $(shell ls ../Set/*.cpp | sed "s/\.cpp/.o/g")

util:
	cd ../Util; make
set:
	cd ../Set; make

draws_2d_gen_q2_theta_correction_show:draws_2d_gen_q2_theta_correction_show.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_gen_q2_theta_correction_table:draws_2d_gen_q2_theta_correction_table.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_q2_theta_eff:draws_2d_q2_theta_eff.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

make_correction_function:make_correction_function.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

.cpp.o:
	$(CXX) $(ROOTFLAGS) -c $< -o $@
draws_2d_gen_q2_theta_correction_show.o  : draws_.h draws_2d_gen_q2_theta_correction_show.cpp   $(other)
draws_2d_gen_q2_theta_correction_table.o : draws_.h draws_2d_gen_q2_theta_correction_table.cpp  $(other)
draws_2d_q2_theta_eff.o                  : draws_.h draws_2d_q2_theta_eff.cpp                   $(other)
make_correction_function.o               : draws_.h make_correction_function.cpp                $(other)

.PHONY: clean
clean:
	rm $(shell ls *.cpp | sed "s/\.cpp/\.o/g")
distclean:
	rm $(shell ls *.cpp | sed "s/\.cpp/\.o/g") $(shell ls ../Util/*.cpp | sed "s/\.cpp/\.o/g") $(shell ls ../Set/*.cpp | sed "s/\.cpp/\.o/g")
