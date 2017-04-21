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

draws_q2_AFB:draws_q2_AFB.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_q2_AFB_kll:draws_q2_AFB_kll.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

theory_curve:theory_curve.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

theory_curve_renormalization_scale:theory_curve_renormalization_scale.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

theory_curve_dev:theory_curve_dev.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

theory_curve_error:theory_curve_error.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

signif_prob:signif_prob.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

weighted_mean_asym_error:weighted_mean_asym_error.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

.cpp.o:
	$(CXX) $(ROOTFLAGS) -c $< -o $@

draws_q2_AFB.o                             : draws_.h draws_q2_AFB.cpp                                $(other)
draws_q2_AFB_kll.o                         : draws_.h draws_q2_AFB_kll.cpp                            $(other)
theory_curve.o                             : draws_.h theory_curve.cpp                                $(other)
theory_curve_renormalization_scale.o       : draws_.h theory_curve_renormalization_scale.cpp          $(other)
theory_curve_dev.o                         : draws_.h theory_curve_dev.cpp                            $(other)
theory_curve_error.o                       : draws_.h theory_curve_error.cpp                          $(other)
signif_prob.o                              : draws_.h signif_prob.cpp                                 $(other)
weighted_mean_asym_error.o                 : draws_.h weighted_mean_asym_error.cpp                    $(other)

.PHONY: clean
clean:
	rm $(shell ls *.cpp | sed "s/\.cpp/\.o/g")
distclean:
	rm $(shell ls *.cpp | sed "s/\.cpp/\.o/g") $(shell ls ../Util/*.cpp | sed "s/\.cpp/\.o/g") $(shell ls ../Set/*.cpp | sed "s/\.cpp/\.o/g")
