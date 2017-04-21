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

draws_1d_gen_M_xsll:draws_1d_gen_M_xsll.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_Mbc_mode:draws_1d_Mbc_mode.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_Mbc_Ks_pi0:draws_1d_Mbc_Ks_pi0.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_gen_q2_theta_correction_table:draws_2d_gen_q2_theta_correction_table.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_gen_q2_theta_correction_table_mb465:draws_2d_gen_q2_theta_correction_table_mb465.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_gen_q2_theta_correction_table_mb495:draws_2d_gen_q2_theta_correction_table_mb495.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_gen_q2_theta_correction_table_fm200:draws_2d_gen_q2_theta_correction_table_fm200.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_gen_q2_theta_correction_table_fm480:draws_2d_gen_q2_theta_correction_table_fm480.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_gen_q2_theta_correction_table_lambdaone429:draws_2d_gen_q2_theta_correction_table_lambdaone429.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_gen_q2_theta_correction_table_lambdaone362:draws_2d_gen_q2_theta_correction_table_lambdaone362.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_gen_q2_theta_correction_table_hadronization:draws_2d_gen_q2_theta_correction_table_hadronization.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_gen_q2_theta_correction_table_transition10:draws_2d_gen_q2_theta_correction_table_transition10.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_gen_q2_theta_correction_table_transition12:draws_2d_gen_q2_theta_correction_table_transition12.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_gen_q2_theta_correction_table_fraction_kll_p:draws_2d_gen_q2_theta_correction_table_fraction_kll_p.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_gen_q2_theta_correction_table_fraction_kll_m:draws_2d_gen_q2_theta_correction_table_fraction_kll_m.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_gen_q2_theta_correction_table_fraction_kstrll_p:draws_2d_gen_q2_theta_correction_table_fraction_kstrll_p.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_gen_q2_theta_correction_table_fraction_kstrll_m:draws_2d_gen_q2_theta_correction_table_fraction_kstrll_m.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_gen_q2_theta_correction_table_fraction_xsll_p:draws_2d_gen_q2_theta_correction_table_fraction_xsll_p.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_gen_q2_theta_correction_table_fraction_xsll_m:draws_2d_gen_q2_theta_correction_table_fraction_xsll_m.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_q2_theta_eff:draws_2d_q2_theta_eff.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_q2_theta_eff_transition_10:draws_2d_q2_theta_eff_transition_10.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_q2_theta_eff_transition_12:draws_2d_q2_theta_eff_transition_12.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_q2_theta_eff_hadronization:draws_2d_q2_theta_eff_hadronization.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_q2_theta_eff_fraction_kll_p:draws_2d_q2_theta_eff_fraction_kll_p.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_q2_theta_eff_fraction_kll_m:draws_2d_q2_theta_eff_fraction_kll_m.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_q2_theta_eff_fraction_kstrll_p:draws_2d_q2_theta_eff_fraction_kstrll_p.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_q2_theta_eff_fraction_kstrll_m:draws_2d_q2_theta_eff_fraction_kstrll_m.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_q2_theta_eff_fraction_xsll_p:draws_2d_q2_theta_eff_fraction_xsll_p.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_q2_theta_eff_fraction_xsll_m:draws_2d_q2_theta_eff_fraction_xsll_m.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

make_correction_function:make_correction_function.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

.cpp.o:
	$(CXX) $(ROOTFLAGS) -c $< -o $@

draws_1d_gen_M_xsll.o                     : draws_.h draws_1d_gen_M_xsll.cpp                      $(other)
draws_1d_Mbc_mode.o                       : draws_.h draws_1d_Mbc_mode.cpp                        $(other)
draws_1d_Mbc_Ks_pi0.o                     : draws_.h draws_1d_Mbc_Ks_pi0.cpp                      $(other)
draws_2d_gen_q2_theta_correction_table.o  : draws_.h draws_2d_gen_q2_theta_correction_table.cpp   $(other)
draws_2d_gen_q2_theta_correction_table_mb465.o  : draws_.h draws_2d_gen_q2_theta_correction_table_mb465.cpp   $(other)
draws_2d_gen_q2_theta_correction_table_mb495.o  : draws_.h draws_2d_gen_q2_theta_correction_table_mb495.cpp   $(other)
draws_2d_gen_q2_theta_correction_table_lambdaone429.o  : draws_.h draws_2d_gen_q2_theta_correction_table_lambdaone429.cpp   $(other)
draws_2d_gen_q2_theta_correction_table_lambdaone362.o  : draws_.h draws_2d_gen_q2_theta_correction_table_lambdaone362.cpp   $(other)
draws_2d_gen_q2_theta_correction_table_fm200.o  : draws_.h draws_2d_gen_q2_theta_correction_table_fm200.cpp   $(other)
draws_2d_gen_q2_theta_correction_table_fm480.o  : draws_.h draws_2d_gen_q2_theta_correction_table_fm480.cpp   $(other)
draws_2d_gen_q2_theta_correction_table_hadronization.o  : draws_.h draws_2d_gen_q2_theta_correction_table_hadronization.cpp   $(other)
draws_2d_gen_q2_theta_correction_table_transition10.o  : draws_.h draws_2d_gen_q2_theta_correction_table_transition10.cpp   $(other)
draws_2d_gen_q2_theta_correction_table_transition12.o  : draws_.h draws_2d_gen_q2_theta_correction_table_transition12.cpp   $(other)
draws_2d_gen_q2_theta_correction_table_fraction_kll_p.o  : draws_.h draws_2d_gen_q2_theta_correction_table_fraction_kll_p.cpp   $(other)
draws_2d_gen_q2_theta_correction_table_fraction_kll_m.o  : draws_.h draws_2d_gen_q2_theta_correction_table_fraction_kll_m.cpp   $(other)
draws_2d_gen_q2_theta_correction_table_fraction_kstrll_p.o  : draws_.h draws_2d_gen_q2_theta_correction_table_fraction_kstrll_p.cpp   $(other)
draws_2d_gen_q2_theta_correction_table_fraction_kstrll_m.o  : draws_.h draws_2d_gen_q2_theta_correction_table_fraction_kstrll_m.cpp   $(other)
draws_2d_gen_q2_theta_correction_table_fraction_xsll_p.o  : draws_.h draws_2d_gen_q2_theta_correction_table_fraction_xsll_p.cpp   $(other)
draws_2d_gen_q2_theta_correction_table_fraction_xsll_m.o  : draws_.h draws_2d_gen_q2_theta_correction_table_fraction_xsll_m.cpp   $(other)
draws_2d_q2_theta_eff.o                   : draws_.h draws_2d_q2_theta_eff.cpp                    $(other)
draws_2d_q2_theta_eff_transition_10.o     : draws_.h draws_2d_q2_theta_eff_transition_10.cpp      $(other)
draws_2d_q2_theta_eff_transition_12.o     : draws_.h draws_2d_q2_theta_eff_transition_12.cpp      $(other)
draws_2d_q2_theta_eff_hadronization.o     : draws_.h draws_2d_q2_theta_eff_hadronization.cpp      $(other)
draws_2d_q2_theta_eff_fraction_kll_p.o    : draws_.h draws_2d_q2_theta_eff_fraction_kll_p.cpp     $(other)
draws_2d_q2_theta_eff_fraction_kll_m.o    : draws_.h draws_2d_q2_theta_eff_fraction_kll_m.cpp     $(other)
draws_2d_q2_theta_eff_fraction_kstrll_p.o : draws_.h draws_2d_q2_theta_eff_fraction_kstrll_p.cpp  $(other)
draws_2d_q2_theta_eff_fraction_kstrll_m.o : draws_.h draws_2d_q2_theta_eff_fraction_kstrll_m.cpp  $(other)
draws_2d_q2_theta_eff_fraction_xsll_p.o   : draws_.h draws_2d_q2_theta_eff_fraction_xsll_p.cpp    $(other)
draws_2d_q2_theta_eff_fraction_xsll_m.o   : draws_.h draws_2d_q2_theta_eff_fraction_xsll_m.cpp    $(other)
make_correction_function.o                : draws_.h make_correction_function.cpp                 $(other)

.PHONY: clean
clean:
	rm $(shell ls *.cpp | sed "s/\.cpp/\.o/g")
distclean:
	rm $(shell ls *.cpp | sed "s/\.cpp/\.o/g") $(shell ls ../Util/*.cpp | sed "s/\.cpp/\.o/g") $(shell ls ../Set/*.cpp | sed "s/\.cpp/\.o/g")
