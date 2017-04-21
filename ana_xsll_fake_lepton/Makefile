CXX = g++

ROOTFLAGS = `root-config --cflags`
ROOTLIBS  = `root-config --libs`
ROOTLIBS_ROOFIT = `root-config --libs`
ROOTLIBS_ROOFIT += -lRooFit -lRooFitCore -lMinuit 


all: util set $(shell ls *.cpp | sed "s/\.cpp//g")
other = $(shell ls ../Util/*.cpp | sed "s/\.cpp/.o/g") $(shell ls ../Set_fake_lepton/*.cpp | sed "s/\.cpp/.o/g")

util:
	cd ../Util; make
set:
	cd ../Set_fake_lepton; make

through_double_cut2:through_double_cut2.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

through_single_cut2:through_single_cut2.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

merge_double_lrnb_cut2:merge_double_lrnb_cut2.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

merge_single_lrnb_cut2:merge_single_lrnb_cut2.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

mk_weight_double:mk_weight_double.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

mk_weight_double_data:mk_weight_double_data.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

mk_weight_single:mk_weight_single.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

mk_weight_single_jpsi:mk_weight_single_jpsi.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

mk_weight_single_jpsi_data:mk_weight_single_jpsi_data.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_Mbc_bb_peak1:draws_1d_Mbc_bb_peak1.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_Mbc_bb_peak2:draws_1d_Mbc_bb_peak2.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_peak_q2_theta:draws_2d_peak_q2_theta.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

invert_coslp:invert_coslp.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

ck_pi_fake:ck_pi_fake.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

.cpp.o:
	$(CXX) $(ROOTFLAGS) -c $< -o $@

through_double_cut2.o               : draws_.h through_double_cut2.cpp                    $(other)
through_single_cut2.o               : draws_.h through_single_cut2.cpp                    $(other)
merge_double_lrnb_cut2.o            : draws_.h merge_double_lrnb_cut2.cpp                 $(other)
merge_single_lrnb_cut2.o            : draws_.h merge_single_lrnb_cut2.cpp                 $(other)
mk_weight_double.o                  : draws_.h mk_weight_double.cpp                       $(other)
mk_weight_double_data.o             : draws_.h mk_weight_double_data.cpp                  $(other)
mk_weight_single.o                  : draws_.h mk_weight_single.cpp                       $(other)
mk_weight_single_jpsi.o             : draws_.h mk_weight_single_jpsi.cpp                  $(other)
mk_weight_single_jpsi_data.o        : draws_.h mk_weight_single_jpsi_data.cpp             $(other)
draws_1d_Mbc_bb_peak1.o             : draws_.h draws_1d_Mbc_bb_peak1.cpp                  $(other)
draws_1d_Mbc_bb_peak2.o             : draws_.h draws_1d_Mbc_bb_peak2.cpp                  $(other)
draws_2d_peak_q2_theta.o            : draws_.h draws_2d_peak_q2_theta.cpp                 $(other)
invert_coslp.o                      : draws_.h invert_coslp.cpp                           $(other)
ck_pi_fake.o                        : draws_.h ck_pi_fake.cpp                             $(other)

.PHONY: clean
clean:
	rm $(shell ls *.cpp | sed "s/\.cpp/\.o/g")
distclean:
	rm $(shell ls *.cpp | sed "s/\.cpp/\.o/g") $(shell ls ../Util/*.cpp | sed "s/\.cpp/\.o/g") $(shell ls ../Set_fake_lepton/*.cpp | sed "s/\.cpp/\.o/g")
