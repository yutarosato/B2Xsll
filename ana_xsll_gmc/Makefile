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

draws_1d_Mbc_type:draws_1d_Mbc_type.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_Mbc_bb_lepself1:draws_1d_Mbc_bb_lepself1.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_Mbc_bb_lepself2:draws_1d_Mbc_bb_lepself2.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_Mbc_bb_smid1:draws_1d_Mbc_bb_smid1.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_Mbc_bb_smid2:draws_1d_Mbc_bb_smid2.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_Mbc_bb_peak1:draws_1d_Mbc_bb_peak1.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_Mbc_bb_peak2:draws_1d_Mbc_bb_peak2.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_Mbc_bb_peak_tf1:draws_1d_Mbc_bb_peak_tf1.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_Mbc_bb_peak_tf2:draws_1d_Mbc_bb_peak_tf2.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_Mbc_bb_rest_peak2:draws_1d_Mbc_bb_rest_peak2.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_Mbc_bb_dmid1:draws_1d_Mbc_bb_dmid1.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_Mbc_bb_dmid2:draws_1d_Mbc_bb_dmid2.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_Mbc_cc_fit:draws_1d_Mbc_cc_fit.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_nopeak_fit:draws_1d_Mbc_nopeak_fit.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_emu_fit:draws_1d_Mbc_emu_fit.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_entry_exp:draws_1d_entry_exp.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_bgsup:draws_1d_bgsup.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_multi-axis_Fmiss:draws_1d_multi-axis_Fmiss.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

merge_cut:merge_cut.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

merge_cut2:merge_cut2.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

merge_cut3:merge_cut3.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

merge_cut4:merge_cut4.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

merge_cut5:merge_cut5.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

merge_cut7:merge_cut7.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

merge_cut8:merge_cut8.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_diagonal_deltaE_true_bcs:draws_1d_diagonal_deltaE_true_bcs.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_Mbc_shape:draws_1d_Mbc_shape.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^



.cpp.o:
	$(CXX) $(ROOTFLAGS) -c $< -o $@

draws_1d_Mbc_type.o                 : draws_.h draws_1d_Mbc_type.cpp                      $(other)
draws_1d_Mbc_bb_lepself1.o          : draws_.h draws_1d_Mbc_bb_lepself1.cpp               $(other)
draws_1d_Mbc_bb_lepself2.o          : draws_.h draws_1d_Mbc_bb_lepself2.cpp               $(other)
draws_1d_Mbc_bb_smid1.o             : draws_.h draws_1d_Mbc_bb_smid1.cpp                  $(other)
draws_1d_Mbc_bb_smid2.o             : draws_.h draws_1d_Mbc_bb_smid2.cpp                  $(other)
draws_1d_Mbc_bb_peak1.o             : draws_.h draws_1d_Mbc_bb_peak1.cpp                  $(other)
draws_1d_Mbc_bb_peak2.o             : draws_.h draws_1d_Mbc_bb_peak2.cpp                  $(other)
draws_1d_Mbc_bb_peak_tf1.o          : draws_.h draws_1d_Mbc_bb_peak_tf1.cpp               $(other)
draws_1d_Mbc_bb_peak_tf2.o          : draws_.h draws_1d_Mbc_bb_peak_tf2.cpp               $(other)
draws_1d_Mbc_bb_rest_peak2.o        : draws_.h draws_1d_Mbc_bb_rest_peak2.cpp             $(other)
draws_1d_Mbc_bb_dmid1.o             : draws_.h draws_1d_Mbc_bb_dmid1.cpp                  $(other)
draws_1d_Mbc_bb_dmid2.o             : draws_.h draws_1d_Mbc_bb_dmid2.cpp                  $(other)
draws_1d_Mbc_cc_fit.o               : draws_.h draws_1d_Mbc_cc_fit.cpp                    $(other)
draws_1d_Mbc_nopeak_fit.o           : draws_.h draws_1d_Mbc_nopeak_fit.cpp                $(other)
draws_1d_Mbc_emu_fit.o              : draws_.h draws_1d_Mbc_emu_fit.cpp                   $(other)
draws_1d_entry_exp.o                : draws_.h draws_1d_entry_exp.cpp                     $(other)
draws_1d_bgsup.o                    : draws_.h draws_1d_bgsup.cpp                         $(other)
draws_1d_multi-axis_Fmiss.o         : draws_.h draws_1d_multi-axis_Fmiss.cpp              $(other)
merge_cut.o                         : draws_.h merge_cut.cpp                              $(other)
merge_cut2.o                        : draws_.h merge_cut2.cpp                             $(other)
merge_cut3.o                        : draws_.h merge_cut3.cpp                             $(other)
merge_cut4.o                        : draws_.h merge_cut4.cpp                             $(other)
merge_cut5.o                        : draws_.h merge_cut5.cpp                             $(other)
merge_cut7.o                        : draws_.h merge_cut7.cpp                             $(other)
merge_cut8.o                        : draws_.h merge_cut8.cpp                             $(other)
draws_1d_diagonal_deltaE_true_bcs.o : draws_.h draws_1d_diagonal_deltaE_true_bcs.cpp      $(other)
draws_1d_Mbc_shape.o                : draws_.h draws_1d_Mbc_shape.cpp                     $(other)

.PHONY: clean
clean:
	rm $(shell ls *.cpp | sed "s/\.cpp/\.o/g")
distclean:
	rm $(shell ls *.cpp | sed "s/\.cpp/\.o/g") $(shell ls ../Util/*.cpp | sed "s/\.cpp/\.o/g") $(shell ls ../Set/*.cpp | sed "s/\.cpp/\.o/g")
