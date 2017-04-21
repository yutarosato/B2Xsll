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

draws_1d_Mbc_fit_bin_dt_mod_PRL:draws_1d_Mbc_fit_bin_dt_mod_PRL.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_fit_bin_dt_mod_PRL_sub:draws_1d_Mbc_fit_bin_dt_mod_PRL_sub.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_fit_bin_dt_mod_PRL_1_6q2:draws_1d_Mbc_fit_bin_dt_mod_PRL_1_6q2.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_fit_bin_dt_mod_PRL_1_6q2_sub:draws_1d_Mbc_fit_bin_dt_mod_PRL_1_6q2_sub.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_self_cf2:draws_1d_Mbc_self_cf2.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_peak_cc2:draws_1d_Mbc_peak_cc2.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_peak_double2:draws_1d_Mbc_peak_double2.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_peak_swap2:draws_1d_Mbc_peak_swap2.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

.cpp.o:
	$(CXX) $(ROOTFLAGS) -c $< -o $@

draws_1d_Mbc_fit_bin_dt_mod_PRL.o           : draws_.h draws_1d_Mbc_fit_bin_dt_mod_PRL.cpp             $(other)
draws_1d_Mbc_fit_bin_dt_mod_PRL_sub.o       : draws_.h draws_1d_Mbc_fit_bin_dt_mod_PRL_sub.cpp         $(other)
draws_1d_Mbc_fit_bin_dt_mod_PRL_1_6q2.o     : draws_.h draws_1d_Mbc_fit_bin_dt_mod_PRL_1_6q2.cpp       $(other)
draws_1d_Mbc_fit_bin_dt_mod_PRL_1_6q2_sub.o : draws_.h draws_1d_Mbc_fit_bin_dt_mod_PRL_1_6q2_sub.cpp   $(other)
draws_1d_Mbc_self_cf2.o                     : draws_.h draws_1d_Mbc_self_cf2.cpp                       $(other)
draws_1d_Mbc_peak_cc2.o                     : draws_.h draws_1d_Mbc_peak_cc2.cpp                       $(other)
draws_1d_Mbc_peak_double2.o                 : draws_.h draws_1d_Mbc_peak_double2.cpp                   $(other)
draws_1d_Mbc_peak_swap2.o                   : draws_.h draws_1d_Mbc_peak_swap2.cpp                     $(other)

.PHONY: clean
clean:
	rm $(shell ls *.cpp | sed "s/\.cpp/\.o/g")
distclean:
	rm $(shell ls *.cpp | sed "s/\.cpp/\.o/g") $(shell ls ../Util/*.cpp | sed "s/\.cpp/\.o/g") $(shell ls ../Set/*.cpp | sed "s/\.cpp/\.o/g")
