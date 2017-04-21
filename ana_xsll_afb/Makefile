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

draws_2d_q2_theta_gmc1:draws_2d_q2_theta_gmc1.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_q2_theta_gmc2:draws_2d_q2_theta_gmc2.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_q2_theta_sig:draws_2d_q2_theta_sig.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_q2_theta_rd1:draws_2d_q2_theta_rd1.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_q2_theta_rd2:draws_2d_q2_theta_rd2.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_Mbc:draws_1d_Mbc.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_sep_peak_scf:draws_1d_Mbc_sep_peak_scf.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_body:draws_1d_Mbc_body.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_jpsi:draws_1d_Mbc_jpsi.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_emu:draws_1d_Mbc_emu.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_2d_nb_lep_rd:draws_2d_nb_lep_rd.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_var:draws_1d_var.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_ksfw:draws_1d_ksfw.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_deltaE:draws_1d_deltaE.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_dzll3d:draws_1d_dzll3d.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_bvtxchi2:draws_1d_bvtxchi2.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_bvtxchi2_scan:draws_1d_bvtxchi2_scan.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

.cpp.o:
	$(CXX) $(ROOTFLAGS) -c $< -o $@

draws_2d_q2_theta_gmc1.o            : draws_.h draws_2d_q2_theta_gmc1.cpp                 $(other)
draws_2d_q2_theta_gmc2.o            : draws_.h draws_2d_q2_theta_gmc2.cpp                 $(other)
draws_2d_q2_theta_sig.o             : draws_.h draws_2d_q2_theta_sig.cpp                  $(other)
draws_2d_q2_theta_rd1.o             : draws_.h draws_2d_q2_theta_rd1.cpp                  $(other)
draws_2d_q2_theta_rd2.o             : draws_.h draws_2d_q2_theta_rd2.cpp                  $(other)
draws_1d_Mbc.o                      : draws_.h draws_1d_Mbc.cpp                           $(other)
draws_1d_Mbc_sep_peak_scf.o         : draws_.h draws_1d_Mbc_sep_peak_scf.cpp              $(other)
draws_1d_Mbc_body.o                 : draws_.h draws_1d_Mbc_body.cpp                      $(other)
draws_1d_Mbc_jpsi.o                 : draws_.h draws_1d_Mbc_jpsi.cpp                      $(other)
draws_1d_Mbc_emu.o                  : draws_.h draws_1d_Mbc_emu.cpp                       $(other)
draws_2d_nb_lep_rd.o                : draws_.h draws_2d_nb_lep_rd.cpp                     $(other)
draws_1d_var.o                      : draws_.h draws_1d_var.cpp                           $(other)
draws_1d_ksfw.o                     : draws_.h draws_1d_ksfw.cpp                          $(other)
draws_1d_deltaE.o                   : draws_.h draws_1d_deltaE.cpp                        $(other)
draws_1d_dzll3d.o                   : draws_.h draws_1d_dzll3d.cpp                        $(other)
draws_1d_bvtxchi2.o                 : draws_.h draws_1d_bvtxchi2.cpp                      $(other)
draws_1d_bvtxchi2_scan.o            : draws_.h draws_1d_bvtxchi2_scan.cpp                 $(other)

.PHONY: clean
clean:
	rm $(shell ls *.cpp | sed "s/\.cpp/\.o/g")
distclean:
	rm $(shell ls *.cpp | sed "s/\.cpp/\.o/g") $(shell ls ../Util/*.cpp | sed "s/\.cpp/\.o/g") $(shell ls ../Set/*.cpp | sed "s/\.cpp/\.o/g")
