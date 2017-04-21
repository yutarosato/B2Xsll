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

draws_1d_likelihood:draws_1d_likelihood.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_likelihood:draws_2d_likelihood.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_likelihood_proj:draws_2d_likelihood_proj.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_likelihood_lep:draws_1d_likelihood_lep.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_likelihood_lep5:draws_1d_likelihood_lep5.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_likelihood_lep:draws_2d_likelihood_lep.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_likelihood_lep5:draws_2d_likelihood_lep5.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_likelihood_lep_proj:draws_2d_likelihood_lep_proj.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_likelihood_lep5_proj:draws_2d_likelihood_lep5_proj.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_diagonal_deltaE_true_bcs:draws_1d_diagonal_deltaE_true_bcs.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_bgsup:draws_1d_bgsup.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_Mbc:draws_1d_Mbc.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_Mbc_fit:draws_1d_Mbc_fit.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_fit_linearity:draws_1d_Mbc_fit_linearity.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

linearity_check:linearity_check.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_Mbc_tmp:draws_1d_Mbc_tmp.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_Mbc_mode:draws_1d_Mbc_mode.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_Mbc_body:draws_1d_Mbc_body.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_ksfw:draws_1d_ksfw.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_ksfw_lr:draws_1d_ksfw_lr.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_ksfw_var:draws_1d_ksfw_var.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_ksfw_var_mode:draws_1d_ksfw_var_mode.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_fmiss:draws_1d_fmiss.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_sideband:draws_1d_sideband.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

mklikelihood:mklikelihood.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

bcs_lr_sig:bcs_lr_sig.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

bcs_lr_bkg:bcs_lr_bkg.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

bcs_lr_bkg_emu:bcs_lr_bkg_emu.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

bcs_lr_lep_sig:bcs_lr_lep_sig.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

bcs_lr_lep_bkg:bcs_lr_lep_bkg.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

bcs_lr_lep_bkg_cc:bcs_lr_lep_bkg_cc.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

bcs_lr_lep_bkg_emu:bcs_lr_lep_bkg_emu.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

bcs_merge:bcs_merge.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

through_cut_lrnb:through_cut_lrnb.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

through_cut_lrnb_emu:through_cut_lrnb_emu.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

calib:calib.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

cal_delh:cal_delh.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

.cpp.o:
	$(CXX) $(ROOTFLAGS) -c $< -o $@

draws_1d_likelihood.o               : draws_.h draws_1d_likelihood.cpp                    $(other)
draws_2d_likelihood.o               : draws_.h draws_2d_likelihood.cpp                    $(other)
draws_2d_likelihood_proj.o          : draws_.h draws_2d_likelihood_proj.cpp               $(other)
draws_1d_likelihood_lep.o           : draws_.h draws_1d_likelihood_lep.cpp                $(other)
draws_1d_likelihood_lep5.o          : draws_.h draws_1d_likelihood_lep5.cpp               $(other)
draws_2d_likelihood_lep.o           : draws_.h draws_2d_likelihood_lep.cpp                $(other)
draws_2d_likelihood_lep5.o          : draws_.h draws_2d_likelihood_lep5.cpp               $(other)
draws_2d_likelihood_lep_proj.o      : draws_.h draws_2d_likelihood_lep_proj.cpp           $(other)
draws_2d_likelihood_lep5_proj.o     : draws_.h draws_2d_likelihood_lep5_proj.cpp          $(other)
draws_1d_Mbc.o                      : draws_.h draws_1d_Mbc.cpp                           $(other)
draws_1d_Mbc_mode.o                 : draws_.h draws_1d_Mbc_mode.cpp                      $(other)
draws_1d_Mbc_body.o                 : draws_.h draws_1d_Mbc_body.cpp                      $(other)
draws_1d_Mbc_fit.o                  : draws_.h draws_1d_Mbc_fit.cpp                       $(other)
draws_1d_Mbc_fit_linearity.o        : draws_.h draws_1d_Mbc_fit_linearity.cpp             $(other)
linearity_check.o                   : draws_.h linearity_check.cpp                        $(other)
draws_1d_ksfw.o                     : draws_.h draws_1d_ksfw.cpp                          $(other)
draws_1d_ksfw_lr.o                  : draws_.h draws_1d_ksfw_lr.cpp                       $(other)
draws_1d_ksfw_var.o                 : draws_.h draws_1d_ksfw_var.cpp                      $(other)
draws_1d_ksfw_var_mode.o            : draws_.h draws_1d_ksfw_var_mode.cpp                 $(other)
draws_1d_fmiss.o                    : draws_.h draws_1d_fmiss.cpp                         $(other)
draws_1d_sideband.o                 : draws_.h draws_1d_sideband.cpp                      $(other)
draws_1d_diagonal_deltaE_true_bcs.o : draws_.h draws_1d_diagonal_deltaE_true_bcs.cpp      $(other)
draws_1d_bgsup.o                    : draws_.h draws_1d_bgsup.cpp                         $(other)
mklikelihood.o                      : draws_.h mklikelihood.h mklikelihood.cpp            $(other)
bcs_lr_sig.o                        : draws_.h bcs_lr_sig.h bcs_lr_sig.cpp                $(other)
bcs_lr_bkg.o                        : draws_.h bcs_lr_bkg.h bcs_lr_bkg.cpp                $(other)
bcs_lr_lep_sig.o                    : draws_.h bcs_lr_sig.h bcs_lr_lep_sig.cpp            $(other)
bcs_lr_lep_bkg.o                    : draws_.h bcs_lr_bkg.h bcs_lr_lep_bkg.cpp            $(other)
bcs_lr_lep_bkg_cc.o                 : draws_.h bcs_lr_bkg.h bcs_lr_lep_bkg_cc.cpp         $(other)
bcs_merge.o                         : draws_.h bcs_merge.cpp                              $(other)
through_cut_lrnb.o                  : draws_.h through_cut_lrnb.cpp                       $(other)
through_cut_lrnb_emu.o              : draws_.h through_cut_lrnb_emu.cpp                   $(other)
calib.o                             : draws_.h calib.cpp                                  $(other)
cal_delh.o                          : draws_.h cal_delh.cpp                               $(other)

.PHONY: clean
clean:
	rm $(shell ls *.cpp | sed "s/\.cpp/\.o/g")
distclean:
	rm $(shell ls *.cpp | sed "s/\.cpp/\.o/g") $(shell ls ../Util/*.cpp | sed "s/\.cpp/\.o/g") $(shell ls ../Set/*.cpp | sed "s/\.cpp/\.o/g")
