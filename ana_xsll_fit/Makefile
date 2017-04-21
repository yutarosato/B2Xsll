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

draws_1d_Mbc_cc_fit:draws_1d_Mbc_cc_fit.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_fit_linearity:draws_1d_Mbc_fit_linearity.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_fit_linearity_bin:draws_1d_Mbc_fit_linearity_bin.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_fit_linearity_bin2:draws_1d_Mbc_fit_linearity_bin2.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_fit_linearity_bin2_toy:draws_1d_Mbc_fit_linearity_bin2_toy.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_fit_linearity_bin2_ensemble:draws_1d_Mbc_fit_linearity_bin2_ensemble.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_fit_linearity_bin3:draws_1d_Mbc_fit_linearity_bin3.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_fit_linearity_bin3_toy:draws_1d_Mbc_fit_linearity_bin3_toy.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_fit_dt:draws_1d_Mbc_fit_dt.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_fit_bin_dt:draws_1d_Mbc_fit_bin_dt.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_fit_bin_dt_mod:draws_1d_Mbc_fit_bin_dt_mod.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_fit_bin_dt_mod_PRL:draws_1d_Mbc_fit_bin_dt_mod_PRL.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_fit_bin_dt_mod_PRL_1_6q2:draws_1d_Mbc_fit_bin_dt_mod_PRL_1_6q2.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_fit_bin_dt_mod_PRL_plus_single_nu:draws_1d_Mbc_fit_bin_dt_mod_PRL_plus_single_nu.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_fit_bin_dt_totFB:draws_1d_Mbc_fit_bin_dt_totFB.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_fit_bin_dt_kll:draws_1d_Mbc_fit_bin_dt_kll.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_fit_bin_dt_kstrll:draws_1d_Mbc_fit_bin_dt_kstrll.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_fit_bin_dt_xsll:draws_1d_Mbc_fit_bin_dt_xsll.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_fit_bin_dt_syst_cfslope:draws_1d_Mbc_fit_bin_dt_syst_cfslope.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_fit_bin_dt_syst_sigshape:draws_1d_Mbc_fit_bin_dt_syst_sigshape.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_fit_bin_dt_syst_scf:draws_1d_Mbc_fit_bin_dt_syst_scf.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_fit_bin_dt_syst_scf_f:draws_1d_Mbc_fit_bin_dt_syst_scf_f.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_fit_bin_dt_syst_cfwidth:draws_1d_Mbc_fit_bin_dt_syst_cfwidth.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_fit_bin_dt_syst_cfwidth_beta:draws_1d_Mbc_fit_bin_dt_syst_cfwidth_beta.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_fit_bin_dt_syst_peak_charmonium_gmc:draws_1d_Mbc_fit_bin_dt_syst_peak_charmonium_gmc.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_fit_bin_dt_syst_peak_charmonium:draws_1d_Mbc_fit_bin_dt_syst_peak_charmonium.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_fit_bin_dt_syst_peak_double:draws_1d_Mbc_fit_bin_dt_syst_peak_double.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_fit_bin_dt_syst_peak_swap:draws_1d_Mbc_fit_bin_dt_syst_peak_swap.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

linearity_check:linearity_check.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

linearity_check_bin:linearity_check_bin.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

linearity_check_bin2:linearity_check_bin2.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

linearity_check_bin2_toy:linearity_check_bin2_toy.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

linearity_check_bin2_ensemble:linearity_check_bin2_ensemble.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

linearity_check_bin3:linearity_check_bin3.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

linearity_check_bin3_toy:linearity_check_bin3_toy.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

fluctuation_check:fluctuation_check.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_Mbc_fit_shape:draws_1d_Mbc_fit_shape.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_fit_shape_for_rd_emu:draws_1d_Mbc_fit_shape_for_rd_emu.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_self_cf:draws_1d_Mbc_self_cf.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_self_cf2:draws_1d_Mbc_self_cf2.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_self_cf2_totFB:draws_1d_Mbc_self_cf2_totFB.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_self_cf_tot:draws_1d_Mbc_self_cf_tot.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_peak_cc:draws_1d_Mbc_peak_cc.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_peak_cc2:draws_1d_Mbc_peak_cc2.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_peak_cc2_totFB:draws_1d_Mbc_peak_cc2_totFB.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_peak_file_fit:draws_1d_Mbc_peak_file_fit.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_peak_double:draws_1d_Mbc_peak_double.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_peak_double2:draws_1d_Mbc_peak_double2.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_peak_double2_totFB:draws_1d_Mbc_peak_double2_totFB.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_peak_swap:draws_1d_Mbc_peak_swap.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_peak_swap2:draws_1d_Mbc_peak_swap2.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_Mbc_peak_swap2_totFB:draws_1d_Mbc_peak_swap2_totFB.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

draws_1d_M_ll:draws_1d_M_ll.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_M_llpi0:draws_1d_M_llpi0.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_entry_gen:draws_1d_entry_gen.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_entry_rec:draws_1d_entry_rec.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_entry_eff:draws_1d_entry_eff.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_var:draws_1d_var.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_ksfw:draws_1d_ksfw.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

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

draws_1d_Mbc_peak_single_nu2:draws_1d_Mbc_peak_single_nu2.o $(other)
	$(CXX) $(ROOTLIBS_ROOFIT) -o $@ $^

.cpp.o:
	$(CXX) $(ROOTFLAGS) -c $< -o $@

draws_1d_Mbc_cc_fit.o                      : draws_.h draws_1d_Mbc_cc_fit.cpp                         $(other)
draws_1d_Mbc_fit_linearity.o               : draws_.h draws_1d_Mbc_fit_linearity.cpp                  $(other)
draws_1d_Mbc_fit_linearity_bin.o           : draws_.h draws_1d_Mbc_fit_linearity_bin.cpp              $(other)
draws_1d_Mbc_fit_linearity_bin2.o          : draws_.h draws_1d_Mbc_fit_linearity_bin2.cpp             $(other)
draws_1d_Mbc_fit_linearity_bin2_toy.o      : draws_.h draws_1d_Mbc_fit_linearity_bin2_toy.cpp         $(other)
draws_1d_Mbc_fit_linearity_bin2_ensemble.o : draws_.h draws_1d_Mbc_fit_linearity_bin2_ensemble.cpp    $(other)
draws_1d_Mbc_fit_linearity_bin3.o          : draws_.h draws_1d_Mbc_fit_linearity_bin3.cpp             $(other)
draws_1d_Mbc_fit_linearity_bin3_toy.o      : draws_.h draws_1d_Mbc_fit_linearity_bin3_toy.cpp         $(other)
draws_1d_Mbc_fit_dt.o                      : draws_.h draws_1d_Mbc_fit_dt.cpp                         $(other)
draws_1d_Mbc_fit_bin_dt.o                  : draws_.h draws_1d_Mbc_fit_bin_dt.cpp                     $(other)
draws_1d_Mbc_fit_bin_dt_mod.o              : draws_.h draws_1d_Mbc_fit_bin_dt_mod.cpp                 $(other)
draws_1d_Mbc_fit_bin_dt_mod_PRL.o          : draws_.h draws_1d_Mbc_fit_bin_dt_mod_PRL.cpp             $(other)
draws_1d_Mbc_fit_bin_dt_mod_PRL_1_6q2.o    : draws_.h draws_1d_Mbc_fit_bin_dt_mod_PRL_1_6q2.cpp       $(other)
draws_1d_Mbc_fit_bin_dt_mod_PRL_plus_single_nu.o : draws_.h draws_1d_Mbc_fit_bin_dt_mod_PRL_plus_single_nu.cpp $(other)
draws_1d_Mbc_fit_bin_dt_totFB.o            : draws_.h draws_1d_Mbc_fit_bin_dt_totFB.cpp               $(other)
draws_1d_Mbc_fit_bin_dt_kll.o              : draws_.h draws_1d_Mbc_fit_bin_dt_kll.cpp                 $(other)
draws_1d_Mbc_fit_bin_dt_kstrll.o           : draws_.h draws_1d_Mbc_fit_bin_dt_kstrll.cpp              $(other)
draws_1d_Mbc_fit_bin_dt_xsll.o             : draws_.h draws_1d_Mbc_fit_bin_dt_xsll.cpp                $(other)
draws_1d_Mbc_fit_bin_dt_syst_cfslope.o     : draws_.h draws_1d_Mbc_fit_bin_dt_syst_cfslope.cpp        $(other)
draws_1d_Mbc_fit_bin_dt_syst_sigshape.o    : draws_.h draws_1d_Mbc_fit_bin_dt_syst_sigshape.cpp       $(other)
draws_1d_Mbc_fit_bin_dt_syst_scf.o         : draws_.h draws_1d_Mbc_fit_bin_dt_syst_scf.cpp            $(other)
draws_1d_Mbc_fit_bin_dt_syst_scf_f.o       : draws_.h draws_1d_Mbc_fit_bin_dt_syst_scf_f.cpp          $(other)
draws_1d_Mbc_fit_bin_dt_syst_cfwidth.o     : draws_.h draws_1d_Mbc_fit_bin_dt_syst_cfwidth.cpp        $(other)
draws_1d_Mbc_fit_bin_dt_syst_cfwidth_beta.o: draws_.h draws_1d_Mbc_fit_bin_dt_syst_cfwidth_beta.cpp   $(other)
draws_1d_Mbc_fit_bin_dt_syst_peak_charmonium_gmc.o     : draws_.h draws_1d_Mbc_fit_bin_dt_syst_peak_charmonium_gmc.cpp    $(other)
draws_1d_Mbc_fit_bin_dt_syst_peak_charmonium.o         : draws_.h draws_1d_Mbc_fit_bin_dt_syst_peak_charmonium.cpp        $(other)
draws_1d_Mbc_fit_bin_dt_syst_peak_double.o : draws_.h draws_1d_Mbc_fit_bin_dt_syst_peak_double.cpp    $(other)
draws_1d_Mbc_fit_bin_dt_syst_peak_swap.o   : draws_.h draws_1d_Mbc_fit_bin_dt_syst_peak_swap.cpp      $(other)
linearity_check.o                          : draws_.h linearity_check.cpp                             $(other)
linearity_check_bin.o                      : draws_.h linearity_check_bin.cpp                         $(other)
linearity_check_bin2.o                     : draws_.h linearity_check_bin2.cpp                        $(other)
linearity_check_bin2_toy.o                 : draws_.h linearity_check_bin2_toy.cpp                    $(other)
linearity_check_bin2_ensemble.o            : draws_.h linearity_check_bin2_ensemble.cpp               $(other)
linearity_check_bin3.o                     : draws_.h linearity_check_bin3.cpp                        $(other)
linearity_check_bin3_toy.o                 : draws_.h linearity_check_bin3_toy.cpp                    $(other)
fluctuation_check.o                        : draws_.h fluctuation_check.cpp                           $(other)
draws_1d_Mbc_fit_shape.o                   : draws_.h draws_1d_Mbc_fit_shape.cpp                      $(other)
draws_1d_Mbc_fit_shape_for_rd_emu.o        : draws_.h draws_1d_Mbc_fit_shape_for_rd_emu.cpp           $(other)
draws_1d_Mbc_self_cf.o                     : draws_.h draws_1d_Mbc_self_cf.cpp                        $(other)
draws_1d_Mbc_self_cf2.o                    : draws_.h draws_1d_Mbc_self_cf2.cpp                       $(other)
draws_1d_Mbc_self_cf2_totFB.o              : draws_.h draws_1d_Mbc_self_cf2_totFB.cpp                 $(other)
draws_1d_Mbc_self_cf_tot.o                 : draws_.h draws_1d_Mbc_self_cf_tot.cpp                    $(other)
draws_1d_Mbc_peak_cc.o                     : draws_.h draws_1d_Mbc_peak_cc.cpp                        $(other)
draws_1d_Mbc_peak_cc2.o                    : draws_.h draws_1d_Mbc_peak_cc2.cpp                       $(other)
draws_1d_Mbc_peak_cc2_totFB.o              : draws_.h draws_1d_Mbc_peak_cc2_totFB.cpp                 $(other)
draws_1d_Mbc_peak_file_fit.o               : draws_.h draws_1d_Mbc_peak_file_fit.cpp                  $(other)
draws_1d_Mbc_peak_double.o                 : draws_.h draws_1d_Mbc_peak_double.cpp                    $(other)
draws_1d_Mbc_peak_double2.o                : draws_.h draws_1d_Mbc_peak_double2.cpp                   $(other)
draws_1d_Mbc_peak_double2_totFB.o          : draws_.h draws_1d_Mbc_peak_double2_totFB.cpp             $(other)
draws_1d_Mbc_peak_swap.o                   : draws_.h draws_1d_Mbc_peak_swap.cpp                      $(other)
draws_1d_Mbc_peak_swap2.o                  : draws_.h draws_1d_Mbc_peak_swap2.cpp                     $(other)
draws_1d_Mbc_peak_swap2_totFB.o            : draws_.h draws_1d_Mbc_peak_swap2_totFB.cpp               $(other)
draws_1d_M_ll.o                            : draws_.h draws_1d_M_ll.cpp                               $(other)
draws_1d_M_llpi0.o                         : draws_.h draws_1d_M_llpi0.cpp                            $(other)
draws_1d_entry_gen.o                       : draws_.h draws_1d_entry_gen.cpp                          $(other)
draws_1d_entry_rec.o                       : draws_.h draws_1d_entry_rec.cpp                          $(other)
draws_1d_entry_eff.o                       : draws_.h draws_1d_entry_eff.cpp                          $(other)
draws_1d_var.o                             : draws_.h draws_1d_var.cpp                                $(other)
draws_1d_ksfw.o                            : draws_.h draws_1d_ksfw.cpp                               $(other)
draws_q2_AFB.o                             : draws_.h draws_q2_AFB.cpp                                $(other)
draws_q2_AFB_kll.o                         : draws_.h draws_q2_AFB_kll.cpp                            $(other)
theory_curve.o                             : draws_.h theory_curve.cpp                                $(other)
theory_curve_renormalization_scale.o       : draws_.h theory_curve_renormalization_scale.cpp          $(other)
theory_curve_dev.o                         : draws_.h theory_curve_dev.cpp                            $(other)
theory_curve_error.o                       : draws_.h theory_curve_error.cpp                          $(other)
signif_prob.o                              : draws_.h signif_prob.cpp                                 $(other)
weighted_mean_asym_error.o                 : draws_.h weighted_mean_asym_error.cpp                    $(other)
draws_1d_Mbc_peak_single_nu2.o             : draws_.h draws_1d_Mbc_peak_single_nu2.cpp                $(other)

.PHONY: clean
clean:
	rm $(shell ls *.cpp | sed "s/\.cpp/\.o/g")
distclean:
	rm $(shell ls *.cpp | sed "s/\.cpp/\.o/g") $(shell ls ../Util/*.cpp | sed "s/\.cpp/\.o/g") $(shell ls ../Set/*.cpp | sed "s/\.cpp/\.o/g")
