CXX = g++

ROOTFLAGS = `root-config --cflags`
ROOTLIBS  = `root-config --libs`
ROOTLIBS_ROOFIT = `root-config --libs`
ROOTLIBS_ROOFIT += -lRooFit -lRooFitCore -lMinuit 

all: util set $(shell ls *.cpp | sed "s/\.cpp//g")
other = $(shell ls ../Util/*.cpp | sed "s/\.cpp/.o/g") $(shell ls ../Set_dstr_fl/*.cpp | sed "s/\.cpp/.o/g")

util:
	cd ../Util; make
set:
	cd ../Set_dstr_fl; make

draws_1d_var:draws_1d_var.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_pid:draws_1d_pid.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_d0_plot:draws_1d_d0_plot.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_d0_plot_data:draws_1d_d0_plot_data.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_d0_fit:draws_1d_d0_fit.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_1d_d0_fit_data:draws_1d_d0_fit_data.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_2d_p_cos:draws_2d_p_cos.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

draws_fake_rate:draws_fake_rate.o $(other)
	$(CXX) $(ROOTLIBS) -o $@ $^

.cpp.o:
	$(CXX) $(ROOTFLAGS) -c $< -o $@

draws_1d_var.o            : draws_.h draws_1d_var.cpp                $(other)
draws_1d_pid.o            : draws_.h draws_1d_pid.cpp                $(other)
draws_1d_d0_plot.o        : draws_.h draws_1d_d0_plot.cpp            $(other)
draws_1d_d0_plot_data.o   : draws_.h draws_1d_d0_plot_data.cpp       $(other)
draws_1d_d0_fit.o         : draws_.h draws_1d_d0_fit.cpp             $(other)
draws_1d_d0_fit_data.o    : draws_.h draws_1d_d0_fit_data.cpp        $(other)
draws_2d_p_cos.o          : draws_.h draws_2d_p_cos.cpp              $(other)
draws_fake_rate.o         : draws_.h draws_fake_rate.cpp             $(other)

.PHONY: clean
clean:
	rm $(shell ls *.cpp | sed "s/\.cpp/\.o/g")
distclean:
	rm $(shell ls *.cpp | sed "s/\.cpp/\.o/g") $(shell ls ../Util/*.cpp | sed "s/\.cpp/\.o/g") $(shell ls ../Set_dstr_fl/*.cpp | sed "s/\.cpp/\.o/g")
