#ifndef Manip_H
#define Manip_H

#include <TMath.h>
#include <TH2D.h>
#include <iostream>
#include <string.h>

// Error is calculated by poisson
TH1D* divide( TH1D* hist1, TH1D* hist2, Char_t* name="", Char_t* title="" ); // [hist1 / hist2]
TH2D* divide( TH2D* hist1, TH2D* hist2, Char_t* name="", Char_t* title="" ); // [hist1 / hist2]

// Error is calculated by binomial
TH1D* divide_binomial( TH1D* hist1, TH1D* hist2, Char_t* name="", Char_t* title="" ); // [hist1 / hist2]
TH2D* divide_binomial( TH2D* hist1, TH2D* hist2, Char_t* name="", Char_t* title="" ); // [hist1 / hist2]


TH1D* multiply( TH1D* hist1, TH1D* hist2, Bool_t fl_wE, Char_t* name="", Char_t* title="" ); // [hist1(data) * hist2(weight,efficiency)]
TH2D* multiply( TH2D* hist1, TH2D* hist2, Bool_t fl_wE, Char_t* name="", Char_t* title="" ); // [hist1(data) * hist2(weight,efficiency)]
// fl_wE : true ( consider stat. and weight error)
//         false( consider    only stat.    error)

#endif
