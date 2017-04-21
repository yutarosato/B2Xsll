#ifndef FitFunc_H
#define FitFunc_H

#include <TROOT.h>
#include <TF1.h>
#include <TH1D.h>
#include <TMath.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>

#include "const.h"

typedef Double_t ( *fit_func )( const Double_t*, const Double_t* );
typedef void     ( *fit_iter )( TF1* func );
static Double_t window_width = 0.08; //[MeV]   for gauss+straight_area (case6)
static Double_t window_mean  = PDGmass::b0;//    for gauss+straight_area (case6)

Int_t iterative_fit( TH1D* hist, TF1* func, Int_t nparam, Double_t* init_var, Char_t* opt_plot,                      Int_t max    );
Int_t iterative_fit( TH1D* hist, TF1* func, Int_t nparam, Double_t* init_var, Char_t* opt_plot, fit_iter manip=NULL, Int_t max=20 );


Int_t n_fitfunc_par( Int_t select );

// < make fit function >
fit_func make_func( Int_t select );

// < parameter setting >
void func_set_parameters( Int_t select, TF1* func,   TH1D* hist,
			  Int_t xbin,   Double_t xmin, Double_t xmax );

// < calculation of signal yields>
// < calculate integral >
void func_get_integral( Int_t select, TF1* func,
			Double_t bin_width,
			Double_t& yield, Double_t& yieldE );
// *************************************************************************************::

// 0. horizontal line
// p[0]=slope, p[1]=offset
Double_t func_horizontal(const Double_t* x, const Double_t* p);

// 1. straight line
// p[0]=slope, p[1]=offset
Double_t func_straight( const Double_t* x, const Double_t* p );


// 2. parabola
// p[0]=ax^2, p[1]=bx, p[2]=c
Double_t func_parabola( const Double_t* x, const Double_t* p );

// 3. cubic
// p[0]=ax^3, p[1]=bx^2, p[2]=cx, p[3]=d
Double_t func_cubic( const Double_t* x, const Double_t* p );

// 4. fourth
// p[0]=ax^4, p[1]=bx^3, p[2]=cx^2, p[3]=dx, p[4]=e
Double_t func_fourth( const Double_t* x, const Double_t* p );

// 5. fourth
// p[0]=ax^5, p[1]=bx^4, p[2]=cx^3, p[3]=dx^2, p[4]=ex, p[5]=f
Double_t func_fifth( const Double_t* x, const Double_t* p );

// 6. straight line with area parametrization
// p[0]=area, p[1]=offset,
Double_t func_straight_area( const Double_t* x, const Double_t* p );

// 10. gaussian
// p[0]=area, p[1]=mean, p[2]=sigma,
Double_t func_gauss( const Double_t* x, const Double_t* p );

// 11. gaussian(p[0-2]) + straight(p[3-4])
Double_t func_gauss_straight( const Double_t* x, const Double_t* p );

// 12. gaussian(p[0-2]) + parabola(p[3-5])
Double_t func_gauss_parabola( const Double_t* x, const Double_t* p );

// 15. gaussian(p[0-2]) + argus(p[3-5])
Double_t func_gauss_argus( const Double_t* x, const Double_t* p );
// 151. gaussian(p[0-2]) + modified-argus(p[3-6])
Double_t func_gauss_modargus( const Double_t* x, const Double_t* p );
// 152. gaussian(p[0-2]) + modified-argus(p[3-5])
Double_t func_gauss_modargus2( const Double_t* x, const Double_t* p );
// 153. gaussian(p[0-2]) + modified-argus(p[3-5])
Double_t func_gauss_modargus3( const Double_t* x, const Double_t* p );

// 16. gaussian(p[0-2]) + straight_area(p[3-4])
Double_t func_gauss_straight_area( const Double_t* x, const Double_t* p );

// 20. double gaussian
// p[0]=area,  p[1]=area_ratio,
// p[2]=mean,
// p[3]=sigma, p[4]=sigma_ratio
Double_t func_2gauss( const Double_t* x, const Double_t* p );

// 21. double gaussian(p[0-4]) + straight(p[5-6])
Double_t func_2gauss_straight( const Double_t* x, const Double_t* p );

// 22. double gaussian(p[0-4]) + parabola(p[5-7])
Double_t func_2gauss_parabola( const Double_t* x, const Double_t* p );

// 25. triple gaussian
// p[0]=area,  p[1]=area_ratio1,   p[2]=area_ratio2,
// p[3]=mean,
// p[4]=sigma, p[5]=sigma_ratio1, p[6]=sigma_ratio2
Double_t func_3gauss( const Double_t* x, const Double_t* p );

// 26. triple gaussian(p[0-4]) + straight(p[5-6])
Double_t func_3gauss_straight( const Double_t* x, const Double_t* p );

// 27. triple gaussian(p[0-4]) + parabola(p[5-7])
Double_t func_3gauss_parabola( const Double_t* x, const Double_t* p );

// 100. bif-gaussian
// p[0]=area, p[1]=mean, p[2]=sigma(low), p[3]=sigma(high)
Double_t func_gauss_bif( const Double_t* x, const Double_t* p );

// 110. bif-gaussian + gaussian
// p[0]=area, p[1]=mean, p[2]=sigma(low), p[3]=sigma(high)
Double_t func_gauss_bif_gauss( const Double_t* x, const Double_t* p );

// 111. bif-gaussian + gaussian + straight
// p[0]=area, p[1]=mean, p[2]=sigma(low), p[3]=sigma(high)
Double_t func_gauss_bif_gauss_straight( const Double_t* x, const Double_t* p );

// 50. argus
Double_t func_argus    (const Double_t* x, const Double_t* p);
Double_t func_roo_argus(const Double_t  x, const Double_t p0, const Double_t p1);

// 51. modified argus
Double_t func_modargus    (const Double_t* x, const Double_t* p);
Double_t func_roo_modargus(const Double_t  x, const Double_t p0, const Double_t p1, const Double_t p2);
//Double_t func_roo_modargus(const Double_t  x, const Double_t p1, const Double_t p2);

// 52. modified argus2
Double_t func_modargus2    (const Double_t* x, const Double_t* p);
Double_t func_roo_modargus2(const Double_t  x, const Double_t p0, const Double_t p1);

// 53. modified argus3
Double_t func_modargus3    (const Double_t* x, const Double_t* p);
Double_t func_roo_modargus3(const Double_t  x, const Double_t p0, const Double_t p1);

// 310. cb function (low-side tail)
Double_t func_cb_low(const Double_t* x, const Double_t*p);

// 311. cb function (low-side tail)(p[0-5]) + straight(p[6-7])
Double_t func_cb_low_straight(const Double_t* x, const Double_t*p);

// 410. cb function (low-side tail,bif)
Double_t func_cb_low_bif(const Double_t* x, const Double_t*p);

// 411. cb function (low-side tail,bif)(p[0-5]) + straight(p[6-7])
Double_t func_cb_low_bif_straight(const Double_t* x, const Double_t*p);

// 20000. cb functione
Double_t func_cb(const Double_t* x, const Double_t*p);

// 20001. cb functione(p[0-7]) + straight(p[8-9])
Double_t func_cb_straight(const Double_t* x, const Double_t*p);

// 30000. cb functione (bifurcated gaussian)
Double_t func_cb_bif(const Double_t* x, const Double_t*p);

// 30001. cb functione (bifurcated gaussian)(p[0-7]) + straight(p[8-9])
Double_t func_cb_bif_straight(const Double_t* x, const Double_t*p);

// 1091
Double_t func_chi( const Double_t *x, const Double_t *p );

// 1092
Double_t func_chi2( const Double_t *x, const Double_t *p );

// 81
Double_t func_threshold1    ( const Double_t *x, const Double_t *p );
Double_t func_roo_threshold1( const Double_t  x, const Double_t p0, const Double_t p1 );

// 82
Double_t func_threshold2    ( const Double_t *x, const Double_t *p );
Double_t func_roo_threshold2( const Double_t  x, const Double_t p0, const Double_t p1, const Double_t p2 );

// 181
Double_t func_gauss_threshold1( const Double_t *x, const Double_t *p );

// 182
Double_t func_gauss_threshold2( const Double_t *x, const Double_t *p );

#endif
