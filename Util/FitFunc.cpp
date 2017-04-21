#include "FitFunc.h"
#include "const.h"

Int_t iterative_fit( TH1D* hist, TF1* func, Int_t nparam, Double_t* init, Char_t* opt_plot, Int_t max ){
  return iterative_fit( hist, func, nparam, init, opt_plot, NULL, max );
}

Int_t iterative_fit( TH1D* hist, TF1* func, Int_t nparam, Double_t* init, Char_t* opt_plot, fit_iter manip, Int_t max ){
  if( manip!=NULL ) manip(func);
  //TFitResultPtr fit_result = hist->Fit(func, "RQSU", opt_plot );
  TFitResultPtr fit_result = hist->Fit(func, "RSU", opt_plot );
  Int_t fit_status = fit_result->Status(); // 0(success), other(false), 4(call limit)
  Int_t cnt_n =  0; // counting variable for iteration
  Int_t cnt_v =  0; // counting variable for fitting parameter
  if( fit_status==0 ){
    fit_result->Print();
    //fit_result->Print("V"); // including covariance matrix
    std::cout << "fit_status = " << fit_status << "(" << (Int_t)fit_result->IsValid() << ")" << std::endl
	      << "FCN        = " << fit_result->MinFcnValue() << std::endl
	      << "chi2/NDF   = " << func->GetChisquare()
	      << " / "           << func->GetNDF()
	      << " = "           << func->GetChisquare()/func->GetNDF()
	      << std::endl;
    return fit_status;
  }
  
  while( fit_status!=0 || cnt_v!=0 ){
    std::cout << "[ Fitting false : " << cnt_n << " : " << cnt_v << " ]" << std::endl;
    
    for(Int_t i=0; i<nparam; i++ ){
      func->ReleaseParameter(i);
      if     ( cnt_v<0  ) func->SetParameter(i, init[i]);
      else if( i!=cnt_v ) func->FixParameter(i, init[i]);
      else                func->SetParameter(i, init[i]);
    }
    if( manip!=NULL ) manip(func);
    
    if( cnt_v<0 ) fit_result = hist->Fit( func,"RSU",  opt_plot );
    else          fit_result = hist->Fit( func,"RQSU", opt_plot );
    fit_status = fit_result->Status();
    
    if( cnt_n%nparam != 0 ){
      for( Int_t i=0; i<nparam; i++ ) init[i] = func->GetParameter(i);
    }
    
    if( cnt_n > max*nparam ){
      std::cout << "[ERROR] Iterative Fitting is false !!!" << std::endl;
      std::cerr << "[ERROR] Iterative Fitting is false !!!" << std::endl;
      return -1;
    }
    
    cnt_n++;
    cnt_v++;
    
    if( cnt_v == nparam ) cnt_v = -1;
  }
  
  fit_result->Print();
  //fit_result->Print("V"); // including covariance matrix
  std::cout << "fit_status = " << fit_status << "(" << (Int_t)fit_result->IsValid() << ")" << std::endl
	    << "FCN        = " << fit_result->MinFcnValue() << std::endl
	    << "chi2/NDF   = " << func->GetChisquare()
	    << " / "           << func->GetNDF()
	    << " = "           << func->GetChisquare()/func->GetNDF()
	    << std::endl;
  return fit_status;
}

Int_t n_fitfunc_par( Int_t select ){
  switch( select ){
  case 0:
    return 1;
    break;
  case 1:
    return 2;
    break;
  case    2:
  case   10:
  case   50:
  case   52:
  case   53:
  case   81:
  case 1091:
    return 3;
    break;
  case   3:
  case   6:
  case  51:
  case  82:
  case 100:
    return 4;
    break;
  case   4:
  case  11:
  case  16:
  case  20:
  case 310:
    return 5;
    break;
  case    5:
  case   12:
  case   15:
  case  152:
  case  153:
  case  110:
  case  181:
  case  410:
  case 1092:
    return 6;
    break;
  case    21:
  case    30:
  case   151:
  case   182:
  case   311:
  case 20000:
    return 7;
    break;
  case    22:
  case   111:
  case   411:
  case 30000:
    return 8;
    break;
  case    31:
  case 20001:
    return 9;
    break;
  case    32:
  case 30001:
    return 10;
    break;
  default:
    std::cerr << "indefinite fit function -> abort()" << std::endl;
    abort();
    break;
  }
}

// < make fit function >
fit_func make_func( Int_t select ){
  fit_func func;
  switch( select ){
  case 0:
    func = func_horizontal;
    break;
  case 1:
    func = func_straight;
    break;
  case 2:
    func = func_parabola;
    break;
  case 3:
    func = func_cubic;
    break;
  case 4:
    func = func_fourth;
    break;
  case 5:
    func = func_fifth;
    break;
  case 10:
    func = func_gauss;
    break;
  case 11:
    func = func_gauss_straight;
    break;
  case 12:
    func = func_gauss_parabola;
    break;
  case 15:
    func = func_gauss_argus;
    break;
  case 151:
    func = func_gauss_modargus;
    break;
  case 152:
    func = func_gauss_modargus2;
    break;
  case 153:
    func = func_gauss_modargus3;
    break;
  case 16:
    func = func_gauss_straight_area;
    break;
  case 20:
    func = func_2gauss;
    break;
  case 21:
    func = func_2gauss_straight;
    break;
  case 22:
    func = func_2gauss_parabola;
    break;
  case 30:
    func = func_3gauss;
    break;
  case 31:
    func = func_3gauss_straight;
    break;
  case 32:
    func = func_3gauss_parabola;
    break;
  case 81:
    func = func_threshold1;
    break;
  case 82:
    func = func_threshold2;
    break;
  case 181:
    func = func_gauss_threshold1;
    break;
  case 182:
    func = func_gauss_threshold2;
    break;
  case 100:
    func = func_gauss_bif;
    break;
  case 110:
    func = func_gauss_bif_gauss;
    break;
  case 111:
    func = func_gauss_bif_gauss_straight;
    break;
  case  50:
    func = func_argus;
    break;
  case  51:
    func = func_modargus;
    break;
  case  52:
    func = func_modargus2;
    break;
  case  53:
    func = func_modargus3;
    break;
  case 1091:
    func = func_chi;
    break;
  case 1092:
    func = func_chi2;
    break;
  case 310:
    func = func_cb_low;
    break;
  case 311:
    func = func_cb_low_straight;
    break;
  case 410:
    func = func_cb_low_bif;
    break;
  case 411:
    func = func_cb_low_bif_straight;
    break;
  case 20000:
    func = func_cb;
    break;
  case 20001:
    func = func_cb_straight;
    break;
  case 30000:
    func = func_cb_bif;
    break;
  case 30001:
    func = func_cb_bif_straight;
    break;
  default:
    std::cerr << Form("Check fit function(%d)",select) << std::endl;
    abort();
  }
  return func;
}


// < parameter setting >
void func_set_parameters( Int_t select, TF1* func,   TH1D* hist,
			  Int_t xbin,   Double_t xmin, Double_t xmax ){

  Double_t a,b,c,d,e,f,sigma,mu,area,mass,gamma,norm,mu1,mu2,sigma1,sigma2,area1,area2,x1,x2,x3;
  switch( select ){
  case 0 : // horizontal line
    a = hist->GetBinContent(1);
    func->SetParNames( "height" );
    func->SetParameter( 0, a );
    break;
  case 1 : // straight
    a = ( hist->GetBinContent(xbin)-hist->GetBinContent(1))/(hist->GetBinCenter(xbin)-hist->GetBinCenter(1) );
    b = hist->GetBinContent(1)-a*hist->GetBinCenter(1);
    func->SetParNames( "slope","offset" );
    func->SetParameters( a,b );
    break;
  case 2: // parabola
    a     = 0;
    b     = ( hist->GetBinContent(xbin)-hist->GetBinContent(1))/(hist->GetBinCenter(xbin)-hist->GetBinCenter(1) );
    c     = hist->GetBinContent(1)-b*hist->GetBinCenter(1);
    func->SetParNames( "a","b","c" );
    func->SetParameters( a,b,c );
    break;
  case 3: // cubic
    a     = 0;
    b     = 0;
    c     = ( hist->GetBinContent(xbin)-hist->GetBinContent(1))/(hist->GetBinCenter(xbin)-hist->GetBinCenter(1) );
    d     = hist->GetBinContent(1)-b*hist->GetBinCenter(1);
    func->SetParNames  ( "a","b","c","d" );
    func->SetParameters(  a,  b,  c,  d  );
    break;
  case 4: // fourth
    a     = 0;
    b     = 0;
    c     = 0;
    d     = ( hist->GetBinContent(xbin)-hist->GetBinContent(1))/(hist->GetBinCenter(xbin)-hist->GetBinCenter(1) );
    e     = hist->GetBinContent(1)-b*hist->GetBinCenter(1);
    func->SetParNames  ( "a","b","c","d","e" );
    func->SetParameters(  a,  b,  c,  d,  e  );
    break;
  case 5: // fifth
    a     = 0;
    b     = 0;
    c     = 0;
    d     = 0;
    e     = ( hist->GetBinContent(xbin)-hist->GetBinContent(1))/(hist->GetBinCenter(xbin)-hist->GetBinCenter(1) );
    f     = hist->GetBinContent(1)-b*hist->GetBinCenter(1);
    func->SetParNames  ( "a","b","c","d","e","f" );
    func->SetParameters(  a,  b,  c,  d,  e,  f  );
    break;
  case 10: // gaussian
    sigma = 2*hist->GetBinWidth(1);
    mu    = hist->GetBinCenter( hist->GetMaximumBin() );
    area  = hist->GetMaximum()*sqrt( TMath::TwoPi() )*sigma;
    func->SetParNames( "area","mean","sigma" );
    func->SetParameters( area,mu,sigma );
    break;
  case 11: // gaussian + straight
    a     = ( hist->GetBinContent(xbin)-hist->GetBinContent(1))/(hist->GetBinCenter(xbin)-hist->GetBinCenter(1) );
    b     = hist->GetBinContent(1)-a*hist->GetBinCenter(1);
    mu    = hist->GetBinCenter(hist->GetMaximumBin());
    sigma = 2.5*hist->GetBinWidth(1);
    area  = ( hist->GetMaximum() - (a*mu+b) )*sqrt( TMath::TwoPi() )*sigma;
    func->SetParNames  ( "area", "mean", "sigma", "slope", "offset" );
    func->SetParameters( area,   mu,     sigma,   a,       b );
    break;
  case 12: // gaussian + parabola
    a     = 0;
    b     = ( hist->GetBinContent(xbin)-hist->GetBinContent(1))/(hist->GetBinCenter(xbin)-hist->GetBinCenter(1) );
    c     = hist->GetBinContent(1)-b*hist->GetBinCenter(1);
    mu    = hist->GetBinCenter(hist->GetMaximumBin());
    sigma = 3*hist->GetBinWidth(1);
    area  = ( hist->GetMaximum() - (a*mu*mu+b*mu+c) )*sqrt( TMath::TwoPi() )*sigma;
    func->SetParNames( "area","mean","sigma","a","b","c" );
    func->SetParameters( area,mu,sigma,a,b,c );
    break;
  case  15:
  case 152:
  case 153:
    a     = ( hist->GetBinContent(xbin)-hist->GetBinContent(1))/(hist->GetBinCenter(xbin)-hist->GetBinCenter(1) );
    b     = hist->GetBinContent(1)-a*hist->GetBinCenter(1);
    mu    = hist->GetBinCenter( hist->GetMaximumBin() );
    sigma = 2*hist->GetBinWidth(1);
    area  = ( hist->GetMaximum() - (a*mu+b) )*sqrt( TMath::TwoPi() )*sigma;
    func->SetParNames  ( "area","mean","sigma","norm","Ebeam","a" );
    func->SetParameters( area,   mu,   sigma,  100*hist->GetBinContent(1),  5.289,    -10 );
    break;
  case 151:
    a     = ( hist->GetBinContent(xbin)-hist->GetBinContent(1))/(hist->GetBinCenter(xbin)-hist->GetBinCenter(1) );
    b     = hist->GetBinContent(1)-a*hist->GetBinCenter(1);
    mu    = hist->GetBinCenter( hist->GetMaximumBin() );
    sigma = 2*hist->GetBinWidth(1);
    area  = ( hist->GetMaximum() - (a*mu+b) )*sqrt( TMath::TwoPi() )*sigma;
    func->SetParNames  ( "area","mean","sigma","norm","Ebeam","a", "new" );
    func->SetParameters( area,   mu,   sigma,  100*hist->GetBinContent(1),  5.289,  -10,  0.5 );
    break;
  case 181: // for gauss + threshold1
    func->SetParNames  ( "area","mean","sigma","norm","Ebeam","k" );
    func->SetParameters( area,   mu,   sigma,  100*hist->GetBinContent(1),  5.289,  -10 );
    break;
  case 182: // for gauss + threshold1
    func->SetParNames  ( "area","mean","sigma","norm","Ebeam","k", "new" );
    func->SetParameters( area,   mu,   sigma,  100*hist->GetBinContent(1),  5.289,  -10, 0.5 );
    break;
  case 16: // gaussian + straight(area parametrization)
    a     = (hist->GetBinContent(xbin)-hist->GetBinContent(1))/(hist->GetBinCenter(xbin)-hist->GetBinCenter(1));
    b     = hist->GetBinContent(1)-a*hist->GetBinCenter(1);
    mu    = hist->GetBinCenter(hist->GetMaximumBin());
    sigma = 3*hist->GetBinWidth(1);
    area  = (hist->GetMaximum() - (a*mu+b))*sqrt(TMath::TwoPi())*sigma;
    func->SetParNames  ("area", "mean", "sigma", "linear_area",       "offset");
    func->SetParameters(area,   mu,     sigma,   (a*window_mean+b)*window_width, b);
    break;
  case 20: // double gaussian
    sigma = 10*hist->GetBinWidth(1);
    //sigma = 2.5*hist->GetBinWidth(1);
    mu    = hist->GetBinCenter(hist->GetMaximumBin());
    area  = hist->GetMaximum()*sqrt(TMath::TwoPi())*sigma;
    func->SetParNames  ("area",  "are_ratio", "mean", "sigma",  "sigma_ratio");
    func->SetParameters(area,     0.6,        mu,     sigma,     0.4*sigma);
    break;
  case 21: // double gaussian + straight line
    a     = (hist->GetBinContent(xbin)-hist->GetBinContent(1))/(hist->GetBinCenter(xbin)-hist->GetBinCenter(1));
    b     = hist->GetBinContent(1)-a*hist->GetBinCenter(1);
    mu    = hist->GetBinCenter(hist->GetMaximumBin());
    sigma = 2*hist->GetBinWidth(1);
    area  = (hist->GetMaximum() - (a*mu+b))*sqrt(TMath::TwoPi())*sigma;
    func->SetParNames  ("area", "area_ratio", "mean", "sigma", "sigma_ratio","slope","offset");
    func->SetParameters(area,    0.6,         mu,     sigma,   0.4,           a,      b);
    func->SetParLimits(0, 0, 10000);
    func->SetParLimits(1, 0, 1.0);
    //func->SetParLimits(3, 0, 0.05);
    //func->SetParLimits(4, 0, 0.5);
    break;
  case 22: // double gaussian + parabola
    a     = 0;
    b     = (hist->GetBinContent(xbin)-hist->GetBinContent(1))/(hist->GetBinCenter(xbin)-hist->GetBinCenter(1));
    c     = hist->GetBinContent(1)-b*hist->GetBinCenter(1);
    mu    = hist->GetBinCenter(hist->GetMaximumBin());
    sigma = 3*hist->GetBinWidth(1);
    area  = (hist->GetMaximum() - (a*mu*mu+b*mu+c))*sqrt(TMath::TwoPi())*sigma;
    func->SetParNames  ("area", "area_ratio", "mean", "sigma", "sigma_ratio", "a", "b","c");
    func->SetParameters(area,   0.6,          mu,     sigma,   0.4           , a,   b,  c);
    func->SetParLimits(0, 0, 10000);
    func->SetParLimits(1, 0, 1.0);
    func->SetParLimits(3, 0, 10000);
    func->SetParLimits(4, 0, 0.5);
    break;
  case 30: // triple gaussian
    sigma = 10*hist->GetBinWidth(1);
    //sigma = 2.5*hist->GetBinWidth(1);
    mu    = hist->GetBinCenter(hist->GetMaximumBin());
    area  = hist->GetMaximum()*sqrt(TMath::TwoPi())*sigma;
    func->SetParNames  ("area",  "are_ratio1", "area_ratio2", "mean", "sigma",  "sigma_ratio1", "sigma_ratio2" );
    func->SetParameters(area,     0.6,               0.3,       mu,     sigma,     0.4,            0.4         );
    func->SetParLimits(1, 0, 1.0);
    func->SetParLimits(2, 0, 1.0);
    func->SetParLimits(5, 0, 1.0);
    func->SetParLimits(6, 0, 1.0);
    break;
  case 31: // triple gaussian + straight line
    a     = (hist->GetBinContent(xbin)-hist->GetBinContent(1))/(hist->GetBinCenter(xbin)-hist->GetBinCenter(1));
    b     = hist->GetBinContent(1)-a*hist->GetBinCenter(1);
    mu    = hist->GetBinCenter(hist->GetMaximumBin());
    sigma = 30*hist->GetBinWidth(1);
    area  = (hist->GetMaximum() - (a*mu+b))*sqrt(TMath::TwoPi())*sigma;
    func->SetParNames  ("area", "area_ratio1", "area_ratio2", "mean", "sigma", "sigma_ratio1", "sigma_ratio2", "slope","offset");
    func->SetParameters(area,    0.6,                0.3,       mu,     sigma,   0.1,             0.3      ,     a,      b);
    func->SetParLimits(1, 0, 1.0);
    func->SetParLimits(2, 0, 1.0);
    func->SetParLimits(5, 0, 1.0);
    func->SetParLimits(6, 0, 1.0);
    break;
  case 32: // triple gaussian + parabola
    a     = 0;
    b     = (hist->GetBinContent(xbin)-hist->GetBinContent(1))/(hist->GetBinCenter(xbin)-hist->GetBinCenter(1));
    c     = hist->GetBinContent(1)-b*hist->GetBinCenter(1);
    mu    = hist->GetBinCenter(hist->GetMaximumBin());
    sigma = 3*hist->GetBinWidth(1);
    area  = (hist->GetMaximum() - (a*mu*mu+b*mu+c))*sqrt(TMath::TwoPi())*sigma;
    func->SetParNames  ("area", "area_ratio1", "area_ratio2", "mean", "sigma", "sigma_ratio1", "sigma_ratio2", "a","b", "c");
    func->SetParameters(area,    0.6,                0.3,       mu,     sigma,   0.4,             0.4,          a,  b,   c );
    func->SetParLimits(1, 0, 1.0);
    func->SetParLimits(2, 0, 1.0);
    func->SetParLimits(5, 0, 1.0);
    func->SetParLimits(6, 0, 1.0);
    break;
  case 100: // bif-gaussian
    sigma = 2*hist->GetBinWidth(1);
    mu    = hist->GetBinCenter( hist->GetMaximumBin() );
    area  = hist->GetMaximum()*sqrt( TMath::TwoPi() )*sigma;
    func->SetParNames( "area","mean","sigmal", "sigmah" );
    func->SetParameters( area,mu,sigma, sigma );
    break;
  case 110: // bif-gaussian + gaussian
    sigma = 2*hist->GetBinWidth(1);
    mu    = hist->GetBinCenter( hist->GetMaximumBin() );
    area  = hist->GetMaximum()*sqrt( TMath::TwoPi() )*sigma;
    func->SetParNames( "area",   "mean","sigmal", "sigmah", "area_ratio", "sigma" );
    func->SetParameters( 0.5*area,mu,    sigma,    sigma,  2,  sigma );
    break;
  case 111: // bif-gaussian + gaussian + straight
    a     = (hist->GetBinContent(xbin)-hist->GetBinContent(1))/(hist->GetBinCenter(xbin)-hist->GetBinCenter(1));
    b     = hist->GetBinContent(1)-a*hist->GetBinCenter(1);
    sigma = 4*hist->GetBinWidth(1);
    mu    = hist->GetBinCenter( hist->GetMaximumBin() );
    area  = hist->GetMaximum()*sqrt( TMath::TwoPi() )*sigma;
    func->SetParNames( "area",   "mean","sigmal", "sigmah", "area_ratio", "sigma", "slope", "offset" );
    func->SetParameters( 0.5*area, mu,   sigma,    sigma,  3,  sigma,    a,       b      );
    break;
  case 50: // for argus
  case 52: // for argus2
  case 53: // for argus2
    func->SetParNames  ("norm","Ebeam","a");
    func->SetParameters( 80,    5.29, -1700);
    //func->FixParameter(1,5.41569);
    break;
  case 51: // for modified argus
    func->SetParNames  ("norm","Ebeam","a", "new");
    func->SetParameters( 100*hist->GetBinContent(1),    5.29, -10, 0.5);
    break;
  case 81: // for threshold1
    func->SetParNames  ("norm","Ebeam","k");
    func->SetParameters( 100*hist->GetBinContent(1),    5.29, -0.2);
    break;
  case 82: // for threshold2
    func->SetParNames  ("norm","Ebeam","k", "new");
    func->SetParameters( 100*hist->GetBinContent(1),    5.29, -0.2, 0.5 );
    break;
  case 1091: // chi
    func->SetParNames  (  "area", "alpha", "n" );
    func->SetParameters(    1,     1,      1   );
    break;
  case 1092: // chi2
    func->SetParNames( "area1","alpha1","n1","area2","alpha2", "n2" );
    func->SetParameters(  1,     1,      1,     1,      1,       1  );
    break;
  case 310: // cb( bifuracated gaussian )
    func->SetParNames( "area","alpha","n","mean","sigma" );
    sigma = 5*hist->GetBinWidth(1);
    mu    = hist->GetBinCenter( hist->GetMaximumBin() );
    area  = hist->GetMaximum()*sqrt( TMath::TwoPi() )*sigma;
    func->SetParameters(area,  1.0,  1.3, mu, sigma );
    break;
  case 311: // cb( bifuracated gaussian ) + straight
    func->SetParNames( "area","alpha","n","mean","sigma", "slope", "offset" );
    sigma = 5*hist->GetBinWidth(1);
    mu    = hist->GetBinCenter( hist->GetMaximumBin() );
    area  = hist->GetMaximum()*sqrt( TMath::TwoPi() )*sigma;
    a = ( hist->GetBinContent(xbin)-hist->GetBinContent(1))/(hist->GetBinCenter(xbin)-hist->GetBinCenter(1) );
    b = hist->GetBinContent(1)-a*hist->GetBinCenter(1);
    func->SetParameters(area,  1.0,  1.3, mu, sigma,  a, b);
    break;
  case 410: // cb( bifuracated gaussian )
    func->SetParNames( "area","alpha","n","mean","sigma_l", "sigma_h" );
    sigma = 5*hist->GetBinWidth(1);
    mu    = hist->GetBinCenter( hist->GetMaximumBin() );
    area  = hist->GetMaximum()*sqrt( TMath::TwoPi() )*sigma;
    func->SetParameters(area,  1.0,  1.3, mu, sigma, sigma);
    break;
  case 411: // cb( bifuracated gaussian ) + straight
    func->SetParNames( "area","alpha","n","mean","sigma_l", "sigma_h", "slope", "offset" );
    sigma = 5*hist->GetBinWidth(1);
    mu    = hist->GetBinCenter( hist->GetMaximumBin() );
    area  = hist->GetMaximum()*sqrt( TMath::TwoPi() )*sigma;
    a = ( hist->GetBinContent(xbin)-hist->GetBinContent(1))/(hist->GetBinCenter(xbin)-hist->GetBinCenter(1) );
    b = hist->GetBinContent(1)-a*hist->GetBinCenter(1);
    func->SetParameters(area,  1.0,  1.3, mu, sigma, sigma, a, b);
    break;
  case 20000: // cb( bifuracated gaussian )
    sigma = 4.0*hist->GetBinWidth(1);
    mu    = hist->GetBinCenter(hist->GetMaximumBin());
    area  = hist->GetMaximum()*sqrt(TMath::TwoPi())*sigma;
    func->SetParNames  ("area","alpha_l","n_l","mean","sigma","alpha_h","n_h" );
    func->SetParameters(area,      1.5,   1,  mu,  sigma,        2.5,     0.5   );
    func->SetParLimits(3, -0.01, 0.01);
    break;
  case 20001: // cb( bifuracated gaussian ) + straight
    sigma = 4.0*hist->GetBinWidth(1);
    mu    = hist->GetBinCenter(hist->GetMaximumBin());
    area  = hist->GetMaximum()*sqrt(TMath::TwoPi())*sigma;
    a = ( hist->GetBinContent(xbin)-hist->GetBinContent(1))/(hist->GetBinCenter(xbin)-hist->GetBinCenter(1) );
    b = hist->GetBinContent(1)-a*hist->GetBinCenter(1);
    func->SetParNames  ("area","alpha_l","n_l","mean","sigma","alpha_h","n_h", "slope", "offset");
    func->SetParameters(area,      1.5,   1,  mu,  sigma,        2.5,     0.5,   a,         b   );
    func->SetParLimits(3, -0.01, 0.01);
    break;
  case 30000: // cb( bifuracated gaussian )
    sigma = 4.0*hist->GetBinWidth(1);
    mu    = hist->GetBinCenter(hist->GetMaximumBin());
    area  = hist->GetMaximum()*sqrt(TMath::TwoPi())*sigma;
    func->SetParNames  ("area","alpha_l","n_l","mean","sigma_l","alpha_h","n_h","sigma_h");
    func->SetParameters(area,      1.5,   1,  mu,  sigma,        2.5,     0.5,    sigma);
    func->SetParLimits(3, -0.01, 0.01);
    break;
  case 30001: // cb( bifuracated gaussian ) + straight
    sigma = 4.0*hist->GetBinWidth(1);
    mu    = hist->GetBinCenter(hist->GetMaximumBin());
    area  = hist->GetMaximum()*sqrt(TMath::TwoPi())*sigma;
    a = ( hist->GetBinContent(xbin)-hist->GetBinContent(1))/(hist->GetBinCenter(xbin)-hist->GetBinCenter(1) );
    b = hist->GetBinContent(1)-a*hist->GetBinCenter(1);
    func->SetParNames  ("area","alpha_l","n_l","mean","sigma_l","alpha_h","n_h","sigma_h", "slope", "offset");
    func->SetParameters(area,      1.5,   1,  mu,  sigma,        2.5,     0.5,    sigma,     a,         b);
    func->SetParLimits(3, -0.01, 0.01);
    break;
  default:
    std::cerr << "indefinite fit function -> abort()" << std::endl;
    abort();
  }
  return;
}

// *************************************************************************************::

void func_get_integral(Int_t select, TF1* func,
		       Double_t bin_width,
		       Double_t& yield, Double_t& yieldE){
  if( select == 21 || select == 31 || select == 111 ){
    yield  = func->GetParameter(0);
    yieldE = func->GetParError(0);
  }else{
    std::cerr << "[ABORT] not supported yet ( sel_fun : " << select << " )" << std::endl, abort();
  }
  
  yield  /= bin_width;
  yieldE /= bin_width;
}

// *************************************************************************************::

// 0. horizontal line
// p[0]=slope, p[1]=offset
Double_t func_horizontal(const Double_t* x, const Double_t* p){
  Double_t f = p[0];
  return f;
}

// 1. straight line
// p[0]=slope, p[1]=offset
Double_t func_straight(const Double_t* x, const Double_t* p){
  Double_t f = p[0]*x[0]+p[1];
  return f;
}

// 2. parabola
// p[0]=ax^2, p[1]=bx, p[2]=c
Double_t func_parabola(const Double_t* x, const Double_t* p){
  Double_t f =  p[0]*x[0]*x[0] + p[1]*x[0] + p[2];
  return f;
}

// 3. cubic
// p[0]=ax^3, p[1]=bx^2, p[2]=cx, p[3]=d
Double_t func_cubic(const Double_t* x, const Double_t* p){
  Double_t f =  p[0]*x[0]*x[0]*x[0] + p[1]*x[0]*x[0] + p[2]*x[0] + p[3];
  return f;
}

// 4. fourth
// p[0]=ax^4, p[1]=bx^3, p[2]=cx^2, p[3]=dx, p[4]=e
Double_t func_fourth(const Double_t* x, const Double_t* p){
  Double_t f =  p[0]*x[0]*x[0]*x[0]*x[0] + p[1]*x[0]*x[0]*x[0] + p[2]*x[0]*x[0] + p[3]*x[0] + p[4];
  return f;
}

// 5. fifth
// p[0]=ax^5, p[1]=bx^4, p[2]=cx^3, p[3]=dx^2, p[4]=ex, p[5]=f
Double_t func_fifth(const Double_t* x, const Double_t* p){
  Double_t f =  p[0]*x[0]*x[0]*x[0]*x[0]*x[0] + p[1]*x[0]*x[0]*x[0]*x[0] + p[2]*x[0]*x[0]*x[0] + p[3]*x[0]*x[0] + p[4]*x[0] + p[5];
  return f;
}

// 6. straight line with area parametrization
// p[0]=area, p[1]=offset,
Double_t func_straight_area(const Double_t* x, const Double_t* p){
  Double_t f = (p[0]/window_width-p[1])/window_mean*x[0]+p[1];
  return f;
}

// 10. gaussian
// p[0]=area, p[1]=mean, p[2]=sigma,
Double_t func_gauss(const Double_t* x, const Double_t* p){
  Double_t f = p[0]/sqrt(TMath::TwoPi())/p[2]*exp(-0.5*((x[0]-p[1])/p[2])*((x[0]-p[1])/p[2]));
  return f;
}

// 11. gaussian(p[0-2]) + straight(p[3-4])
Double_t func_gauss_straight(const Double_t* x, const Double_t* p){
  Double_t f = func_gauss(x,&p[0]) + func_straight(x,&p[3]);
  return f;
}

// 12. gaussian(p[0-2]) + parabola(p[3-5])
Double_t func_gauss_parabola(const Double_t* x, const Double_t* p){
  Double_t f = func_gauss(x,&p[0]) + func_parabola(x,&p[3]);
  return f;
}

// 15. gaussian(p[0-2]) + argus(p[3-5])
Double_t func_gauss_argus(const Double_t* x, const Double_t* p){
  Double_t f = func_gauss(x,&p[0]) + func_argus(x,&p[3]);
  return f;
}

// 151. gaussian(p[0-2]) + modified-argus(p[3-6])
Double_t func_gauss_modargus(const Double_t* x, const Double_t* p){
  Double_t f = func_gauss(x,&p[0]) + func_modargus(x,&p[3]);
  return f;
}

// 152. gaussian(p[0-2]) + modified-argus2(p[3-5])
Double_t func_gauss_modargus2(const Double_t* x, const Double_t* p){
  Double_t f = func_gauss(x,&p[0]) + func_modargus2(x,&p[3]);
  return f;
}

// 153. gaussian(p[0-2]) + modified-argus3(p[3-5])
Double_t func_gauss_modargus3(const Double_t* x, const Double_t* p){
  Double_t f = func_gauss(x,&p[0]) + func_modargus3(x,&p[3]);
  return f;
}

// 16. gaussian(p[0-2]) + straight_area(p[3-4])
Double_t func_gauss_straight_area(const Double_t* x, const Double_t* p){
  Double_t f = func_gauss(x,&p[0]) + func_straight_area(x,&p[3]);
  return f;
}

// 20. double gaussian
// p[0]=area,  p[1]=area_ratio,
// p[2]=mean,
// p[3]=sigma, p[4]=sigma_ratio
Double_t func_2gauss(const Double_t* x, const Double_t* p){
  Double_t f = p[0]*  p[1]  /sqrt(TMath::TwoPi())/(p[3]     ) * exp( -0.5 * ((x[0]-p[2])/(p[3]     )) * ((x[0]-p[2])/(p[3]     )) )
    +          p[0]*(1-p[1])/sqrt(TMath::TwoPi())/(p[3]*p[4]) * exp( -0.5 * ((x[0]-p[2])/(p[3]*p[4])) * ((x[0]-p[2])/(p[3]*p[4])) );
  return f;
}

// 21. double gaussian(p[0-4]) + straight(p[5-6])
Double_t func_2gauss_straight(const Double_t* x, const Double_t* p){
  Double_t f = func_2gauss(x,&p[0]) + func_straight(x,&p[5]);
  return f;
}

// 22. double gaussian(p[0-4]) + parabola(p[5-7])
Double_t func_2gauss_parabola(const Double_t* x, const Double_t* p){
  Double_t f = func_2gauss(x,&p[0]) + func_parabola(x,&p[5]);
  return f;
}

// 30. triple gaussian
// p[0]=area,  p[1]=area_ratio1,   p[2]=area_ratio2,
// p[3]=mean,
// p[4]=sigma, p[5]=sigma_ratio1, p[6]=sigma_ratio2
Double_t func_3gauss(const Double_t* x, const Double_t* p){
  Double_t f = p[0]*  p[1]        /sqrt(TMath::TwoPi())/(p[4]      ) * exp(-0.5*((x[0]-p[3])/(p[4]      )) * ((x[0]-p[3])/(p[4]      )))
    +          p[0]*  p[2]        /sqrt(TMath::TwoPi())/(p[4]*p[5] ) * exp(-0.5*((x[0]-p[3])/(p[4]*p[5] )) * ((x[0]-p[3])/(p[4]*p[5] )))
    +          p[0]*(1-p[1]-p[2]) /sqrt(TMath::TwoPi())/(p[4]*p[6] ) * exp(-0.5*((x[0]-p[3])/(p[4]*p[6] )) * ((x[0]-p[3])/(p[4]*p[6] )));
  return f;
}

// 31. triple gaussian(p[0-6]) + straight(p[7-8])
Double_t func_3gauss_straight(const Double_t* x, const Double_t* p){
  Double_t f = func_3gauss(x,&p[0]) + func_straight(x,&p[7]);
  return f;
}

// 32. triple gaussian(p[0-6]) + parabola(p[7-9])
Double_t func_3gauss_parabola(const Double_t* x, const Double_t* p){
  Double_t f = func_3gauss(x,&p[0]) + func_parabola(x,&p[7]);
  return f;
}

// 100. bif-gaussian
// p[0]=area, p[1]=mean, p[2]=sigma(low), p[3]=sigma(high)
Double_t func_gauss_bif(const Double_t* x, const Double_t* p){
  Double_t f;
  if( x[0]<p[1] ) f = 2*p[0]/(1+p[3]/p[2])/sqrt(TMath::TwoPi())/p[2]*exp(-0.5*((x[0]-p[1])/p[2])*((x[0]-p[1])/p[2]));
  else            f = 2*p[0]/(1+p[3]/p[2])/sqrt(TMath::TwoPi())/p[2]*exp(-0.5*((x[0]-p[1])/p[3])*((x[0]-p[1])/p[3]));
  return f;
}

// 110. bif-gaussian + gaussian
// p[0]=area, p[1]=mean, p[2]=sigma(low), p[3]=sigma(high), p[4]=area_ratio, p[5]=sigma
Double_t func_gauss_bif_gauss(const Double_t* x, const Double_t* p){
  Double_t f;
  if( x[0]<p[1] ) f = 2*p[0]*(1-p[4])/(1+p[3]/p[2])/sqrt(TMath::TwoPi())/p[2]*exp(-0.5*((x[0]-p[1])/p[2])*((x[0]-p[1])/p[2]));
  else            f = 2*p[0]*(1-p[4])/(1+p[3]/p[2])/sqrt(TMath::TwoPi())/p[2]*exp(-0.5*((x[0]-p[1])/p[3])*((x[0]-p[1])/p[3]));
  f += p[0]*p[4]/sqrt(TMath::TwoPi())/p[5]*exp(-0.5*((x[0]-p[1])/p[5])*((x[0]-p[1])/p[5]));
  return f;
}

// 111. bif-gaussian+gaussian(p[0-5]) + straight(p[6-7])
Double_t func_gauss_bif_gauss_straight(const Double_t* x, const Double_t* p){
  Double_t f = func_gauss_bif_gauss(x,&p[0]) + func_straight(x,&p[6]);
  return f;
}

// 50. argus
// p[0]=normalization, p[1]=Ebeam, p[2]=a
Double_t func_argus(const Double_t* x, const Double_t* p){
  Double_t f;
  if(x[0]<p[1]){
    f = p[0]*x[0]*sqrt(1-x[0]*x[0]/(p[1]*p[1]))*exp(p[2]*(1-x[0]*x[0]/p[1]/p[1]));
  }else{
    f=0;
  }
  return f;
}

// 51. modified argus
// p[0]=normalization, p[1]=Ebeam, p[2]=a, p[3]=new parameter
Double_t func_modargus(const Double_t* x, const Double_t* p){
  Double_t f;
  if(x[0]<p[1]){
    f = p[0]*x[0]*pow(1-x[0]*x[0]/(p[1]*p[1]), p[3])*exp(p[2]*(1-x[0]*x[0]/p[1]/p[1]));
  }else{
    f=0;
  }
  return f;
}

// 52. modified argus2
// p[0]=normalization, p[1]=Ebeam, p[2]=a
Double_t func_modargus2(const Double_t* x, const Double_t* p){
  Double_t f;
  if(x[0]<p[1]){
    f = p[0]*x[0]*sqrt(1-x[0]*x[0]/(p[1]*p[1]))*exp(p[2]*(1-x[0]*x[0]/p[1]/p[1])*(1-x[0]*x[0]/p[1]/p[1]));
  }else{
    f=0;
  }
  return f;
}

// 53. modified argus3
// p[0]=normalization, p[1]=Ebeam, p[2]=a
Double_t func_modargus3(const Double_t* x, const Double_t* p){
  Double_t f;
  if(x[0]<p[1]){
    f = p[0]*x[0]*sqrt(1-x[0]*x[0]/(p[1]*p[1]))*exp(p[2]*(1-x[0]*x[0]/p[1]/p[1])*(1-x[0]*x[0]/p[1]/p[1])*(1-x[0]*x[0]/p[1]/p[1]));
  }else{
    f=0;
  }
  return f;
}

//  1091
Double_t func_chi( const Double_t *x, const Double_t *p ){
  Double_t f = p[0]*(TMath::Power(p[1]*x[0],p[2]-1.)*
		     TMath::Exp(-1.* p[1]*x[0])/
		     TMath::Gamma(p[2]));
  return f;
}

//  1092
Double_t func_chi2( const Double_t *x, const Double_t *p ){
  Double_t f = func_chi(x,&p[0]) + func_chi(x,&p[3]);  
  return f;
}

// 310. cb function (low-side tail)
// p[0]=area(gaussian)
// p[1]=alpha(range of gauss), p[2]=n(tail shape),
// p[3]=xbar(mean), p[4]=sigma
Double_t func_cb_low(const Double_t* x, const Double_t*p){
  Double_t z   = x[0]-p[3];
  Double_t dim = 1.0-p[1]/p[2]*z/p[4]-p[1]*p[1]/p[2];
  Double_t f;

  if( z<-p[4]*p[1]){
    f = p[0]/sqrt(TMath::TwoPi())/p[4]*exp(-0.5*p[1]*p[1])/pow(dim,p[2]); // high-side tail
  }else{
    f = p[0]/sqrt(TMath::TwoPi())/p[4]*exp(-0.5*z*z/(p[4]*p[4])); // gauss
  }
  return f;
}

// 311. cb function (low-side tail)(p[0-4]) + straight(p[5-6])
Double_t func_cb_low_straight(const Double_t* x, const Double_t*p){
  Double_t f = func_cb_low(x,&p[0]) + func_straight(x,&p[5]);
  return f;
}


// 410. cb function (low-side tail,bif)
// p[0]=area(gaussian)
// p[1]=alpha(range of gauss), p[2]=n(tail shape),
// p[3]=xbar(mean), p[4]=sigma_low, p[5]=sigma_high
Double_t func_cb_low_bif(const Double_t* x, const Double_t*p){
  Double_t z   = x[0]-p[3];
  Double_t dim = 1.0-p[1]/p[2]*z/p[4]-p[1]*p[1]/p[2];
  Double_t f;

  if( z<-p[4]*p[1]){
    f = p[0]/sqrt(TMath::TwoPi())/p[4]*exp(-0.5*p[1]*p[1])/pow(dim,p[2]); // high-side tail
  }else if( z<0 ){
    f = p[0]/sqrt(TMath::TwoPi())/p[4]*exp(-0.5*z*z/(p[4]*p[4])); // gauss(low-side)
  }else{
    f = p[0]/sqrt(TMath::TwoPi())/p[4]*exp(-0.5*z*z/(p[5]*p[5])); // gauss(high-side)
  }
  return f;
}

// 411. cb function (low-side tail,bif)(p[0-5]) + straight(p[6-7])
Double_t func_cb_low_bif_straight(const Double_t* x, const Double_t*p){
  Double_t f = func_cb_low_bif(x,&p[0]) + func_straight(x,&p[6]);
  return f;
}

// 20000. cb function (bifurcated gaussian)
// p[0]=area(gaussian)
// p[1]=alpha_low(range of gauss), p[2]=n_low(tail shape),
// p[3]=xbar(mean), p[4]=sigma
// p[5]=alpha_high(range of gauss), p[6]=n_high(tail shape),
Double_t func_cb(const Double_t* x, const Double_t*p){
  Double_t z   = x[0]-p[3];
  Double_t dim_l = 1.0+p[1]/p[2]*fabs(z)/p[4]-p[1]*p[1]/p[2];
  Double_t dim_h = 1.0+p[5]/p[6]*fabs(z)/p[4]-p[5]*p[5]/p[6];
  Double_t f;
  if(      z <-p[4]*p[1] ) f = p[0]/sqrt(TMath::TwoPi())/p[4]*exp(-0.5*p[1]*p[1])/pow(dim_l,p[2]); // low-side tail
  else if( z > p[4]*p[5] ) f = p[0]/sqrt(TMath::TwoPi())/p[4]*exp(-0.5*p[5]*p[5])/pow(dim_h,p[6]); // high-side tail
  else                     f = p[0]/sqrt(TMath::TwoPi())/p[4]*exp(-0.5*z*z/(p[4]*p[4])); // gaussian (high-side)
  return f;
}

// 20001. cb functione (bifurcated gaussian)(p[0-6]) + straight(p[7-8])
Double_t func_cb_straight(const Double_t* x, const Double_t*p){
  Double_t f = func_cb(x,&p[0]) + func_straight(x,&p[7]);
  return f;
}

// 30000. cb function (bifurcated gaussian)
// p[0]=area(gaussian)
// p[1]=alpha_low(range of gauss), p[2]=n_low(tail shape),
// p[3]=xbar(mean), p[4]=sigma_low
// p[5]=alpha_high(range of gauss), p[6]=n_high(tail shape),
// p[7]=sigma_high
Double_t func_cb_bif(const Double_t* x, const Double_t*p){
  Double_t z   = x[0]-p[3];
  Double_t dim_l = 1.0+p[1]/p[2]*fabs(z)/p[4]-p[1]*p[1]/p[2];
  Double_t dim_h = 1.0+p[5]/p[6]*fabs(z)/p[7]-p[5]*p[5]/p[6];
  Double_t f;
  if(      z <-p[4]*p[1] ) f = p[0]/sqrt(TMath::TwoPi())/p[4]*exp(-0.5*p[1]*p[1])/pow(dim_l,p[2]); // low-side tail
  else if( z > p[7]*p[5] ) f = p[0]/sqrt(TMath::TwoPi())/p[4]*exp(-0.5*p[5]*p[5])/pow(dim_h,p[6]); // high-side tail
  else if( z < 0         ) f = p[0]/sqrt(TMath::TwoPi())/p[4]*exp(-0.5*z*z/(p[4]*p[4])); // gaussian (low-side)
  else                     f = p[0]/sqrt(TMath::TwoPi())/p[4]*exp(-0.5*z*z/(p[7]*p[7])); // gaussian (high-side)
  return f;
}

// 30001. cb functione (bifurcated gaussian)(p[0-7]) + straight(p[8-9])
Double_t func_cb_bif_straight(const Double_t* x, const Double_t*p){
  Double_t f = func_cb_bif(x,&p[0]) + func_straight(x,&p[8]);
  return f;
}

// 81. threshold function
// p[0]=norm, p[1]=Ebeam, p[2]=shape
Double_t func_threshold1( const Double_t* x, const Double_t* p ){
  Double_t f;
  if( x[0]<p[1] ){
    //f = p[0]*x[0]*sqrt(1-x[0]*x[0]/(p[1]*p[1]))*exp(p[2]*(1-x[0]/p[1])*(1-x[0]/p[1]));
    f = p[0]*sqrt(1-x[0]*x[0]/(p[1]*p[1]))*exp((p[1]-x[0])/p[2]); // threshold function in delta-M fitting
  }else{
    f=0;
  }
  return f;
}

// 82. threshold function
// p[0]=norm, p[1]=Ebeam, p[2]=shape, p[3]=new parameter
Double_t func_threshold2( const Double_t* x, const Double_t* p ){
  Double_t f;
  if( x[0]<p[1] ){
    //f = p[0]*x[0]*pow(1-x[0]*x[0]/(p[1]*p[1]), p[3])*exp(p[2]*(1-x[0]/p[1])*(1-x[0]/p[1]));
    f = p[0]*pow(1-x[0]*x[0]/(p[1]*p[1]), p[3])*exp((p[1]-x[0])/p[2]); // threshold function in delta-M fitting
  }else{
    f=0;
  }
  return f;
}

// 181. gaussian(p[0-2]) + threshold1(p[3-5])
Double_t func_gauss_threshold1(const Double_t* x, const Double_t* p){
  Double_t f = func_gauss(x,&p[0]) + func_threshold1(x,&p[3]);
  return f;
}

// 182. gaussian(p[0-2]) + threshold2(p[3-5])
Double_t func_gauss_threshold2(const Double_t* x, const Double_t* p){
  Double_t f = func_gauss(x,&p[0]) + func_threshold2(x,&p[3]);
  return f;
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++           PDF FOR ROOFIT           ++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// -50. argus
// p0=Ebeam, p1=a
Double_t func_roo_argus(const Double_t x, const Double_t p0, const Double_t p1){
  Double_t f;
  if(x<p0){
    f = x*sqrt(1-x*x/(p0*p0))*exp(p1*(1-x*x/p0/p0));
  }else{
    f=0;
  }
  return f;
}

// -51. modified argus
// p0=Ebeam, p1=a
Double_t func_roo_modargus(const Double_t x, const Double_t p0, const Double_t p1, const Double_t p2){
  //Double_t func_roo_modargus(const Double_t x, const Double_t p1, const Double_t p2){
  //Double_t p0 = 5.289;
  Double_t f;
  if(x<p0){
    f = x*pow(1-x*x/(p0*p0),p2)*exp(p1*(1-x*x/p0/p0));
  }else{
    f=0;
  }
  return f;
}

// -52. modified-argus2
// p0=Ebeam, p1=a
Double_t func_roo_modargus2(const Double_t x, const Double_t p0, const Double_t p1){
  Double_t f;
  if(x<p0){
    f = x*sqrt(1-x*x/(p0*p0))*exp(p1*(1-x*x/p0/p0)*(1-x*x/p0/p0));
  }else{
    f=0;
  }
  return f;
}

// -53. modified-argus3
// p0=Ebeam, p1=a
Double_t func_roo_modargus3(const Double_t x, const Double_t p0, const Double_t p1){
  Double_t f;
  if(x<p0){
    f = x*sqrt(1-x*x/(p0*p0))*exp(p1*(1-x*x/p0/p0)*(1-x*x/p0/p0)*(1-x*x/p0/p0));
  }else{
    f=0;
  }
  return f;
}

// -81. threshold function
// p0=Ebeam, p1=shape
Double_t func_roo_threshold1( const Double_t x, const Double_t p0, const Double_t p1 ){
  Double_t f;
  if( x<p0 ){
    //f = x*sqrt(1-x/p0)*exp(p1*(1-x/p0)*(1-x/p0)); // mod-sqrt*gauss
    //f = x*sqrt(1-x*x/(p0*p0))*exp(p1*(1-x/p0)*(1-x/p0)); // sqrt*gauss
    f = sqrt(1-x*x/(p0*p0))*exp((p0-x)/p1); // threshold function in delta-M fitting
  }else{
    f=0;
  }
  return f;
}

// -82. threshold function
// p0=Ebeam, p1=shape, p2=new
Double_t func_roo_threshold2( const Double_t x, const Double_t p0, const Double_t p1, const Double_t p2 ){
  Double_t f;
  if( x<p0 ){
    //f = x*pow(1-x/p0,p2)*exp(p1*(1-x/p0)*(1-x/p0)); // mod-sqrt*gauss
    //f = x*pow(1-x*x/(p0*p0),p2)*exp(p1*(1-x/p0)*(1-x/p0)); // sqrt*gauss
    f = pow(1-x*x/(p0*p0),p2)*exp((p0-x)/p1); // threshold function in delta-M fitting
  }else{
    f=0;
  }
  return f;
}
