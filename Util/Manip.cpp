#include "Manip.h"

TH1D* divide( TH1D* hist1, TH1D* hist2, Char_t* name, Char_t* title ){
  TH1D* hist = (TH1D*)hist1->Clone();
  hist->Reset();
  if( strlen(name ) ) hist->SetName ( name );
  else                hist->SetName ( Form("%s_%s",hist1->GetName(),hist2->GetName()) );
  if( strlen(title) ) hist->SetTitle( title );
  else                hist->SetTitle( Form("%s_%s",hist1->GetTitle(),hist2->GetTitle()) );
  for( Int_t m=1; m<=hist->GetNbinsX(); m++ ){
    // set value
    if( hist2->GetBinContent(m) ) hist->SetBinContent( m,hist1->GetBinContent(m)/hist2->GetBinContent(m) );
    else                          hist->SetBinContent( m,0 );
    // set error
    if( hist1->GetBinContent(m) && hist2->GetBinContent(m) ) hist->SetBinError( m,hist1->GetBinContent(m)/hist2->GetBinContent(m)
										*sqrt(
										      hist1->GetBinError(m)*hist1->GetBinError(m)/hist1->GetBinContent(m)/hist1->GetBinContent(m) +
										      hist2->GetBinError(m)*hist2->GetBinError(m)/hist2->GetBinContent(m)/hist2->GetBinContent(m)
										      )
										);
    else                                                         hist->SetBinError( m,0 );
  }
  
  return hist;
}


TH2D* divide( TH2D* hist1, TH2D* hist2, Char_t* name, Char_t* title ){
  TH2D* hist = (TH2D*)hist1->Clone();
  hist->Reset();
  if( strlen(name ) ) hist->SetName ( name );
  else                hist->SetName ( Form("%s_%s",hist1->GetName(),hist2->GetName()) );
  if( strlen(title) ) hist->SetTitle( title );
  else                hist->SetTitle( Form("%s_%s",hist1->GetTitle(),hist2->GetTitle()) );
  for( Int_t m=1; m<=hist->GetNbinsX(); m++ ){
    for( Int_t n=1; n<=hist->GetNbinsY(); n++ ){
      // set value
      if( hist2->GetBinContent(m,n) ) hist->SetBinContent( m,n,hist1->GetBinContent(m,n)/hist2->GetBinContent(m,n) );
      else                              hist->SetBinContent( m,n,0 );
      // set error
      if( hist1->GetBinContent(m,n) && hist2->GetBinContent(m,n) ) hist->SetBinError( m,n,hist1->GetBinContent(m,n)/hist2->GetBinContent(m,n)
										      *sqrt(
											    hist1->GetBinError(m,n)*hist1->GetBinError(m,n)/hist1->GetBinContent(m,n)/hist1->GetBinContent(m,n) +
											    hist2->GetBinError(m,n)*hist2->GetBinError(m,n)/hist2->GetBinContent(m,n)/hist2->GetBinContent(m,n)
											    )
										      );
      else                                                         hist->SetBinError( m,n,0 );
    }
  }
  
  return hist;
}


TH1D* divide_binomial( TH1D* hist1, TH1D* hist2, Char_t* name, Char_t* title ){
  TH1D* hist = (TH1D*)hist1->Clone();
  hist->Reset();
  if( strlen(name ) ) hist->SetName ( name );
  else                hist->SetName ( Form("%s_%s",hist1->GetName(),hist2->GetName()) );
  if( strlen(title) ) hist->SetTitle( title );
  else                hist->SetTitle( Form("%s_%s",hist1->GetTitle(),hist2->GetTitle()) );
  for( Int_t m=1; m<=hist->GetNbinsX(); m++ ){
    // set value
    if( hist2->GetBinContent(m) ) hist->SetBinContent( m,hist1->GetBinContent(m)/hist2->GetBinContent(m) );
    else                              hist->SetBinContent( m,0 );
    // set error
    if( hist1->GetBinContent(m) && hist2->GetBinContent(m) ) hist->SetBinError( m,
										(hist1->GetBinContent(m)/hist2->GetBinContent(m))
										*sqrt(
										      1/hist1->GetBinContent(m)*
										      ( 1 - hist1->GetBinContent(m)/hist2->GetBinContent(m))
										      )
										);
    else                                                         hist->SetBinError( m,0 );
  }
  
  return hist;
}


TH2D* divide_binomial( TH2D* hist1, TH2D* hist2, Char_t* name, Char_t* title ){
  TH2D* hist = (TH2D*)hist1->Clone();
  hist->Reset();
  if( strlen(name ) ) hist->SetName ( name );
  else                hist->SetName ( Form("%s_%s",hist1->GetName(),hist2->GetName()) );
  if( strlen(title) ) hist->SetTitle( title );
  else                hist->SetTitle( Form("%s_%s",hist1->GetTitle(),hist2->GetTitle()) );
  for( Int_t m=1; m<=hist->GetNbinsX(); m++ ){
    for( Int_t n=1; n<=hist->GetNbinsY(); n++ ){
      // set value
      if( hist2->GetBinContent(m,n) ) hist->SetBinContent( m,n,hist1->GetBinContent(m,n)/hist2->GetBinContent(m,n) );
      else                              hist->SetBinContent( m,n,0 );
      // set error
      if( hist1->GetBinContent(m,n) && hist2->GetBinContent(m,n) ) hist->SetBinError( m,n,
										      (hist1->GetBinContent(m,n)/hist2->GetBinContent(m,n))
										      *sqrt(
											    1/hist1->GetBinContent(m,n)*
											    ( 1 - hist1->GetBinContent(m,n)/hist2->GetBinContent(m,n))
											    )
										      );
      
      else                                                         hist->SetBinError( m,n,0 );
    }
  }
  
  return hist;
}


TH1D* multiply( TH1D* hist1, TH1D* hist2, Bool_t fl_wE, Char_t* name, Char_t* title ){
  TH1D* hist = (TH1D*)hist1->Clone();
  hist->Reset();
  if( strlen(name ) ) hist->SetName ( name );
  else                hist->SetName ( Form("%s_x_%s",hist1->GetName(),hist2->GetName()) );
  if( strlen(title) ) hist->SetTitle( title );
  else                hist->SetTitle( Form("%s_x_%s",hist1->GetTitle(),hist2->GetTitle()) );
  for( Int_t m=1; m<=hist->GetNbinsX(); m++ ){
    // set value
    hist->SetBinContent( m, hist1->GetBinContent(m)*hist2->GetBinContent(m) );
    // set error
    Double_t  error2  = hist1->GetBinError  (m) * hist1->GetBinError  (m) * hist2->GetBinContent(m) * hist2->GetBinContent(m);
    if(fl_wE) error2 += hist1->GetBinContent(m) * hist1->GetBinContent(m) * hist2->GetBinError  (m) * hist2->GetBinError  (m);
    hist->SetBinError( m, sqrt(error2) );
  }
  
  return hist;
}

TH2D* multiply( TH2D* hist1, TH2D* hist2, Bool_t fl_wE, Char_t* name, Char_t* title ){
  TH2D* hist = (TH2D*)hist1->Clone();
  hist->Reset();
  if( strlen(name ) ) hist->SetName ( name );
  else                hist->SetName ( Form("%s_x_%s",hist1->GetName(),hist2->GetName()) );
  if( strlen(title) ) hist->SetTitle( title );
  else                hist->SetTitle( Form("%s_x_%s",hist1->GetTitle(),hist2->GetTitle()) );
  for( Int_t m=1; m<=hist->GetNbinsX(); m++ ){
    for( Int_t n=1; n<=hist->GetNbinsY(); n++ ){
      // set value
      hist->SetBinContent( m,n,hist1->GetBinContent(m,n)*hist2->GetBinContent(m,n) );
      // set error
      Double_t  error2  = hist1->GetBinError  (m,n) * hist1->GetBinError  (m,n) * hist2->GetBinContent(m,n) * hist2->GetBinContent(m,n);
      if(fl_wE) error2 += hist1->GetBinContent(m,n) * hist1->GetBinContent(m,n) * hist2->GetBinError  (m,n) * hist2->GetBinError  (m,n);
      hist->SetBinError( m,n,sqrt(error2) );
    }
  }
  
  return hist;
}
