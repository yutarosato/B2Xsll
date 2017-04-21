#include "Canvas.h"
#include <iostream>
TCanvas* Canvas( const Char_t* name, const Char_t* title, Int_t n1, Int_t n2, Int_t width, Int_t height ){
  int nx=1;
  int ny=1;
  if( n2==0 ){
    while( nx*ny < n1 ){
      if(nx>ny) ny++;
      else nx++;
    }
  }else{
    nx = n1;
    ny = n2;
  }
  //if( width*nx  > 1800 ) width = 1800/nx;
  //if( height*ny > 1200 ) width = 1200/ny;
  TCanvas* c =  new TCanvas( name, title, width*nx, height*ny );
  //gROOT->SetStyle("Plain");
  c->Divide( nx, ny );
  return c;
}

void Deco( TH1D* hist,    Int_t mode, Int_t color_l, Int_t color_m ){
  switch( mode ){
  case 0: // Marker with error bar
    hist->SetLineColor  ( color_l );
    hist->SetMarkerColor( color_m );
    break;
  case 1: // Marker with error bar (slim width)
    hist->SetLineWidth  (1);
    //hist->SetMarkerStyle(1);
    hist->SetLineColor  ( color_l );
    hist->SetMarkerColor( color_m );
    break;
  case 2: // Histogram
    hist->SetLineWidth  (1);
    hist->SetMarkerStyle(1);
    hist->SetLineColor  ( color_l );
    hist->SetMarkerColor( color_m );
    break;
  case 3: // Fill
    hist->SetLineWidth  (1);
    hist->SetFillColor  ( color_l );
    hist->SetMarkerColor( color_m );
    break;
  defalut:
    std::cerr << Form("Check Deco-mode(%d)",mode) << std::endl;
    abort();
    break;
  }
  return;
}

void Deco( TF1* func, Int_t mode, Int_t color_l, Int_t color_m ){
  switch( mode ){
  case 0: //
    func->SetLineColor  ( color_l );
    func->SetMarkerColor( color_m );
    break;
  case 1: //
    func->SetLineColor  ( color_l );
    func->SetMarkerColor( color_m );
    func->SetLineWidth  ( 2 );
    break;
  defalut:
    std::cerr << Form("Check Deco-mode(%d)",mode) << std::endl;
    abort();
    break;
  }
  return;
}

void Deco( TGraph* graph, Int_t mode, Int_t color_l, Int_t color_m ){
  switch( mode ){
  case 0: // point plot (closed marker)
    graph->SetLineColor  ( color_l );
    graph->SetMarkerColor( color_m );
    graph->SetMarkerStyle(  20 );
    graph->SetMarkerSize ( 0.4 );
    break;
  case 1: // point plot (open marker)
    graph->SetLineColor  ( color_l );
    graph->SetMarkerColor( color_m );
    graph->SetMarkerStyle(  24 );
    graph->SetMarkerSize ( 0.4 );
    break;
  defalut:
    std::cerr << Form("Check Deco-mode(%d)",mode) << std::endl;
    abort();
    break;
  }
  return;
}

void Deco( TGraphErrors* graph, Int_t mode, Int_t color_l, Int_t color_m ){
  switch( mode ){
  case 0: // point plot (closed marker)
    graph->SetLineColor  ( color_l );
    graph->SetMarkerColor( color_m );
    graph->SetMarkerStyle(  20 );
    graph->SetMarkerSize ( 0.4 );
    break;
  case 1: // point plot (open marker)
    graph->SetLineColor  ( color_l );
    graph->SetMarkerColor( color_m );
    graph->SetMarkerStyle(  24 );
    graph->SetMarkerSize ( 0.4 );
    break;
  defalut:
    std::cerr << Form("Check Deco-mode(%d)",mode) << std::endl;
    abort();
    break;
  }
  return;
}

TH2D* Waku( Int_t n, TH1D** hist, const char* xlabel, const char* name, const char* title, Int_t ymin_flag, Double_t space ){
  Double_t range[4]; // {xmin,xmax,ymin,ymax}
  range[0] = hist[0]->GetXaxis()->GetXmin();
  range[1] = hist[0]->GetXaxis()->GetXmax();
  range[2] = hist[0]->GetMinimum();
  range[3] = hist[0]->GetMaximum();
  Int_t xbin = hist[0]->GetNbinsX();
  for( Int_t i=1; i<n; i++){
    if( range[0] > hist[i]->GetXaxis()->GetXmin() ) range[0] = hist[i]->GetXaxis()->GetXmin();
    if( range[1] < hist[i]->GetXaxis()->GetXmax() ) range[1] = hist[i]->GetXaxis()->GetXmax();
    if( range[2] > hist[i]->GetMinimum() ) range[2] = hist[i]->GetMinimum();
    if( range[3] < hist[i]->GetMaximum() ) range[3] = hist[i]->GetMaximum();
  }

  if( ymin_flag==1 ) range[2] = 0;
  Double_t dy = range[3] - range[2];
  range[3] += dy*space;
  //range[3] = 250; // kkkkkkk
  if     ( ymin_flag==2 && range[2] <=0 ) range[2] = 0.0001*range[3];
  else if( ymin_flag!=1 && ymin_flag!=2 ) range[2] -= dy*space;
    
  TH2D* waku = new TH2D(name,title,
			2,range[0],range[1],
			2,range[2],range[3]
			);
  waku->GetXaxis()->CenterTitle();
  waku->GetYaxis()->CenterTitle();
  waku->SetXTitle( xlabel );
  std::stringstream sTmp;
  sTmp << "Events/" <<  Double_t((range[1]-range[0])/xbin);
  waku->SetYTitle( sTmp.str().c_str() );
  //((TGaxis*)waku->GetXaxis())->SetMaxDigits(3); // move to the "Style.cpp"
  //((TGaxis*)waku->GetYaxis())->SetMaxDigits(3); // move to the "Style.cpp"

  return waku;
}

TH2D* Waku( Int_t n, TGraph** graph, const char* xlabel, const char* name, const char* title, Int_t ymin_flag, Double_t space ){
  Double_t range[4]; // {xmin,xmax,ymin,ymax}
  range[0] = graph[0]->GetXaxis()->GetXmin();
  range[1] = graph[0]->GetXaxis()->GetXmax();
  range[2] = graph[0]->GetYaxis()->GetXmin();
  range[3] = graph[0]->GetYaxis()->GetXmax();

  for( Int_t i=1; i<n; i++){
    if( range[0] > graph[i]->GetXaxis()->GetXmin() ) range[0] = graph[i]->GetXaxis()->GetXmin();
    if( range[1] < graph[i]->GetXaxis()->GetXmax() ) range[1] = graph[i]->GetXaxis()->GetXmax();
    if( range[2] > graph[i]->GetYaxis()->GetXmin() ) range[2] = graph[i]->GetYaxis()->GetXmin();
    if( range[3] < graph[i]->GetYaxis()->GetXmax() ) range[3] = graph[i]->GetYaxis()->GetXmax();
  }

  if( ymin_flag ) range[2] = 0;
  Double_t dy = range[3] - range[2];
  range[3] += dy*space;
  //range[3] = 4.5; // kkkkkkk
  if     ( ymin_flag==2 && range[2] <=0 ) range[2] = 0.0001*range[3];
  else if( ymin_flag!=1 && ymin_flag!=2 ) range[2] -= dy*space;
    
  TH2D* waku = new TH2D(name,title,
			2,range[0],range[1],
			2,range[2],range[3]
			);

  waku->GetXaxis()->CenterTitle();
  waku->GetYaxis()->CenterTitle();
  waku->SetXTitle( xlabel );
  
  waku->SetYTitle( "Events" );
  //((TGaxis*)waku->GetXaxis())->SetMaxDigits(3); // move to the "Style.cpp"
  //((TGaxis*)waku->GetYaxis())->SetMaxDigits(3); // move to the "Style.cpp"

  return waku;
}

TH2D* Waku( Int_t n, TGraphErrors** graph, const char* xlabel, const char* name, const char* title, Int_t ymin_flag, Double_t space ){
  Double_t range[4]; // {xmin,xmax,ymin,ymax}
  range[0] = graph[0]->GetXaxis()->GetXmin();
  range[1] = graph[0]->GetXaxis()->GetXmax();
  range[2] = graph[0]->GetYaxis()->GetXmin();
  range[3] = graph[0]->GetYaxis()->GetXmax();

  for( Int_t i=1; i<n; i++){
    if( range[0] > graph[i]->GetXaxis()->GetXmin() ) range[0] = graph[i]->GetXaxis()->GetXmin();
    if( range[1] < graph[i]->GetXaxis()->GetXmax() ) range[1] = graph[i]->GetXaxis()->GetXmax();
    if( range[2] > graph[i]->GetYaxis()->GetXmin() ) range[2] = graph[i]->GetYaxis()->GetXmin();
    if( range[3] < graph[i]->GetYaxis()->GetXmax() ) range[3] = graph[i]->GetYaxis()->GetXmax();
  }

  if( ymin_flag ) range[2] = 0;
  Double_t dy = range[3] - range[2];
  range[3] += dy*space;
  //range[3] = 4.5; // kkkkkkk
  if     ( ymin_flag==2 && range[2] <=0 ) range[2] = 0.0001*range[3];
  else if( ymin_flag!=1 && ymin_flag!=2 ) range[2] -= dy*space;
    
  TH2D* waku = new TH2D(name,title,
			2,range[0],range[1],
			2,range[2],range[3]
			);

  waku->GetXaxis()->CenterTitle();
  waku->GetYaxis()->CenterTitle();
  waku->SetXTitle( xlabel );
  
  waku->SetYTitle( "Events" );
  //((TGaxis*)waku->GetXaxis())->SetMaxDigits(3); // move to the "Style.cpp"
  //((TGaxis*)waku->GetYaxis())->SetMaxDigits(3); // move to the "Style.cpp"

  return waku;
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Double_t GetYMin( Int_t n, TH1D** hist ){
  Double_t ymin = hist[0]->GetMinimum();
  for( Int_t i=1; i<n; i++){
    if( ymin > hist[i]->GetMinimum() ) ymin = hist[i]->GetMinimum();
  }
  return ymin;
}

Double_t GetYMin( Int_t n, TGraph** graph ){
  Double_t ymin = graph[0]->GetYaxis()->GetXmin();
  for( Int_t i=1; i<n; i++){
    if( ymin > graph[i]->GetYaxis()->GetXmin() ) ymin = graph[i]->GetYaxis()->GetXmin();
  }
  return ymin;
}

Double_t GetYMin( Int_t n, TGraphErrors** graph ){
  Double_t ymin = graph[0]->GetYaxis()->GetXmin();
  for( Int_t i=1; i<n; i++){
    if( ymin > graph[i]->GetYaxis()->GetXmin() ) ymin = graph[i]->GetYaxis()->GetXmin();
  }
  return ymin;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Double_t GetYMax( Int_t n, TH1D** hist ){
  Double_t ymax = hist[0]->GetMaximum();
  for( Int_t i=1; i<n; i++){
    if( ymax < hist[i]->GetMaximum() ) ymax = hist[i]->GetMaximum();
  }
  return ymax;
}

Double_t GetYMax( Int_t n, TGraph** graph ){
  Double_t ymax = graph[0]->GetYaxis()->GetXmax();
  for( Int_t i=1; i<n; i++){
    if( ymax < graph[i]->GetYaxis()->GetXmax() ) ymax = graph[i]->GetYaxis()->GetXmax();
  }
  return ymax;
}

Double_t GetYMax( Int_t n, TGraphErrors** graph ){
  Double_t ymax = graph[0]->GetYaxis()->GetXmax();
  for( Int_t i=1; i<n; i++){
    if( ymax < graph[i]->GetYaxis()->GetXmax() ) ymax = graph[i]->GetYaxis()->GetXmax();
  }
  return ymax;
}


TH2D* Waku( Int_t n, TH1D** hist, const char* xlabel, Int_t ymin_flag, Double_t space ){
  return Waku(n, hist, xlabel, xlabel, xlabel, ymin_flag, space );
}
TH2D* Waku( Int_t n, TGraph** graph, const char* xlabel, Int_t ymin_flag, Double_t space ){
  return Waku(n, graph, xlabel, xlabel, xlabel, ymin_flag, space );
}
TH2D* Waku( Int_t n, TGraphErrors** graph, const char* xlabel, Int_t ymin_flag, Double_t space ){
  return Waku(n, graph, xlabel, xlabel, xlabel, ymin_flag, space );
}
