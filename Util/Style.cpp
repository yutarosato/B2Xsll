#include "Style.h"

TStyle* Style( Int_t fl ){
  

  //my TStyle
  TStyle* myStyle = new  TStyle("myStyle", "MY Style");
  //set the background color to white
  myStyle->SetFillColor(10);
  myStyle->SetFrameFillColor(10);
  myStyle->SetCanvasColor(10);
  myStyle->SetPadColor(10);
  myStyle->SetTitleFillColor(10);
  myStyle->SetStatColor(10);

  //dont put a colored frame around the plots
  myStyle->SetFrameBorderMode(0);
  myStyle->SetCanvasBorderMode(0);
  myStyle->SetPadBorderMode(0);
  myStyle->SetLegendBorderSize(1);

  //use the primary color palette
  myStyle->SetPalette(1,0);
  
  //set the default line color for a histogram to be black
  myStyle->SetHistLineColor(kBlack);
  
  //set the default line color for a fit function to be red
  myStyle->SetFuncColor(kRed);

  //make the axis labels black
  myStyle->SetLabelColor(kBlack,"xyz");

  // set the length of error bar
  myStyle->SetEndErrorSize(0);

  //set the default title color to be black
  myStyle->SetTitleColor(kBlack);
  
  //set the margins
  myStyle->SetPadBottomMargin(0.15);
  myStyle->SetPadTopMargin   (0.30);
  myStyle->SetPadRightMargin (0.15);
  myStyle->SetPadLeftMargin  (0.15);

  //set axis label and title text sizes
  //myStyle->SetLabelFont(42,"xyz");
  myStyle->SetLabelSize(0.05,"xyz");
  myStyle->SetLabelOffset(0.01,"xyz");
  //myStyle->SetTitleFont(42,"xyz");
  myStyle->SetTitleSize(0.06,"xyz");
  myStyle->SetTitleOffset(1.0,"x");
  myStyle->SetTitleOffset(1.3,"y");
  myStyle->SetTitleOffset(0.4,"z");
  //myStyle->SetStatFont(42);
  //myStyle->SetStatFontSize(0.07);
  myStyle->SetTitleBorderSize(0);
  myStyle->SetStatBorderSize(0);
  //myStyle->SetTextFont(42);

  // set label-style
  myStyle->SetStripDecimals(false);
  //TGaxis::SetMaxDigits(3);

  //set line widths
  myStyle->SetFrameLineWidth(1);
  myStyle->SetFuncWidth(1);
  myStyle->SetHistLineWidth(1);
  myStyle->SetFuncStyle(2);
  
  //set the number of divisions to show
  myStyle->SetNdivisions(506, "xy");

  //turn off xy grids
  if( fl!=1 ){
    myStyle->SetPadGridX(1);
    myStyle->SetPadGridY(1);
  }
  //set the tick mark style
  myStyle->SetPadTickX(0);
  myStyle->SetPadTickY(0);
  if( fl==1) myStyle->SetTickLength( 0, "XY" );

  //turn off stats
  myStyle->SetOptStat(0);
  myStyle->SetOptFit(1111);
  
  //marker settings
  myStyle->SetMarkerStyle(20);
  myStyle->SetMarkerSize(0.7);
  myStyle->SetLineWidth(1); 

  //done
  myStyle->cd();

  //gStyle->ls();
  gROOT->GetColor( 3)->SetRGB(0.0, 0.0, 1.0); // blue
  gROOT->GetColor( 4)->SetRGB(0.0, 0.5, 0.0); // green
  gROOT->GetColor( 5)->SetRGB(0.8, 0.0, 0.8); // purple
  gROOT->GetColor( 6)->SetRGB(1.0, 0.4, 0.1); // orange
  gROOT->GetColor( 7)->SetRGB(0.2, 0.6, 0.6); // light blue
  gROOT->GetColor( 8)->SetRGB(0.6, 0.3, 0.3); // brown
  gROOT->GetColor( 9)->SetRGB(1.0, 1.0, 0.0); // yellow
  gROOT->GetColor(10)->SetRGB(1.0, 1.0, 1.0); // white (canvas color)
  gROOT->GetColor(11)->SetRGB(0.4, 0.4, 0.4); // gray
  gROOT->GetColor(12)->SetRGB(1.0, 0.0, 1.0); // light purple
  gROOT->GetColor(13)->SetRGB(1.0, 0.5, 1.0); // light red
  gROOT->GetColor(14)->SetRGB(1.0, 0.6, 0.0); // light orange
  gROOT->GetColor(15)->SetRGB(0.4, 1.0, 1.0); // very light blue
  gROOT->GetColor(16)->SetRGB(0.0, 1.0, 0.3); // light green

  gROOT->ForceStyle();

  return myStyle;
}
