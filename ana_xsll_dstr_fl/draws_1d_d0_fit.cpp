#include <iostream>
#include <math.h>
#include <string.h>

#include "../Util/Style.h"
#include "../Util/Canvas.h"
#include "../Util/const.h"
#include "../Util/MChain.h"
#include "../Util/MCut.h"
#include "../Util/MCut_array.h"
#include "../Util/FitFunc.h"
#include "../Set_dstr_fl/Nominal_cut_selection.h"
#include "../Set_dstr_fl/Branch.h"

#include "draws_.h"

#include <vector>
#include <stdlib.h>
#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TString.h>
#include <TF1.h>
#include <TFile.h>
#include <TArrow.h>

Int_t main( Int_t argc, Char_t** argv ){
  using namespace gmc;
  //using namespace rd;
  using namespace d0mass;
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==5 || argc==6) ) std::cerr << "wrong input" << std::endl
					<< " Usage : ./draws_1d (int)fl_par (int)fl_mom (int)fl_cos (int)fl_adef [(int)fl_appRun]" << std::endl
					<< "[fl_par ] ; 0(K), 1(pi)"            << std::endl
					<< "[fl_adef] : 0(0.4GeV), 1(0.6GeV) " << std::endl, abort();

  const Int_t   fl_par    = atoi( argv[1] );
  const Int_t   fl_mom    = atoi( argv[2] );
  const Int_t   fl_cos    = atoi( argv[3] );
  const Int_t   fl_adef   = atoi( argv[4] );
  Int_t         fl_appRun = 1;
  if( argc==6 ) fl_appRun = atoi( argv[5] );
  if( fl_mom > n_mom[fl_adef] ) std::cerr << "[ABORT] Invalid fl_mom : " << fl_mom << std::endl, abort();
  if( fl_cos > n_cos          ) std::cerr << "[ABORT] Invalid fl_cos : " << fl_cos << std::endl, abort();

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TFile infile0 ( Form("pic_d0_plot_xdst050/d0_plot_pid0_par%d_mom%d_cos%d_adef%d.root", fl_par, fl_mom, fl_cos, fl_adef) ); // change the xdst cut value !
  TFile infile1 ( Form("pic_d0_plot_xdst050/d0_plot_pid1_par%d_mom%d_cos%d_adef%d.root", fl_par, fl_mom, fl_cos, fl_adef) ); // change the xdst cut value !
  TFile infile2 ( Form("pic_d0_plot_xdst050/d0_plot_pid2_par%d_mom%d_cos%d_adef%d.root", fl_par, fl_mom, fl_cos, fl_adef) ); // change the xdst cut value !

  const Int_t Ncategory = 5;
  const Int_t Nhist     = 3*Ncategory; // [mu-id, e-id, no-pid]
  
  const Int_t  sel_fun     = 111;
  const Int_t  sel_fun_bkg = 1;
  const Bool_t flag_save   = true; // outfile.eps and outfile.root
  Int_t        fl_err      = 0;
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TH1D**    hist     = new TH1D*  [Nhist];
  TF1**     func     = new TF1*   [Nhist/Ncategory];
  TCanvas*  c1       = Canvas( "c1","c1", 3, 3 );
  Double_t* init_var = new Double_t[n_fitfunc_par(sel_fun)];
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  for(Int_t i=0; i<Ncategory; i++ ){
    hist[i+0*Ncategory] = (TH1D*)infile0.Get( Form("hist%d",i) );
    hist[i+1*Ncategory] = (TH1D*)infile1.Get( Form("hist%d",i) );
    //hist[i+1*Ncategory] = (TH1D*)infile0.Get( Form("hist%d",i) ); // tmppppppppppp
    hist[i+2*Ncategory] = (TH1D*)infile2.Get( Form("hist%d",i) );
    std::cout << i+0*Ncategory << " : " << hist[i+0*Ncategory]->GetEntries() << std::endl;
    std::cout << i+1*Ncategory << " : " << hist[i+1*Ncategory]->GetEntries() << std::endl;
    std::cout << i+2*Ncategory << " : " << hist[i+2*Ncategory]->GetEntries() << std::endl;
  }
  if( hist[14]->Integral() < 10000 ) std::cout << "[WARNING] too small sample : " << hist[14]->Integral() << std::endl;

  for(Int_t i=0; i<Nhist; i++ ) Deco( hist[i], 3, col_fil[i%Ncategory], 1 );

  // +++++++ make fit-func ++++++++++++++++++++++++++++++++++
  Double_t slope = (hist[14]->Integral(xbin-5,xbin) - hist[14]->Integral(1,5) )/ (xmax-xmin);
  for( Int_t j=0; j<Nhist/Ncategory; j++ ){
    func[j] = new TF1( Form("pdf_func_d0_pid%d", j), make_func(sel_fun),offset+xmin_fit,offset+xmax_fit,n_fitfunc_par(sel_fun) );
    func[j]->SetLineColor(3);
    func[j]->SetLineWidth(2);
    if( sel_fun==111 ){
      Double_t area = 0;
      if     ( j==0 ) area = 0.50*hist[4+j*Ncategory]->Integral()*hist[4+j*Ncategory]->GetBinWidth(0);
      else if( j==1 ) area = 0.08*hist[4+j*Ncategory]->Integral()*hist[4+j*Ncategory]->GetBinWidth(0);
      else if( j==2 ) area = 0.85*hist[4+j*Ncategory]->Integral()*hist[4+j*Ncategory]->GetBinWidth(0);

      Double_t a = (hist[4+j*Ncategory]->Integral(xbin-5,xbin) - hist[4+j*Ncategory]->Integral(1,5) )/ (xmax-xmin);
      Double_t b = hist[4+j*Ncategory]->Integral(1,5)/5 - a*xmin;
      func[j]->SetParNames( "area",  "mean","sigmal", "sigmah",  "area_ratio", "sigma", "slope", "offset" );
      func[j]->SetParameters( area,    1.865,  0.014,     0.010,      0.86,    0.0046,   a,     b  );
      func[j]->SetParLimits( 1, 1.861,  1.868 );
      func[j]->SetParLimits( 2, 0.008,  0.020 );
      func[j]->SetParLimits( 3, 0.008,  0.020 );
      func[j]->SetParLimits( 4, 0.700,  1.000 );
      func[j]->SetParLimits( 5, 0.002,  0.008 );
      func[j]->SetParLimits( 6, -slope-300, 0     );
      //func[j]->SetParameters( area,    1.865,  0.014,     0.010,      0.86,    0.0046,   -300,     600  );
    }else if( sel_fun==110 ){
      func[j]->SetParNames( "area",  "mean","sigmal", "sigmah",  "area_ratio", "sigma" );
      func[j]->SetParameters( 7,    1.865,   0.014,     0.011,      5.5,       0.0046 );
    }else if( sel_fun==21 ){
      func[j]->SetParNames  ("area", "area_ratio", "mean", "sigma", "sigma_ratio","slope","offset");
      func[j]->SetParameters( 50,    0.85,       1.865,   0.005,    2.0,       -280,    550   );
    }else if( sel_fun==31 ){
      func[j]->SetParNames  ("area", "area_ratio1", "area_ratio2", "mean", "sigma", "sigma_ratio1", "sigma_ratio2", "slope","offset");
      func[j]->SetParameters( 50,       0.85,             0.10,    1.865,   0.005,       2.0,            3.0,        -280,    550   );
    }else{
      func_set_parameters(sel_fun, func[j], hist[Nhist-2], xbin, offset+xmin, offset+xmax);
    }

    if( j==2 ){
      for(Int_t i=0; i<n_fitfunc_par(sel_fun); i++ ) init_var[i] = func[j]->GetParameter(i);
    }
  }
  
  // +++++++ tlegend1 ++++++++++++++++++++++++++++++++++
  TLegend* legend1 = new TLegend( 0.75,0.75,0.99,0.99 );
  legend1->AddEntry( hist[3],"true-D*",                "F" );
  legend1->AddEntry( hist[2],"true-D + fake-slow-#pi", "F" );
  legend1->AddEntry( hist[1],"false-D(switch K-#pi)",  "F" );
  legend1->AddEntry( hist[0],"false-D",                "F" );
  
  
  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();
  c1->cd(2); gPad->SetLogy();
  c1->cd(5); gPad->SetLogy();
  c1->cd(8); gPad->SetLogy();

  for(Int_t i=Ncategory-2; i>=0; i-- ){
    if( i==Ncategory-2 ){
      c1->cd(1); hist[i+0*Ncategory]->Draw();
      c1->cd(2); hist[i+0*Ncategory]->Draw();
      c1->cd(4); hist[i+1*Ncategory]->Draw();
      c1->cd(5); hist[i+1*Ncategory]->Draw();
      c1->cd(7); hist[i+2*Ncategory]->Draw();
      c1->cd(8); hist[i+2*Ncategory]->Draw();
    }else{
      c1->cd(1); hist[i+0*Ncategory]->Draw( "same" );
      c1->cd(2); hist[i+0*Ncategory]->Draw( "same" );
      c1->cd(4); hist[i+1*Ncategory]->Draw( "same" );
      c1->cd(5); hist[i+1*Ncategory]->Draw( "same" );
      c1->cd(7); hist[i+2*Ncategory]->Draw( "same" );
      c1->cd(8); hist[i+2*Ncategory]->Draw( "same" );
    }
  }

  for( Int_t j=Nhist/Ncategory-1; j>=0; j-- ){
    c1->cd(3+3*j);
    hist[Ncategory*(j+1)-1]->SetMarkerStyle(1);
    Int_t fit_result = hist[Ncategory*(j+1)-1]->Fit(func[j],"R","PE0");
    Int_t cnt_n =  0; // counting variable for iteration
    Int_t cnt_v =  0; // counting variable for fitting parameter

    if( j==Nhist/Ncategory-1 ){
      while( fit_result!=0 || cnt_n==0 || cnt_v!=0 ){
	std::cout << "[ Fitting false : " << cnt_n << " : " << cnt_v << " ]" << std::endl;
	
	for(Int_t i=0; i<n_fitfunc_par(sel_fun); i++ ){
	  func[2]->ReleaseParameter(i);
	  if     ( cnt_v<0  ) func[2]->SetParameter(i, init_var[i]);
	  else if( i!=cnt_v ) func[2]->FixParameter(i, init_var[i]);
	  else                func[2]->SetParameter(i, init_var[i]);
	}
	func[2]->SetParLimits( 1,  1.861, 1.868 );
	func[2]->SetParLimits( 2,  0.008, 0.020 );
	func[2]->SetParLimits( 3,  0.008, 0.020 );
	func[2]->SetParLimits( 4,  0.700, 1.000 );
	func[2]->SetParLimits( 5,  0.002, 0.008 );
	func[2]->SetParLimits( 6, -slope-300, 0     );
	
	if( cnt_v<0 ) fit_result = hist[Ncategory*(j+1)-1]->Fit(func[j],"R", "PE0");
	else          fit_result = hist[Ncategory*(j+1)-1]->Fit(func[j],"RQ","PE0");

	if( cnt_n%n_fitfunc_par(sel_fun) != 0 ){
	  for( Int_t i=0; i<n_fitfunc_par(sel_fun); i++ ) init_var[i] = func[2]->GetParameter(i);
	}

	if( cnt_n > 20*n_fitfunc_par(sel_fun) ){
	  std::cout << "[ERROR] Fitting is false !!!" << std::endl;
	  std::cerr << "[ERROR] Fitting is false !!! : ";
	  for(Int_t i=0; i<argc; i++ ) std::cerr << argv[i] << " ";
	  std::cerr << std::endl;
	  fl_err = 1;
	  //continue;
	}

	cnt_n++;
	cnt_v++;

	if( cnt_v == n_fitfunc_par(sel_fun) ) cnt_v = -1;

      }
    }
    
    if( j==Nhist/Ncategory-1 ){
      if( sel_fun==111 ){
	func[0]->FixParameter(1, func[2]->GetParameter(1) );
	func[0]->FixParameter(2, func[2]->GetParameter(2) );
	func[0]->FixParameter(3, func[2]->GetParameter(3) );
	func[0]->FixParameter(4, func[2]->GetParameter(4) );
	func[0]->FixParameter(5, func[2]->GetParameter(5) );
	func[1]->FixParameter(1, func[2]->GetParameter(1) );
	func[1]->FixParameter(2, func[2]->GetParameter(2) );
	func[1]->FixParameter(3, func[2]->GetParameter(3) );
	func[1]->FixParameter(4, func[2]->GetParameter(4) );
	func[1]->FixParameter(5, func[2]->GetParameter(5) );
      }
    }
    std::cout << "chi2/NDF = " << func[j]->GetChisquare() << " / " << func[j]->GetNDF()
    << " = " << func[j]->GetChisquare()/func[j]->GetNDF()
    << std::endl;

    c1->cd(3*j+1); func[j]->Draw("same");
    c1->cd(3*j+2); func[j]->Draw("same");
    c1->cd(3*j+1); legend1->Draw();
  }
  // +++++++ calculation of fake rate  ++++++++++++++++++++++++++++++++++

  Double_t yields [Nhist/Ncategory];
  Double_t yieldsE[Nhist/Ncategory];
  for( Int_t j=0; j<Nhist/Ncategory; j++ ){
    func_get_integral( sel_fun, func[j], hist[0]->GetBinWidth(0), yields[j], yieldsE[j]);
    std::cout << "lep" << j << " : " << yields[j] << " +- " << yieldsE[j] << std::endl;
  }
  Double_t fake_rate [2]; // mu-id, e-id
  Double_t fake_rateE[2]; // mu-id, e-id
  for( Int_t j=0; j<2; j++ ){
    fake_rate[j]  = yields[j]/yields[2];
    fake_rateE[j] = fake_rate[j] * sqrt(yieldsE[j]*yieldsE[j]/yields[j]/yields[j] + yieldsE[2]*yieldsE[2]/yields[2]/yields[2] );
    std::cout << "fake rate : " << fake_rate[j] << " +- " << fake_rateE[j] << std::endl;
  }



  // ++++++ save fake-rate ++++++++++++++++++++
  if( !fl_err ){
    std::ofstream fout( Form("log/log_d0_fit_par%d_mom%d_cos%d_adef%d.dat",  fl_par, fl_mom, fl_cos, fl_adef) );
    for( Int_t j=0; j<3; j++ ){
      // Error check
      if( fake_rate[j] < 0 ){
	std::cout << "[ERROR] negative yields !!!!! : " << yields[j] << std::endl;
	std::cerr << "[ERROR] negative yields !!!!! : ";
	for(Int_t i=0; i<argc; i++ ) std::cerr << argv[i] << " ";
	std::cerr << std::endl;
	continue;
      }else if( fake_rateE[j] > 1 ){
	std::cout << "[ERROR] too large fake rate error !!!!! : " << fake_rateE[j] << std::endl;
	std::cerr << "[ERROR] too large fake rate error !!!!! : ";
	for(Int_t i=0; i<argc; i++ ) std::cerr << argv[i] << " ";
	std::cerr << std::endl;
	continue;
      }
      
      fout << std::setw( 3) << std::right << j
	   << std::setw( 3) << std::right << fl_par
	   << std::setw( 3) << std::right << fl_mom
	   << std::setw( 3) << std::right << fl_cos
	   << std::setw(20) << std::right << func[j]->GetChisquare()/func[j]->GetNDF();
      if( j<2 ) fout << std::setw(20) << std::right << fake_rate[j]
		     << std::setw(20) << std::right << fake_rateE[j];
      fout << std::endl;
    }
    fout.close();

    /////////===========================================
    /*
    std::ofstream fout_text( Form("log_test/log_d0_fit_par%d_mom%d_cos%d_adef%d.dat",  fl_par, fl_mom, fl_cos, fl_adef) );
    for( Int_t j=0; j<3; j++ ){
      // Error check
      if( fake_rate[j] < 0 ){
	std::cout << "[ERROR] negative yields !!!!! : " << yields[j] << std::endl;
	std::cerr << "[ERROR] negative yields !!!!! : ";
	for(Int_t i=0; i<argc; i++ ) std::cerr << argv[i] << " ";
	std::cerr << std::endl;
	continue;
      }else if( fake_rateE[j] > 1 ){
	std::cout << "[ERROR] too large fake rate error !!!!! : " << fake_rateE[j] << std::endl;
	std::cerr << "[ERROR] too large fake rate error !!!!! : ";
	for(Int_t i=0; i<argc; i++ ) std::cerr << argv[i] << " ";
	std::cerr << std::endl;
	continue;
      }
      
      fout_text << std::setw( 3) << std::right << j
		<< std::setw( 3) << std::right << fl_par
		<< std::setw( 3) << std::right << fl_mom
		<< std::setw( 3) << std::right << fl_cos
		<< std::setw(20) << std::right << func[j]->GetChisquare()/func[j]->GetNDF();
      if( j<2 ) fout_text << std::setw(20) << std::right << (hist[5*j+4]->GetEntries() - hist[5*j+2]->GetEntries())/(hist[5*2+4]->GetEntries() - hist[5*2+2]->GetEntries())
			  << std::setw(20) << std::right << fake_rateE[j];
      //if( j<2 ) fout_text << std::setw(20) << std::right << fake_rate[j]
      //<< std::setw(20) << std::right << fake_rateE[j];
      fout_text << std::endl;
    }
    fout_text.close();
    */
    /////////===========================================
  }
  
  // +++++++ bkg function ++++++++++++++++++++++++++++++++++
  ///*
  TF1** func_bkg = new TF1*[Nhist/Ncategory];
  for( Int_t j=0; j<Nhist/Ncategory; j++ ){
    func_bkg[j] = new TF1( "pdf_func_d0_bkg", make_func(sel_fun_bkg),offset+xmin_fit,offset+xmax_fit,n_fitfunc_par(sel_fun_bkg) );
    func_bkg[j]->SetLineColor(1);
    if( sel_fun == 111 ){
      func_bkg[j]->SetParameter( 0, func[j]->GetParameter(6) );
      func_bkg[j]->SetParameter( 1, func[j]->GetParameter(7) );
    }else if( sel_fun == 31 ){
      func_bkg[j]->SetParameter( 0, func[j]->GetParameter(7) );
      func_bkg[j]->SetParameter( 1, func[j]->GetParameter(8) );
    }
    c1->cd(3*j+1); func_bkg[j]->Draw("same");
    c1->cd(3*j+2); func_bkg[j]->Draw("same");
  }
  //*/
  // ++++++++ save +++++++++++++++++++++++++++++++++++++++++++++  
  c1->Update();
  if( flag_save ){
    c1->Print(     Form("pic/d0_fit_par%d_mom%d_cos%d_adef%d.eps",  fl_par, fl_mom, fl_cos, fl_adef) );
    TFile outfile( Form("pic/d0_fit_par%d_mom%d_cos%d_adef%d.root", fl_par, fl_mom, fl_cos, fl_adef), "RECREATE" );
    for( Int_t i=0; i<Nhist;           i++ ) hist[i]->Write();
    for( Int_t j=0; j<Nhist/Ncategory; j++ ) func[j]->Write();
    c1->Write();
    
    outfile.Close();
  }
  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();

  delete[] hist;
  delete[] func;
  delete   c1;
  delete   legend1;
    
  return 0;
}

