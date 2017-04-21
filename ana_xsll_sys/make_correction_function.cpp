#include <iostream>
#include <math.h>

#include "../Util/Style.h"
#include "../Util/Canvas.h"
#include "../Util/const.h"
#include "../Util/MChain.h"
#include "../Util/MCut.h"
#include "../Util/MCut_array.h"
#include "../Util/FitFunc.h"
#include "../Set/Nominal_cut_selection.h"
#include "../Set/Branch.h"
#include "../Set/makeCut.h"

#include "draws_.h"

#include <TApplication.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1D.h>
#include <TFile.h>

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++;

Int_t main( Int_t argc, Char_t** argv ){
  TApplication app( "app", &argc, argv );
  Style();
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++;
  Int_t fl_appRun  = 1;
  if( !(argc==2 || argc==3) ) std::cerr << "Wrong input" << std::endl
					<< " Usage : ./draw (int)fl_dat [(int)fl_appRun]"
					<< std::endl, abort();
  const Int_t   fl_dat   = atoi( argv[1] );
  if( argc==3 ) fl_appRun = atoi( argv[2] );
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++;
  const Int_t Ndat = 19;
  const Char_t* infile[Ndat]   = {
    "correction_table",
    "correction_table_fm200",
    "correction_table_fm480",
    "correction_table_mb465",
    "correction_table_mb495",
    "correction_table_transition10",
    "correction_table_transition12",
    "correction_table_hadronization",
    "correction_table_fraction_kll_p",
    "correction_table_fraction_kll_m",
    "correction_table_fraction_kstrll_p",
    "correction_table_fraction_kstrll_m",
    "correction_table_fraction_xsll_p",
    "correction_table_fraction_xsll_m",
    "correction_table_test_kll",
    "correction_table_test_kstrll",
    "correction_table_test_xsll",
    "correction_table_lambdaone362",
    "correction_table_lambdaone429",
  };
  const Int_t   nbin       = 6;
  const Bool_t  fl_message = !true;
  const Bool_t  fl_save    = true;
  const Bool_t  fl_waku    = !true; // default(false)
  const Bool_t  fl_error   = true; // default(true )
  const Bool_t  fl_fit     = true; // default(true )

  //++++++++++++++++++++++++++++++
  TGraphErrors**** g  = new TGraphErrors***[2]; // [ee,mm][norm,flip][q2-bin]
  TF1****          f1 = new TF1***         [2]; // [ee,mm][norm,flip][q2-bin] // ideal line (straight)
  TF1****          f2 = new TF1***         [2]; // [ee,mm][norm,flip][q2-bin] // correction function (ax)
  TF1****          f3 = new TF1***         [2]; // [ee,mm][norm,flip][q2-bin] // correction function (ax^3+bx )
  for( Int_t i=0; i<2; i++ ){
    g[i]  = new TGraphErrors**[2];
    f1[i] = new TF1**         [2];
    f2[i] = new TF1**         [2];
    f3[i] = new TF1**         [2];
    for( Int_t j=0; j<2; j++ ){
      g [i][j] = new TGraphErrors*[nbin];
      f1[i][j] = new TF1*         [nbin];
      f2[i][j] = new TF1*         [nbin];
      f3[i][j] = new TF1*         [nbin];
      for( Int_t k=0; k<nbin; k++ ){
	g [i][j][k] = new TGraphErrors();
	f1[i][j][k] = new TF1( Form("%d_%d_%d",i,j,k), "x", -0.5, 0.5 );
	f1[i][j][k]->SetLineColor(2);	
	f2[i][j][k] = new TF1( Form("correction_line_%d_%d_%d",i,j,k+1), "[0]*x", -0.5, 0.5 );
	f2[i][j][k]->SetParameter(0, 1.0);
	f2[i][j][k]->SetParameter(1, 0.0);
	f2[i][j][k]->SetLineColor(3);	
	f3[i][j][k] = new TF1( Form("correction_cubic_%d_%d_%d",i,j,k+1), "[0]*x+[1]*x*x*x", -0.5, 0.5 );
	f3[i][j][k]->SetParameter(0, 1.0);
	f3[i][j][k]->SetParameter(1, 0.0);
	f3[i][j][k]->SetLineColor(4);	
      }
    }
  }
  Int_t cnt[2][2][nbin] = {0};
  //++++++++++++++++++++++++++++++
  Char_t   buf[1024]; 
  ifstream fin;

  fin.open( Form("%s.dat", infile[fl_dat]) );
  while(!fin.eof()){
    fin.getline(buf,1024);

    if(buf[0]=='#') continue; // comment
    if(buf[0]=='*') break;    // finish
    std::istringstream sTmp(buf);
    Int_t    rm_l, fl_A7, fl_A9, fl_A10, fl_q2;
    Double_t afb_meas, afb_measE, afb_gen, afb_genE;
    sTmp >> rm_l >> fl_A7 >> fl_A9 >> fl_A10 >> fl_q2 >> afb_meas >> afb_measE >> afb_gen >> afb_genE;
    if( fl_message ) std::cout << std::setw( 5) << std::right << rm_l
			       << std::setw( 5) << std::right << fl_A7
			       << std::setw( 5) << std::right << fl_A9
			       << std::setw( 5) << std::right << fl_A10
			       << std::setw( 5) << std::right << fl_q2
			       << std::setw(15) << std::right << afb_meas
			       << std::setw(15) << std::right << afb_measE
			       << std::setw(15) << std::right << afb_gen
			       << std::setw(15) << std::right << afb_genE
			       << std::endl;
    Int_t tmp_rm_l = 0;
    if     ( rm_l==1 ) tmp_rm_l = 0; // ee
    else if( rm_l==0 ) tmp_rm_l = 1; // mm
    else               std::cerr << "[ABORT] Wrong rm_l : " << rm_l << std::endl, abort();

    Int_t tmp_fl_A7 = 0;
    if     ( fl_A7== 100 ) tmp_fl_A7 = 0; // norm
    else if( fl_A7==-100 ) tmp_fl_A7 = 1; // flip
    else                   std::cerr << "[ABORT] Wrong fl_A7 : " << fl_A7 << std::endl, abort();

    g[tmp_rm_l][tmp_fl_A7][fl_q2]->SetPoint     ( cnt[tmp_rm_l][tmp_fl_A7][fl_q2], afb_meas,  afb_gen  );
    if( fl_error ) g[tmp_rm_l][tmp_fl_A7][fl_q2]->SetPointError( cnt[tmp_rm_l][tmp_fl_A7][fl_q2], afb_measE, afb_genE );
    cnt[tmp_rm_l][tmp_fl_A7][fl_q2]++;
  }
  fin.close();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  for( Int_t i=0; i<2; i++ ){
    for( Int_t j=0; j<2; j++ ){
      for( Int_t k=0; k<nbin; k++ ){
	//if( cnt[i][j][k] != 21*21) std::cerr << "Wrong number of input " << i << " " << j << " " << k << " " << cnt[i][j][k] << std::endl, abort();
	g[i][j][k]->Sort();
	g[i][j][k]->SetTitle( Form("%d q^{2} bin", k+1) );
	g[i][j][k]->GetXaxis()->SetTitle("A_{FB}^{meas.}");
	g[i][j][k]->GetYaxis()->SetTitle("A_{FB}^{gen}");
	g[i][j][k]->GetXaxis()->CenterTitle();
	g[i][j][k]->GetYaxis()->CenterTitle();
	g[i][j][k]->SetMarkerStyle(6);
      }
    }
  }

  TCanvas*** c1 = new TCanvas**[2];
  TCanvas*** c2 = new TCanvas**[2];
  TH2D****    w = new TH2D***  [2];
  for( Int_t i=0; i<2; i++ ){
    c1[i] = new TCanvas*[2];
    c2[i] = new TCanvas*[2];
    w[i]  = new TH2D**  [2];
    for( Int_t j=0; j<2; j++ ){
      if     ( i==0 && j==0 ) c1[i][j] = Canvas( "ee_norm", "ee norm", 3, 3 ), c2[i][j] = Canvas( "ee_norm_proj", "ee norm proj", 3, 3 );
      else if( i==0 && j==1 ) c1[i][j] = Canvas( "ee_flip", "ee flip", 3, 3 ), c2[i][j] = Canvas( "ee_flip_proj", "ee flip proj", 3, 3 );
      if     ( i==1 && j==0 ) c1[i][j] = Canvas( "mm_norm", "mm norm", 3, 3 ), c2[i][j] = Canvas( "mm_norm_proj", "mm norm proj", 3, 3 );
      else if( i==1 && j==1 ) c1[i][j] = Canvas( "mm_flip", "mm flip", 3, 3 ), c2[i][j] = Canvas( "mm_flip_proj", "mm flip proj", 3, 3 );
      c1[i][j]->Draw();
      w[i][j] = new TH2D*[nbin];
      for( Int_t k=0; k<nbin; k++ ){
	if     ( i==0 && j==0 ) w[i][j][k] = new TH2D( Form("ee_norm %d bin",k+1), Form("ee norm, %d bin", k+1), 2, -0.40, 0.40, 2, -0.40, 0.40 );
	else if( i==0 && j==1 ) w[i][j][k] = new TH2D( Form("ee_flip %d bin",k+1), Form("ee flip, %d bin", k+1), 2, -0.40, 0.40, 2, -0.40, 0.40 );
	if     ( i==1 && j==0 ) w[i][j][k] = new TH2D( Form("mm_norm %d bin",k+1), Form("mm norm, %d bin", k+1), 2, -0.40, 0.40, 2, -0.40, 0.40 );
	else if( i==1 && j==1 ) w[i][j][k] = new TH2D( Form("mm_flip %d bin",k+1), Form("mm flip, %d bin", k+1), 2, -0.40, 0.40, 2, -0.40, 0.40 );
	w[i][j][k]->GetXaxis()->CenterTitle();
	w[i][j][k]->GetYaxis()->CenterTitle();
	w[i][j][k]->SetXTitle( "A_{FB}^{meas.}" );
	w[i][j][k]->SetYTitle( "A_{FB}^{gen}"   );
      }
    }
  }

  for( Int_t i=0; i<2; i++ ){
    for( Int_t j=0; j<2; j++ ){
      for( Int_t k=0; k<nbin; k++ ){
	if     ( nbin==6 && (k==2 || k==4) ) continue;
	else if( nbin==9 && (k==3 || k==5) ) continue;
	std::cout << std::endl
		  << Form( "*****************************[i=%d, j=%d, k=%d]*************************************", i,j,k )
		  << std::endl;
	c1[i][j]->cd(k+1);
	if( fl_waku ){
	  w[i][j][k]->Draw();
	  if( fl_fit ){
	    g [i][j][k]->Fit(f2[i][j][k]);
	    //g [i][j][k]->Fit(f3[i][j][k]);
	    g [i][j][k]->Draw("Psame");
	    f2[i][j][k]->Draw("same");
	    //f3[i][j][k]->Draw("same");
	  }else g[i][j][k]->Draw("Psame");
	}else{
	  if( fl_fit ){
	    g [i][j][k]->Fit(f2[i][j][k]);
	    //g [i][j][k]->Fit(f3[i][j][k]);
	    g [i][j][k]->Draw("AP");
	    f2[i][j][k]->Draw("same");
	    //f3[i][j][k]->Draw("same");
	  }else g[i][j][k]->Draw("AP");
	}
	f1[i][j][k]->Draw("same");
      }
    }
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Estimate line width
  TH1D**** h  = new TH1D***[2];
  TF1****  fg = new TF1*** [2]; // gaussian
  const Int_t sel_fun = 10;
  for( Int_t i=0; i<2; i++ ){
    h [i] = new TH1D**[2];
    fg[i] = new TF1** [2];
    for( Int_t j=0; j<2; j++ ){
      h [i][j] = new TH1D*[nbin];
      fg[i][j] = new TF1* [nbin];
      c2[i][j]->Draw();
      for( Int_t k=0; k<nbin; k++ ){
	if     ( nbin==6 && (k==2 || k==4) ) continue;
	else if( nbin==9 && (k==3 || k==5) ) continue;
	c2[i][j]->cd(k+1);
	if     ( nbin==6 && k==0 ) h [i][j][k] = new TH1D( Form("proj_%d_%d_%d",i,j,k), Form("proj_%d_%d_%d",i,j,k), 60, -0.030, 0.030 );
	else if( nbin==6 && k==1 ) h [i][j][k] = new TH1D( Form("proj_%d_%d_%d",i,j,k), Form("proj_%d_%d_%d",i,j,k), 60, -0.020, 0.020 );
	else if( nbin==6 && k==3 ) h [i][j][k] = new TH1D( Form("proj_%d_%d_%d",i,j,k), Form("proj_%d_%d_%d",i,j,k), 60, -0.010, 0.010 );
	else if( nbin==6 && k==5 ) h [i][j][k] = new TH1D( Form("proj_%d_%d_%d",i,j,k), Form("proj_%d_%d_%d",i,j,k), 60, -0.010, 0.010 );
	Double_t* tmpX = g[i][j][k]->GetX();
	Double_t* tmpY = g[i][j][k]->GetY();
	Double_t slope = f2[i][j][k]->GetParameter(0);
	  if( fl_message ) std::cout << "[i="      << std::setw(2) << std::right << i                  << ", "
				     <<  "j="      << std::setw(2) << std::right << j                  << ", "
				     <<  "k="      << std::setw(2) << std::right << k                  << "] "
				     <<  "N="      << std::setw(3) << std::right << g[i][j][k]->GetN() << ", "
				     <<  "slope="  << std::setw(9) << std::right << slope              << std::endl;
	for( Int_t cnt=0; cnt<g[i][j][k]->GetN(); cnt++ ){
	  Double_t distance = (slope*tmpX[cnt]-tmpY[cnt])/sqrt(slope*slope+1);
	  h[i][j][k]->Fill( distance );
	  if( fl_message ) std::cout << "["            << std::setw(3) << std::right << cnt       << "] "
				     << "(x,y) = "     << std::setw(9) << std::right << tmpX[cnt] << ", "
				     <<                   std::setw(9) << std::right << tmpY[cnt] << ") -> "
				     << " distance = " << std::setw(9) << std::right << distance  << std::endl;
	}
	if     ( nbin==6 && k==0 ) fg[i][j][k] = new TF1 ( Form("gaus_%d_%d_%d",i,j,k), make_func(sel_fun),-0.03, 0.03, n_fitfunc_par(sel_fun) );
	else if( nbin==6 && k==1 ) fg[i][j][k] = new TF1 ( Form("gaus_%d_%d_%d",i,j,k), make_func(sel_fun),-0.02, 0.02, n_fitfunc_par(sel_fun) );
	else if( nbin==6 && k==3 ) fg[i][j][k] = new TF1 ( Form("gaus_%d_%d_%d",i,j,k), make_func(sel_fun),-0.01, 0.01, n_fitfunc_par(sel_fun) );
	else if( nbin==6 && k==5 ) fg[i][j][k] = new TF1 ( Form("gaus_%d_%d_%d",i,j,k), make_func(sel_fun),-0.01, 0.01, n_fitfunc_par(sel_fun) );

	Double_t mean  = h[i][j][k]->GetBinCenter( h[i][j][k]->GetMaximumBin() );
	Double_t sigma = h[i][j][k]->GetRMS();
	fg[i][j][k]->SetParNames  ( "area","mean","sigma" );
	fg[i][j][k]->SetParameters(   1,   mean,   sigma  );
	h[i][j][k]->Fit( fg[i][j][k] );
      }
    }
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  
  // SAVE
  if( fl_save ){
    c1[0][0]->Print( Form("pic/%s_ee_norm.eps",      infile[fl_dat]) );
    c1[0][1]->Print( Form("pic/%s_ee_flip.eps",      infile[fl_dat]) );
    c1[1][0]->Print( Form("pic/%s_mm_norm.eps",      infile[fl_dat]) );
    c1[1][1]->Print( Form("pic/%s_mm_flip.eps",      infile[fl_dat]) );
    c2[0][0]->Print( Form("pic/%s_ee_norm_proj.eps", infile[fl_dat]) );
    c2[0][1]->Print( Form("pic/%s_ee_flip_proj.eps", infile[fl_dat]) );
    c2[1][0]->Print( Form("pic/%s_mm_norm_proj.eps", infile[fl_dat]) );
    c2[1][1]->Print( Form("pic/%s_mm_flip_proj.eps", infile[fl_dat]) );
    TFile outfile  ( Form("pic/%s.root",             infile[fl_dat]), "RECREATE" );
    
    for( Int_t i=0; i<2; i++ ){
      for( Int_t j=0; j<2; j++ ){
	for( Int_t k=0; k<nbin; k++ ){
	  f2[i][j][k]->Write();
	  //f3[i][j][k]->Write();
	}
      }
    }
    outfile.Close();
  }

  // LOG
  std::cout << "  RooRealVar*** cf_slope = new RooRealVar**[Nroohist];"                        << std::endl
	    << "  for(Int_t i=0; i<Nroohist; i++ ) cf_slope[i] = new RooRealVar*[Nbin_afb/2];" << std::endl;
  for( Int_t i=0; i<2; i++ ){
    Int_t j=0; // norm
    Int_t cnt = 0; // bin count
    for( Int_t k=0; k<nbin; k++ ){
      if     ( nbin==6 && (k==2 || k==4) ) continue;
      else if( nbin==9 && (k==3 || k==5) ) continue;
      std::cout << "    cf_slope["
		<< std::setw( 2) << std::right << i   << "]["
		<< std::setw( 2) << std::right << cnt << "] = new RooRealVar( "
		<< std::setw(10) << std::right << Form( "\"cf_%d_%d\", ", i, cnt )
		<< std::setw(10) << std::right << Form( "\"cf_%d_%d\", ", i, cnt )
		<< std::setw(10) << std::right << f2[i][j][k]->GetParameter(0) << " );"
		<< std::endl;
      cnt++;
    }
  }
  
  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();
  
  return 0;
}
