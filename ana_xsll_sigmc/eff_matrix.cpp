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

#include "draws_.h"

#include <vector>
#include <map>
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
#include <TRandom.h>

const Bool_t  flag_scale   = !true;
const Double_t scale_event = sigmc_amount;

Int_t main( Int_t argc, Char_t** argv ){
  TApplication app( "app", &argc, argv );
  Style(1);
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==6 || argc==7) ) std::cerr << "wrong input" << std::endl
					<< " Usage : ./eff_matrix (int)fl_mode_ll (int)fl_xsid (int)fl_low (int)setname_i (int)setname_f [(int)fl_appRun]" << std::endl
					<< " [fl_xsid] : 0(all), 1(K), 2(K*), 3(Xs)" << std::endl
					<< " [fl_low] : 1(efficiecy), 2(BCS false ratio)" << std::endl, abort();
  Int_t   fl_mode_ll = atoi( argv[1] ); // 1(e), 0(mu)
  Int_t   fl_xsid    = atoi( argv[2] ); // 0(all), 1(K), 2(K*), 3(Xs)
  Int_t   fl_low     = atoi( argv[3] ); // 1(efficiecy), 2(BCS false ratio)
  Int_t   setname_i  = atoi( argv[4] );
  Int_t   setname_f  = atoi( argv[5] );
  Int_t   fl_appRun  = 1;
  if( argc==7 ) fl_appRun = atoi( argv[6] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Char_t* indir_rec  = "../data/sigmc/hbk6/right/eff_matrix/log_entry_rec_calib_190dzll_ccpi0_20121108/";
  //Char_t* indir_rec  = "../data/sigmc/hbk6/right/eff_matrix/log_entry_rec_calib_190dzll_ccpi0_after_bgsup_20121214/";
  //Char_t* indir_rec  = "eff_matrix_6th_q2bin/";
  //Char_t* indir_rec  = "eff_matrix_20131106/";
  Char_t* indir_gen  = "../data/sigmc/hbk5/right/log_entry_gen_20120822/";
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  // Make maps
  std::map< Int_t, const Char_t* > mode_rec_map;
  std::map< Int_t, const Char_t* > mode_gen_map;
  for( Int_t k=0; k<nmode; k++ ){
    mode_rec_map.insert( std::pair<Int_t, const Char_t*> (mode[k], mode_name[k]) ); // 4 digits -> mode_name
    mode_gen_map.insert( std::pair<Int_t, const Char_t*> (      k, mode_name[k]) ); // id(0,1,...) -> mode_name
  }
  mode_gen_map.insert( std::pair<Int_t, const Char_t*> (nmode,    "Others"    ) );
  mode_gen_map.insert( std::pair<Int_t, const Char_t*> (nmode+1,  "Total"     ) );
  mode_gen_map.insert( std::pair<Int_t, const Char_t*> (nmode+2,  "True"      ) );
  mode_gen_map.insert( std::pair<Int_t, const Char_t*> (nmode+3,  "q2=t, fl=t") );
  mode_gen_map.insert( std::pair<Int_t, const Char_t*> (nmode+4,  "q2=t, fl=u") );
  mode_gen_map.insert( std::pair<Int_t, const Char_t*> (nmode+5,  "q2=t, fl=f") );
  mode_gen_map.insert( std::pair<Int_t, const Char_t*> (nmode+6,  "q2=f, fl=t") );
  mode_gen_map.insert( std::pair<Int_t, const Char_t*> (nmode+7,  "q2=f, fl=u") );
  mode_gen_map.insert( std::pair<Int_t, const Char_t*> (nmode+8,  "q2=f, fl=f") );
  mode_gen_map.insert( std::pair<Int_t, const Char_t*> (nmode+9,  "S-CF"      ) );
  mode_gen_map.insert( std::pair<Int_t, const Char_t*> (nmode+10, "O-CF"      ) );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  // Make hists
  TH2D** hist = new TH2D* [3];
  hist[0] = new TH2D("hist_black","hist_black",                                200,0,1,200,0,1);
  hist[1] = new TH2D("hist_red",  "hist_red",                                  200,0,1,200,0,1);
  hist[2] = new TH2D("hist_back", Form("lepton_%d_gmxs%d_set%s-%s",fl_mode_ll,fl_xsid,set_name[setname_i],set_name[setname_f]), 200,0,1,200,0,1);
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  // Set Labels
  // X-axis
  for( Int_t n=0; n<3; n++ ){
    for( Int_t k=0; k<2*nmode; k++ ){
      ((TAxis*)hist[n]->GetXaxis())->SetBinLabel( k+1, mode_name[Int_t(k/2)] );
      k++;
      ((TAxis*)hist[n]->GetXaxis())->SetBinLabel( k+1, Form("%s_true",mode_name[Int_t(k/2)]) );
    }
    ((TAxis*)hist[n]->GetXaxis())->SetBinLabel( 2*nmode+1,"TotalL" );
    ((TAxis*)hist[n]->GetXaxis())->SetBinLabel( 2*nmode+2,"Total"  );
    ((TAxis*)hist[n]->GetXaxis())->SetBinLabel( 2*nmode+3,"TotalR" );
    
    // Y-axis
    for( Int_t k=0; k<2*nmode; k++ ){
      ((TAxis*)hist[n]->GetYaxis())->SetBinLabel( 2*nmode-k+13,mode_name[Int_t(k/2)] );
      k++;
      ((TAxis*)hist[n]->GetYaxis())->SetBinLabel( 2*nmode-k+13,Form("%s_low",mode_name[Int_t(k/2)]) );
    }
    
    ((TAxis*)hist[n]->GetYaxis())->SetBinLabel(  1, "O-CF"       );
    ((TAxis*)hist[n]->GetYaxis())->SetBinLabel(  2, "S-CF"       );
    ((TAxis*)hist[n]->GetYaxis())->SetBinLabel(  3, "q2=f, fl=f" );
    ((TAxis*)hist[n]->GetYaxis())->SetBinLabel(  4, "q2=f, fl=u" );
    ((TAxis*)hist[n]->GetYaxis())->SetBinLabel(  5, "q2=f, fl=t" );
    ((TAxis*)hist[n]->GetYaxis())->SetBinLabel(  6, "q2=t, fl=f" );
    ((TAxis*)hist[n]->GetYaxis())->SetBinLabel(  7, "q2=t, fl=u" );
    ((TAxis*)hist[n]->GetYaxis())->SetBinLabel(  8, "q2=t, fl=t" );
    ((TAxis*)hist[n]->GetYaxis())->SetBinLabel(  9, "True"       );
    ((TAxis*)hist[n]->GetYaxis())->SetBinLabel( 10, "Total"      );
    ((TAxis*)hist[n]->GetYaxis())->SetBinLabel( 11, "Total_up"   );
    ((TAxis*)hist[n]->GetYaxis())->SetBinLabel( 12, "Others_low" );
    ((TAxis*)hist[n]->GetYaxis())->SetBinLabel( 13, "Others"     );
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // background color (hist[2])
  hist[2]->SetAxisRange(0,10,"Z");
  for( Int_t k1=0; k1<nmode; k1++ ){
    for( Int_t k2=0; k2<nmode; k2++ ){ 
      if( k1==k2 ){
	if( fl_low == 1 ){
	  hist[2]->Fill( mode_name[k1],                 Form("%s_low",mode_name[k2]), 9.0 ); // Red (diagonal)
	  hist[2]->Fill( Form("%s_true",mode_name[k1]), Form("%s_low",mode_name[k2]), 9.0 ); // Red (diagonal)
	}else if( fl_low == 2 ){
	  hist[2]->Fill( mode_name[k1],                 Form("%s_low",mode_name[k2]), 8.0 ); // Orange (diagonal)
	  hist[2]->Fill( Form("%s_true",mode_name[k1]), Form("%s_low",mode_name[k2]), 8.0 ); // Orange (diagonal)
	}
      }else{
	hist[2]->Fill( mode_name[k1],                 Form("%s_low",mode_name[k2]), k1%2 ? 7.5 : 7.0 ); // Yellow or Light Green (off-diagonal)
	hist[2]->Fill( Form("%s_true",mode_name[k1]), Form("%s_low",mode_name[k2]), k1%2 ? 7.5 : 7.0 ); // Yellow or Light Green (off-diagonal)
      }
    }
    hist[2]->Fill( mode_name[k1],                 "Others_low", k1%2 ? 7.5 : 7.0 ); // Yellow or Light Green (off-diagonal)
    hist[2]->Fill( Form("%s_true",mode_name[k1]), "Others_low", k1%2 ? 7.5 : 7.0 ); // Yellow or Light Green (off-diagonal)
    


    // Total Line
    hist[2]->Fill( mode_name[k1],                 "Total_up",   6.5 ); // Dark Green (Total of rec-events)
    hist[2]->Fill( Form("%s_true",mode_name[k1]), "Total_up",   6.5 ); // Dark Green (Total of rec-events)
    hist[2]->Fill( mode_name[k1],                 "Total",      6.5 ); // Dark Green (Total of rec-events)
    hist[2]->Fill( Form("%s_true",mode_name[k1]), "Total",      6.5 ); // Dark Green (Total of rec-events)
    
    hist[2]->Fill( "TotalL", mode_name[k1],                      6.5 ); // Dark Green (Total of gen-events)
    hist[2]->Fill( "TotalL", Form("%s_low",mode_name[k1]),       6.5 ); // Dark Green (Total of gen-events)
    hist[2]->Fill( "Total",  mode_name[k1],                      6.5 ); // Dark Green (Total of gen-events)
    hist[2]->Fill( "Total",  Form("%s_low",mode_name[k1]),       6.5 ); // Dark Green (Total of gen-events)
    hist[2]->Fill( "TotalR", mode_name[k1],                      6.5 ); // Dark Green (Total of gen-events)
    hist[2]->Fill( "TotalR", Form("%s_low",mode_name[k1]),       6.5 ); // Dark Green (Total of gen-events)

    // Calculation Line
    hist[2]->Fill( mode_name[k1],                 "True",       0.0              ); //
    hist[2]->Fill( Form("%s_true",mode_name[k1]), "True",       9.0              ); // Red
    hist[2]->Fill( mode_name[k1],                 "q2=t, fl=t", k1%2 ? 7.5 : 7.0 ); // Yellow or Light Green
    hist[2]->Fill( Form("%s_true",mode_name[k1]), "q2=t, fl=t", k1%2 ? 7.5 : 7.0 ); // Yellow or Light Green
    hist[2]->Fill( mode_name[k1],                 "q2=t, fl=u", 0.0              ); //
    hist[2]->Fill( Form("%s_true",mode_name[k1]), "q2=t, fl=u", 0.0              ); //
    hist[2]->Fill( mode_name[k1],                 "q2=t, fl=f", k1%2 ? 7.5 : 7.0 ); // Yellow or Light Green
    hist[2]->Fill( Form("%s_true",mode_name[k1]), "q2=t, fl=f", k1%2 ? 7.5 : 7.0 ); // Yellow or Light Green
    hist[2]->Fill( mode_name[k1],                 "q2=f, fl=t", 0.0              ); //
    hist[2]->Fill( Form("%s_true",mode_name[k1]), "q2=f, fl=t", 0.0              ); //
    hist[2]->Fill( mode_name[k1],                 "q2=f, fl=u", k1%2 ? 7.5 : 7.0 ); // Yellow or Light Green
    hist[2]->Fill( Form("%s_true",mode_name[k1]), "q2=f, fl=u", k1%2 ? 7.5 : 7.0 ); // Yellow or Light Green
    hist[2]->Fill( mode_name[k1],                 "q2=f, fl=f", 0.0              ); //
    hist[2]->Fill( Form("%s_true",mode_name[k1]), "q2=f, fl=f", 0.0              ); //
    hist[2]->Fill( mode_name[k1],                 "S-CF",       k1%2 ? 7.5 : 7.0 ); // Yellow or Light Green
    hist[2]->Fill( Form("%s_true",mode_name[k1]), "S-CF",       k1%2 ? 7.5 : 7.0 ); // Yellow or Light Green
    hist[2]->Fill( mode_name[k1],                 "O-CF",       0.0              ); //
    hist[2]->Fill( Form("%s_true",mode_name[k1]), "O-CF",       0.0              ); //
  }
  // [MIGISHITA]
  hist[2]->Fill( "TotalL", "Others",     6.5 ); // Dark Green  (Total of gen-events)
  hist[2]->Fill( "TotalL", "Others_low", 6.5 ); // Dark Green  (Total of gen-events)
  hist[2]->Fill( "TotalL", "Total_up",   6.5 ); // Dark Green  (Total of gen-events)
  hist[2]->Fill( "TotalL", "Total",      6.5 ); // Dark Green  (Total of gen-events)
  hist[2]->Fill( "TotalL", "True",       6.5 ); // Dark Green  (Total of gen-events)
  hist[2]->Fill( "TotalL", "q2=t, fl=t", 6.5 ); // Dark Green  (Total of gen-events)
  hist[2]->Fill( "TotalL", "q2=t, fl=u", 6.5 ); // Dark Green  (Total of gen-events)
  hist[2]->Fill( "TotalL", "q2=t, fl=f", 6.5 ); // Dark Green  (Total of gen-events)
  hist[2]->Fill( "TotalL", "q2=f, fl=t", 6.5 ); // Dark Green  (Total of gen-events)
  hist[2]->Fill( "TotalL", "q2=f, fl=u", 6.5 ); // Dark Green  (Total of gen-events)
  hist[2]->Fill( "TotalL", "q2=f, fl=f", 6.5 ); // Dark Green  (Total of gen-events)
  hist[2]->Fill( "TotalL", "S-CF",       6.5 ); // Dark Green  (Total of gen-events)
  hist[2]->Fill( "TotalL", "O-CF",       6.5 ); // Dark Green  (Total of gen-events)
  
  hist[2]->Fill( "Total",  "Others",     6.5 ); // Dark Green  (Total of gen-events)
  hist[2]->Fill( "Total",  "Others_low", 6.5 ); // Dark Green  (Total of gen-events)
  hist[2]->Fill( "Total",  "Total_up",   6.5 ); // Dark Green  (Total of gen-events)
  hist[2]->Fill( "Total",  "Total",      6.5 ); // Dark Green  (Total of gen-events)
  hist[2]->Fill( "Total",  "True",       6.5 ); // Dark Green  (Total of gen-events)
  hist[2]->Fill( "Total",  "q2=t, fl=t", 6.5 ); // Dark Green  (Total of gen-events)
  hist[2]->Fill( "Total",  "q2=t, fl=u", 6.5 ); // Dark Green  (Total of gen-events)
  hist[2]->Fill( "Total",  "q2=t, fl=f", 6.5 ); // Dark Green  (Total of gen-events)
  hist[2]->Fill( "Total",  "q2=f, fl=t", 6.5 ); // Dark Green  (Total of gen-events)
  hist[2]->Fill( "Total",  "q2=f, fl=u", 6.5 ); // Dark Green  (Total of gen-events)
  hist[2]->Fill( "Total",  "q2=f, fl=f", 6.5 ); // Dark Green  (Total of gen-events)
  hist[2]->Fill( "Total",  "S-CF",       6.5 ); // Dark Green  (Total of gen-events)
  hist[2]->Fill( "Total",  "O-CF",       6.5 ); // Dark Green  (Total of gen-events)
  
  hist[2]->Fill( "TotalR", "Others",     6.5 ); // Dark Green  (Total of gen-events)
  hist[2]->Fill( "TotalR", "Others_low", 6.5 ); // Dark Green  (Total of gen-events)
  hist[2]->Fill( "TotalR", "Total_up",   6.5 ); // Dark Green  (Total of gen-events)
  hist[2]->Fill( "TotalR", "Total",      9.0 ); // Red         (Total of gen-events)
  hist[2]->Fill( "TotalR", "True",       9.0 ); // Red         (Total of gen-events)
  hist[2]->Fill( "TotalR", "q2=t, fl=t", 7.0 ); // Light Green (Total of gen-events)
  hist[2]->Fill( "TotalR", "q2=t, fl=u", 0.0 ); //             (Total of gen-events)
  hist[2]->Fill( "TotalR", "q2=t, fl=f", 7.0 ); // Light Green (Total of gen-events)
  hist[2]->Fill( "TotalR", "q2=f, fl=t", 0.0 ); //             (Total of gen-events)
  hist[2]->Fill( "TotalR", "q2=f, fl=u", 7.0 ); // Light Green (Total of gen-events)
  hist[2]->Fill( "TotalR", "q2=f, fl=f", 0.0 ); //             (Total of gen-events)
  hist[2]->Fill( "TotalR", "S-CF",       7.0 ); // Light Green (Total of gen-events)
  hist[2]->Fill( "TotalR", "O-CF",       0.0 ); //             (Total of gen-events)

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Read & Fill
  ifstream fin;
  Int_t    xs, ll, id;
  Double_t entry;
  Char_t   buf[1024];
  Double_t gen_event[nmode+2]         ={0}; // [mode + other + total]
  Double_t rec_event[nmode+1][nmode+11] ={{0},{0}}; // [mode + total][mode + other + total + true + 6(q2fl) + S-CF + O-CF]
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // # of generated events
  std::cout << " [ # of generated events ]" << std::endl;
  for( Int_t n=setname_i; n<=setname_f; n++ ){
    fin.open( Form("%s/log_gen_entry_lep%d_gmxs%d_set%s.dat",indir_gen,fl_mode_ll,fl_xsid,set_name[n]) );
    std::cout << std::setw(2) << std::right << n+1 << ": "
	      << Form( "%s/log_gen_entry_lep%d_gmxs%d_set%s.dat",indir_gen,fl_mode_ll,fl_xsid,set_name[n] )
	      << std::endl;
    if(!fin){
      std::cerr << Form( "%s/log_gen_entry_lep%d_gmxs%d_set%s.dat",indir_gen,fl_mode_ll,fl_xsid,set_name[n] )
		<< "-> [ABORT] can not open data file"
		<< std::endl;
      abort();
    }
    while(!fin.eof()){
      fin.getline(buf,1024);
      if(buf[0]=='#') continue; // comment
      if(buf[0]=='*') break;    // finish
      std::istringstream sTmp(buf);
      sTmp >> ll >> id >> entry;
      if( flag_scale ) entry /= scale_event;
      Double_t tmp = 100*entry;
      tmp += 0.5;
      tmp  = (Int_t)tmp;
      tmp /= 100;
      std::map<Int_t,const Char_t*>::iterator gen = mode_gen_map.find(id);
      if( gen==mode_gen_map.end() ) std::cerr << "[ABORT] gen==mode_gen_map.end()" << std::endl, abort();
      if( gen->second=="Total" ) hist[0]->Fill( "Total", "Total_up", tmp );
      else hist[0]->Fill( "Total", gen->second, tmp );
      gen_event[id] += entry;
    }
    fin.close();
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // each mode
  std::cout << " [ # of reconstructed events ]" << std::endl;

  for( Int_t n=setname_i; n<=setname_f; n++ ){
    for( Int_t k=0; k<nmode; k++ ){
      fin.open( Form("%s/log_rec_entry_%d_lep%d_gmxs%d_set%s.dat",indir_rec, mode[k], fl_mode_ll, fl_xsid, set_name[n]) );
      std::cout << std::setw(2) << std::right << k+1 << " : "
		<< Form( "%s/log_rec_entry_%d_lep%d_gmxs%d_set%s.dat", indir_rec, mode[k],fl_mode_ll, fl_xsid, set_name[n] )
		<< std::endl;
      if(!fin){
	std::cerr << Form( "%s/log_rec_entry_%d_lep%d_gmxs%d_set%s.dat", indir_rec, mode[k],fl_mode_ll, fl_xsid, set_name[n] )
		  << " -> [WARNING] can not open data file"
		  << std::endl;
	continue;
      }

      while(!fin.eof()){
	fin.getline(buf,1024);
	if(buf[0]=='#') continue; // comment
	if(buf[0]=='*') break;    // finish
	std::istringstream sTmp(buf);
	sTmp >> xs >> ll >> id >> entry;
	if( flag_scale ) entry /= scale_event;
	Double_t tmp = 100*entry;
	tmp += 0.5;
	tmp  = (Int_t)tmp;
	tmp /= 100;
	std::map<Int_t,const Char_t*>::iterator rec = mode_rec_map.find(xs);
	std::map<Int_t,const Char_t*>::iterator gen = mode_gen_map.find(id);
	if( rec==mode_rec_map.end() ) abort();
	if( gen==mode_gen_map.end() ) abort();

	hist[0]->Fill( rec->second, gen->second, tmp );
	rec_event[k][id]     += entry;
	rec_event[nmode][id] += entry;
	
	if( gen->second == "True" ){ // True component
	  hist[0]->Fill( Form("%s_true",rec->second), rec->second, tmp );
	  hist[0]->Fill( rec->second,                 rec->second, tmp );
	}
	if( id>nmode ) hist[0]->Fill( "TotalL", gen->second, tmp ); // Total, True, q2fl component
	if( id<nmode && rec->second==gen->second ){ // S-CF component
	  hist[0]->Fill( rec->second, "S-CF", tmp );
	  hist[0]->Fill( "TotalL",    "S-CF", tmp );
	  rec_event[k][nmode+9] += entry;
	  rec_event[nmode][nmode+9] += entry;
	}
	if( (id<nmode && rec->second!=gen->second) || id==nmode ){ // O-CF component
	  hist[0]->Fill( rec->second, "O-CF", tmp );
	  hist[0]->Fill( "TotalL",    "O-CF", tmp );
	  rec_event[k][nmode+10] += entry;
	  rec_event[nmode][nmode+10] += entry;
	}
      }
      fin.close();
    }
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  // calculation
  for( Int_t k1=0; k1<nmode; k1++ ){
    Double_t tmp=0;
    std::map<Int_t,const Char_t*>::iterator gen = mode_gen_map.find(k1);
    for( Int_t m=0; m<9; m++ ){ // total ratio of true, q2fl and CF
      tmp = rec_event[k1][nmode+1]==0 ? 0 : 10000*Double_t(rec_event[k1][nmode+2+m])/rec_event[k1][nmode+1];
      tmp += 0.5;
      tmp  = (Int_t)tmp;
      tmp /= 100;
      if( m==0 ) hist[1]->Fill( Form("%s_true",gen->second), mode_gen_map.find(nmode+2+m)->second, tmp );
      else       hist[0]->Fill( Form("%s_true",gen->second), mode_gen_map.find(nmode+2+m)->second, tmp );
    }
    if( fl_low == 1 ){ // Efficiency
      // (Only true)
      tmp  = gen_event[k1]==0 ? 0 : 10000*Double_t(rec_event[k1][nmode+2])/gen_event[k1];
      tmp += 0.5;
      tmp  = (Int_t)tmp;
      tmp /= 100;
      hist[1]->Fill( Form("%s_true",gen->second), Form("%s_low",gen->second), tmp );
    }else if( fl_low == 2 ){ // BCS False Ratio
      for( Int_t k2=0; k2<nmode; k2++ ){
	std::map<Int_t,const Char_t*>::iterator rec = mode_gen_map.find(k2);
	if( k1==k2 ){
	  tmp  = rec_event[k1][nmode+2]==0 ? 0: 10000*(rec_event[k1][k2]-rec_event[k1][nmode+2])/rec_event[k1][nmode+2];
	}else{
	  tmp  = rec_event[k1][nmode+2]==0 ? 0: 10000*rec_event[k1][k2]/rec_event[k1][nmode+2];
	}
	tmp += 0.5;
	tmp  = (Int_t)tmp;
	tmp /= 100;
	hist[1]->Fill( rec->second, Form("%s_low",gen->second), tmp );
      }
    }
  }
  
  // Total efficiency
  Double_t tmp =  gen_event[nmode+1]==0 ? 0 : 10000*rec_event[nmode][nmode+2]/(gen_event[nmode+1]);
  tmp += 0.5;
  tmp  = (Int_t)tmp;
  tmp /= 100;
  hist[1]->Fill( "TotalR", "Total", tmp );
  for( Int_t m=0; m<9; m++ ){ // Total ratio of true, q2fl and CF
    Double_t tmp = rec_event[nmode][nmode+1]==0 ? 0 : 10000*Double_t(rec_event[nmode][nmode+2+m])/rec_event[nmode][nmode+1];
    tmp += 0.5;
    tmp  = (Int_t)tmp;
    tmp /= 100;
    if( m==0 ) hist[1]->Fill( "TotalR", mode_gen_map.find(nmode+2+m)->second, tmp );
    else       hist[0]->Fill( "TotalR", mode_gen_map.find(nmode+2+m)->second, tmp );
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Decorate
  for( Int_t n=0; n<3; n++ ){
    hist[n]->LabelsDeflate("X");
    hist[n]->LabelsDeflate("Y");
    hist[n]->LabelsOption ("v");
    hist[n]->SetMarkerSize(0.55);
    hist[n]->SetLabelSize(0.03,"XY");
    
    // Correct Labels
    // X-axis
    for( Int_t k=0; k<2*nmode; k++ ){
      ((TAxis*)hist[n]->GetXaxis())->SetBinLabel( k+1, mode_name[Int_t(k/2)] );
      k++;
      ((TAxis*)hist[n]->GetXaxis())->SetBinLabel( k+1, "" );
    }
    ((TAxis*)hist[n]->GetXaxis())->SetBinLabel( 2*nmode+1,"" );
    ((TAxis*)hist[n]->GetXaxis())->SetBinLabel( 2*nmode+3,"" );
    
    // Y-axis
    
    for( Int_t k=0; k<2*nmode; k++ ){
      ((TAxis*)hist[n]->GetYaxis())->SetBinLabel( 2*nmode-k+13,mode_name[Int_t(k/2)] );
      k++;
      ((TAxis*)hist[n]->GetYaxis())->SetBinLabel( 2*nmode-k+13,"" );
    }
    
    ((TAxis*)hist[n]->GetYaxis())->SetBinLabel(1,  "O-CF"       );
    ((TAxis*)hist[n]->GetYaxis())->SetBinLabel(2,  "S-CF"       );
    ((TAxis*)hist[n]->GetYaxis())->SetBinLabel(3,  "q2=f, fl=f" );
    ((TAxis*)hist[n]->GetYaxis())->SetBinLabel(4,  "q2=f, fl=u" );
    ((TAxis*)hist[n]->GetYaxis())->SetBinLabel(5,  "q2=f, fl=t" );
    ((TAxis*)hist[n]->GetYaxis())->SetBinLabel(6,  "q2=t, fl=f" );
    ((TAxis*)hist[n]->GetYaxis())->SetBinLabel(7,  "q2=t, fl=u" );
    ((TAxis*)hist[n]->GetYaxis())->SetBinLabel(8,  "q2=t, fl=t" );
    ((TAxis*)hist[n]->GetYaxis())->SetBinLabel(9,  "True"       );
    ((TAxis*)hist[n]->GetYaxis())->SetBinLabel(10, "Total"      );
    ((TAxis*)hist[n]->GetYaxis())->SetBinLabel(11, ""           );
    ((TAxis*)hist[n]->GetYaxis())->SetBinLabel(12, ""           );
    ((TAxis*)hist[n]->GetYaxis())->SetBinLabel(13, "Others"     );
  }
  if( fl_low == 1 ) hist[1]->SetMarkerColor(0);
  else if( fl_low == 2 ) hist[1]->SetMarkerColor(1);
  hist[2]->SetLabelOffset(0.01,"Z");
  
  // Draw
  TCanvas* c1 = Canvas( "c1","c1",1,1,1600,800 );
  c1->Draw();
  
  hist[2]->Draw("COLZ");
  hist[0]->Draw("textsame");
  hist[1]->Draw("textsame");

  // Display
  std::cout << "-----------------------------------------------------------------------------" << std::endl
	    << std::setw(20) << std::right << "[mode name]"
	    << " : "
	    << std::setw(10) << std::right << "[# of true]"
	    << std::setw(10) << std::right << "[# of total]"
	    << std::setw(10) << std::right << "[efficiency]"
	    << std::setw(10) << std::right << "[purity]"
	    << std::endl
	    << "-----------------------------------------------------------------------------" << std::endl;
  for( Int_t k=0; k<nmode; k++ ){
    std::cout << std::setw(20) << std::right << mode_name[k]
	      << " : "
	      << std::setw(10) << std::right << rec_event[k][nmode+2]
	      << " "
	      << std::setw(10) << std::right << rec_event[k][nmode+1]
	      << " ";
    if( gen_event[k] ) std::cout << std::setw(9) << std::right << Double_t(Int_t(10000*rec_event[k][nmode+2]/gen_event[k]))/100;
    else std::cout << std::setw(9) << std::right << 0;
    std::cout << "% ";
    if( rec_event[k][nmode+1] ) std::cout << std::setw(9) << std::right << Double_t(Int_t(10000*rec_event[k][nmode+2]/rec_event[k][nmode+1]))/100;
    else std::cout << std::setw(9) << std::right << 0;
    std::cout << "%"
	      << std::endl;
  }
  std::cout << "-----------------------------------------------------------------------------" << std::endl
	    << std::setw(20) << std::right << "total"
	    << " : "
	    << std::setw(10) << std::right << rec_event[nmode][nmode+2]
	    << " "
	    << std::setw(10) << std::right << rec_event[nmode][nmode+1]
	    << " "
	    << std::setw(9)  << std::right << Double_t( Int_t(10000*rec_event[nmode][nmode+2]/(gen_event[nmode+1])) )/100
	    << "% "
	    << std::setw(9)  << std::right << Double_t(Int_t(10000*rec_event[nmode][nmode+2]/rec_event[nmode][nmode+1]))/100
	    << "%"
	    << std::endl;

  c1->Update();
  c1->Print( Form("pic/matrix_eff_lep%d_gmxs%d_low%d_set%s-%s.eps", fl_mode_ll, fl_xsid, fl_low, set_name[setname_i], set_name[setname_f] ) );
  c1->Print( Form("pic/matrix_eff_lep%d_gmxs%d_low%d_set%s-%s.png", fl_mode_ll, fl_xsid, fl_low, set_name[setname_i], set_name[setname_f] ) );
  c1->Print( Form("pic/matrix_eff_lep%d_gmxs%d_low%d_set%s-%s.pdf", fl_mode_ll, fl_xsid, fl_low, set_name[setname_i], set_name[setname_f] ) );

  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();

  delete[] hist;
    
  return 0;
}
