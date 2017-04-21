#include "bcs_lr_bkg.h"

Int_t main( Int_t argc, Char_t** argv ){
  // ========================================
  //               BCS criteria
  // ========================================
  // 1. Select maximum likelihood ratio value
  if( argc != 4 ){
    std::cerr << "wrong input" << std::endl
	      << "Usage : ./bcs_lr_bkg (char*)infile (char*)outdir (char*)brname" << std::endl
	      << std::endl;
    abort();
  }
  const Bool_t  fl_message = !true;
  //const Char_t* type_chg = argv[1];
  //const Char_t* type_psi = argv[2];
  //const Char_t* indir   = argv[1];
  // --------------------------------------------------------------------  
  const Char_t* tname_B = "h511";
  const Char_t* infile  = argv[1];
  const Char_t* outdir  = argv[2];
  const Char_t* brname  = argv[3];
  std::string basename  = gSystem->BaseName( infile );  
  if( basename.find(".root") == std::string::npos ){
    std::cerr << "[infile] " << std::setw(80) << infile
	      << " -> wrong file name" << std::endl;
    abort();
  }
  basename.erase( basename.rfind(".root") );
  // --------------------------------------------------------------------    
  TChain* chain_B = new TChain( tname_B );
  Int_t nfile = chain_B->Add( infile );
  // --------------------------------------------------------------------
  MCut_array* cut1 = new MCut_array( branch_table() );
  make_cut_ee_dembc( cut1, tname_B ); // only dE&Mbc cut
  cut1->Set( "Mbc", 1,  5.22, 0.0, 5.30 );
  //cut1->Set( "Mbc", 1,  5.20, 0.0, 5.30 ); // wide Mbc region
  TCut cut_ee = cut1->Output();
  MCut_array* cut2 = new MCut_array( branch_table() );
  make_cut_mm_dembc( cut2, tname_B ); // only dE&Mbc cut
  cut2->Set( "Mbc", 1,  5.22, 0.0, 5.30 );
  //cut2->Set( "Mbc", 1,  5.20, 0.0, 5.30 ); // wide Mbc region
  TCut cut_mm = cut2->Output();

  //TCut cut_5body_veto = "( rm_xs!=401 && rm_xs!=410 && rm_xs!=1301 && rm_xs!=1310 )"; // remove 5-body modes
  //TCut cut_dzll3d = "dzll3d<0.0150";
  // --------------------------------------------------------------------
  TTree* tree_B = new TTree();
  tree_B = chain_B->CopyTree( cut_ee || cut_mm );
  //tree_B = chain_B->CopyTree( cut_5body_veto && (cut_ee || cut_mm) ); // remove 5-body modes
  //tree_B = chain_B->CopyTree( cut_dzll3d && (cut_ee || cut_mm) ); // apply tight dzll3d selection
  Int_t dupli=0, sum_dupli=0;
  tree_B->Branch("dupli", &dupli, "dupli/I");
  Int_t tr_exist=0, tr_rmxs=0;
  tree_B->Branch("tr_exist", &tr_exist, "tr_exist/I" );
  tree_B->Branch("tr_rmxs",  &tr_rmxs,  "tr_rmxs/I"  );
  // --------------------------------------------------------------------
  TTree* newtree_B = tree_B->CloneTree( 0, "newtree" );
  Float_t  exprun, event, bcs;
  //Double_t lr;
  Float_t lr_ee;
  Float_t lr_mm;
  Float_t  rm_xs, rm_l;

  tree_B->SetBranchAddress( "exprun",       &exprun  );
  tree_B->SetBranchAddress( "event",        &event   );
  tree_B->SetBranchAddress( Form(brname,1), &lr_ee   );
  tree_B->SetBranchAddress( Form(brname,0), &lr_mm   );
  tree_B->SetBranchAddress( "rm_xs",        &rm_xs   );
  tree_B->SetBranchAddress( "rm_l",         &rm_l    );

  long long int n = 0;
  std::multimap<long long int, long long int> m;
  while( tree_B->GetEntry(n, 0) ) {
    m.insert( std::pair<long long int, long long int>((long long int)exprun*10000000000+(long long int)event, (long long int)n) ); // Xs0
    n++;
  }
  if( fl_message ) std::cout <<  n << " events (total) -> " << m.size() << " events" << std::endl;
  if( fl_message ) std::cout << " ************************************************************************************************************************" << std::endl;
  if( m.size()==0 ){
    std::cout << std::setw(7) << std::right << chain_B->GetEntries() << " " // total
	      << std::setw(7) << std::right << tree_B->GetEntries()  << " " // bcs-region cut
	      << std::setw(5) << std::right << 0                     << " " // select mode(xs and lepton)
	      << std::setw(5) << std::right << 0                     << " " // after BCS
	      << std::setw(5) << std::right << nfile                 << "  "
	      << infile       << "   "
	      << int( 100*((double)tree_B->GetEntries()/chain_B->GetEntries()) ) << " %(bcs-region cut) "
	      << 0 << " %(BCS)"
	      << std::endl;
    
    
    TFile* rootf = new TFile( Form("%s/%s_bcs.root",outdir,basename.c_str()), "RECREATE" );
    newtree_B->Write();
    rootf->Close();
    
    delete chain_B;
    delete cut1;
    delete cut2;
    delete tree_B;
    delete newtree_B;
    delete rootf;
    return 0;
  }

  std::multimap<double, int> m_event;
  double prev_event = 0;
  int cnt = 0;
  std::multimap<long long int, long long int>::iterator it_last = m.end();
  it_last--;

  for( std::multimap<long long int, long long int>::iterator i = m.begin(); i != m.end(); i++ ){
    if (i == m.begin()) {
      tree_B->GetEntry( i->second, 0 );
      if(      rm_l==1 ) bcs = lr_ee;
      else if( rm_l==0 ) bcs = lr_mm;
      else               bcs = ( lr_ee > lr_mm ? lr_ee : lr_mm );
      tr_exist++;
      tr_rmxs = (Int_t)rm_xs;

      m_event.insert( std::pair<double, int>(bcs, i->second) );
      prev_event = i->first;
      if( fl_message ) std::cout << "[" << std::setw(5) << std::right << i->second << "] "
				 << "exprun = "  << std::setw(8)  << std::right << exprun
				 << ", event = " << std::setw(6)  << std::right << event
				 << ", bcs = "   << std::setw(15) << std::right << bcs
				 << ", lr_ee = " << std::setw(15) << std::right << lr_ee
				 << ", lr_mm = " << std::setw(15) << std::right << lr_mm
				 << ", rm_xs = " << std::setw(5)  << std::right << rm_xs
				 << ", rm_l = "  << std::setw(5)  << std::right << rm_l
				 << std::endl;
    } else if (i->first != prev_event) {
      //tree_B->GetEntry( (m_event.begin()->second), 0 ); // RRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
      tree_B->GetEntry( select(m_event), 0 );
      if(      rm_l==1 ) bcs = lr_ee;
      else if( rm_l==0 ) bcs = lr_mm;
      else               bcs = ( lr_ee > lr_mm ? lr_ee : lr_mm );
      dupli      = m_event.size();
      sum_dupli += dupli;
      if( fl_message ) std::cout << " -------> [FILL!!!!!!] "
				 << "exprun = "     << exprun   << ", event = "   << event << ", bcs = " << bcs
				 << ", lr_ee = "    << lr_ee    << ", lr_mm = "   << lr_mm
				 << ", rm_xs = "    << rm_xs    << ", rm_l = "     << rm_l     
				 << ", dupli = "    << dupli   
				 << ", tr_exist = " << tr_exist << ", tr_rmxs = " << tr_rmxs
				 << std::endl;

      newtree_B->Fill();
      tr_exist=0;
      tr_rmxs=0;
      cnt++;
      m_event.clear();
      if( fl_message ) std::cout << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
				 << std::endl << std::endl;
      tree_B->GetEntry( i->second, 0 );
      if(      rm_l==1 ) bcs = lr_ee;
      else if( rm_l==0 ) bcs = lr_mm;
      else               bcs = ( lr_ee > lr_mm ? lr_ee : lr_mm );
      tr_exist++;
      tr_rmxs = (Int_t)rm_xs;
      
      m_event.insert( std::pair<double, int>(bcs, i->second) );
      prev_event = i->first;
      if( fl_message ) std::cout << "[" << std::setw(5) << std::right << i->second << "] "
				 << "exprun = "  << std::setw(8)  << std::right << exprun
				 << ", event = " << std::setw(6)  << std::right << event
				 << ", bcs = "   << std::setw(15) << std::right << bcs
				 << ", lr_ee = " << std::setw(15) << std::right << lr_ee
				 << ", lr_mm = " << std::setw(15) << std::right << lr_mm
				 << ", rm_xs = " << std::setw(5)  << std::right << rm_xs
				 << ", rm_l = "  << std::setw(5)  << std::right << rm_l
				 << std::endl;
    } else {
      tree_B->GetEntry(i->second, 0);
      if(      rm_l==1 ) bcs = lr_ee;
      else if( rm_l==0 ) bcs = lr_mm;
      else               bcs = ( lr_ee > lr_mm ? lr_ee : lr_mm );
      tr_exist++;
      tr_rmxs = (Int_t)rm_xs;

      m_event.insert( std::pair<double, int>(bcs, i->second) );
      prev_event = i->first;
      if( fl_message ) std::cout << "[" << std::setw(5) << std::right << i->second << "] "
				 << "exprun = "  << std::setw(8)  << std::right << exprun
				 << ", event = " << std::setw(6)  << std::right << event
				 << ", bcs = "   << std::setw(15) << std::right << bcs
				 << ", lr_ee = " << std::setw(15) << std::right << lr_ee
				 << ", lr_mm = " << std::setw(15) << std::right << lr_mm
				 << ", rm_xs = " << std::setw(5)  << std::right << rm_xs
				 << ", rm_l = "  << std::setw(5)  << std::right << rm_l
				 << std::endl;
    }

    if ( i == it_last ) {
      //tree_B->GetEntry( (m_event.begin()->second), 0 ); // RRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
      tree_B->GetEntry( select(m_event), 0 );
      if(      rm_l==1 ) bcs = lr_ee;
      else if( rm_l==0 ) bcs = lr_mm;
      else               bcs = ( lr_ee > lr_mm ? lr_ee : lr_mm );
      dupli      = m_event.size();
      sum_dupli += dupli;
      if( fl_message ) std::cout << " -------> [FILL!!!!!!] "
				 << "exprun = "     << exprun   << ", event = "   << event << ", bcs = " << bcs
				 << ", lr_ee = "    << lr_ee    << ", lr_mm = "   << lr_mm
				 << ", rm_xs = "    << rm_xs    << ", rm_l = "    << rm_l
				 << ", dupli = "    << dupli   
				 << ", tr_exist = " << tr_exist << ", tr_rmxs = " << tr_rmxs
				 << std::endl;
      newtree_B->Fill();
      tr_exist=0;
      tr_rmxs=0;
      cnt++;
      if( fl_message ) std::cout << " ************************************************************************************************************************" << std::endl;
    }
  }

  if( fl_message ){
    std::cout <<  n << " events (total) -> " << m.size() << " events" << std::endl;
    for( std::multimap<long long int, long long int>::iterator i = m.begin(); i != m.end(); i++ ){
      tree_B->GetEntry( i->second, 0 );
      if(      rm_l==1 ) bcs = lr_ee;
      else if( rm_l==0 ) bcs = lr_mm;
      else               bcs = ( lr_ee > lr_mm ? lr_ee : lr_mm );
      
      std::cout << "[" << std::setw(5) << std::right << i->second << "] "
		<< "exprun = "  << std::setw(8)  << std::right << exprun
		<< ", event = " << std::setw(6)  << std::right << event
		<< ", bcs = "   << std::setw(15) << std::right << bcs
		<< ", lr_ee = " << std::setw(15) << std::right << lr_ee
		<< ", lr_mm = " << std::setw(15) << std::right << lr_mm
		<< ", rm_xs = " << std::setw(5)  << std::right << rm_xs	
		<< ", rm_l = "  << std::setw(5)  << std::right << rm_l
		<< std::endl;
      
    }
  }

  // -------------< CHECK > ---------------------------------------------------------------------------
  if( m.size() != sum_dupli )          std::cerr << infile << " -> [ABORT] Event is missing : "    << m.size() << " ?=? " << sum_dupli << std::endl, abort();
  if( newtree_B->GetEntries() != cnt ) std::cerr << infile << " -> [ABORT] Event is not filled : " << newtree_B->GetEntries() << " ?=? " << cnt << std::endl, abort();
  
  // -------------< DISPLAY > ---------------------------------------------------------------------------  
  std::cout << std::setw(7) << std::right << chain_B->GetEntries() << " " // total
	    << std::setw(7) << std::right << tree_B->GetEntries()  << " " // bcs-region cut
	    << std::setw(5) << std::right << sum_dupli             << " " // select mode(xs and lepton)
	    << std::setw(5) << std::right << cnt                   << " " // after BCS
	    << std::setw(5) << std::right << nfile                 << "  "
    	    << infile       << "   "
    	    << int( 100*((double)tree_B->GetEntries()/chain_B->GetEntries()) ) << " %(bcs-region cut) "
	    << int( 100*((double)                 cnt/            sum_dupli) ) << " %(BCS)"
	    << std::endl;

  // --------------------------------------------------------------------------------------------

  TFile* rootf = new TFile( Form("%s/%s_bcs.root",outdir,basename.c_str()), "RECREATE" );
  newtree_B->Write();
  rootf->Close();

  delete chain_B;
  delete cut1;
  delete cut2;
  delete tree_B;
  delete newtree_B;
  delete rootf;
  return 0;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Int_t select( std::multimap<double, int>& event ){
  const Bool_t fl_message  = !true;
  if( fl_message ){
    std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "event.size() = " << event.size() << std::endl;
    for( std::multimap<double, int>::iterator i = event.begin(); i != event.end(); i++ ){
      std::cout << "[" << i->second << "] bcs = " << i->first << std::endl;
    }
  }
  
  std::multimap<double, int>::iterator it_last = event.end();
  it_last--;
  double  min = event.begin()->first;
  double  max = it_last->first;
  double& var = max;
  return it_last->second;
}
