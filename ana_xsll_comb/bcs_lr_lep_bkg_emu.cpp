#include "bcs_lr_bkg.h"

Int_t main( Int_t argc, Char_t** argv ){
  // ========================================
  //               BCS criteria
  // ========================================
  if( argc != 7 ){
    std::cerr << "wrong input" << std::endl
	      << "Usage : ./bcs_lr_lep_bkg_emu (char*)expno (int)stream (char*)type (char*)indir (char*)outdir (char*)brname" << std::endl
	      << std::endl;
    abort();
  }

  const Bool_t  fl_message = !true;
  const Char_t* expno  = argv[1];
  const Int_t   stream = atoi( argv[2] );
  const Char_t* type   = argv[3];
  // --------------------------------------------------------------------  
  const Char_t* tname_B = "h511";
  const Char_t* indir   = argv[4];
  const Char_t* outdir  = argv[5];
  const Char_t* brname  = argv[6];
  // --------------------------------------------------------------------    
  Int_t fl_rd  = -1; // 1(rd)  0(gmc)
  std::stringstream sTmp;
  if( strstr(type, "rd")!=NULL ){
    fl_rd = 1; // RD
    sTmp << indir << "/RD_e0" << expno << "_*.root";
  }else{
    fl_rd = 0; // gMC
    sTmp << indir << "/gMC_" << type << "_e0" << expno << "_s0" << stream << "_*.root";
  }
  Char_t tmp_infile[255];
  strcpy( tmp_infile, (Char_t*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();
  const Char_t* infile = (const Char_t*)tmp_infile;

  TChain* chain_B = new TChain( tname_B );
  Int_t nfile = chain_B->Add( infile );

  // --------------------------------------------------------------------
  MCut_array* cut = new MCut_array( branch_table() );
  make_cut_emu( cut, tname_B );
  cut->Set( "Mbc", 1,  5.22, 0.0, 5.30 );
  //cut->Set( "Mbc", 1,  5.20, 0.0, 5.30 ); // wide Mbc region
  TCut cut_emu = cut->Output();
  // --------------------------------------------------------------------
  TTree* tree_B = new TTree();
  tree_B = chain_B->CopyTree( cut_emu );
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

    if( fl_rd ) sTmp << Form( "%s/RD_e0%s_caseB_emu_bcs.root",             outdir,       expno         ); // [RD]
    else        sTmp << Form( "%s/gMC_%s_e0%s_s0%d_caseB_emu_bcs.root",    outdir, type, expno, stream ); // [gMC]
    Char_t outfile[1024];
    strcpy( outfile, (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();

    TFile* rootf = new TFile( outfile, "RECREATE" );
    newtree_B->Write();
    rootf->Close();
    delete chain_B;
    delete cut;
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
  if( fl_rd ) sTmp << Form( "%s/RD_e0%s_caseB_emu_bcs.root",             outdir,       expno         ); // [RD]
  else        sTmp << Form( "%s/gMC_%s_e0%s_s0%d_caseB_emu_bcs.root",    outdir, type, expno, stream ); // [gMC]
  
  Char_t outfile[1024];
  strcpy( outfile, (char*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();
  
  TFile* rootf = new TFile( outfile, "RECREATE" );
  newtree_B->Write();
  rootf->Close();
  delete chain_B;
  delete cut;
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
