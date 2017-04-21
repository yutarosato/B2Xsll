#include "MChain.h"

MChain::MChain()
  :fFile(""),fMode(""),fChange(""),fNfile(0)
{
  fChain = new TChain();
  fTree  = new TTree();
  fCut   = new MCut_array();
}
MChain::MChain( const MChain& a ){
  fFile   = a.GetFileName();
  fChain  = a.GetChain();
  fTree   = a.GetTree();
  fCut    = new MCut_array( *a.GetCut() );
  fMode   = a.GetModeName();
  fChange = a.GetChange();
  fNfile  = a.GetNfile();
}
MChain::~MChain(){
  delete fChain;
  delete fTree;
  delete fCut;
}

MChain::MChain( const Char_t* file, const Char_t* tname, MCut_array cut, const Int_t nfile, const Char_t* tail, const Int_t file_i )
  :fFile(file),fMode(""),fChange("")
{
  fChain = new TChain(tname);
  fTree  = new TTree();
  fCut   = new MCut_array(cut);
  Int_t count = 0;
  if( nfile==0 ){
    if( !(strchr(tail,'*') || strchr(tail,'?') || strchr(tail,'[') || strchr(tail,']')) ){ // for not use of wild card
      std::ifstream f( Form("%s%s",file,tail) );
      if( !f ){
	std::cerr << Form("%s%s",file,tail) << " does not exist" << std::endl;
	abort();
      }
    }
    count = fChain->Add( Form("%s%s",file,tail) );
    std::cout << file << tail
	      << " ( " << count  << " files )"
	      << std::endl;
    if( count==0 ) std::cerr << Form("%s%s",file,tail) << " does not exist" << std::endl, abort();
  }else{
    for( Int_t i=file_i; i<=file_i+nfile-1; i++){
      if( !(strchr(tail,'*') || strchr(tail,'?') || strchr(tail,'[') || strchr(tail,']')) ){
	std::ifstream f( Form("%s%d%s",file,i,tail) );
	if( !f ){
	  std::cout << Form("%s%d%s",file,i,tail) << " does not exist" << std::endl;
	  continue;
	}
      }
      count += fChain->Add( Form("%s%d%s",file,i,tail) );
    }
    std::cout << file << Form("[%d->%d]%s", file_i, file_i+nfile-1, tail)
	      << " ( " << count << " files )"
	      << std::endl;
  }
  fNfile = count;
}

void MChain::MakeTree(){
  fTree = fChain->CopyTree( fCut->Output() );
  TString sTmp = "";
  for( Int_t i=0; i<fCut->GetN(); i++ ){
    if( fCut->GetChange(i) ) sTmp += fCut->GetName(i);
  }
  fChange = sTmp;
  return;
}

MChain& MChain::operator = ( const MChain& a ){
  if( this == &a ) return *this;
  delete[] fChain;
  delete[] fTree;
  delete[] fCut;
  fFile   = a.GetFileName();
  fChain  = a.GetChain();
  fTree   = a.GetTree();
  fCut    = new MCut_array( *a.GetCut() );
  fMode   = a.GetModeName();
  fChange = a.GetChange();
  fNfile  = a.GetNfile();
  return *this;
}
