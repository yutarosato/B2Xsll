#include "MCut_array.h"

MCut_array::MCut_array()
{
}
MCut_array::MCut_array( const MCut_array& a ){
  fN      = a.GetN();
  fCut    = new MCut[a.GetN()];
  fChange = new Int_t[fN];
  for(Int_t i=0; i<fN; i++){
    fCut[i]    = a.GetCut(i);
    fChange[i] = a.GetChange(i);
  }
}

MCut_array::~MCut_array(){
  delete[] fCut;
  delete[] fChange;
}

MCut_array::MCut_array( Int_t n ){
  fN      = n;
  fCut    = new MCut[fN];
  fChange = new Int_t[fN];
}


void MCut_array::Init( Int_t i, Char_t* name, Int_t tag, Int_t type ){
  if( i>=fN ){ std::cerr << "Change the number of branches" << std::endl; abort(); }
  MCut tmp( name, tag, type );
  fCut[i] = tmp;
  return;
}

void MCut_array::ClearChange(){
  for(Int_t i=0; i<fN; i++) fChange[i]=0;
  return;
}

void MCut_array::Set( Int_t tag, Int_t flag, Double_t val1, Double_t offset, Double_t val2 ){
  for( Int_t i=0; i<fN; i++){
    if( fCut[i].GetTag()==tag ){
      fCut[i].Set( flag, val1, offset, val2 );
      fChange[i]++;
    }
  }
  return;
}

void MCut_array::Set( Char_t* name, Int_t flag, Double_t val1, Double_t offset, Double_t val2 ){
  Int_t tmp=0;
  for(Int_t i=0; i<fN; i++){
    if( !strcmp(name,fCut[i].GetName()) ){
      fCut[i].Set( flag, val1, offset, val2 );
      fChange[i]++;
      tmp++;
    }
  }
  if( !tmp ){ std::cerr << Form("Check the branch name(%s)", name) << std::endl; abort(); }
  return;
}

void MCut_array::Set_second( Char_t* name, Int_t flag, Int_t combine, Double_t val1, Double_t offset, Double_t val2 ){
  Int_t tmp=0;
  for(Int_t i=0; i<fN; i++){
    if( !strcmp(name,fCut[i].GetName_second()) ){
      fCut[i].Set_second( flag, val1, offset, val2 );
      fCut[i].SetCombine( combine );
      fChange[i]++;
      tmp++;
    }
  }
  if( !tmp ){ std::cerr << Form("Check the branch name(%s)", name) << std::endl; abort(); }
  return;
}

void MCut_array::SetCombine( Int_t tag, Int_t combine ){
  for( Int_t i=0; i<fN; i++){
    if( fCut[i].GetTag()==tag ){
      fCut[i].SetCombine( combine );
      fChange[i]++;
    }
  }
  return;
}

void MCut_array::SetCombine( Char_t* name, Int_t combine ){
  Int_t tmp=0;
  for(Int_t i=0; i<fN; i++){
    if( !strcmp(name,fCut[i].GetName()) ){
      fCut[i].SetCombine( combine );
      fChange[i]++;
      tmp++;
    }
  }
  if( !tmp ){ std::cerr << Form("Check the branch name(%s)", name) << std::endl; abort(); }
  return;
}

void MCut_array::SetFunc( Int_t tag, Int_t type ){
  for( Int_t i=0; i<fN; i++){
    if( fCut[i].GetTag()==tag ){
      fCut[i].SetFunc( type );
      fChange[i]++;
    }
  }
  return;
}

void MCut_array::SetFunc( Char_t* name, Int_t type ){
  Int_t tmp=0;
  for(Int_t i=0; i<fN; i++){
    if( !strcmp(name,fCut[i].GetName()) ){
      fCut[i].SetFunc( type );
      fChange[i]++;
      tmp++;
    }
  }
  if( !tmp ){ std::cerr << Form("Check the branch name(%s)", name) << std::endl; abort(); }
  return;
}

void MCut_array::SetFunc_second( Int_t tag, Int_t type ){
  for( Int_t i=0; i<fN; i++){
    if( fCut[i].GetTag()==tag ){
      fCut[i].SetFunc_second( type );
      fChange[i]++;
    }
  }
  return;
}

void MCut_array::SetFunc_second( Char_t* name, Int_t type ){
  Int_t tmp=0;
  for(Int_t i=0; i<fN; i++){
    if( !strcmp(name,fCut[i].GetName_second()) ){
      fCut[i].SetFunc_second( type );
      fChange[i]++;
      tmp++;
    }
  }
  if( !tmp ){ std::cerr << Form("Check the branch name(%s)", name) << std::endl; abort(); }
  return;
}

void MCut_array::AllOff(){
  TCut cut;
  for(Int_t i=0; i<fN; i++){
    fCut[i].SetFlag(0);
    fChange[i]++;
  }
  return;
}


TCut MCut_array::Output()const{
  TCut cut;
  for(Int_t i=0; i<fN; i++){
    if( fCut[i].GetFlag() ) cut += fCut[i].Output();
  }
  return cut;
}

const Char_t* MCut_array::Output2()const{
  std::stringstream cut;
  int fl=0;
  for(Int_t i=0; i<fN; i++){
    if( fCut[i].GetFlag() ){
      if( fl ) cut << " && ";
      cut << fCut[i].Output2();
      fl++;
    }
  }
  return cut.str().c_str();
}

void MCut_array::Display( Int_t all )const{
  std::cout << "     <No.>               <Name>"
	    <<"<Tag>      <Type><Offset><Val1><Val2>"
	    << " |<Combine>| "
	    << "<Name> <Type><Offset><Val1><Val2>"
	    << std::endl;
  
  for(Int_t i=0; i<fN; i++){
    if( fChange[i] ) std::cout << std::setw(5) << "*";
    else if( all ) std::cout << std::setw(5) << "";
    if( all || fChange[i] ){
      std::cout << std::setw(3) << i << " : "
		<< std::setw(20) << fCut[i].GetName() << "("
		<< std::setw(6) << fCut[i].GetTag()  << ") -> ";
      
      if( fCut[i].GetFlag() ){
	fCut[i].Display(0);
      }else{
	std::cout << "not applied" << std::endl;
      }
    }
  }
  return;
}

MCut_array& MCut_array::operator = ( const MCut_array& a ){
  if( this == &a ) return *this;
  delete[] fCut;
  delete[] fChange;
  fN      = a.GetN();
  fCut    = new MCut[a.GetN()];
  fChange = new Int_t[fN];
  for(Int_t i=0; i<fN; i++){
    fCut[i]    = a.GetCut(i);
    fChange[i] = a.GetChange(i);
  }
  return *this;
}

