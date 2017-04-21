#ifndef MChain_H
#define MChain_H

#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TCut.h>
#include <TString.h>
#include "MCut_array.h"

#include <iostream>
#include <fstream>

class MChain{
 public:
  MChain();
  MChain( const MChain& a );
  ~MChain();
  MChain( const Char_t* file, const Char_t* tname, MCut_array cut, const Int_t nfile, const Char_t* tail=".root", const Int_t file_i=1 );
  
 private:
  TString     fFile;  // file name
  TChain*     fChain;
  TTree*      fTree;
  MCut_array* fCut;
  TString     fMode;
  TString     fChange;
  Int_t       fNfile; // # of files
  

 public:
  const char* GetFileName()const{ return fFile.Data();   };
  TChain*     GetChain()   const{return fChain;          };
  TTree*      GetTree()    const{ return fTree;          };
  MCut_array* GetCut()     const{return  fCut;           };
  const char* GetModeName()const{ return fMode.Data();   };
  const char* GetChange()  const{ return fChange.Data(); };
  Int_t       GetNfile()   const{ return fNfile;         };
  
  void SetTreeName( Char_t* tname ){ fChain->SetName(tname); };
  void SetModeName( Char_t* name  ){ fMode = name;           };
  void MakeTree   ();

  MChain& operator = ( const MChain& a );
  
};

#endif
