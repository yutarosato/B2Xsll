#ifndef MCut_array_H
#define MCut_array_H

#include <TROOT.h>
#include <TCut.h>
#include <TString.h>

#include "MCut.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <string.h>
#include <vector>

typedef Char_t*(*cut_func)(TString, Double_t, Double_t, Double_t);

class MCut_array
{
 public:
  MCut_array();
  MCut_array( const MCut_array& a );
  ~MCut_array();
  MCut_array( Int_t n );

 private:
  
  Int_t  fN;
  MCut*  fCut;
  Int_t* fChange;
  
 public:
  MCut&        Get      ( Int_t i )const{ return fCut[i];             };
  Int_t        GetN     (         )const{ return fN;                  };
  MCut*        GetCut   (         )const{ return fCut;                };
  MCut         GetCut   ( Int_t i )const{ return fCut[i];             };
  Int_t        GetChange( Int_t i )const{ return fChange[i];          };
  
  Int_t        GetTag   ( Int_t i )const{ return fCut[i].GetTag();    };
  
  const char*  GetName  ( Int_t i )const{ return fCut[i].GetName();   };
  cut_func     GetFunc  ( Int_t i )const{ return fCut[i].GetFunc();   };
  Int_t        GetType  ( Int_t i )const{ return fCut[i].GetType();   };
  Int_t        GetFlag  ( Int_t i )const{ return fCut[i].GetFlag();   };
  Double_t     GetOffset( Int_t i )const{ return fCut[i].GetOffset(); };
  Double_t     GetVal1  ( Int_t i )const{ return fCut[i].GetVal1();   };
  Double_t     GetVal2  ( Int_t i )const{ return fCut[i].GetVal2();   };

  Int_t        GetCombine ( Int_t i )const{ return fCut[i].GetCombine(); };

  const char*  GetName_second  ( Int_t i )const{ return fCut[i].GetName_second();   };
  cut_func     GetFunc_second  ( Int_t i )const{ return fCut[i].GetFunc_second();   };
  Int_t        GetType_second  ( Int_t i )const{ return fCut[i].GetType_second();   };
  Int_t        GetFlag_second  ( Int_t i )const{ return fCut[i].GetFlag_second();   };
  Double_t     GetOffset_second( Int_t i )const{ return fCut[i].GetOffset_second(); };
  Double_t     GetVal1_second  ( Int_t i )const{ return fCut[i].GetVal1_second();   };
  Double_t     GetVal2_second  ( Int_t i )const{ return fCut[i].GetVal2_second();   };

  void Init( Int_t   i,    Char_t* name, Int_t tag,       Int_t    type );
  void Init_second( Int_t i, Char_t* name, Int_t type ){ fCut[i].SetName_second(name); fCut[i].SetFunc_second(type); };
  void ClearChange();
  void Set ( Int_t   tag,  Int_t flag,   Double_t val1=1, Double_t offset=0, Double_t val2=0 );
  void Set ( Char_t* name, Int_t flag,   Double_t val1=1, Double_t offset=0, Double_t val2=0 );
  void Set_second ( Char_t* name, Int_t flag, Int_t combine, Double_t val1=1, Double_t offset=0, Double_t val2=0 );

  void SetCombine( Int_t   tag,  Int_t combine );
  void SetCombine( Char_t* name, Int_t combine );
  
  void SetFunc( Int_t   tag,  Int_t type );
  void SetFunc( Char_t* name, Int_t type );
  void SetFunc_second( Int_t   tag,  Int_t type );
  void SetFunc_second( Char_t* name, Int_t type );
  void AllOff();
  
  TCut          Output()  const;
  const Char_t* Output2() const;
  void          Display( Int_t all=1 ) const;
  
  MCut_array& operator = ( const MCut_array& a );
};

#endif
