#ifndef MCut_H
#define MCut_H

#include <TROOT.h>
#include <TString.h>
#include <TCut.h>

#include <iostream>
#include <iomanip>

typedef Char_t*(*cut_func)(TString, Double_t, Double_t, Double_t);

class MCut
{
 public:
  MCut();
  MCut( const MCut& a );
  ~MCut();
  MCut( Char_t* name, Int_t tag=0,     Int_t type=4,
	Int_t flag=0, Double_t val1=0, Double_t offset=0, Double_t val2=0 );
  
 private:
  Int_t    fTag;    // cut-type() ,ex : k-id, pi-id...

  TString  fName;   // branch name
  cut_func fFunc;   // function to make cut
  Int_t    fType;   // function-type(1: window, 2:low, 3:high, 4:esc, 5:bool(0or1), 6:bool(-1or+1), 7:window + minus value, (0:optional)
  Int_t    fFlag;   // enable flag
  Double_t fVal[3]; // cut value[val1,offset,val2]

  Int_t    fCombine; // combine-type (1 -> ||, other -> &&)

  TString  fName_second;
  cut_func fFunc_second;   
  Int_t    fType_second;   
  Int_t    fFlag_second;   
  Double_t fVal_second[3];
  
 public:
  Int_t       GetTag()   const{ return fTag;         };
  
  const char* GetName()  const{ return fName.Data(); };
  cut_func    GetFunc()  const{ return fFunc;        };
  Int_t       GetType()  const{ return fType;        };
  Int_t       GetFlag()  const{ return fFlag;        };
  Double_t    GetVal1()  const{ return fVal[0];      };
  Double_t    GetOffset()const{ return fVal[1];      };
  Double_t    GetVal2()  const{ return fVal[2];      };

  Int_t       GetCombine()const{ return fCombine; };
  
  const char* GetName_second()  const{ return fName_second.Data(); };
  cut_func    GetFunc_second()  const{ return fFunc_second;        };
  Int_t       GetType_second()  const{ return fType_second;        };
  Int_t       GetFlag_second()  const{ return fFlag_second;        };
  Double_t    GetVal1_second()  const{ return fVal_second[0];      };
  Double_t    GetOffset_second()const{ return fVal_second[1];      };
  Double_t    GetVal2_second()  const{ return fVal_second[2];      };

  
  void SetTag   ( Int_t    tag      ){ fTag    = tag;             };
  
  void SetName  ( char*    name     ){ fName   = name;            };
  void SetFunc  ( cut_func func     ){ fFunc   = func; fType = 0; };
  void SetFunc  ( Int_t    type     );
  void SetFlag  ( Int_t    flag=1   ){ fFlag   = flag;            };
  void SetVal1  ( Double_t val1     ){ fVal[0] = val1;            };
  void SetOffset( Double_t offset=0 ){ fVal[1] = offset;          };
  void SetVal2  ( Double_t val2     ){ fVal[2] = val2;            };

  void SetCombine ( Int_t combine=1 ){ fCombine = combine; };
  
  void SetName_second  ( char*    name     ){ fName_second   = name;                   };
  void SetFunc_second  ( cut_func func     ){ fFunc_second   = func; fType_second = 0; };
  void SetFunc_second  ( Int_t    type     );
  void SetFlag_second  ( Int_t    flag=1   ){ fFlag_second   = flag;                   };
  void SetVal1_second  ( Double_t val1     ){ fVal_second[0] = val1;                   };
  void SetOffset_second( Double_t offset=0 ){ fVal_second[1] = offset;                 };
  void SetVal2_second  ( Double_t val2     ){ fVal_second[2] = val2;                   };

  
  void Set        ( Int_t flag=1, Double_t val1=1, Double_t offset=0, Double_t val2=0 );
  void Set_second ( Int_t flag=1, Double_t val1=1, Double_t offset=0, Double_t val2=0 );
  
  TCut    Output  ()const;
  Char_t* Output2 ()const;
  void    Display (Int_t n=1) const ; // "n" is flag for display name ?

  MCut& operator = ( const MCut& a );
};

Char_t* make_cut_window      ( TString p, Double_t low,  Double_t offset,     Double_t high      ); // [1] L <= p-o  < H
Char_t* make_cut_low         ( TString p, Double_t low,  Double_t offset=0.0, Double_t dust=0.0  ); // [2] L <= p-o
Char_t* make_cut_high        ( TString p, Double_t high, Double_t offset=0.0, Double_t dust=0.0  ); // [3]      p-o  < H
Char_t* make_cut_esc         ( TString p, Double_t low,  Double_t offset=0.0, Double_t high=-0.5 ); // [4] ( (L>p-o) || (p-o>=H) )
Char_t* make_cut_tf          ( TString p, Double_t tr=2, Double_t offset=0.0, Double_t dust=0.0  ); // [5] if (int)tr!=0, p-o!= 0. else, p-o== 0.
Char_t* make_cut_pm          ( TString p, Double_t pm=1, Double_t offset=0.0, Double_t dust=0.0  ); // [6] if (int)tr>=0, p-o==+1. else, p-o==-1.
Char_t* make_cut_window_plus ( TString p, Double_t low,  Double_t offset,     Double_t high      ); // [7] L <= p-o  < H || p < 0

#endif
