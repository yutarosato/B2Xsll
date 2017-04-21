#include "MCut.h"

MCut::MCut()
{
}
MCut::MCut( const MCut& a ){
  fTag    = a.GetTag();
  
  fName   = a.GetName();
  fFunc   = a.GetFunc();
  fType   = a.GetType();
  fFlag   = a.GetFlag();
  fVal[0] = a.GetVal1();
  fVal[1] = a.GetOffset();
  fVal[2] = a.GetVal2();
  
  fCombine = a.GetCombine();
  
  fName_second   = a.GetName_second();
  fFunc_second   = a.GetFunc_second();
  fType_second   = a.GetType_second();
  fFlag_second   = a.GetFlag_second();
  fVal_second[0] = a.GetVal1_second();
  fVal_second[1] = a.GetOffset_second();
  fVal_second[2] = a.GetVal2_second();
}
MCut::~MCut(){

}

MCut::MCut( Char_t* name,  Int_t tag,    Int_t type,
	    Int_t flag,    Double_t val1, Double_t offset, Double_t val2 )
  :fName(name),fTag(tag),fType(type),fFlag(flag)
{
  SetFunc( type );
  fVal[0] = val1;
  fVal[1] = offset;
  fVal[2] = val2;

  fCombine     = 0;
  
  fName_second   = "";
  fType_second   = 0;
  fFlag_second   = 0;
  fVal_second[0] = 0;
  fVal_second[1] = 0;
  fVal_second[2] = 0;
}

void MCut::SetFunc( Int_t type ){
  switch( type ){
  case 0:
    fType = 0;
    break;
  case 1:
    fType = 1;
    fFunc = make_cut_window;
    break;
  case 2:
    fType = 2;
    fFunc = make_cut_low;
    break;
  case 3:
    fType = 3;
    fFunc = make_cut_high;
    break;
  case 4:
    fType = 4;
    fFunc = make_cut_esc;
    break;
  case 5:
    fType = 5;
    fFunc = make_cut_tf;
    break;
  case 6:
    fType = 6;
    fFunc = make_cut_pm;
    break;
  case 7:
    fType = 7;
    fFunc = make_cut_window_plus;
    break;
  default:
    std::cerr << "not defined cut-funciton type(" << fName << ") --> abort!!" << std::endl;
    abort();
    break;
  }
}

void MCut::SetFunc_second( Int_t type ){
  switch( type ){
  case 0:
    fType_second = 0;
    break;
  case 1:
    fType_second = 1;
    fFunc_second = make_cut_window;
    break;
  case 2:
    fType_second = 2;
    fFunc_second = make_cut_low;
    break;
  case 3:
    fType_second = 3;
    fFunc_second = make_cut_high;
    break;
  case 4:
    fType_second = 4;
    fFunc_second = make_cut_esc;
    break;
  case 5:
    fType_second = 5;
    fFunc_second = make_cut_tf;
    break;
  case 6:
    fType_second = 6;
    fFunc_second = make_cut_pm;
    break;
  case 7:
    fType_second = 7;
    fFunc_second = make_cut_window_plus;
    break;
  default:
    std::cerr << "not defined cut-funciton type(" << fName << ") --> abort!!" << std::endl;
    abort();
    break;
  }
}

void MCut::Set( Int_t flag, Double_t val1, Double_t offset, Double_t val2 ){
  fFlag   = flag;
  fVal[0] = val1;
  fVal[1] = offset;
  fVal[2] = val2;
}

void MCut::Set_second( Int_t flag, Double_t val1, Double_t offset, Double_t val2 ){
  fFlag_second   = flag;
  fVal_second[0] = val1;
  fVal_second[1] = offset;
  fVal_second[2] = val2;
}

TCut MCut::Output ()const{
  if( fFlag_second ){
    if( fCombine == 1 ){
      //std::cout << Form( "(%s || %s)",fFunc(fName,fVal[0],fVal[1],fVal[2]),fFunc(fName_second, fVal_second[0], fVal_second[1], fVal_second[2]) ) << std::endl;
      return TCut( Form("(%s || %s)",
			fFunc       (fName,        fVal[0],        fVal[1],        fVal[2]       ),
			fFunc_second(fName_second, fVal_second[0], fVal_second[1], fVal_second[2])
			) );
    }else{
      //std::cout << Form( "(%s && %s)",fFunc(fName,fVal[0],fVal[1],fVal[2]),fFunc(fName_second, fVal_second[0], fVal_second[1], fVal_second[2]) ) << std::endl;
      return TCut( Form("(%s && %s)",
			fFunc       (fName,        fVal[0],        fVal[1],        fVal[2]       ),
			fFunc_second(fName_second, fVal_second[0], fVal_second[1], fVal_second[2])
			) );
    }
  }else{
    //std::cout << fFunc(fName, fVal[0], fVal[1], fVal[2]) << std::endl;
    return TCut( fFunc(fName, fVal[0], fVal[1], fVal[2]) );
  }
};

Char_t* MCut::Output2 ()const{
  if( fFlag_second ){
    if( fCombine == 1 ){
      return Form("(%s || %s)",
		  fFunc_second(fName,        fVal[0],        fVal[1],        fVal[2]       ),
		  fFunc_second(fName_second, fVal_second[0], fVal_second[1], fVal_second[2])
		  );
    }else{
      return Form("(%s && %s)",
		  fFunc_second(fName,        fVal[0],        fVal[1],        fVal[2]       ),
		  fFunc_second(fName_second, fVal_second[0], fVal_second[1], fVal_second[2])
		  );
    }
  }else{
    return fFunc(fName, fVal[0], fVal[1], fVal[2]);
  }
};

void MCut::Display( Int_t n ) const {
  if( n ){
    std::cout<< "               <Name> <Tag,Flag>  <Type><Offset><Val1><Val2>"
	     << " |<Combine>| "
	     << "<Name> <Type><Offset><Val1><Val2>"
	     << std::endl;
    std::cout << std::setw(20) << fName << " ("
	      << std::setw(3) << fTag  << ", "
	      << std::setw(4) << fFlag << ") : ";
  }
  if( fType ){
    std::cout << std::setw(3) << fType
	      << std::setw(7) << std::setprecision(4) << fVal[1]
	      << std::setw(8) << fVal[0]
	      << std::setw(6) << fVal[2];
  }else{
    std::cout << "(optional) " << Output();
  }
  if( fFlag_second ){
    std::cout << std::setw(8)  << fCombine;
    std::cout << std::setw(10) << fName_second;
    if( fType_second ){
      std::cout << std::setw(7) << fType_second
		<< std::setw(7) << std::setprecision(4) << fVal_second[1]
		<< std::setw(8) << fVal_second[0]
		<< std::setw(6) << fVal_second[2]
		<< std::endl;
    }else{
      std::cout << "(optional) " << Output() << std::endl;
    }
  }else{
    std::cout << "     none" << std::endl;
  }
}

Char_t* make_cut_window( TString p, Double_t low, Double_t offset, Double_t high ){
  //return Form( "(%s && %s )",make_cut_low(p.Data(),offset,low),make_cut_high(p.Data(),offset,high) );
  return Form( "(%s && %s )",make_cut_low(p.Data(),low,offset),make_cut_high(p.Data(),high,offset) );
}
Char_t* make_cut_low( TString p, Double_t low, Double_t offset, Double_t dust ){
  return Form( "( %f <= %s-%f )",low,p.Data(),offset );
}
Char_t* make_cut_high( TString p, Double_t high, Double_t offset, Double_t dust ){
  return Form( "( %f > %s-%f )",high,p.Data(),offset );
}
Char_t* make_cut_esc( TString p, Double_t low, Double_t offset, Double_t high ){
  return Form( "( %s || %s )", make_cut_high(p.Data(),low,offset), make_cut_low(p.Data(),high,offset));
}

Char_t* make_cut_tf( TString p, Double_t tr, Double_t offset, Double_t dust ){
  if( (Int_t)tr!=0 ) return Form( "(%s-%d)",p.Data(),(Int_t)offset );
  else return Form( "(%s-%d==0)",p.Data(),(Int_t)offset );
}
Char_t* make_cut_pm( TString p, Double_t pm, Double_t offset, Double_t dust ){
  if( (Int_t)pm>=0 ) return Form( "(%s-%d==1)",p.Data(),(Int_t)offset );
  else return Form( "(%s-%d==-1)",p.Data(),(Int_t)offset );
}
Char_t* make_cut_window_plus( TString p, Double_t low, Double_t offset, Double_t high ){
  //return Form( "( (%s && %s ) || (%s) )",make_cut_low(p.Data(),offset,low),make_cut_high(p.Data(),offset,high), make_cut_high(p.Data(),0,0) );
  return Form( "( (%s && %s ) || (%s) )",make_cut_low(p.Data(),low,offset),make_cut_high(p.Data(),high,offset), make_cut_high(p.Data(),0,0) );
}

MCut& MCut::operator = ( const MCut& a ){
  if( this == &a ) return *this;
  fTag    = a.GetTag();
  
  fName   = a.GetName();
  fFunc   = a.GetFunc();
  fType   = a.GetType();
  fFlag   = a.GetFlag();
  fVal[0] = a.GetVal1();
  fVal[1] = a.GetOffset();
  fVal[2] = a.GetVal2();

  fCombine = a.GetCombine();
  
  fName_second   = a.GetName_second();
  fFunc_second   = a.GetFunc_second();
  fType_second   = a.GetType_second();
  fFlag_second   = a.GetFlag_second();
  fVal_second[0] = a.GetVal1_second();
  fVal_second[1] = a.GetOffset_second();
  fVal_second[2] = a.GetVal2_second();
  return *this;
}
