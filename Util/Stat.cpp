#include "Stat.h"

Double_t chisq2cl( const Int_t n, const Double_t chi2 ){
  Double_t SRTOPI = 0.7978846;
  Double_t UPL    = 170.0;
  Double_t ROOT2I = 0.70710678;

  Double_t prob = 0.0;
  Double_t sum,term;
  Int_t m;
  Int_t i,k;
  Double_t temp_i,temp_n;
  Double_t srty;

  if((n <= 0)||(chi2 < 0.0)){
    return prob;
  }
  if(n > 60){
    temp_n = (Double_t)n;
    srty = std::sqrt(chi2) - std::sqrt(temp_n-0.5);
    if (srty < 12.0){
      prob = 0.5*TMath::Erfc(srty);
      return prob;
    }
    return prob;
  }
  if(chi2 > UPL){
    return prob;
  }
  sum = exp( -0.5 * chi2 );
  term = sum;
  m = (Int_t)floor(n/2);

  if( 2*m == n ){
    if( m == 1 ){
      prob = sum;
      return prob;
    }else{
      for(i=2;i<m+1;i++){
	temp_i = (Double_t)i;
	term = 0.5*chi2*term/(temp_i-1.0);
	sum = sum + term;
      }
      prob = sum;
      return prob;
    }
  }else{
    srty = std::sqrt(chi2);
    prob = TMath::Erfc(ROOT2I*srty);
    if(n == 1){
      return prob;
    }
    if(n == 3){
      prob = SRTOPI*srty*sum + prob;
      return prob;
    }else{
      k = m - 1;
      for(i=1;i<k+1;i++){
	temp_i = (Double_t)i;
	term = term*chi2/(2.0*temp_i + 1.0);
	sum = sum + term;
      }
      prob = SRTOPI*srty*sum + prob;
      return prob;
    }
  }
}
