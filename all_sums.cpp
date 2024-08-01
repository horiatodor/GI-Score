#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector all_sums(NumericVector x, NumericVector y) {
  
  //initialize the correct size vector
  NumericVector results(x.length()*y.length());
  
  //nested for loop
  for(int i=0; i < x.length(); i++){
    for(int j=0; j < y.length(); j++){
      int index = i*y.length()+j;
      results[index] = x[i]+y[j];
    }
  }
  
  //return statement
  return results;
}

