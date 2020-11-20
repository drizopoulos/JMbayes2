//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
//Define field with chosen dimensions
int field_3d (const int &x) {
  field<mat> Test_field(2,2,2);
  
  //Fill it with 2x2 random matrices in each position
  for (int i=0;i<2;i++){
    for (int j=0;j<2;j++){
      for (int k=0;k<2;k++){
        Test_field(i,j,k)=randu(2,2);
      }}}
  
  //Display results
  Rcout<<Test_field<<endl;
  Rcout<<Test_field.slice(0)<<endl;
  Rcout<<Test_field.slice(1)<<endl;
  
  return 1;
}
