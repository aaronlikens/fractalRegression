#include <Rcpp.h>
using namespace Rcpp;

// count the unique elements in a vector
int countUnique2(NumericVector y) {
  NumericVector x = clone(y);
  std::sort(x.begin(), x.end());
  int count = 0;
  if (x[0] == x[1]){
    count += 1;
  }
  
  double diffx = 0.0;
  for ( int i = 1; i < x.size(); ++i){
    diffx = x[i] - x[i-1];
    if (diffx != 0){
      count +=1;
    }
  }
  
  return count;
}