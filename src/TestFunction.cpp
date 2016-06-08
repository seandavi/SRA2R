#include <Rcpp.h>
#include "read.h"
#include <Rdefines.h>
#include <string>
#include <vector>
#include <iostream>
using namespace Rcpp;
using namespace std;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


// http://stackoverflow.com/questions/12975341/to-string-is-not-a-member-of-std-says-so-g
// namespace patch
// {
// template < typename T > std::string to_string( const T& n )
// {
//   std::ostringstream stm ;
//   stm << n ;
//   return stm.str() ;
// }
// }

std::string to_string(int n) {
  std::ostringstream stm;
  stm << n;
  return stm.str();
}

// [[Rcpp::export]]
CharacterVector getAccessions(std::string prefix = "SRR") {
  vector<string> accs;
  long tester;
  for (int i = 100000; i < 100010; i++) {
    try {
      string acc = prefix;
      acc += to_string(i);
      if (getFastqCount(acc, false) > 0) {
        accs.push_back(acc);
      }
    } catch (...) {
      continue;
    }
  }
  int n = accs.size();
  CharacterVector returnVector(n);
  for (int i = 0; i < n; i++) {
    returnVector[i] = accs.at(i);
  }
  return returnVector;
  return CharacterVector();
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//


