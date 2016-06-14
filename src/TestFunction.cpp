#include <Rcpp.h>
#include "read.h"
#include <Rdefines.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <ncbi-vdb/NGS.hpp>
#include <ngs-bam/ngs-bam.hpp>
#include <ngs/ErrorMsg.hpp>
#include <ngs/ReadCollection.hpp>
#include <ngs/ReadIterator.hpp>
#include <ngs/Read.hpp>

#include <ngs/Reference.hpp>
#include <ngs/Alignment.hpp>
#include <ngs/PileupIterator.hpp>

#include <math.h>
#include <iostream>
using namespace Rcpp;
using namespace std;
using namespace ngs;

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

std::string to_string(long int n) {
  std::ostringstream stm;
  stm << n;
  return stm.str();
}

// [[Rcpp::export]]
CharacterVector getAccessions(long int start = 100000, long int stop = 100010, std::string prefix = "SRR", std::string file = "SRR.txt") {
  vector<string> accs;
  long tester;
  ofstream ofs;
  ofs.open(file.c_str());
  for (long int i = start; i < stop; i++) {
    try {
      string acc = prefix;
      acc += to_string(i);
      cout << i;
      if (getFastqCount(acc, false) > 0) {
        accs.push_back(acc);
        ofs << acc << " ";
        cout << " ✓";
      }
      cout << endl;
    } catch (...) {
      continue;
    }
  }
  ofs.close();
  int n = accs.size();
  CharacterVector returnVector(n);
  for (int i = 0; i < n; i++) {
    returnVector[i] = accs.at(i);
  }
  return returnVector;
  return CharacterVector();
}

// [[Rcpp::export]]
DataFrame getReferenceCount(CharacterVector accs) {
  vector<string> refs;
  vector<int> refCount;
  int n = accs.size();
  for (int i = 0; i < n; i++) {
    try {
      cout << i << endl;
      CharacterVector c = getReference(as<string>(accs[i]))["Length"];
      if (c.size() > 0) {
        refs.push_back(as<string>(accs[i]));
        refCount.push_back(c.size());
      }
    } catch(...) {
    }
  }
  return DataFrame::create(_["Run"]= refs, _["ReferenceCount"]= refCount);
}

