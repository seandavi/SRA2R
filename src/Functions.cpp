
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

//[[Rcpp::plugins(cpp11)]]
int cigarToQueryWidth(std::string cigar) {
  int arr[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  int sum = 0;
  if (cigar.length() == 0) {
    return 0;
  }
  
  while(cigar.length() > 0) {
    string num = "";
    while (cigar[0] >= 48 && cigar[0] <= 57) {
      num = num + cigar[0];
      cigar = cigar.substr(1);
    }
    if (cigar[0] == 'S' || cigar[0] == 'M' || cigar[0] == 'I') {
      sum = sum + std::stoi(num);
    }
    cigar = cigar.substr(1);
  }
  return sum;
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
      acc += std::to_string(i);
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
// [[Rcpp::export]]
Rcpp::DataFrame getGAlignmentsData(std::string acc) {
  ReadCollection run = ncbi::NGS::openReadCollection(acc);
  CharacterVector seqnames, strand, cigar;
  IntegerVector qwidth, start, end, width, njunc;
  
  AlignmentIterator ai = run.getAlignments(Alignment::primaryAlignment);
  std::vector<std::vector<std::string>> references = getRefs(acc);
  std::map<std::string,std::string> refs;
  for (int i = 0; i < references.at(0).size(); i++) {
    refs[references.at(0).at(i)] = references.at(1).at(i);
  }
  while(ai.nextAlignment()) {
    if(refs.find(ai.getReferenceSpec()) != refs.end()) {
      seqnames.push_back(refs[ai.getReferenceSpec()]);
    } else {
      seqnames.push_back(ai.getReferenceSpec());
    }
    if (ai.getIsReversedOrientation()) {
      strand.push_back("-");
    } else {
      strand.push_back("+");
    }
    cigar.push_back(ai.getShortCigar(false).toString());
    qwidth.push_back(cigarToQueryWidth(ai.getShortCigar(false).toString()));
    start.push_back(ai.getAlignmentPosition() + 1);
    end.push_back(ai.getAlignmentPosition() + ai.getAlignmentLength());
    width.push_back(ai.getAlignmentLength());
    njunc.push_back(0); //idk what this does
  }
  DataFrame d = DataFrame::create(
    _["seqnames"] = seqnames,
    _["strand"] = strand,
    _["cigar"] = cigar,
    _["qwidth"] = qwidth,
    _["start"] = start,
    _["end"] = end,
    _["width"] = width,
    _["njunc"] = njunc
  );
  
  Environment myEnv = Environment::global_env();
  Function sort = myEnv["sortAlignments"];
  return sort(d);
}


