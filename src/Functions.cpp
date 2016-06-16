
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
  if (cigar.length() == 0) {
    return 0;
  }
  string num;
  while (cigar[0] >= 48 && cigar[0] <= 57) {
    num = num + cigar[0];
    cigar = cigar.substr(1);
  }
  if (cigar[0] == 'S' || cigar[0] == 'M' || cigar[0] == 'I') {
    return std::stoi(num) + cigarToQueryWidth(cigar.substr(1));
  } else { //'D'
    return cigarToQueryWidth(cigar);
  }
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
// [[Rcpp::export]]
Rcpp::DataFrame getGAlignmentsData(std::string acc) {
  ReadCollection run = ncbi::NGS::openReadCollection(acc);
  CharacterVector seqnames, strand, cigar;
  IntegerVector qwidth, start, end, width, njunc;
  
  AlignmentIterator ai = run.getAlignments(Alignment::primaryAlignment);
  int counter = 1;
  while(ai.nextAlignment()) {
    cout << counter << endl;
    seqnames.push_back(ai.getReferenceSpec());
    if (ai.getIsReversedOrientation()) {
      strand.push_back("-");
    } else {
      strand.push_back("+");
    }
    cigar.push_back(ai.getShortCigar(false).toString()); //test true or false; short or long
    qwidth.push_back(cigarToQueryWidth(ai.getShortCigar(false).toString()));
    start.push_back(ai.getAlignmentPosition());
    end.push_back(ai.getAlignmentPosition() + ai.getAlignmentLength() - 1);
    width.push_back(ai.getAlignmentLength());
    njunc.push_back(0); //idk what this does
    counter = counter + 1;
  }
  return DataFrame::create(
    _["seqnames"] = seqnames,
    _["strand"] = strand,
    _["cigar"] = cigar,
    _["qwidth"] = qwidth,
    _["start"] = start,
    _["end"] = end,
    _["width"] = width,
    _["njunc"] = njunc
  );
}

