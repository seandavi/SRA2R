
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
int cigarAnalysis(std::string cigar, bool qwidth = true) {
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
    if (qwidth && (cigar[0] == 'S' || cigar[0] == 'M' || cigar[0] == 'I')) {
      sum = sum + std::stoi(num);
    } else if (!qwidth && cigar[0] == 'N') {
      sum = sum + std::stoi(num);
    }
    cigar = cigar.substr(1);
  }
  return sum;
}

// [[Rcpp::export]]
DataFrame getReferenceCount(CharacterVector accs, bool track=false) {
  vector<string> refs;
  vector<int> refCount;
  int n = accs.size();
  for (int i = 0; i < n; i++) {
    try {
      if (track) {
        cout << i << endl;
      }
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

//'
//' This returns data relavent to GAlignments objects
//' @author Rosa Choe
//' @param acc An accession or a path to an actual SRA file (with .sra suffix)
//' @param seqname Reference name for alignments
//' @param low_bound A lower bound for the locations of interest
//' @param up_bound An upper bound for the locations of interest
//' @param track A boolean value that represents whether or not to print indices while running
//' @return a Rcpp::DataFrame containing data relavent to GAlignments objects
// [[Rcpp::export]]
DataFrame cpp_getGAlignments(std::string acc, std::string seqname, int low_bound = 1, int up_bound = 0, bool track = false) {
  ReadCollection run = ncbi::NGS::openReadCollection(acc);
  CharacterVector seqnames, strand, cigar;
  IntegerVector qwidth, start, end, width, njunc;
  
  bool all = false;
  bool any = false;
  if(low_bound == 1 && up_bound == 0) {
    all = true;
  }
  if (seqname.compare("") == 0 || seqname.compare("*") == 0) {
    any = true;
  }
  AlignmentIterator ai = run.getAlignments(Alignment::primaryAlignment);
  std::vector<std::vector<std::string>> references = getRefs(acc);
  std::map<std::string,std::string> refs;
  int count = 1;
  for (int i = 0; i < references.at(0).size(); i++) {
    refs[references.at(0).at(i)] = references.at(1).at(i);
  }
  while(ai.nextAlignment()) {
    if (track) {
      cout << count << endl;
    }
    
    int curr_start = ai.getAlignmentPosition() + 1;
    int curr_end = ai.getAlignmentPosition() + ai.getAlignmentLength();
    string canon = ai.getReferenceSpec();
    string common = refs[ai.getReferenceSpec()];
    if ((any || seqname.compare(canon) == 0 || seqname.compare(common) == 0) && (all || (curr_start <= up_bound && curr_end >= low_bound))) {
      start.push_back(curr_start);
      end.push_back(curr_end);
      if (refs.find(ai.getReferenceSpec()) != refs.end()) {
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
      qwidth.push_back(cigarAnalysis(ai.getShortCigar(false).toString()));
      width.push_back(ai.getAlignmentLength());
      njunc.push_back(cigarAnalysis(ai.getShortCigar(false).toString(), false));
      
      count = count + 1;
    }
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