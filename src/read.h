#ifndef READ_H
#define READ_H
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

long getFastqCount(Rcpp::String acc, bool forward_to_r);
DataFrame getReference(Rcpp::String acc);
std::vector<std::vector<std::string>> getRefs(std::string acc);
#endif
